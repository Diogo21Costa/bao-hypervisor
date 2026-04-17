#include <arch/pmu.h>
#include <cpu.h>
#include <bitmap.h>
#include <platform.h>

static size_t pmu_event_counters_num[CPU_MAX];
static BITMAP_ALLOC_ARRAY(pmu_events_bitmap, PMU_CNTR_MAX_NUM, CPU_MAX);
static irq_handler_t pmu_event_handlers[CPU_MAX][PMU_CNTR_MAX_NUM];
static irqid_t pmu_irq_ids[CPU_MAX];
static bool pmu_irq_reserved[CPU_MAX];

static inline void pmu_restore_el2_control(void)
{
    uint64_t mdcr = sysreg_mdcr_el2_read();
    mdcr &= ~((uint64_t)MDCR_EL2_HPMN_MASK);
    mdcr |= (uint64_t)MDCR_EL2_HPME | (uint64_t)PMU_N_CNTR_GIVEN;
    sysreg_mdcr_el2_write(mdcr);
}

static inline size_t pmu_cpu_index(void)
{
    if (cpu()->id >= CPU_MAX) {
        ERROR("Invalid CPU id = %lu\n", cpu()->id);
    }

    return (size_t)cpu()->id;
}

ssize_t pmu_cntr_alloc(void)
{
    size_t cpu_index = pmu_cpu_index();

    for (size_t index = PMU_N_CNTR_GIVEN; index < pmu_event_counters_num[cpu_index]; index++) {
        if (bitmap_get(pmu_events_bitmap[cpu_index], index) == BITMAP_NOT_SET) {
            bitmap_set(pmu_events_bitmap[cpu_index], index);
            return (ssize_t)index;
        }
    }

    return ERROR_NO_MORE_EVENT_COUNTERS;
}

void pmu_cntr_free(size_t counter)
{
    size_t cpu_index = pmu_cpu_index();

    if (counter >= PMU_CNTR_MAX_NUM) {
        ERROR("Invalid PMU counter index = %zu\n", counter);
    }

    bitmap_clear(pmu_events_bitmap[cpu_index], counter);
}

static void pmu_interrupt_handler(irqid_t int_id)
{
    size_t cpu_index = pmu_cpu_index();
    uint64_t overflow_status = sysreg_pmovsclr_el0_read();
    sysreg_pmovsclr_el0_write(overflow_status);

    for (size_t counter = PMU_N_CNTR_GIVEN; counter < pmu_event_counters_num[cpu_index];
         counter++) {
        if (bit64_get(overflow_status, counter) != 0U) {
            irq_handler_t handler = pmu_event_handlers[cpu_index][counter];
            if (handler != NULL) {
                handler(int_id);
            }
        }
    }
}

void pmu_enable(void)
{
    size_t cpu_index = pmu_cpu_index();

    uint64_t pmcr = sysreg_pmcr_el0_read();
    size_t counters_num = (size_t)((pmcr & PMCR_EL0_N_MASK) >> PMCR_EL0_N_POS);

    if (counters_num > PMU_CNTR_MAX_NUM) {
        counters_num = PMU_CNTR_MAX_NUM;
    }

    pmu_event_counters_num[cpu_index] = counters_num;
    bitmap_clear_consecutive(pmu_events_bitmap[cpu_index], 0, PMU_CNTR_MAX_NUM);
    for (size_t i = 0; i < PMU_CNTR_MAX_NUM; i++) {
        pmu_event_handlers[cpu_index][i] = NULL;
    }

    pmu_restore_el2_control();

    pmcr = sysreg_pmcr_el0_read();
    pmcr = bit_set(pmcr, 0);   // PMCR_EL0.E
    sysreg_pmcr_el0_write(pmcr);
    __asm__ volatile("isb");
}

void pmu_interrupt_enable(cpuid_t cpu_id)
{
    UNUSED_ARG(cpu_id);
    size_t cpu_index = pmu_cpu_index();
    pmu_restore_el2_control();

    irqid_t irq_offset = platform.arch.events.events_irq_offset;
    if (irq_offset == INVALID_IRQID) {
        ERROR("Invalid PMU interrupt offset\n");
    }

    irqid_t irq_id = (irqid_t)(irq_offset + (irqid_t)cpu_index);
    if (!pmu_irq_reserved[cpu_index]) {
        pmu_irq_ids[cpu_index] = interrupts_reserve(irq_id, pmu_interrupt_handler);
        if (pmu_irq_ids[cpu_index] == INVALID_IRQID) {
            ERROR("Failed to assign PMU interrupt id = %u\n", irq_id);
        }
        pmu_irq_reserved[cpu_index] = true;
    }

    interrupts_arch_enable(pmu_irq_ids[cpu_index], true);
}

void pmu_interrupt_disable(cpuid_t cpu_id)
{
    UNUSED_ARG(cpu_id);
    size_t cpu_index = pmu_cpu_index();

    if (pmu_irq_reserved[cpu_index]) {
        interrupts_arch_enable(pmu_irq_ids[cpu_index], false);
    }
}

void pmu_define_event_cntr_irq_callback(irq_handler_t handler, size_t counter)
{
    size_t cpu_index = pmu_cpu_index();

    if (counter >= PMU_CNTR_MAX_NUM) {
        ERROR("Invalid PMU counter index = %zu\n", counter);
    }

    pmu_event_handlers[cpu_index][counter] = handler;
}
