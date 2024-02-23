/**
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) Bao Project and Contributors. All rights reserved.
 */

#include <arch/smmuv2.h>
#include <arch/spinlock.h>
#include <bitmap.h>
#include <bit.h>
#include <arch/sysregs.h>
#include <platform.h>
#include <cpu.h>
#include <mem.h>

#define SME_MAX_NUM 128
#define CTX_MAX_NUM 128

#define SMMU_PMEVTYPER_P_OFF        (31)
#define SMMU_PMEVTYPER_U_OFF        (30)
#define SMMU_PMEVTYPER_NSP_OFF      (29)
#define SMMU_PMEVTYPER_NSU_OFF      (28)

#define SMMU_PMEVTYPER_EVENT_OFF    (0)
#define SMMU_PMEVTYPER_EVENT_LEN    (16)

#define SMMU_IRQ_ID         187

struct smmu_hw {
    volatile struct smmu_glbl_rs0_hw* glbl_rs0;
    volatile struct smmu_glbl_rs1_hw* glbl_rs1;
    volatile struct smmu_cntxt_hw* cntxt;
    volatile struct smmu_pmu_hw* pmu;
};

struct smmu_priv {
    struct smmu_hw hw;

    /* For easier book keeping */
    spinlock_t sme_lock;
    size_t sme_num;
    BITMAP_ALLOC(sme_bitmap, SME_MAX_NUM);
    BITMAP_ALLOC(grp_bitmap, SME_MAX_NUM);

    spinlock_t ctx_lock;
    size_t ctx_num;
    BITMAP_ALLOC(ctxbank_bitmap, CTX_MAX_NUM);

    spinlock_t pmu_lock;
    size_t pmu_cntrs_num;
    size_t cntr_groups_num;
};

struct smmu_priv smmu;

/**
 * Iterate stream match entries.
 *
 * @sme: starting point of the loop cursor
 */
#define smmu_for_each_sme(sme)                                                \
    for (size_t __bit = bitmap_get(smmu.sme_bitmap, sme); sme < smmu.sme_num; \
         __bit = bitmap_get(smmu.sme_bitmap, ++sme))                          \
        if (!__bit)                                                           \
            continue;                                                         \
        else

/**
 * Accessors inline functions.
 */
inline bool smmu_sme_is_group(size_t sme)
{
    return bitmap_get(smmu.grp_bitmap, sme);
}

inline size_t smmu_sme_get_ctx(size_t sme)
{
    return S2CR_CBNDX(smmu.hw.glbl_rs0->S2CR[sme]);
}

inline streamid_t smmu_sme_get_id(size_t sme)
{
    return SMMU_SMR_ID(smmu.hw.glbl_rs0->SMR[sme]);
}

inline streamid_t smmu_sme_get_mask(size_t sme)
{
    return SMMU_SMR_MASK(smmu.hw.glbl_rs0->SMR[sme]);
}

static void smmu_check_features()
{
    unsigned version =
        bit32_extract(smmu.hw.glbl_rs0->IDR7, SMMUV2_IDR7_MAJOR_OFF, SMMUV2_IDR7_MAJOR_LEN);
    if (version != 2) {
        ERROR("smmu unsupported version: %d", version);
    }

    if (!(smmu.hw.glbl_rs0->IDR0 & SMMUV2_IDR0_S2TS_BIT)) {
        ERROR("smmuv2 does not support 2nd stage translation");
    }

    if (!(smmu.hw.glbl_rs0->IDR0 & SMMUV2_IDR0_SMS_BIT)) {
        ERROR("smmuv2 does not support stream match");
    }

    /**
     * TODO: the most common smmuv2 implementation (mmu-500) does not provide ptw coherency. So we
     * must add some mechanism software-managed coherency mechanism for the vms using the smmu
     * according to the result of this feature test.
     */
    if (!(smmu.hw.glbl_rs0->IDR0 & SMMUV2_IDR0_CTTW_BIT)) {
        WARNING("smmuv2 does not support coherent page table walks");
    }

    if (!(smmu.hw.glbl_rs0->IDR0 & SMMUV2_IDR0_BTM_BIT)) {
        ERROR("smmuv2 does not support tlb maintenance broadcast");
    }

    if (!(smmu.hw.glbl_rs0->IDR2 & SMMUV2_IDR2_PTFSv8_4kB_BIT)) {
        ERROR("smmuv2 does not support 4kb page granule");
    }

    size_t pasize = bit32_extract(smmu.hw.glbl_rs0->IDR2, SMMUV2_IDR2_OAS_OFF, SMMUV2_IDR2_OAS_LEN);
    size_t ipasize =
        bit32_extract(smmu.hw.glbl_rs0->IDR2, SMMUV2_IDR2_IAS_OFF, SMMUV2_IDR2_IAS_LEN);

    if (pasize < parange) {
        ERROR("smmuv2 does not support the full available pa range");
    } else if (ipasize < parange) {
        ERROR("smmuv2 does not support the full available ipa range");
    }
}

void smmu_init()
{
    /*
     * Alloc pages for global address space.
     *
     * Map the first 4k so we can read all the info we need to further allocate smmu registers.
     */
    vaddr_t smmu_glbl_rs0 = mem_alloc_map_dev(&cpu()->as, SEC_HYP_GLOBAL, INVALID_VA,
        platform.arch.smmu.base, NUM_PAGES(sizeof(struct smmu_glbl_rs0_hw)));

    smmu.hw.glbl_rs0 = (struct smmu_glbl_rs0_hw*)smmu_glbl_rs0;

    size_t pg_size = smmu.hw.glbl_rs0->IDR1 & SMMUV2_IDR1_PAGESIZE_BIT ? 0x10000 : 0x1000;
    size_t num_page = 1ULL << (bit32_extract(smmu.hw.glbl_rs0->IDR1, SMMUV2_IDR1_NUMPAGEDXB_OFF,
                                   SMMUV2_IDR1_NUMPAGEDXB_LEN) +
                          1);
    size_t ctx_bank_num =
        bit32_extract(smmu.hw.glbl_rs0->IDR1, SMMUV2_IDR1_NUMCB_OFF, SMMUV2_IDR1_NUMCB_LEN);

    vaddr_t smmu_glbl_rs1 = mem_alloc_map_dev(&cpu()->as, SEC_HYP_GLOBAL, INVALID_VA,
        platform.arch.smmu.base + pg_size, NUM_PAGES(sizeof(struct smmu_glbl_rs1_hw)));

    vaddr_t smmu_cntxt = mem_alloc_map_dev(&cpu()->as, SEC_HYP_GLOBAL, INVALID_VA,
        platform.arch.smmu.base + (num_page * pg_size), NUM_PAGES(pg_size * ctx_bank_num));


    vaddr_t smmu_pmu = mem_alloc_map_dev(&cpu()->as, SEC_HYP_GLOBAL, INVALID_VA,
        platform.arch.smmu.base + (3 * pg_size), NUM_PAGES(pg_size * ctx_bank_num));

    smmu.hw.glbl_rs1 = (struct smmu_glbl_rs1_hw*)smmu_glbl_rs1;
    smmu.hw.cntxt = (struct smmu_cntxt_hw*)smmu_cntxt;
    // TODO: check if it possible to map the PMU before mapping the 16 pages
    smmu.hw.pmu = (struct smmu_pmu_hw*)smmu_pmu;

    /* Everything is mapped. Initialize book-keeping data. */

    smmu_check_features();

    smmu.ctx_lock = SPINLOCK_INITVAL;
    smmu.ctx_num = ctx_bank_num;
    bitmap_clear_consecutive(smmu.ctxbank_bitmap, 0, smmu.ctx_num);

    smmu.sme_lock = SPINLOCK_INITVAL;
    smmu.sme_num = smmu.hw.glbl_rs0->IDR0 & SMMUV2_IDR0_MASK;
    bitmap_clear_consecutive(smmu.sme_bitmap, 0, smmu.sme_num);
    bitmap_clear_consecutive(smmu.grp_bitmap, 0, smmu.sme_num);

    /* Clear random reset state. */
    smmu.hw.glbl_rs0->GFSR = smmu.hw.glbl_rs0->GFSR;
    smmu.hw.glbl_rs0->NSGFSR = smmu.hw.glbl_rs0->NSGFSR;

    for (size_t i = 0; i < smmu.sme_num; i++) {
        smmu.hw.glbl_rs0->SMR[i] = 0;
    }

    for (size_t i = 0; i < smmu.ctx_num; i++) {
        smmu.hw.cntxt[i].SCTLR = 0;
        smmu.hw.cntxt[i].FSR = -1;
    }

    /* Enable IOMMU. */
    uint32_t cr0 = smmu.hw.glbl_rs0->CR0;
    cr0 = SMMUV2_CR0_CLEAR(cr0);
    cr0 |= SMMUV2_CR0_USFCFG | SMMUV2_CR0_SMCFCFG;
    cr0 &= ~SMMUV2_CR0_CLIENTPD;
    smmu.hw.glbl_rs0->CR0 = cr0;
}

ssize_t smmu_alloc_ctxbnk()
{
    spin_lock(&smmu.ctx_lock);
    /* Find a free context bank. */
    ssize_t nth = bitmap_find_nth(smmu.ctxbank_bitmap, smmu.ctx_num, 1, 0, false);
    if (nth >= 0) {
        bitmap_set(smmu.ctxbank_bitmap, nth);
    }
    spin_unlock(&smmu.ctx_lock);

    return nth;
}

static size_t smmu_cb_ttba_offset(size_t t0sz)
{
    size_t offset = 12;

    if (parange_table[parange] < 44) {
        /* SMMUV2_TCR_SL0_1 */
        if (t0sz >= 21 && t0sz <= 33) {
            offset = 37 - t0sz;
        }
    } else {
        /* SMMUV2_TCR_SL0_0 */
        if (t0sz >= 16 && t0sz <= 24) {
            offset = 28 - t0sz;
        }
    }

    return offset;
}

void smmu_write_ctxbnk(size_t ctx_id, paddr_t root_pt, asid_t vm_id)
{
    spin_lock(&smmu.ctx_lock);
    if (!bitmap_get(smmu.ctxbank_bitmap, ctx_id)) {
        ERROR("smmu ctx %d is already allocated", ctx_id);
    } else {
        /* Set type as stage 2 only. */
        smmu.hw.glbl_rs1->CBAR[ctx_id] = SMMUV2_CBAR_VMID(vm_id);
        smmu.hw.glbl_rs1->CBA2R[ctx_id] = SMMUV2_CBAR_VA64;

        /**
         * This should closely match to the VTCR configuration set up in vmm_arch_init as we're
         * sharing page table between the VM and its smmu context.
         */
        uint32_t tcr = ((parange << SMMUV2_TCR_PS_OFF) & SMMUV2_TCR_PS_MSK);
        size_t t0sz = 64 - parange_table[parange];
        tcr |= SMMUV2_TCR_TG0_4K;
        tcr |= SMMUV2_TCR_ORGN0_WB_RA_WA;
        tcr |= SMMUV2_TCR_IRGN0_WB_RA_WA;
        tcr |= SMMUV2_TCR_T0SZ(t0sz);
        tcr |= SMMUV2_TCR_SH0_IS;
        tcr |= ((parange_table[parange] < 44) ? SMMUV2_TCR_SL0_1 : SMMUV2_TCR_SL0_0);
        smmu.hw.cntxt[ctx_id].TCR = tcr;
        smmu.hw.cntxt[ctx_id].TTBR0 = root_pt & SMMUV2_CB_TTBA(smmu_cb_ttba_offset(t0sz));

        uint32_t sctlr = smmu.hw.cntxt[ctx_id].SCTLR;
        sctlr = SMMUV2_SCTLR_CLEAR(sctlr);
        sctlr |= SMMUV2_SCTLR_DEFAULT;
        smmu.hw.cntxt[ctx_id].SCTLR |= sctlr;
    }
    spin_unlock(&smmu.ctx_lock);
}

ssize_t smmu_alloc_sme()
{
    spin_lock(&smmu.sme_lock);
    /* Find a free sme. */
    ssize_t nth = bitmap_find_nth(smmu.sme_bitmap, smmu.sme_num, 1, 0, false);
    if (nth >= 0) {
        bitmap_set(smmu.sme_bitmap, nth);
    }
    spin_unlock(&smmu.sme_lock);

    return nth;
}

/*
 * When adding a new stream match entry there are a two different cases:
 *
 *      1. sme is a group;
 *      2. sme is a device.
 *
 * Groups can be merged together if one is found to be inclusive or equal of the other.
 *
 * Devices can be added (i.e. merged into) a group, but not together.
 *
 * This function searches for existing smes that are compatible for merging with the new sme,
 * raising an ERROR when conflicting attributes are found.
 */
bool smmu_compatible_sme_exists(streamid_t mask, streamid_t id, size_t ctx, bool group)
{
    bool included = false;
    size_t sme = 0;

    spin_lock(&smmu.sme_lock);
    smmu_for_each_sme(sme)
    {
        streamid_t sme_mask = smmu_sme_get_mask(sme);
        streamid_t mask_r = mask & sme_mask;
        streamid_t diff_id = (smmu_sme_get_id(sme) ^ id) & ~(mask | sme_mask);

        if (!diff_id) {
            /* Only group-to-group or device-to-group can be merged */
            if (((group || smmu_sme_is_group(sme)) && (mask_r == mask || mask_r == sme_mask)) &&
                ctx == smmu_sme_get_ctx(sme)) {
                /* Compatible entry found.
                 *
                 * If the new entry includes an existing one, there is the possibility that it will
                 * include other existing entries, it is therefore necessary to remove the existing
                 * entry and keep searching.
                 */
                if (mask > sme_mask) {
                    bitmap_clear(smmu.sme_bitmap, sme);
                } else {
                    included = true;
                    break;
                }

            } else {
                ERROR("SMMU sme conflict");
            }
        }
    }
    spin_unlock(&smmu.sme_lock);

    return included;
}

void smmu_write_sme(size_t sme, streamid_t mask, streamid_t id, bool group)
{
    spin_lock(&smmu.sme_lock);
    if (!bitmap_get(smmu.sme_bitmap, sme)) {
        ERROR("smmu: trying to write unallocated sme %d", sme);
    } else {
        smmu.hw.glbl_rs0->SMR[sme] = mask << SMMU_SMR_MASK_OFF;
        smmu.hw.glbl_rs0->SMR[sme] |= id & SMMU_ID_MSK;
        smmu.hw.glbl_rs0->SMR[sme] |= SMMUV2_SMR_VALID;

        if (group) {
            bitmap_set(smmu.grp_bitmap, sme);
        }
    }
    spin_unlock(&smmu.sme_lock);
}

void smmu_write_s2c(size_t sme, size_t ctx_id)
{
    spin_lock(&smmu.sme_lock);
    if (!bitmap_get(smmu.ctxbank_bitmap, ctx_id)) {
        ERROR("smmu: trying to write unallocated s2c %d", ctx_id);
    } else if (!bitmap_get(smmu.sme_bitmap, sme)) {
        ERROR("smmu: trying to bind unallocated sme %d", sme);
    } else {
        /* Initial contex is a translation context. */
        uint32_t s2cr = smmu.hw.glbl_rs0->S2CR[sme];

        s2cr = S2CR_CLEAR(s2cr);
        s2cr |= S2CR_DFLT;
        s2cr |= ctx_id & S2CR_CBNDX_MASK;

        smmu.hw.glbl_rs0->S2CR[sme] = s2cr;
    }
    spin_unlock(&smmu.sme_lock);
}

size_t smmu_implemented_event_cntrs() {
    size_t pmcfgr = smmu.hw.cntxt[0].PMCFGR;
    pmcfgr = pmcfgr & SMMU_PMCFGR_N_MASK;
    return pmcfgr;
}

size_t implemented_cntr_groups() {
    size_t pmcfgr = smmu.hw.cntxt[0].PMCFGR;
    pmcfgr = pmcfgr & SMMU_PMCFGR_NCG_MASK;
    return pmcfgr;
}

bool smmu_is_valid_event(uint32_t smmu_event) {
    for (enum smmuv2_pmu_events event = SMMU_PME_CYCLE_COUNT; event <= SMMU_PME_ACC_WRITE; event++) {
        if (smmu_event == event) {
            return true;
        }
    }

    return false;
}

// Dummy IRQ Handler
void smmu_irq_handler(){
    return;
}

bool smmu_setup_counter(size_t ctx_id, size_t counter_id, uint32_t smmu_event, bool en_irq) {
    counter_id = smmu_implemented_event_cntrs() + 1;
    if (counter_id >= SMMU_PMU_MAX_COUNTERS) {
        return false;
    }

    if (!smmu_is_valid_event(smmu_event)) {
        return false;
    }
    /*
    Procedure:
        1. Set evente type
        2. Set counter
        3. Clear counter overflows
        4. Set IRQ callback (optional)
        5. Enable the IRQ (only if 4. is performed)
    */

    smmu_set_event_type(counter_id, smmu_event);

    // smmu_set_event_cntr(counter_id, UINT32_MAX);
    smmu_set_event_cntr(counter_id, 10);
    
    smmu_event_ctr_ovf_clr(counter_id);

    if(en_irq) {
        smmu_pmu_interrupt_enable(counter_id, smmu_irq_handler);
    }

    return true;

}

void smmu_set_event_type(size_t counter, size_t event){

    // PMEVTYPERn, Performance Monitors Event Type Registers
    // uint32_t pmevtyper = smmu.hw.cntxt->PMEVTYPERm[counter];
    uint32_t pmevtyper = smmu.hw.cntxt[0].PMEVTYPERm[counter];
    // console_printk("pmevtyper: %d\n", pmevtyper);

    /*  P, bit[31] Privileged transactions filtering bit. Controls the counting
    *   of Secure privileged transactions.
    *       0 Count events relating to Secure privileged transactions.
    *       1 Do not count events relating to Secure privileged transactions.
    */
    pmevtyper = bit32_clear(pmevtyper, SMMU_PMEVTYPER_P_OFF);
    // console_printk("pmevtyper: %d\n", pmevtyper);

    /*  U, bit[30] Unprivileged transactions filtering bit. Controls the
    *    counting of Secure unprivileged transactions.
    *       0 Count events relating to Secure unprivileged transactions.
    *       1 Do not count events relating to Secure unprivileged transactions.
    */
    pmevtyper = bit32_clear(pmevtyper, SMMU_PMEVTYPER_U_OFF);
    // console_printk("pmevtyper: %d\n", pmevtyper);

    /*  NSP, bit[29] Non-secure Privileged transactions filtering bit. Controls
    *   the counting of Non-secure privileged transactions.
    *       P==NSP Count events relating to Non-secure privileged transactions.
    *       P!=NSP Do not count events relating to Non-secure privileged 
    *           transactions.
    */
    pmevtyper = bit32_set(pmevtyper, SMMU_PMEVTYPER_NSP_OFF);
    // console_printk("pmevtyper: %d\n", pmevtyper);

    /*   NSU, bit[28] Non-secure unprivileged transactions filtering bit.
    *    Controls counting of Non-secure unprivileged transactions.
    *        U==NSU Count events relating to Non-secure unprivileged 
    *           transactions.
    *        U!=NSU Do not count events relating to Non-secure unprivileged 
    *            transactions.Non-secure unprivileged transactions filtering
    *            bit. Controls counting of Non-secure unprivileged transactions.
    */
    pmevtyper = bit32_set(pmevtyper, SMMU_PMEVTYPER_NSU_OFF);
    console_printk("pmevtyper: %d\n", pmevtyper);

    pmevtyper = bit32_insert(pmevtyper, event, SMMU_PMEVTYPER_EVENT_OFF, SMMU_PMEVTYPER_EVENT_LEN);
    
    smmu.hw.cntxt[0].PMEVTYPERm[counter] = pmevtyper;
    console_printk("pmevtyper: %d\n", pmevtyper);
    console_printk("smmu.hw.cntxt[0].PMEVTYPERm[counter]: %d\n", smmu.hw.cntxt[0].PMEVTYPERm[counter]);

}

void smmu_set_event_cntr(size_t counter, size_t value){
    smmu.hw.cntxt[0].PMEVCNTRm[counter] = value;
}

void smmu_event_ctr_ovf_clr(size_t counter){
    uint32_t pmovsclr = smmu.hw.cntxt[0].PMOVSCLR;
    pmovsclr = bit32_clear(pmovsclr, counter);
    smmu.hw.cntxt[0].PMOVSCLR = pmovsclr;
}

void smmu_pmu_interrupt_enable(size_t counter, irq_handler_t handler){
    interrupts_reserve(SMMU_IRQ_ID, handler);
    interrupts_arch_enable(SMMU_IRQ_ID, true);

    volatile uint32_t pmintenset = smmu.hw.cntxt[0].PMINTENSET;
    pmintenset = bit_set(pmintenset, counter);
    smmu.hw.cntxt[0].PMINTENSET = pmintenset;
}

void smmu_events_init(size_t ctx_id) {

    // smmu_enable_pmc(ctx_id);

    size_t counter_id = 0;
    uint32_t event = SMMU_PME_TLB_ENTRY_WRITE;

    console_printk("vm ctx_id: %d\n", ctx_id);

    smmu_setup_counter(ctx_id, counter_id, event, true);
    
    uint32_t counter_val = 0;
    counter_val = smmu_read_counter(counter_id);
    counter_val += 1;
    console_printk("counter_val: %d\n", counter_val);
}

uint32_t smmu_read_counter(size_t counter){
    return smmu.hw.cntxt[0].PMEVCNTRm[counter];
}

// void smmu_enable_pmc(size_t ctx_id/*, size_t counter_id*/) {

//     size_t counter_id = 0;
//     spin_lock(&smmu.ctx_lock);

//     uint32_t pmcgcr = smmu.hw.pmu->PMCFGR[counter_id];
//     pmcgcr = bit32_set(pmcgcr, SMMU_PMCGCR_CBAEN_OFF);

//     pmcgcr = bit32_insert(pmcgcr, ctx_id, SMMU_PMCGCR_NDX_OFF, SMMU_PMCGCR_NDX_LEN);

//     uint32_t pmcr = smmu.hw.cntxt[ctx_id].PMCR;

//     // enable export of events
//     pmcr = bit32_set(pmcr, SMMU_PMCR_X_OFF);

//     // reset all event counters to zero
//     pmcr = bit32_set(pmcr, SMMU_PMCR_P_OFF);

//     // enable all counters
//     pmcr = bit32_set(pmcr, SMMU_PMCR_E_OFF);

//     smmu.hw.cntxt[ctx_id].PMCR = pmcr;
//     spin_unlock(&smmu.ctx_lock);

//     console_printk("pmcr: %d\n", pmcr);
//     console_printk("smmu.hw.cntxt[ctx_id].PMCR: %d\n", smmu.hw.cntxt[ctx_id].PMCR);
// }


// uint32_t smmu_alloc_counter(size_t ctx_id) {
//     spin_lock(&smmu.ctx_lock);

//     // Get number of implemented counters
//     size_t implement_cntrs = smmu_implemented_event_cntrs();

//     uint32_t pmcgcr = smmu.hw.pmu->PMCFGR[implement_cntrs];

//     spin_unlock(&smmu.ctx_lock);

// }

/******************************************************************************* Clean Code frome here*/
/******************************************************************************/
/******************************************************************************/

// #define SMMU_PMU_MAX_EVENT_CNTRS    (256)
// #define SMMU_PMU_MAX_CNTR_GROUPS    (256)
// #define SMMU_PMUX_MAX_CNTR_PER_GRP  (15)


void smmu_pmu_init();
void smmu_pmu_filtering(size_t filter_type, size_t cntr_group_id);
uint32_t smmu_pmu_alloc_cntr_grp();
uint32_t smmu_pmu_alloc_cntr();

void smmu_pmu_config_cntr_group(size_t counter_group_id);
void smmu_pmu_filtering(size_t filter_type, size_t cntr_group_id);
void smmu_en_pmc_event_export(size_t cntr_group_id);
void smmu_en_pmc(size_t cntr_group_id);
void smmu_en_ctxbnk_assignment(size_t cntr_group_id);

enum smmu_pmu_filter {
    smmu_pmu_filter_glb = 0b00,             // Count events on a global basis
    smmu_pmu_filter_rest_pmcgsmr = 0b01,    // Counter events are restricted to matches in the corresponding PMCGSMRn register
    smmu_pmu_filter_rest_trans_bnk = 0b10,  // Counter events are restricted to the translation context bank indicated by PMCGCRn.NDX
    smmu_pmu_filter_reserved = 0b11         // Reserved
};

void smmu_pmu_init() {
    size_t counter_group_id = smmu_pmu_alloc_cntr_grp();

    // Setup SMMU PMU filtering
    smmu_pmu_config_cntr_group(counter_group_id);
    size_t ctxt_bank_id = smmu_cntr_grp_ctx_bank(counter_group_id);
}

/*************************************************************************************************** [Begin] Configure Counter Group*/
/**************************************************************************************************/

#define SMMU_PMCGCR_CGNC_OFF            (24)
#define SMMU_PMCGCR_CGNC_LEN            (4)
#define SMMU_PMCGCR_CGNC_MASK           BIT32_MASK(SMMU_PMCGCR_CGNC_OFF, SMMU_PMCGCR_CGNC_LEN)

#define SMMU_PMCGCR_SIDG_OFF             (16)
#define SMMU_PMCGCR_SIDG_LEN             (7)
#define SMMU_PMCGCR_SIDG_MASK            BIT32_MASK(SMMU_PMCGCR_SIDG_OFF, SMMU_PMCGCR_SIDG_LEN)

#define SMMU_PMCGCR_X_OFF               (12)
#define SMMU_PMCGCR_X_LEN               (1)

#define SMMU_PMCGCR_E_OFF               (11)
#define SMMU_PMCGCR_E_LEN               (1)

#define SMMU_PMCGCR_CBAEN_OFF           (10)
#define SMMU_PMCGCR_CBAEN_LEN           (1)

#define SMMU_PMCGCR_TCEFCFG_OFF         (8)
#define SMMU_PMCGCR_TCEFCFG_LEN         (2)

#define SMMU_PMCGCR_NDX_OFF             (0)
#define SMMU_PMCGCR_NDX_LEN             (8)
#define SMMU_PMCGCR_NDX_MASK            BIT32_MASK(SMMU_PMCGCR_NDX_OFF, SMMU_PMCGCR_NDX_LEN)

void smmu_pmu_config_cntr_group(size_t counter_group_id) {
    smmu_pmu_filtering(smmu_pmu_filter_glb, counter_group_id);
    smmu_en_pmc_event_export(counter_group_id);
    smmu_en_pmc(counter_group_id);
    smmu_en_ctxbnk_assignment(counter_group_id);
}

void smmu_pmu_filtering(size_t filter_type, size_t cntr_group_id) {
    // Performance Monitors Counter Group Configuration Registers
    uint32_t pmcgcr = smmu.hw.pmu->PMCGCRn[cntr_group_id];
    pmcgcr = bit32_insert(pmcgcr, filter_type, SMMU_PMCGCR_TCEFCFG_OFF, SMMU_PMCGCR_TCEFCFG_LEN);
    smmu.hw.pmu->PMCGCRn[cntr_group_id] = pmcgcr;
}

void smmu_en_pmc_event_export(size_t cntr_group_id){
    uint32_t pmcgcr = smmu.hw.pmu->PMCGCRn[cntr_group_id];
    pmcgcr = bit32_set(pmcgcr, SMMU_PMCGCR_X_OFF);
    smmu.hw.pmu->PMCGCRn[cntr_group_id] = pmcgcr;
}

void smmu_en_pmc(size_t cntr_group_id) {
    uint32_t pmcgcr = smmu.hw.pmu->PMCGCRn[cntr_group_id];
    pmcgcr = bit32_set(pmcgcr, SMMU_PMCGCR_E_OFF);
    smmu.hw.pmu->PMCGCRn[cntr_group_id] = pmcgcr;
}

void smmu_en_ctxbnk_assignment(size_t cntr_group_id) {
    uint32_t pmcgcr = smmu.hw.pmu->PMCGCRn[cntr_group_id];
    pmcgcr = bit32_set(pmcgcr, SMMU_PMCGCR_CBAEN_OFF);
    smmu.hw.pmu->PMCGCRn[cntr_group_id] = pmcgcr;
}

size_t smmu_cntr_grp_ctx_bank(size_t cntr_group_id) {
    uint32_t pmcgcr = smmu.hw.pmu->PMCGCRn[cntr_group_id];
    size_t cntx_bank_id = pmcgcr & SMMU_PMCGCR_NDX_MASK;
    return cntx_bank_id;
}

size_t smmu_cntr_grp_implemented_cntrs(size_t cntr_group_id) {
    uint32_t pmcgcr = smmu.hw.pmu->PMCGCRn[cntr_group_id];
    size_t num_cntrs = pmcgcr & SMMU_PMCGCR_CGNC_MASK;
    return num_cntrs;
}

size_t smmu_cntr_grp_stream_id(size_t cntr_group_id) {
    uint32_t pmcgcr = smmu.hw.pmu->PMCGCRn[cntr_group_id];
    size_t stream_id = pmcgcr & SMMU_PMCGCR_CGNC_MASK;
    return stream_id;
}

/*************************************************************************************************** [End] Configure Counter Group*/
/**************************************************************************************************/


/*************************************************************************************************** [Begin] Configure Counter*/
/**************************************************************************************************/

#define SMMU_PMCR_E_OFF               (0)
#define SMMU_PMCR_E_LEN               (1)

void smmu_cb_pmc_enable(size_t ctxt_id) {
    uint32_t pmcr = smmu.hw.cntxt[ctxt_id].PMCR;
    pmcr = bit32_set(pmcr, SMMU_PMCR_E_OFF);
    smmu.hw.cntxt[ctxt_id].PMCR = pmcr;
}

#define SMMU_PMCR_P_OFF               (1)
#define SMMU_PMCR_P_LEN               (1)

void smmu_cb_pmc_reset(size_t ctxt_id) {
    uint32_t pmcr = smmu.hw.cntxt[ctxt_id].PMCR;
    pmcr = bit32_set(pmcr, SMMU_PMCR_P_OFF);
    smmu.hw.cntxt[ctxt_id].PMCR = pmcr;
}

#define SMMU_PMCR_X_OFF               (4)
#define SMMU_PMCR_X_LEN               (1)

void smmu_cb_pmc_enable_export(size_t ctxt_id) {
    uint32_t pmcr = smmu.hw.cntxt[ctxt_id].PMCR;
    pmcr = bit32_set(pmcr, SMMU_PMCR_X_OFF);
    smmu.hw.cntxt[ctxt_id].PMCR = pmcr;
}


/*************************************************************************************************** [End] Configure Counter*/
/**************************************************************************************************/

uint32_t smmu_pmu_alloc_cntr_grp() {
    // Performance Monitors Configuration Register
    uint32_t pmcfgr = smmu.hw.pmu->PMCFGR;
    uint32_t num_cntr_groups = pmcfgr & SMMU_PMCFGR_NCG_MASK;
    
    pmcfgr = bit32_insert(pmcfgr, num_cntr_groups+1, SMMU_PMCFGR_NCG_OFF, SMMU_PMCFGR_NCG_LEN);
    smmu.hw.pmu->PMCFGR = pmcfgr;

    return num_cntr_groups;
}

uint32_t smmu_pmu_alloc_cntr() {
    // Performance Monitors Configuration Register
    uint32_t pmcfgr = smmu.hw.pmu->PMCFGR;
    uint32_t num_cntr_groups = pmcfgr & SMMU_PMCFGR_N_MASK;
    
    pmcfgr = bit32_insert(pmcfgr, num_cntr_groups+1, SMMU_PMCFGR_N_OFF, SMMU_PMCFGR_N_LEN);
    smmu.hw.pmu->PMCFGR = pmcfgr;

    return num_cntr_groups;
}