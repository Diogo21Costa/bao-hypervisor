/**
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) Bao Project and Contributors. All rights reserved.
 */

#include <hypercall.h>
#include <cpu.h>
#include <vm.h>
#include <ipc.h>
#include <arch/smmuv2.h>

long int hypercall(unsigned long id)
{
    long int ret = -HC_E_INVAL_ID;

    unsigned long ipc_id = vcpu_readreg(cpu()->vcpu, HYPCALL_ARG_REG(0));
    unsigned long arg1 = vcpu_readreg(cpu()->vcpu, HYPCALL_ARG_REG(1));
    unsigned long arg2 = vcpu_readreg(cpu()->vcpu, HYPCALL_ARG_REG(2));

    unsigned long event, counter;
    ssize_t ctx_id;

    switch (id) {
        case HC_IPC:
            ret = ipc_hypercall(ipc_id, arg1, arg2);
            break;
        case HC_SMMU_PMU_SETUP_CNTR:
            event = vcpu_readreg(cpu()->vcpu, HYPCALL_ARG_REG(0));
            counter = vcpu_readreg(cpu()->vcpu, HYPCALL_ARG_REG(1));
            ctx_id = cpu()->vcpu->vm->io.prot.mmu.ctx_id;
            smmu_cb_setup_counter(ctx_id, event, counter);
            break;
        case HC_SMMU_PMU_EN_CNTRS:
            ctx_id = cpu()->vcpu->vm->io.prot.mmu.ctx_id;
            smmu_cb_pmc_enable(ctx_id);
            break;
        case HC_SMMU_PMU_RST_CNTRS:
            ctx_id = cpu()->vcpu->vm->io.prot.mmu.ctx_id;
            smmu_cb_pmc_reset(ctx_id);
            break;
        default:
            WARNING("Unknown hypercall id %d", id);
    }

    return ret;
}
