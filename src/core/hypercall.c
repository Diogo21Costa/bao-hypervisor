/**
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) Bao Project and Contributors. All rights reserved.
 */

#include <hypercall.h>
#include <cpu.h>
#include <vm.h>
#include <ipc.h>

long int hypercall(unsigned long id)
{
    long int ret = -HC_E_INVAL_ID;

    unsigned long ipc_id = vcpu_readreg(cpu()->vcpu, HYPCALL_ARG_REG(0));
    unsigned long arg1 = vcpu_readreg(cpu()->vcpu, HYPCALL_ARG_REG(1));
    unsigned long arg2 = vcpu_readreg(cpu()->vcpu, HYPCALL_ARG_REG(2));

    unsigned long event, count;

    switch (id) {
        case HC_IPC:
            ret = ipc_hypercall(ipc_id, arg1, arg2);
            break;
        case HC_SMMU_PMU:
            event = vcpu_readreg(cpu()->vcpu, HYPCALL_ARG_REG(0));
            count = vcpu_readreg(cpu()->vcpu, HYPCALL_ARG_REG(1));
            WARNING("SMMU PMU to be configured here with the following conditions:");
            WARNING("event: %d", event);
            WARNING("counter reset to: %d", count);
            break;
        default:
            WARNING("Unknown hypercall id %d", id);
    }

    return ret;
}
