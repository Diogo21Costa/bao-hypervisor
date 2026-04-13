/**
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) Bao Project and Contributors. All rights reserved.
 */

#include <hypercall.h>

#ifdef CONFIG_PROFILER
#include <profiler_hypercall.h>
#endif

long int hypercall(unsigned long id)
{
    long int ret = -HC_E_INVAL_ID;

    switch (id) {
        case HC_IPC:
            ret = ipc_hypercall();
            break;
        case HC_REMIO:
            ret = remio_hypercall();
            break;

        default:
            #ifdef CONFIG_PROFILER
            if(id >= PROFILER_HYPERCALL_BASE && id < PROFILER_HYPERCALL_MAX) {
                ret = profiler_hypercall(id);
                break;
            }
            #endif
            WARNING("Unknown hypercall id %d\n", id);
    }

    return ret;
}
