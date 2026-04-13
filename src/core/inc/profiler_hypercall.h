/**
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) Bao Project and Contributors. All rights reserved.
 */

#ifndef __PROFILER_HYPERCALL_H__
#define __PROFILER_HYPERCALL_H__

#include <profiler.h>

#define HC_PROFILER_OS_READY  8
#define HC_PROFILE_CPU_READY  9
#define HC_PROFILE_INVALID_ID 15

int profiler_hypercall(unsigned long id);

#endif /* __PROFILER_HYPERCALL_H__ */
