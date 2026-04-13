/**
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) Bao Project and Contributors. All rights reserved.
 */

#ifndef __PROFILER_CONFIG_H__
#define __PROFILER_CONFIG_H__

#include <bao.h>

struct out_of_core_profiler_config {
    size_t num_target_cpus;
    size_t num_profile_events;
    size_t* list_events;
};

#endif /* __PROFILER_CONFIG_H__ */
