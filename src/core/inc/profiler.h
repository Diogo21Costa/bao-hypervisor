/**
 * SPDX-License-Identifier: Apache-2.0
 * Copyright (c) Bao Project and Contributors. All rights reserved.
 */

#ifndef __PROFILER_H__
#define __PROFILER_H__

#include <bao.h>
#include <cpu.h>
#include <spinlock.h>

#define MAX_PROFILING_EVENTS 6
#define MAX_PROFILING_CPUS   4

typedef struct {
    size_t num_target_cpus;
    size_t num_profiling_events;
    size_t* list_events;
} profiler_t;

typedef struct {
    volatile uint32_t num_enabled_cpus;
    volatile uint32_t num_cpus_ready;
    volatile uint32_t num_events_per_cpu;
    volatile uint32_t list_events[MAX_PROFILING_EVENTS];
    volatile uint32_t apu_init_ready;
    volatile uint32_t profiler_ready;
} profiler_manager_t;

extern profiler_t profiler;
extern spinlock_t profiler_lock;

void profiler_manager_init(size_t num_target_cpus, size_t num_profiling_events, size_t* list_events,
    size_t cpu_id);
vaddr_t profiler_init(void);
void profiler_config_core(vaddr_t coresight_base_addr, size_t cpu_id);

#endif /* __PROFILER_H__ */
