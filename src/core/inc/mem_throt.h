/**
 * SPDX-License-Identifier: Apache-2.0 
 * Copyright (c) Bao Project and Contributors. All rights reserved
 */

#ifndef __MEM_THROT_H__
#define __MEM_THROT_H__

#include <timer.h>
#include <events.h>

typedef struct mem_throt_info {
    unsigned long budget;
    bool throttled;
    bool initialized;
    size_t counter_id;
    uint64_t period_us;
    uint64_t period_counts;
} mem_throt_t;


void mem_throt_init(unsigned long budget, uint64_t period_us);
void mem_throt_period_timer_callback(irqid_t int_id);

void mem_throt_event_overflow_callback(irqid_t int_id);
void mem_throt_resume(void);

void mem_throt_timer_init(irq_handler_t handler);
void mem_throt_events_init(events_enum event, unsigned long budget, irq_handler_t handler);

#endif /* __MEM_THROT_H__ */
