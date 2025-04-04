/*
BSD 3-Clause License

Copyright (c) Alliance for Sustainable Energy, LLC. See also https://github.com/NREL/ssc/blob/develop/LICENSE
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/



#ifndef _FUELCELL_COMMON_DATA_H_
#define _FUELCELL_COMMON_DATA_H_

#include <stdio.h>
#include "code_generator_utilities.h"

namespace {
	char load_profile_path_fc[256];
	char ac_watt_output_path[256];
	char ac_watt_lifetime_output_path[256];
	int nfc1 = sprintf(ac_watt_output_path, "%s/test/input_cases/general_data/ac.csv",SSCDIR);
	int nfc2 = sprintf(load_profile_path_fc, "%s/test/input_cases/general_data/commercial_load.csv",SSCDIR);
	int nfc3 = sprintf(ac_watt_lifetime_output_path, "%s/test/input_cases/general_data/ac_lifetime.csv",SSCDIR);
}

/**
*  Default data for no-financial pvsamv1 run that can be further modified
*/
void fuelcell_nofinancial_default(ssc_data_t &data)
{
	set_array(data, "ac",ac_watt_output_path, 8760);
	ssc_data_set_number(data, "system_use_lifetime_output", 0);
	ssc_data_set_number(data, "analysis_period", 30);
	set_array(data, "load",load_profile_path_fc, 8760);
	ssc_data_set_number(data, "fuelcell_degradation", 0.0099999997764825821);
	ssc_data_set_number(data, "fuelcell_degradation_restart", 1);
	ssc_data_set_number(data, "fuelcell_dispatch_choice", 0);
	ssc_data_set_number(data, "fuelcell_fixed_pct", 50);
	ssc_data_set_number(data, "fuelcell_dynamic_response_up", 20);
	ssc_data_set_number(data, "fuelcell_dynamic_response_down", 20);
	ssc_data_set_number(data, "fuelcell_efficiency_choice", 0);
	ssc_number_t p_fuelcell_efficiency[33] = { 0, 0, 50, 16, 21, 50, 25, 25, 50, 34, 32, 50, 44, 37, 50, 53, 42, 50, 62, 47, 49, 72, 50, 48, 82, 52, 47, 90, 52, 46, 100, 51, 45 };
	ssc_data_set_matrix(data, "fuelcell_efficiency", p_fuelcell_efficiency, 11, 3);
	ssc_number_t p_fuelcell_shutdown[2] = { 0, 0};
	ssc_data_set_matrix(data, "fuelcell_availability_schedule", p_fuelcell_shutdown, 1, 2);
	ssc_data_set_number(data, "fuelcell_fuel_available", 28316846);
	ssc_data_set_number(data, "fuelcell_fuel_price", 7.5);
	ssc_data_set_number(data, "fuelcell_fuel_type", 0);
	ssc_data_set_number(data, "fuelcell_lhv", 8917.48046875);
	ssc_data_set_number(data, "fuelcell_number_of_units", 1);
	ssc_data_set_number(data, "fuelcell_operation_options", 1);
	ssc_data_set_number(data, "fuelcell_replacement_option", 0);
	ssc_data_set_number(data, "fuelcell_replacement_percent", 50);
	ssc_number_t p_fuelcell_replacement_schedule[1] = { 0 };
	ssc_data_set_array(data, "fuelcell_replacement_schedule", p_fuelcell_replacement_schedule, 1);
	ssc_data_set_number(data, "fuelcell_startup_time", 24);
	ssc_data_set_number(data, "fuelcell_is_started", 0);
	ssc_data_set_number(data, "fuelcell_shutdown_time", 24);
	ssc_data_set_number(data, "fuelcell_type", 0);
	ssc_data_set_number(data, "fuelcell_unit_max_power", 100);
	ssc_data_set_number(data, "fuelcell_unit_min_power", 20);
	ssc_number_t p_dispatch_manual_fuelcellcharge[6] = { 0, 0, 0, 0, 0, 0 };
	ssc_data_set_array(data, "dispatch_manual_fuelcellcharge", p_dispatch_manual_fuelcellcharge, 6);
	ssc_number_t p_dispatch_manual_fuelcelldischarge[6] = { 0, 0, 0, 0, 0, 0 };
	ssc_data_set_array(data, "dispatch_manual_fuelcelldischarge", p_dispatch_manual_fuelcelldischarge, 6);
	ssc_number_t p_dispatch_manual_percent_fc_discharge[1] = { 0 };
	ssc_data_set_array(data, "dispatch_manual_percent_fc_discharge", p_dispatch_manual_percent_fc_discharge, 1);
	ssc_number_t p_dispatch_manual_units_fc_discharge[1] = { 1 };
	ssc_data_set_array(data, "dispatch_manual_units_fc_discharge", p_dispatch_manual_units_fc_discharge, 1);
	ssc_number_t p_dispatch_manual_fuelcell_sched[288] = { 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 1, 1, 1, 1, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 3, 3, 1, 1, 1, 1 };
	ssc_data_set_matrix(data, "dispatch_manual_fuelcell_sched", p_dispatch_manual_fuelcell_sched, 12, 24);
	ssc_number_t p_dispatch_manual_fuelcell_sched_weekend[288] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
	ssc_data_set_matrix(data, "dispatch_manual_fuelcell_sched_weekend", p_dispatch_manual_fuelcell_sched_weekend, 12, 24);
	ssc_number_t p_dispatch[1] = { 0 };
	ssc_data_set_array(data, "fuelcell_dispatch", p_dispatch, 1);
}

#endif
