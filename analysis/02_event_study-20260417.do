*==============================================================================* 
* Project: DrunkDriving-TrafficDecline-RiskAnalysis         				   *
* Version:    									    				           * 
* Date:  				                                                       *
* Author: Hui Liu  					                                           * 
*==============================================================================*


* set directories 
clear all
set maxvar  30000
set more off


********************************************************************************
* 
* Analysis: Event Study Analysis for Parallel Trends Testing
* 
********************************************************************************

* Data versioning
global mydate    "20260415"
global myversion "v1"


* ==================== Data preparation & file setup ====================

global base "I:\Academic\demo\TrafficDecline-RiskAnalysis\github"

global result_base    "$base\\data\\reg_results"
global level1         "$result_base\\${mydate}-${myversion}"
global result         "$level1\02_analysis_event_study"

display "Result path: $result"

capture mkdir "$level1"
capture mkdir "$result"

capture cd "$result"
if _rc == 0 {
    display as result "Success: $result"
    cd "$result"
}
else {
    display as error "Failed to create: $result"
}

global datapath "$base\\data\\working_data"


* ==================== Function definition ====================

* Function 1: Prepare data and generate event study variables
capture program drop prepare_event_study
program define prepare_event_study
    xtset id year
    
    gen treated = x
    gen post_ = month_01
    gen acc = acc_num
    gen drink = drink_num
    gen acc_ratio__ = acc / drink
    drop if drink_num == 0
    gen did = post_ * treated
    
    capture drop city_id
    encode city, g(city_id)
    
    * Year fixed effects
    gen nyear = 1 if year_whole >= 2001 & year_whole <= 2012
    replace nyear = 2 if year_whole >= 2013 & year_whole <= 2024
    replace nyear = 3 if year_whole >= 2025 & year_whole <= 2036
    replace nyear = 4 if year_whole >= 2037 & year_whole <= 2048
    
    * Generate relative period variable
    gen period = year - 2004
    
    * Generate event study dummy variables
    forvalues i = 3(-1)1 {
        gen pre_`i' = (period == -`i' & treated == 1)
    }
    gen current = (period == 0 & treated == 1)
    forvalues j = 1(1)8 {
        gen time_`j' = (period == `j' & treated == 1)
    }
    
    * Exclude baseline period (t-1)
    drop pre_1
end

* Function 2: Run event study regression
capture program drop run_event_study
program define run_event_study
    args output_name
    
    * ── Descriptive statistics for the baseline period (post_ == 0) ───────────
    qui sum acc       if treated == 1 & post_ == 0
    local t_acc   = r(mean)
    qui sum acc       if treated == 0 & post_ == 0
    local c_acc   = r(mean)

    qui sum drink     if treated == 1 & post_ == 0
    local t_drink = r(mean)
    qui sum drink     if treated == 0 & post_ == 0
    local c_drink = r(mean)

    qui sum acc_ratio__ if treated == 1 & post_ == 0
    local t_ratio = r(mean)
    qui sum acc_ratio__ if treated == 0 & post_ == 0
    local c_ratio = r(mean)

    qui tab city_id if treated == 1
    local n_treat = r(r)
    qui tab city_id if treated == 0
    local n_ctrl  = r(r)

	
	* ── Model 1: drunk driving cases(PPML) ────────────────────────────────────
	ppmlhdfe drink pre_* current time_*, ///
        absorb(month#city_id nyear#city_id) cluster(city_id) irr nolog
    local eff_drink = (exp(_b[current]) - 1) * 100
	outreg2 using "$result/`output_name'.xls", replace se nocons lab dec(3) ///
        keep(pre_* current time_*) eform ci level(95) noaster ///
        addtext(Month×City FE, YES, Year×City FE, YES, SE Clustered, City Level) ///
        addstat("Observations",               e(N),          ///
                "Treated Cities",             `n_treat',     ///
                "Control Cities",             `n_ctrl',      ///
                "Baseline Mean (Treated)",    `t_drink',     ///
                "Baseline Mean (Control)",    `c_drink',     ///
                "Policy Effect (%) at t=0",   `eff_drink',   ///
                "Pseudo R-squared",           e(r2_p))

	* ── Model 2: drunk-driving-related crashes ────────────────────────────────
    ppmlhdfe acc pre_* current time_*, ///
        absorb(month#city_id nyear#city_id) cluster(city_id) irr nolog
    local eff_acc = (exp(_b[current]) - 1) * 100
    outreg2 using "$result/`output_name'.xls", append se nocons lab dec(3) ///
        keep(pre_* current time_*) eform ci level(95) noaster ///
        addtext(Month×City FE, YES, Year×City FE, YES, SE Clustered, City Level) ///
        addstat("Observations",               e(N),          ///
                "Treated Cities",             `n_treat',     ///
                "Control Cities",             `n_ctrl',      ///
                "Baseline Mean (Treated)",    `t_acc',       ///
                "Baseline Mean (Control)",    `c_acc',       ///
                "Policy Effect (%) at t=0",   `eff_acc',     ///
                "Pseudo R-squared",           e(r2_p))

    * ── Model 3: crash incidence ──────────────────────────────────────────────
    reghdfe acc_ratio__ pre_* current time_*, ///
        absorb(month#city_id nyear#city_id) cluster(city_id)
    local eff_ratio = (_b[current] / `t_ratio') * 100
    outreg2 using "$result/`output_name'.xls", append se nocons lab dec(3) ///
        keep(pre_* current time_*) ci level(95) noaster ///
        addtext(Month×City FE, YES, Year×City FE, YES, SE Clustered, City Level) ///
        addstat("Observations",               e(N),          ///
                "Treated Cities",             `n_treat',     ///
                "Control Cities",             `n_ctrl',      ///
                "Baseline Mean (Treated)",    `t_ratio',     ///
                "Baseline Mean (Control)",    `c_ratio',     ///
                "Policy Effect (%) at t=0",   `eff_ratio',   ///
                "R-squared",                  e(r2),         ///
				"R-squared-Within",           e(r2_within))
end

* Function 3: Execute single analysis
capture program drop do_analysis
program define do_analysis
    args output_name data_name
    
    di as result _n ">>> Processing: `output_name'"
    
    use "$datapath/`data_name'.dta", clear
    prepare_event_study
    run_event_study "`output_name'"
end


* ==================== Execute analysis ====================

* Group 1: Baseline
do_analysis "pop_All Pool Samples" "all_vs_saps_pop_pop"

* Group 2: Urban-Rural Classification
do_analysis "crl_Urbanization areas" "all_vs_saps_crl_city"
do_analysis "crl_Suburban areas" "all_vs_saps_crl_city_rural"
do_analysis "crl_Rural areas" "all_vs_saps_crl_rural"

* Group 3: Region
do_analysis "reg_East" "all_vs_saps_reg_East"
do_analysis "reg_Middle" "all_vs_saps_reg_Middle"
do_analysis "reg_West" "all_vs_saps_reg_West"

* Group 4: Road Type
do_analysis "rak_motorway" "all_vs_saps_rak_motorway"
do_analysis "rak_service" "all_vs_saps_rak_service"
do_analysis "rak_residential" "all_vs_saps_rak_residential"
do_analysis "rak_mainRoads" "all_vs_saps_rak_mainRoads"

* Group 5: Vehicle Type
do_analysis "car_C" "all_vs_saps_car_C"
do_analysis "car_D" "all_vs_saps_car_D"

* Group 6: Time Period
do_analysis "tim_afternoon" "all_vs_saps_tim_afternoon"
do_analysis "tim_midnight" "all_vs_saps_tim_midnight"
do_analysis "tim_morning" "all_vs_saps_tim_morning"
do_analysis "tim_night" "all_vs_saps_tim_night"

* Group 7: Prior DUI Crash Incidence
do_analysis "rat_lowerRate" "all_vs_saps_rat_lowerRate"
do_analysis "rat_higherRate" "all_vs_saps_rat_higherRate"


* ==================== Done ====================
display _n "{text:-}" _continu
di as result "All event study analyses done."
display as result "{hline 50}"
