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
* Analysis: Placebo Test
* 
********************************************************************************

* Data versioning
global mydate    "20260415"
global myversion "v1"


* ==================== Data preparation & file setup ====================

global base "I:\Academic\demo\TrafficDecline-RiskAnalysis\github"

global result_base    "$base\\data\\reg_results"
global level1         "$result_base\\${mydate}-${myversion}"
global result         "$level1\03_placebo_test"

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


/*==============================================================================
Step 1: Prepare Data and generate DID variables
==============================================================================*/

use "$datapath/all_vs_saps_pop_pop.dta", clear

xtset id year

gen acc         = acc_num
gen drink       = drink_num
gen acc_ratio__ = acc / drink
drop if drink_num == 0

capture drop city_id
encode city, g(city_id)

gen nyear = 1 if year_whole >= 2001 & year_whole <= 2012
replace nyear = 2 if year_whole >= 2013 & year_whole <= 2024
replace nyear = 3 if year_whole >= 2025 & year_whole <= 2036
replace nyear = 4 if year_whole >= 2037 & year_whole <= 2048


* Scenario II
gen saps_short = 1 if year_whole >= 2040 & year_whole <= 2042
replace saps_short = 0 if saps_short == .
* Scenario III
gen saps_medium = 1 if year_whole >= 2043 & year_whole <= 2048
replace saps_medium = 0 if saps_medium == .


gen ym = ym(year_whole, month)
format ym %tm
xtset city_id ym

bysort city_id: gen obs_n    = _n
bysort city_id: egen max_obs = max(obs_n)


/*==============================================================================
Step 2: Record the actual coefficients (baseline)
==============================================================================*/

quietly ppmlhdfe drink saps_short saps_medium, ///
    absorb(month#city_id nyear#city_id) cluster(city_id) irr nolog
scalar real_s_drink = _b[saps_short]
scalar real_m_drink = _b[saps_medium]

quietly ppmlhdfe acc saps_short saps_medium, ///
    absorb(month#city_id nyear#city_id) cluster(city_id) irr nolog
scalar real_s_acc = _b[saps_short]
scalar real_m_acc = _b[saps_medium]

quietly reghdfe acc_ratio__ saps_short saps_medium, ///
    absorb(month#city_id nyear#city_id) cluster(city_id)
scalar real_s_ratio = _b[saps_short]
scalar real_m_ratio = _b[saps_medium]

di as text "  drink short=" %8.5f real_s_drink "  medium=" %8.5f real_m_drink
di as text "  acc   short=" %8.5f real_s_acc   "  medium=" %8.5f real_m_acc
di as text "  ratio short=" %8.5f real_s_ratio "  medium=" %8.5f real_m_ratio _n


/*==============================================================================
Step 3: Trim the sample to Scenario I
==============================================================================*/

keep if year_whole <= 2039


/*==============================================================================
Step 4: Run 1,000 iterations, simultaneously performing DID + event study
==============================================================================*/

* ---- Initialize DID result storage ----
tempfile did_results
postfile did_handle                         ///
    iteration                               ///
    p_s_drink p_m_drink                     ///
    p_s_acc   p_m_acc                       ///
    p_s_ratio p_m_ratio                     ///
    using `did_results', replace

* ---- Initialize event study result storage ----
* Each iteration × 12 time points (-3 -2 -1 0 1 2 3 4 5 6 7 8)
tempfile es_results
postfile es_handle                          ///
    iteration time_point                    ///
    es_acc es_drink es_ratio                ///
    using `es_results', replace

* Event study time point labels
local tp_list "-3 -2 -1 0 1 2 3 4 5 6 7 8"

forvalues i = 1/1000 {

    if mod(`i', 50) == 0 {
        di as text "  Processing: `i'/1000 (" %3.0f `i'/10 "%)"
    }

    quietly {

        cap restore, not
        preserve

        * ----------------------------------------------------------
        * Generate pseudo-treatment: independent random fake_start for each city
        * ----------------------------------------------------------
        bysort city_id: gen temp_fs = floor(runiform() * 28) + 2004 if _n == 1
        bysort city_id: egen fake_start = max(temp_fs)
        drop temp_fs

        * DID variables
        gen fake_short  = (year_whole >= fake_start     & year_whole <= fake_start + 2)
        gen fake_medium = (year_whole >= fake_start + 3 & year_whole <= fake_start + 8)

        * Event study variables (baseline period = fake_start - 1, omitted)
        gen fake_pre3 = (year_whole == fake_start - 3)
        gen fake_pre2 = (year_whole == fake_start - 2)
        gen fake_d0   = (year_whole == fake_start)
        gen fake_d1   = (year_whole == fake_start + 1)
        gen fake_d2   = (year_whole == fake_start + 2)
        gen fake_d3   = (year_whole == fake_start + 3)
        gen fake_d4   = (year_whole == fake_start + 4)
        gen fake_d5   = (year_whole == fake_start + 5)
        gen fake_d6   = (year_whole == fake_start + 6)
        gen fake_d7   = (year_whole == fake_start + 7)
        gen fake_d8   = (year_whole == fake_start + 8)

        local evvars "fake_pre3 fake_pre2 fake_d0 fake_d1 fake_d2 fake_d3 fake_d4 fake_d5 fake_d6 fake_d7 fake_d8"

        * ----------------------------------------------------------
        * DID regression
        * ----------------------------------------------------------

        * drink
        cap ppmlhdfe drink fake_short fake_medium, ///
            absorb(month#city_id nyear#city_id) cluster(city_id) irr nolog
        local cs_drink = cond(_rc==0, _b[fake_short],  .)
        local cm_drink = cond(_rc==0, _b[fake_medium], .)

        * acc
        cap ppmlhdfe acc fake_short fake_medium, ///
            absorb(month#city_id nyear#city_id) cluster(city_id) irr nolog
        local cs_acc = cond(_rc==0, _b[fake_short],  .)
        local cm_acc = cond(_rc==0, _b[fake_medium], .)

        * ratio
        cap reghdfe acc_ratio__ fake_short fake_medium, ///
            absorb(month#city_id nyear#city_id) cluster(city_id)
        local cs_ratio = cond(_rc==0, _b[fake_short],  .)
        local cm_ratio = cond(_rc==0, _b[fake_medium], .)

        post did_handle                                         ///
            (`i')                                               ///
            (`cs_drink') (`cm_drink')                          ///
            (`cs_acc')   (`cm_acc')                            ///
            (`cs_ratio') (`cm_ratio')

        * ----------------------------------------------------------
        * Event-study approach
        * ----------------------------------------------------------
        
		* drink
        cap ppmlhdfe drink `evvars', ///
            absorb(month#city_id nyear#city_id) cluster(city_id) irr nolog
        if _rc == 0 {
            foreach v of local evvars {
                local ce_drink_`v' = _b[`v']
            }
        }
        else {
            foreach v of local evvars { local ce_drink_`v' = . }
        }
		
        * acc
        cap ppmlhdfe acc `evvars', ///
            absorb(month#city_id nyear#city_id) cluster(city_id) irr nolog
        if _rc == 0 {
            foreach v of local evvars {
                local ce_acc_`v' = _b[`v']
            }
        }
        else {
            foreach v of local evvars { local ce_acc_`v' = . }
        }


        * ratio
        cap reghdfe acc_ratio__ `evvars', ///
            absorb(month#city_id nyear#city_id) cluster(city_id)
        if _rc == 0 {
            foreach v of local evvars {
                local ce_ratio_`v' = _b[`v']
            }
        }
        else {
            foreach v of local evvars { local ce_ratio_`v' = . }
        }

        * Write event study results (by time point)
        * Mapping between variable names and time points
        local vl_list "fake_pre3 fake_pre2 ref fake_d0 fake_d1 fake_d2 fake_d3 fake_d4 fake_d5 fake_d6 fake_d7 fake_d8"

        forvalues t = 1/12 {
            local tp : word `t' of `tp_list'
            local vl : word `t' of `vl_list'

            if "`vl'" == "ref" {
                * Baseline period: log coefficient = 0 (IRR = 1, OLS coef = 0)
                post es_handle (`i') (`tp') (0) (0) (0)
            }
            else {
                post es_handle (`i') (`tp') ///
                    (`ce_acc_`vl'') (`ce_drink_`vl'') (`ce_ratio_`vl'')
            }
        }

        restore
    }
}

postclose did_handle
postclose es_handle

di as text _n "✓ 1000 iterations completed" _n


/*==============================================================================
Step 5: Summarize DID results
==============================================================================*/

use `did_results', clear
save "$result/placebo_estimates_short_medium.dta", replace

* Calculate pseudo p-values
foreach var in drink acc ratio {
    foreach term in s m {
        quietly summarize p_`term'_`var', detail
        if real_`term'_`var' > 0 {
            quietly count if p_`term'_`var' >= real_`term'_`var' & !missing(p_`term'_`var')
        }
        else {
            quietly count if p_`term'_`var' <= real_`term'_`var' & !missing(p_`term'_`var')
        }
        local extreme = r(N)
        quietly count if !missing(p_`term'_`var')
        local pval_`term'_`var' = `extreme' / r(N)
    }
}

di as text "Pseudo p-value summary: "
di as text "  drink short=" %5.3f `pval_s_drink' "  medium=" %5.3f `pval_m_drink'
di as text "  acc   short=" %5.3f `pval_s_acc'   "  medium=" %5.3f `pval_m_acc'
di as text "  ratio short=" %5.3f `pval_s_ratio' "  medium=" %5.3f `pval_m_ratio'

/*------------------------------------------------------------------------------
  Output SI Table 1: Placebo DID summary table
  Content: Actual coefficients (IRR / OLS coef) + pseudo p-values, 3 outcomes × 2 scenarios
------------------------------------------------------------------------------*/

* Convert actual coefficients to IRR (PPML stores log coefficients)
local real_s_acc_irr   = exp(real_s_acc)
local real_m_acc_irr   = exp(real_m_acc)
local real_s_drink_irr = exp(real_s_drink)
local real_m_drink_irr = exp(real_m_drink)
* Ratio is OLS, no conversion needed
local real_s_ratio_v   = real_s_ratio
local real_m_ratio_v   = real_m_ratio

* Calculate mean and SD of the placebo distribution (as reference for distribution characteristics)
foreach var in acc drink ratio {
    foreach term in s m {
        if "`var'" != "ratio" {
            quietly gen irr_p_`term'_`var' = exp(p_`term'_`var')
            quietly sum irr_p_`term'_`var', detail
        }
        else {
            quietly sum p_`term'_`var', detail
        }
        local mean_`term'_`var' = r(mean)
        local sd_`term'_`var'   = r(sd)
        local p5_`term'_`var'   = r(p5)
        local p95_`term'_`var'  = r(p95)
    }
}

* Write txt table
local outf "$result/placebo_SI_table1_DID.txt"
cap file close fh
file open fh using "`outf'", write replace text

file write fh "Supplementary Table: Placebo DID Estimates" _n
file write fh "Sample: October 2016 to December 2019 (pre-treatment period); 1,000 randomizations" _n
file write fh _n

* Table header
file write fh _tab "Drunk Driving Cases (IRR)" _tab _tab ///
                   "Traffic Crashes (IRR)"     _tab _tab ///
                   "Crash Incidence (OLS)"     _n
file write fh _tab "Scenario II" _tab "Scenario III" _tab ///
                   "Scenario II" _tab "Scenario III" _tab ///
                   "Scenario II" _tab "Scenario III" _n
file write fh "----------------------------------------------------------------------" _n

* Actual coefficient row
file write fh "Real estimate" _tab 
file write fh %6.3f (`real_s_acc_irr') _tab
file write fh %6.3f (`real_m_acc_irr') _tab
file write fh %6.3f (`real_s_drink_irr') _tab
file write fh %6.3f (`real_m_drink_irr') _tab
file write fh %6.3f (`real_s_ratio_v') _tab
file write fh %6.3f (`real_m_ratio_v') _n

* Pseudo p-value row
file write fh "Placebo p-value" _tab
file write fh %5.3f (`pval_s_acc') _tab
file write fh %5.3f (`pval_m_acc') _tab
file write fh %5.3f (`pval_s_drink') _tab
file write fh %5.3f (`pval_m_drink') _tab
file write fh %5.3f (`pval_s_ratio') _tab
file write fh %5.3f (`pval_m_ratio') _n

* Placebo distribution mean row
file write fh "Placebo mean" _tab
file write fh %6.3f (`mean_s_acc') _tab
file write fh %6.3f (`mean_m_acc') _tab
file write fh %6.3f (`mean_s_drink') _tab
file write fh %6.3f (`mean_m_drink') _tab
file write fh %6.3f (`mean_s_ratio') _tab
file write fh %6.3f (`mean_m_ratio') _n

* Placebo distribution SD row
file write fh "Placebo SD" _tab
file write fh %6.3f (`sd_s_acc') _tab
file write fh %6.3f (`sd_m_acc') _tab
file write fh %6.3f (`sd_s_drink') _tab
file write fh %6.3f (`sd_m_drink') _tab
file write fh %6.3f (`sd_s_ratio') _tab
file write fh %6.3f (`sd_m_ratio') _n

* Placebo distribution [p5, p95] row
file write fh "Placebo [p5, p95]" _tab
file write fh "[" 
file write fh %5.3f (`p5_s_acc') 
file write fh ", " 
file write fh %5.3f (`p95_s_acc') 
file write fh "]" _tab
file write fh "[" 
file write fh %5.3f (`p5_m_acc') 
file write fh ", " 
file write fh %5.3f (`p95_m_acc') 
file write fh "]" _tab

file write fh "----------------------------------------------------------------------" _n
file write fh "Note: IRR reported for PPML models; OLS coefficient for crash incidence." _n
file write fh "Placebo p-value = proportion of 1,000 randomizations with estimate" _n
file write fh "at least as extreme as the real estimate (one-sided)." _n

file close fh

* Export to Excel
putexcel set "$result/placebo_SI_table1_DID.xlsx", replace
putexcel A1 = "Supplementary Table: Placebo DID Estimates (1,000 randomizations)"
putexcel A2 = ""
putexcel B3 = "Drunk Driving Cases (IRR)" C3 = "" ///
         D3 = "Traffic Crashes (IRR)"     E3 = "" ///
         F3 = "Crash Incidence (OLS)"     G3 = ""
putexcel B4 = "Scenario II" C4 = "Scenario III" ///
         D4 = "Scenario II" E4 = "Scenario III" ///
         F4 = "Scenario II" G4 = "Scenario III"
putexcel A5  = "Real estimate" ///
         B5  = `real_s_acc_irr'   C5  = `real_m_acc_irr'   ///
         D5  = `real_s_drink_irr' E5  = `real_m_drink_irr' ///
         F5  = `real_s_ratio_v'   G5  = `real_m_ratio_v'
putexcel A6  = "Placebo p-value" ///
         B6  = `pval_s_acc'   C6  = `pval_m_acc'   ///
         D6  = `pval_s_drink' E6  = `pval_m_drink' ///
         F6  = `pval_s_ratio' G6  = `pval_m_ratio'
putexcel A7  = "Placebo mean" ///
         B7  = `mean_s_acc'   C7  = `mean_m_acc'   ///
         D7  = `mean_s_drink' E7  = `mean_m_drink' ///
         F7  = `mean_s_ratio' G7  = `mean_m_ratio'
putexcel A8  = "Placebo SD" ///
         B8  = `sd_s_acc'   C8  = `sd_m_acc'   ///
         D8  = `sd_s_drink' E8  = `sd_m_drink' ///
         F8  = `sd_s_ratio' G8  = `sd_m_ratio'
putexcel A9  = "Placebo [p5]" ///
         B9  = `p5_s_acc'   C9  = `p5_m_acc'   ///
         D9  = `p5_s_drink' E9  = `p5_m_drink' ///
         F9  = `p5_s_ratio' G9  = `p5_m_ratio'
putexcel A10 = "Placebo [p95]" ///
         B10 = `p95_s_acc'   C10 = `p95_m_acc'   ///
         D10 = `p95_s_drink' E10 = `p95_m_drink' ///
         F10 = `p95_s_ratio' G10 = `p95_m_ratio'
putexcel A12 = "Note: IRR for PPML; OLS coef for crash incidence. Placebo p-value = proportion of 1,000 randomizations at least as extreme as real estimate."


/*==============================================================================
Step 6: Summarize event study results, export to Excel
==============================================================================*/

use `es_results', clear

* Convert PPML log coefficients to IRR
gen irr_acc   = exp(es_acc)
gen irr_drink = exp(es_drink)
* Ratio uses OLS, no conversion needed

* Summarize by time point: mean + 5th/95th percentiles
collapse                                        ///
    (mean)  mean_irr_acc   = irr_acc            ///
    (mean)  mean_irr_drink = irr_drink          ///
    (mean)  mean_ratio     = es_ratio           ///
    (p5)    p5_irr_acc     = irr_acc            ///
    (p5)    p5_irr_drink   = irr_drink          ///
    (p5)    p5_ratio       = es_ratio           ///
    (p95)   p95_irr_acc    = irr_acc            ///
    (p95)   p95_irr_drink  = irr_drink          ///
    (p95)   p95_ratio      = es_ratio           ///
    , by(time_point)

sort time_point

save "$result/placebo_event_study.dta", replace
export excel using "$result/placebo_event_study.xlsx", ///
    firstrow(varlabels) replace

/*------------------------------------------------------------------------------
  Output SI Table 2: Placebo event study summary table
  Content: 12 time points × 3 outcomes, each cell: mean [p5, p95]
  Note: PPML converted to IRR, OLS kept as original values
------------------------------------------------------------------------------*/

* Data is now in a collapsed 12×1 structure
* Variables: time_point  mean_irr_acc  mean_irr_drink  mean_ratio
*        p5/p95_irr_acc  p5/p95_irr_drink  p5/p95_ratio

* Export Excel
putexcel set "$result/placebo_SI_table2_EventStudy.xlsx", replace

* Table header
putexcel A1 = "Supplementary Table: Placebo Event Study Estimates (1,000 randomizations, ±90% CI = [p5, p95])"
putexcel A3 = "Month" ///
    B3 = "Drunk Driving Cases: Mean (IRR)" C3 = "Drunk Driving Cases: [p5]" D3 = "Drunk Driving Cases: [p95]" ///
    E3 = "Traffic Crashes: Mean (IRR)"     F3 = "Traffic Crashes: [p5]"     G3 = "Traffic Crashes: [p95]" ///
    H3 = "Crash Incidence: Mean (OLS)"     I3 = "Crash Incidence: [p5]"     J3 = "Crash Incidence: [p95]"

* Write row by row (12 time points)
local row = 4
forvalues t = 1/12 {
    local tp = time_point[`t']
    putexcel A`row' = (`tp') ///
        B`row' = (mean_irr_acc[`t'])   C`row' = (p5_irr_acc[`t'])   D`row' = (p95_irr_acc[`t'])   ///
        E`row' = (mean_irr_drink[`t']) F`row' = (p5_irr_drink[`t']) G`row' = (p95_irr_drink[`t']) ///
        H`row' = (mean_ratio[`t'])     I`row' = (p5_ratio[`t'])     J`row' = (p95_ratio[`t'])
    local row = `row' + 1
}

local note_row = `row' + 1
putexcel A`note_row' = "Note: Month = months relative to fake treatment onset. Reference period (month −1) omitted; coefficient set to 1 (IRR) or 0 (OLS)."
putexcel A`=`note_row'+1' = "IRR reported for PPML models (drunk driving cases, traffic crashes); OLS coefficient for crash incidence."
putexcel A`=`note_row'+2' = "90% CI constructed from the 5th and 95th percentiles of the 1,000 placebo distributions."

di as text "✓ SI Table 2 (Event Study) Excel 已写出: $result/placebo_SI_table2_EventStudy.xlsx"

* Also write txt version (for easy verification)
local outf2 "$result/placebo_SI_table2_EventStudy.txt"
cap file close fh2
file open fh2 using "`outf2'", write replace text

file write fh2 "Supplementary Table: Placebo Event Study (1,000 randomizations, ±90% CI = [p5, p95])" _n
file write fh2 "Month" _tab ///
    "acc_mean" _tab "acc_p5"  _tab "acc_p95"  _tab ///
    "drink_mean" _tab "drink_p5" _tab "drink_p95" _tab ///
    "ratio_mean" _tab "ratio_p5" _tab "ratio_p95" _n
file write fh2 "-------------------------------------------------------------------------------------" _n

forvalues t = 1/12 {
    local tp = time_point[`t']
    
    * Write month
    file write fh2 %3.0f (`tp')
	
    * Write the three values for drink
    file write fh2 _tab
    file write fh2 %6.3f (mean_irr_drink[`t'])
    file write fh2 _tab
    file write fh2 %6.3f (p5_irr_drink[`t'])
    file write fh2 _tab
    file write fh2 %6.3f (p95_irr_drink[`t'])
    
    * Write the three values for acc
    file write fh2 _tab
    file write fh2 %6.3f (mean_irr_acc[`t'])
    file write fh2 _tab
    file write fh2 %6.3f (p5_irr_acc[`t'])
    file write fh2 _tab
    file write fh2 %6.3f (p95_irr_acc[`t'])
    
    * Write the three values for ratio
    file write fh2 _tab
    file write fh2 %6.3f (mean_ratio[`t'])
    file write fh2 _tab
    file write fh2 %6.3f (p5_ratio[`t'])
    file write fh2 _tab
    file write fh2 %6.3f (p95_ratio[`t'])
    
    * New line
    file write fh2 _n
}

file write fh2 "-------------------------------------------------------------------------------------" _n
file write fh2 "Note: Reference period (month -1) set to 1 (IRR) / 0 (OLS)." _n
file close fh2

/*==============================================================================
Step 7: Export DID distribution to Excel (for Python plotting)
==============================================================================*/

use "$result/placebo_estimates_short_medium.dta", clear

gen p_s_acc_irr   = round(exp(p_s_acc),   0.001)
gen p_m_acc_irr   = round(exp(p_m_acc),   0.001)
gen p_s_drink_irr = round(exp(p_s_drink), 0.001)
gen p_m_drink_irr = round(exp(p_m_drink), 0.001)
gen p_s_ratio_r   = round(p_s_ratio,      0.001)
gen p_m_ratio_r   = round(p_m_ratio,      0.001)

keep iteration                              ///
     p_s_acc_irr   p_m_acc_irr             ///
     p_s_drink_irr p_m_drink_irr           ///
     p_s_ratio_r   p_m_ratio_r

label var p_s_drink_irr "Placebo short-term: drink (IRR)"
label var p_m_drink_irr "Placebo medium-term: drink (IRR)"
label var p_s_acc_irr   "Placebo short-term: acc (IRR)"
label var p_m_acc_irr   "Placebo medium-term: acc (IRR)"
label var p_s_ratio_r   "Placebo short-term: acc ratio (OLS)"
label var p_m_ratio_r   "Placebo medium-term: acc ratio (OLS)"

export excel using "$result/placebo_for_python.xlsx", ///
    firstrow(varlabels) replace
