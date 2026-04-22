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
* Analysis: Robusrness Checks
* 
********************************************************************************

* Data versioning
global mydate    "20260415"
global myversion "v1"


* ==================== Data preparation & file setup ====================
global base "I:\Academic\demo\TrafficDecline-RiskAnalysis\github"

global result_base    "$base\\data\\reg_results"
global level1         "$result_base\\${mydate}-${myversion}"
global result         "$level1\04_robustness_checks"

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

global datapath "$base\\data"
global maindata "working_data"


* ==================== Auxiliary Macro: Descriptive Statistics =================
capture program drop calc_desc
program define calc_desc
    qui sum acc         if treated==1 & saps_short==0 & saps_medium==0
    scalar t_acc   = r(mean)
    qui sum acc         if treated==0 & saps_short==0 & saps_medium==0
    scalar c_acc   = r(mean)
    qui sum drink       if treated==1 & saps_short==0 & saps_medium==0
    scalar t_drink = r(mean)
    qui sum drink       if treated==0 & saps_short==0 & saps_medium==0
    scalar c_drink = r(mean)
    qui sum acc_ratio__ if treated==1 & saps_short==0 & saps_medium==0
    scalar t_ratio = r(mean)
    qui sum acc_ratio__ if treated==0 & saps_short==0 & saps_medium==0
    scalar c_ratio = r(mean)
    qui tab city_id if treated==1
    scalar n_treat = r(r)
    qui tab city_id if treated==0
    scalar n_ctrl  = r(r)
end

* ==================== Data preparation routine ====================
capture program drop data_prepare
program define data_prepare
    foreach v in gdp car ///
        AQI_max CO_max NO2_max O3_max PM10_max PM2_5_max SO2_max ///
        AQI_min CO_min NO2_min O3_min PM10_min PM2_5_min SO2_min ///
        AQI_mean CO_mean NO2_mean O3_mean PM10_mean PM2_5_mean SO2_mean ///
        joblossing kouzhao pop {
        cap replace `v' = . if `v' == -1
    }
    xtset id year
    gen treated     = x
    gen post_       = month_01
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
    gen saps_short  = (year_whole >= 2040 & year_whole <= 2042)
    gen saps_medium = (year_whole >= 2043 & year_whole <= 2048)
end

* ==================== regression function ====================
capture program drop run_three
program define run_three
    args outfile is_first se_label

    * drink
    ppmlhdfe drink saps_short saps_medium, ///
        absorb(month#city_id nyear#city_id) cluster(city_id) irr nolog
    if "`is_first'" == "replace" {
        outreg2 using "$result/`outfile'.xls", replace se nocons lab dec(3) ///
            keep(saps_short saps_medium) eform ci level(95) noaster ///
            addtext(Month×City FE, YES, Year×City FE, YES, SE Clustered, `se_label') ///
            addstat("Observations", e(N), "Treated Cities", n_treat, ///
                    "Control Cities", n_ctrl, ///
                    "Baseline Mean (Treated)", t_drink, "Baseline Mean (Control)", c_drink, ///
                    "Pseudo R-squared", e(r2_p))
    }
    else {
        outreg2 using "$result/`outfile'.xls", append se nocons lab dec(3) ///
            keep(saps_short saps_medium) eform ci level(95) noaster ///
            addtext(Month×City FE, YES, Year×City FE, YES, SE Clustered, `se_label') ///
            addstat("Observations", e(N), "Treated Cities", n_treat, ///
                    "Control Cities", n_ctrl, ///
                    "Baseline Mean (Treated)", t_drink, "Baseline Mean (Control)", c_drink, ///
                    "Pseudo R-squared", e(r2_p))
    }

    * acc
    ppmlhdfe acc saps_short saps_medium, ///
        absorb(month#city_id nyear#city_id) cluster(city_id) irr nolog
    outreg2 using "$result/`outfile'.xls", append se nocons lab dec(3) ///
        keep(saps_short saps_medium) eform ci level(95) noaster ///
        addtext(Month×City FE, YES, Year×City FE, YES, SE Clustered, `se_label') ///
        addstat("Observations", e(N), "Treated Cities", n_treat, ///
                "Control Cities", n_ctrl, ///
                "Baseline Mean (Treated)", t_acc, "Baseline Mean (Control)", c_acc, ///
                "Pseudo R-squared", e(r2_p))

    * ratio
    reghdfe acc_ratio__ saps_short saps_medium, ///
        absorb(month#city_id nyear#city_id) cluster(city_id)
    outreg2 using "$result/`outfile'.xls", append se nocons lab dec(3) ///
        keep(saps_short saps_medium) ci level(95) noaster ///
        addtext(Month×City FE, YES, Year×City FE, YES, SE Clustered, `se_label') ///
        addstat("Observations", e(N), "Treated Cities", n_treat, ///
                "Control Cities", n_ctrl, ///
                "Baseline Mean (Treated)", t_ratio, "Baseline Mean (Control)", c_ratio, ///
                "R-squared", e(r2), "R-squared Within", e(r2_within))
end


* ============================================================
* 1. Population weighting
* ============================================================
use "$datapath/$maindata/all_vs_saps_pop_pop.dta", clear
data_prepare
calc_desc

bysort city_id: egen pop_mean = mean(pop)
replace pop = pop_mean if pop == .

* drink
ppmlhdfe drink saps_short saps_medium [pweight=pop], ///
    absorb(month#city_id nyear#city_id) cluster(city_id) irr nolog
outreg2 using "$result/robustness_1_weighted.xls", replace se nocons lab dec(3) ///
    keep(saps_short saps_medium) eform ci level(95) noaster ///
    addtext(Month×City FE, YES, Year×City FE, YES, SE Clustered, City Level, Weight, Population) ///
    addstat("Observations", e(N), "Treated Cities", n_treat, "Control Cities", n_ctrl, ///
            "Baseline Mean (Treated)", t_drink, "Baseline Mean (Control)", c_drink, ///
            "Pseudo R-squared", e(r2_p))

* acc
ppmlhdfe acc saps_short saps_medium [pweight=pop], ///
    absorb(month#city_id nyear#city_id) cluster(city_id) irr nolog
outreg2 using "$result/robustness_1_weighted.xls", append se nocons lab dec(3) ///
    keep(saps_short saps_medium) eform ci level(95) noaster ///
    addtext(Month×City FE, YES, Year×City FE, YES, SE Clustered, City Level, Weight, Population) ///
    addstat("Observations", e(N), "Treated Cities", n_treat, "Control Cities", n_ctrl, ///
            "Baseline Mean (Treated)", t_acc, "Baseline Mean (Control)", c_acc, ///
            "Pseudo R-squared", e(r2_p))

* ratio
reghdfe acc_ratio__ saps_short saps_medium [aweight=pop], ///
    absorb(month#city_id nyear#city_id) cluster(city_id)
outreg2 using "$result/robustness_1_weighted.xls", append se nocons lab dec(3) ///
    keep(saps_short saps_medium) ci level(95) noaster ///
    addtext(Month×City FE, YES, Year×City FE, YES, SE Clustered, City Level, Weight, Population) ///
    addstat("Observations", e(N), "Treated Cities", n_treat, "Control Cities", n_ctrl, ///
            "Baseline Mean (Treated)", t_ratio, "Baseline Mean (Control)", c_ratio, ///
            "R-squared", e(r2), "R-squared Within", e(r2_within))

di as text "✓ Check 1 completed"


* ============================================================
* 2. Including control variables
* ============================================================
use "$datapath/$maindata/all_vs_saps_pop_pop.dta", clear
data_prepare
calc_desc

foreach v in wind proc temp suns humi PM2_5_mean AQI_mean SO2_mean CO_mean O3_mean NO2_mean kouzhao newconfirm newcure newdeath ggcz driver joblossing gdp car holiday {
    cap bysort city_id: egen `v'_m = mean(`v')
    cap replace `v' = `v'_m if `v' == .
}

global  ctrl_social    "gdp car holiday ggcz driver joblossing"
global  ctrl_weather   "wind proc temp suns humi"
global  ctrl_air       "AQI_mean SO2_mean CO_mean O3_mean NO2_mean"
global  ctrl_covid     "kouzhao newconfirm"
global  ctrl_full      "gdp car holiday ggcz driver joblossing wind proc temp suns humi AQI_mean SO2_mean CO_mean O3_mean NO2_mean kouzhao newconfirm"

cap mkdir "$result/colldiag"

foreach grp in social weather air covid full {

    di as result _n "  >> ctrl_`grp'"
    local cvars ${ctrl_`grp'}
    local var_list "saps_short saps_medium `cvars'"
    
    * ===== Multicollinearity diagnosis: output VIF for all variables =====
    local allvars "saps_short saps_medium `cvars'"
    local vif_names ""
    local vif_values ""

    foreach var of local allvars {
        local others ""
        foreach v of local allvars {
            if "`v'" != "`var'" {
                local others `others' `v'
            }
        }
        if "`others'" == "" {
            local vif_val = 1
        }
        else {
            capture qui regress `var' `others'
            if _rc == 0 & e(r2) < 1 & !missing(e(r2)) {
                local vif_val = 1 / (1 - e(r2))
            }
            else {
                local vif_val = .
            }
        }
        local vif_names `vif_names' `var'
        local vif_values `vif_values' `vif_val'
    }

    di as text "=== VIF for each variables ==="
    local nvar : word count `vif_names'
    forvalues i = 1/`nvar' {
        local varname : word `i' of `vif_names'
        local vifval : word `i' of `vif_values'
        di "  `varname' : `vifval'"
    }

    * Extract VIF for Scenario II and Scenario III
    local vif_ss = .
    local vif_sm = .
    forvalues i = 1/`nvar' {
        local varname : word `i' of `vif_names'
        if "`varname'" == "saps_short" local vif_ss : word `i' of `vif_values'
        if "`varname'" == "saps_medium" local vif_sm : word `i' of `vif_values'
    }

    * Export VIF to Excel
    putexcel set "$result/colldiag/vif_all_`grp'.xlsx", replace
    putexcel A1 = "VIF for each variable"
    putexcel A2 = "Group: `grp'"
    putexcel A3 = "Date: $S_DATE $S_TIME"
    putexcel A5 = "Variable"
    putexcel B5 = "VIF"
    local row = 6
    forvalues i = 1/`nvar' {
        local varname : word `i' of `vif_names'
        local vifval : word `i' of `vif_values'
        putexcel A`row' = "`varname'"
        putexcel B`row' = `vifval'
        local row = `row' + 1
    }

    * Condition index
    qui corr `allvars'
    matrix C = r(C)
    local cond_idx = .
    capture {
        matrix eigensystem C lambda V_eig
        local lmax = lambda[1,1]
        local lmin = lambda[1, colsof(lambda)]
        if `lmin' > 0 {
            local cond_idx = sqrt(`lmax' / `lmin')
        }
    }
    putexcel D5 = "Condition Index"
    putexcel D6 = `cond_idx'
    di as text "Condition Index: `cond_idx'"

    * ===== regression =====
    * drink (PPML)
	di as text "here"
    ppmlhdfe drink saps_short saps_medium `cvars', ///
        absorb(month#city_id nyear#city_id) cluster(city_id) irr nolog
    outreg2 using "$result/robustness_2_ctrl_`grp'.xls", replace ///
        se nocons lab dec(3) keep(`var_list') eform ci level(95) noaster ///
        addtext(Month×City FE, YES, Year×City FE, YES, SE Clustered, City Level) ///
        addstat("Observations", e(N), "Treated Cities", n_treat, ///
                "Control Cities", n_ctrl, "Baseline Mean (Treated)", t_drink, ///
				"Baseline Mean (Control)", c_drink, ///
				"Pseudo R-squared", e(r2_p))

    * acc (PPML)
	di as text "here"
    ppmlhdfe acc saps_short saps_medium `cvars', ///
        absorb(month#city_id nyear#city_id) cluster(city_id) irr nolog
    outreg2 using "$result/robustness_2_ctrl_`grp'.xls", append ///
        se nocons lab dec(3) keep(`var_list') eform ci level(95) noaster ///
        addtext(Month×City FE, YES, Year×City FE, YES, SE Clustered, City Level) ///
        addstat("Observations", e(N), "Treated Cities", n_treat, ///
                "Control Cities", n_ctrl, "Baseline Mean (Treated)", t_acc, ///
                "Baseline Mean (Control)", c_drink, "Pseudo R-squared", e(r2_p))

    * ratio (OLS)
    reghdfe acc_ratio__ saps_short saps_medium `cvars', ///
        absorb(month#city_id nyear#city_id) cluster(city_id)
    outreg2 using "$result/robustness_2_ctrl_`grp'.xls", append ///
        se nocons lab dec(3) keep(`var_list') ci level(95) noaster ///
        addtext(Month×City FE, YES, Year×City FE, YES, SE Clustered, City Level) ///
        addstat("Observations", e(N), "Treated Cities", n_treat, ///
                "Control Cities", n_ctrl, "Baseline Mean (Treated)", t_ratio, ///
                "Baseline Mean (Control)", c_ratio, "R-squared", e(r2), ///
                "R-squared Within", e(r2_within))
}


di as text _n "✓ Check 2 completed"


* ============================================================
* 3. Sample selection robustness (top 1% / top 5% / provincial capitals)
* ============================================================
use "$datapath/$maindata/all_vs_saps_pop_pop.dta", clear
data_prepare

bysort city_id: egen mean_acc = mean(acc)
egen p99_acc = pctile(mean_acc), p(99)
egen p95_acc = pctile(mean_acc), p(95)

* 3a. Drop top 1%
preserve
keep if mean_acc < p99_acc
calc_desc
run_three "robustness_3_drop_top1pct" "replace" "City Level"
restore

* 3b. Drop top 5%
preserve
keep if mean_acc < p95_acc
calc_desc
run_three "robustness_3_drop_top5pct" "replace" "City Level"
restore

* 3c. Drop provincial capitals
preserve
gen is_capital = 0
bysort prov_code: egen min_id = min(id)
replace is_capital = 1 if id == min_id
drop if is_capital == 1
calc_desc
run_three "robustness_3_drop_capitals" "replace" "City Level"
restore

di as text "✓ Check 3 completed"


* ============================================================
* 4. Thresholds for missing values (5th / 10th /20th)
* ============================================================
local thpath "$datapath/thresholds"
local first_th = 1

foreach th in 5 10 20 {
    local thfile "`thpath'/`th'/all_vs_saps_pop_pop.dta"
    cap use "`thfile'", clear
    if _rc != 0 {
        di as error "File not found: th`th', skipped"
        continue
    }
    data_prepare
    calc_desc
    if `first_th' == 1 {
        run_three "robustness_4_threshold_th`th'" "replace" "City Level"
        local first_th = 0
    }
    else {
        run_three "robustness_4_threshold_th`th'" "replace" "City Level"
    }
}

di as text "✓ Check 4 completed"


* ============================================================
* 5. Alter the time window
* ============================================================
use "$datapath/$maindata/all_vs_saps_pop_pop.dta", clear
data_prepare

* Scenario I, latter period (2025–2039)
preserve
keep if year_whole >= 2025
calc_desc
run_three "robustness_5_window" "replace" "City Level"
restore

di as text "✓ Check 5 completed"


* ============================================================
* 6. Change the fixed effects specification
* ============================================================
use "$datapath/$maindata/all_vs_saps_pop_pop.dta", clear
data_prepare
calc_desc

* 6a. City FE + Year FE + Month FE (non-interacted)
ppmlhdfe drink saps_short saps_medium, ///
    absorb(city_id nyear month) cluster(city_id) irr nolog
outreg2 using "$result/robustness_6_fe_simple.xls", replace se nocons lab dec(3) ///
    keep(saps_short saps_medium) eform ci level(95) noaster ///
    addtext(City FE, YES, Year FE, YES, Month FE, YES, SE Clustered, City Level) ///
    addstat("Observations", e(N), "Treated Cities", n_treat, "Control Cities", n_ctrl, ///
            "Baseline Mean (Treated)", t_drink, "Baseline Mean (Control)", c_drink, ///
            "Pseudo R-squared", e(r2_p))

ppmlhdfe acc saps_short saps_medium, ///
    absorb(city_id nyear month) cluster(city_id) irr nolog
outreg2 using "$result/robustness_6_fe_simple.xls", append se nocons lab dec(3) ///
    keep(saps_short saps_medium) eform ci level(95) noaster ///
    addtext(City FE, YES, Year FE, YES, Month FE, YES, SE Clustered, City Level) ///
    addstat("Observations", e(N), "Treated Cities", n_treat, "Control Cities", n_ctrl, ///
            "Baseline Mean (Treated)", t_acc, "Baseline Mean (Control)", c_acc, ///
            "Pseudo R-squared", e(r2_p))

reghdfe acc_ratio__ saps_short saps_medium, ///
    absorb(city_id nyear month) cluster(city_id)
outreg2 using "$result/robustness_6_fe_simple.xls", append se nocons lab dec(3) ///
    keep(saps_short saps_medium) ci level(95) noaster ///
    addtext(City FE, YES, Year FE, YES, Month FE, YES, SE Clustered, City Level) ///
    addstat("Observations", e(N), "Treated Cities", n_treat, "Control Cities", n_ctrl, ///
            "Baseline Mean (Treated)", t_ratio, "Baseline Mean (Control)", c_ratio, ///
            "R-squared", e(r2), "R-squared Within", e(r2_within))

* 6b. Province–year fixed effects
ppmlhdfe drink saps_short saps_medium, ///
    absorb(month#city_id prov_code#nyear) cluster(city_id) irr nolog
outreg2 using "$result/robustness_6_fe_prov.xls", append se nocons lab dec(3) ///
    keep(saps_short saps_medium) eform ci level(95) noaster ///
    addtext(Month×City FE, YES, Prov×Year FE, YES, SE Clustered, City Level) ///
    addstat("Observations", e(N), "Treated Cities", n_treat, "Control Cities", n_ctrl, ///
            "Baseline Mean (Treated)", t_drink, "Baseline Mean (Control)", c_drink, ///
            "Pseudo R-squared", e(r2_p))

ppmlhdfe acc saps_short saps_medium, ///
    absorb(month#city_id prov_code#nyear) cluster(city_id) irr nolog
outreg2 using "$result/robustness_6_fe_prov.xls", append se nocons lab dec(3) ///
    keep(saps_short saps_medium) eform ci level(95) noaster ///
    addtext(Month×City FE, YES, Prov×Year FE, YES, SE Clustered, City Level) ///
    addstat("Observations", e(N), "Treated Cities", n_treat, "Control Cities", n_ctrl, ///
            "Baseline Mean (Treated)", t_acc, "Baseline Mean (Control)", c_acc, ///
            "Pseudo R-squared", e(r2_p))

reghdfe acc_ratio__ saps_short saps_medium, ///
    absorb(month#city_id prov_code#nyear) cluster(city_id)
outreg2 using "$result/robustness_6_fe_prov.xls", append se nocons lab dec(3) ///
    keep(saps_short saps_medium) ci level(95) noaster ///
    addtext(Month×City FE, YES, Prov×Year FE, YES, SE Clustered, City Level) ///
    addstat("Observations", e(N), "Treated Cities", n_treat, "Control Cities", n_ctrl, ///
            "Baseline Mean (Treated)", t_ratio, "Baseline Mean (Control)", c_ratio, ///
            "R-squared", e(r2), "R-squared Within", e(r2_within))

di as text "✓ Check 6 completed"


* ============================================================
* 7. Bootstrap
* ============================================================
use "$datapath/$maindata/all_vs_saps_pop_pop.dta", clear
data_prepare

* ── Main regression, cluster-robust SE ──
ppmlhdfe drink saps_short saps_medium, ///
    absorb(month#city_id nyear#city_id) cluster(city_id) nolog
scalar b_dr_s       = _b[saps_short]
scalar se_cl_dr_s   = _se[saps_short]
scalar b_dr_m       = _b[saps_medium]
scalar se_cl_dr_m   = _se[saps_medium]

ppmlhdfe acc saps_short saps_medium, ///
    absorb(month#city_id nyear#city_id) cluster(city_id) nolog
scalar b_acc_s      = _b[saps_short]
scalar se_cl_acc_s  = _se[saps_short]
scalar b_acc_m      = _b[saps_medium]
scalar se_cl_acc_m  = _se[saps_medium]

reghdfe acc_ratio__ saps_short saps_medium, ///
    absorb(month#city_id nyear#city_id) cluster(city_id)
scalar b_rt_s       = _b[saps_short]
scalar se_cl_rt_s   = _se[saps_short]
scalar b_rt_m       = _b[saps_medium]
scalar se_cl_rt_m   = _se[saps_medium]

* ── Bootstrap, compute SE ──
bootstrap _b[saps_short] _b[saps_medium], reps(200) seed(20240101) notable: ///
    ppmlhdfe drink saps_short saps_medium, absorb(month#city_id nyear#city_id) nolog
scalar se_bs_dr_s  = _se[_bs_1]
scalar se_bs_dr_m  = _se[_bs_2]

bootstrap _b[saps_short] _b[saps_medium], reps(200) seed(20240101) notable: ///
    ppmlhdfe acc saps_short saps_medium, absorb(month#city_id nyear#city_id) nolog
scalar se_bs_acc_s = _se[_bs_1]
scalar se_bs_acc_m = _se[_bs_2]

bootstrap _b[saps_short] _b[saps_medium], reps(200) seed(20240101) notable: ///
    reghdfe acc_ratio__ saps_short saps_medium, absorb(month#city_id nyear#city_id)
scalar se_bs_rt_s  = _se[_bs_1]
scalar se_bs_rt_m  = _se[_bs_2]

* ── Output comparison table as plain text (copy-paste ready for manuscript) ──
local outf "$result/robustness_7_bootstrap.txt"
cap file close fh
file open fh using "`outf'", write replace text

file write fh "Bootstrap SE vs Cluster SE Comparison (200 replications)" _n
file write fh "==========================================================" _n
file write fh _n

* Table header
file write fh "Outcome    Variable      IRR/Coef   Cluster SE   p(cluster)   Bootstrap SE   p(bootstrap)" _n
file write fh "--------------------------------------------------------------------------" _n

* drink - Scenario II
file write fh "drink      saps_short    " %6.3f (exp(b_dr_s)) "      " ///
    %7.4f (se_cl_dr_s) "      " ///
    %7.4f (2*(1-normal(abs(b_dr_s/se_cl_dr_s)))) "        " ///
    %7.4f (se_bs_dr_s) "          " ///
    %7.4f (2*(1-normal(abs(b_dr_s/se_bs_dr_s)))) _n

* drink - Scenario III
file write fh "drink      saps_medium   " %6.3f (exp(b_dr_m)) "      " ///
    %7.4f (se_cl_dr_m) "      " ///
    %7.4f (2*(1-normal(abs(b_dr_m/se_cl_dr_m)))) "        " ///
    %7.4f (se_bs_dr_m) "          " ///
    %7.4f (2*(1-normal(abs(b_dr_m/se_bs_dr_m)))) _n

file write fh "--------------------------------------------------------------------------" _n

* acc - Scenario II
file write fh "acc        saps_short    " %6.3f (exp(b_acc_s)) "      " ///
    %7.4f (se_cl_acc_s) "      " ///
    %7.4f (2*(1-normal(abs(b_acc_s/se_cl_acc_s)))) "        " ///
    %7.4f (se_bs_acc_s) "          " ///
    %7.4f (2*(1-normal(abs(b_acc_s/se_bs_acc_s)))) _n

* acc - Scenario III
file write fh "acc        saps_medium   " %6.3f (exp(b_acc_m)) "      " ///
    %7.4f (se_cl_acc_m) "      " ///
    %7.4f (2*(1-normal(abs(b_acc_m/se_cl_acc_m)))) "        " ///
    %7.4f (se_bs_acc_m) "          " ///
    %7.4f (2*(1-normal(abs(b_acc_m/se_bs_acc_m)))) _n

file write fh "--------------------------------------------------------------------------" _n

* acc_ratio__ - Scenario II
file write fh "acc_ratio  saps_short    " %6.3f (b_rt_s) "      " ///
    %7.4f (se_cl_rt_s) "      " ///
    %7.4f (2*(1-normal(abs(b_rt_s/se_cl_rt_s)))) "        " ///
    %7.4f (se_bs_rt_s) "          " ///
    %7.4f (2*(1-normal(abs(b_rt_s/se_bs_rt_s)))) _n

* acc_ratio__ - Scenario III
file write fh "acc_ratio  saps_medium   " %6.3f (b_rt_m) "      " ///
    %7.4f (se_cl_rt_m) "      " ///
    %7.4f (2*(1-normal(abs(b_rt_m/se_cl_rt_m)))) "        " ///
    %7.4f (se_bs_rt_m) "          " ///
    %7.4f (2*(1-normal(abs(b_rt_m/se_bs_rt_m)))) _n

file write fh "--------------------------------------------------------------------------" _n
file write fh "Note: IRR reported for PPML models; coefficient for OLS (acc_ratio)." _n
file write fh "      SE in log scale for PPML. p-values two-sided." _n

file close fh
di as result "✓ Check 7 completed"


* ============================================================
* 8. Change clustering level (province level)
* ============================================================
use "$datapath/$maindata/all_vs_saps_pop_pop.dta", clear
data_prepare
calc_desc
run_three "robustness_8_cluster_prov" "replace" "Province Level"
di as text "✓ Check 8 completed"

