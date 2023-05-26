********************************************************************************
*
*	Do-file:		data_preparation_and_descriptives.do
*
*	Project:		sotrovimab-and-Paxlovid
*
*	Programmed by:	Bang Zheng
*
*	Data used:		output/input.csv
*
*	Data created:	output/main.dta  (main analysis dataset)
*
*	Other output:	logs/data_preparation.log
*
********************************************************************************
*
*	Purpose: This do-file creates the variables required for the 
*			 main analysis and saves into Stata dataset, and describes 
*            variables by drug groups.
*  
********************************************************************************

* Open a log file
cap log close
log using ./logs/data_preparation_feasibility, replace t
clear

* import dataset
import delimited ./output/input_feasibility.csv, delimiter(comma) varnames(1) case(preserve) 
describe
codebook

*  Convert strings to dates  *
foreach var of varlist  sotrovimab_covid_therapeutics molnupiravir_covid_therapeutics paxlovid_covid_therapeutics remdesivir_covid_therapeutics	///
        casirivimab_covid_therapeutics tocilizumab_covid_therapeutics sarilumab_covid_therapeutics date_treated  ///
        sotrovimab_covid_hosp paxlovid_covid_hosp molnupiravir_covid_hosp remdesivir_covid_hosp casirivimab_covid_hosp tocilizumab_covid_hosp sarilumab_covid_hosp ///
		date_treated_hosp start_date death_with_covid_date death_with_covid_underly_date death_date {
  capture confirm string variable `var'
  if _rc==0 {
  rename `var' a
  gen `var' = date(a, "YMD")
  drop a
  format %td `var'
  sum `var',f
  sum `var' if covid_therapeutics!="" ,f
  sum `var' if covid_therapeutics_hosp!="" ,f
  }
}
*the following date variables had no observation*
*hiv_aids_nhsd_icd10
*transplant_all_y_codes_opcs4
*transplant_thymus_opcs4
*transplant_conjunctiva_y_code_op
*transplant_conjunctiva_opcs4
*transplant_stomach_opcs4
*transplant_ileum_1_Y_codes_opcs4
*transplant_ileum_2_Y_codes_opcs4
*transplant_ileum_1_opcs4
*transplant_ileum_2_opcs4

tab covid_therapeutics
tab registered_treated
tab covid_therapeutics_hosp
tab registered_treated_hosp
tab high_risk_cohort_covid_therapeut
tab high_risk_cohort_covid_therapeut if covid_therapeutics!=""
tab high_risk_cohort_covid_therapeut if covid_therapeutics_hosp!=""
count  if tocilizumab_covid_hosp!=.&death_with_covid_date!=.
count  if sarilumab_covid_hosp!=.&death_with_covid_date!=.
count  if tocilizumab_covid_hosp!=.&death_date!=.
count  if sarilumab_covid_hosp!=.&death_date!=.

*check hosp/death event date range*
*codebook covid_hosp_outcome_date2 hospitalisation_outcome_date2 death_date

*exclusion criteria*

log close




