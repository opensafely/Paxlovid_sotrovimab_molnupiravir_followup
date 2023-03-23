********************************************************************************
*
*	Do-file:		data_preparation_and_descriptives_raw_rate.do
*
*	Project:		sotrovimab-and-Paxlovid
*
*	Programmed by:	Bang Zheng
*
*	Data used:		output/input_raw_rate.csv
*
*	Output:	        logs/data_preparation_raw_rate.log
*
****************************************************************************************************************************************************************
*
*	Purpose: The first section creates the variables required for the untreated group.
*  
****************************************************************************************************************************************************************

*set start date of the study period*
local start_MDY= "12,16,2021"
local start_DMY: display %td mdy(`start_MDY')

*create result table*
postfile mytab str25 (Cohorts N COVID_hosp_30_day COVID_mortality_30_day all_mortality_30_day Age_mean Female_prop Vaccination_3_or_more) using "./output/table.dta", replace


* Open a log file
cap log close
log using ./logs/data_preparation_raw_rate, replace t
clear

* import dataset
import delimited ./output/input_raw_rate.csv, delimiter(comma) varnames(1) case(preserve) 
describe
keep if registered_eligible==1
drop if cancer_opensafely_snomed_new==""&immunosuppresant_drugs_nhsd==""&oral_steroid_drugs_nhsd==""&immunosupression_nhsd_new==""&solid_organ_transplant_nhsd_new==""&downs_syndrome_nhsd==""&haematological_disease_nhsd==""&ckd_stage_5_nhsd==""&liver_disease_nhsd==""&hiv_aids_nhsd==""&multiple_sclerosis_nhsd==""&motor_neurone_disease_nhsd==""&myasthenia_gravis_nhsd==""&huntingtons_disease_nhsd=="" 



*  Convert strings to dates  *
foreach var of varlist  sotrovimab_covid_therapeutics molnupiravir_covid_therapeutics paxlovid_covid_therapeutics remdesivir_covid_therapeutics	///
        casirivimab_covid_therapeutics date_treated start_date ///
        covid_test_positive_date  primary_covid_hospital_discharge primary_covid_hospital_admission ///
	   any_covid_hospital_discharge_dat any_covid_hospital_admission_dat death_date dereg_date  ///
	   cancer_opensafely_snomed_new   immunosuppresant_drugs_nhsd ///
	   oral_steroid_drugs_nhsd  immunosupression_nhsd_new   solid_organ_transplant_nhsd_new  ///
	   covid_hosp_outcome_date0 covid_hosp_outcome_date1 covid_hosp_outcome_date2 covid_hosp_discharge_date0 covid_hosp_discharge_date1 covid_hosp_discharge_date2 ///
	   death_with_covid_date death_with_covid_underly_date covid_hosp_outcome_date0T covid_hosp_outcome_date1T covid_hosp_outcome_date2T covid_hosp_discharge_date0T ///
	   covid_hosp_discharge_date1T covid_hosp_discharge_date2T death_with_covid_dateT  death_with_covid_underly_dateT death_dateT dereg_dateT  ///
	   downs_syndrome_nhsd haematological_disease_nhsd ckd_stage_5_nhsd liver_disease_nhsd hiv_aids_nhsd  ///
	   multiple_sclerosis_nhsd motor_neurone_disease_nhsd myasthenia_gravis_nhsd huntingtons_disease_nhsd advanced_decompensated_cirrhosis decompensated_cirrhosis_icd10 ///
	   ascitic_drainage_snomed  ckd_stages_3_5 ckd_primis_stage_date ckd3_icd10 ckd4_icd10 ckd5_icd10 dialysis dialysis_icd10 dialysis_procedure kidney_transplant kidney_transplant_icd10 ///
	   kidney_transplant_procedure creatinine_ctv3_date creatinine_snomed_date creatinine_short_snomed_date eGFR_record_date eGFR_short_record_date liver_disease_nhsd_icd10 ///
	   solid_organ_transplant_snomed drugs_do_not_use drugs_consider_risk  {
  capture confirm string variable `var'
  if _rc==0 {
  rename `var' a
  gen `var' = date(a, "YMD")
  drop a
  format %td `var'
  }
}


*check hosp/death event date range*
codebook  covid_hosp_outcome_date2 death_date
sum covid_hosp_outcome_date2  
local end_DMY: disp %td r(max)-31


*exclusion criteria*
count if date_treated!=.
count if date_treated!=.&date_treated<=(covid_test_positive_date+30)
count if date_treated!=.&date_treated<(covid_test_positive_date-30)
drop if sotrovimab_covid_therapeutics!=. | paxlovid_covid_therapeutics!=. | molnupiravir_covid_therapeutics!=. |remdesivir_covid_therapeutics!=. |casirivimab_covid_therapeutics!=.
sum age,de
keep if age>=18 & age<110
tab sex,m
keep if sex=="F"|sex=="M"
keep if has_died==0
tab covid_test_positive covid_positive_previous_30_days,m
keep if covid_test_positive==1 & covid_positive_previous_30_days==0
drop if primary_covid_hospital_discharge!=.|primary_covid_hospital_admission!=.
drop if any_covid_hospital_admission_dat!=.|any_covid_hospital_discharge_dat!=.
*restrict study period! NOTE: end of study period should be 30 days earlier than the latest hosp event date in current extraction*
keep if covid_test_positive_date>=mdy(`start_MDY')&covid_test_positive_date<=date("`end_DMY'", "DMY")
sum covid_test_positive_date
*drop if stp==""

assert start_date==covid_test_positive_date
replace start_date=covid_test_positive_date if start_date!=covid_test_positive_date

*correcting COVID hosp events: further ignore any day cases or sotro initiators who had COVID hosp record with mab procedure codes on Day 0 or 1 *
count if covid_hosp_outcome_date0!=start_date&covid_hosp_outcome_date0!=.
count if covid_hosp_outcome_date1!=(start_date+1)&covid_hosp_outcome_date1!=.
*check if any patient discharged in AM and admitted in PM*
count if primary_covid_hospital_admission!=.&primary_covid_hospital_discharge==.
count if primary_covid_hospital_admission==.&primary_covid_hospital_discharge!=.
count if primary_covid_hospital_admission!=.&primary_covid_hospital_discharge!=.&primary_covid_hospital_admission==primary_covid_hospital_discharge
count if primary_covid_hospital_admission!=.&primary_covid_hospital_discharge!=.&primary_covid_hospital_admission<primary_covid_hospital_discharge
count if primary_covid_hospital_admission!=.&primary_covid_hospital_discharge!=.&primary_covid_hospital_admission>primary_covid_hospital_discharge
*ignore day cases and mab procedures in day 0/1*
replace covid_hosp_outcome_date0=. if covid_hosp_outcome_date0==covid_hosp_discharge_date0&covid_hosp_outcome_date0!=.
replace covid_hosp_outcome_date1=. if covid_hosp_outcome_date1==covid_hosp_discharge_date1&covid_hosp_outcome_date1!=.
*replace covid_hosp_outcome_date0=. if covid_hosp_outcome_date0==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure!=.&drug==1
*replace covid_hosp_outcome_date1=. if covid_hosp_outcome_date1==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure!=.&drug==1

gen covid_hospitalisation_outcome_da=covid_hosp_outcome_date2
replace covid_hospitalisation_outcome_da=covid_hosp_outcome_date1 if covid_hosp_outcome_date1!=.
replace covid_hospitalisation_outcome_da=covid_hosp_outcome_date0 if covid_hosp_outcome_date0!=.

gen days_to_covid_admission=covid_hospitalisation_outcome_da-covid_test_positive_date if covid_hospitalisation_outcome_da!=.

*ignore and censor day cases on or after day 2 from this analysis*
*ignore and censor admissions for mab procedure >= day 2 and with same-day or 1-day discharge*
*gen covid_hosp_date_day_cases_mab=covid_hospitalisation_outcome_da if covid_hosp_outcome_date2==covid_hosp_discharge_date2&covid_hosp_outcome_date2!=.&days_to_covid_admission>=2
*replace covid_hosp_date_day_cases_mab=covid_hospitalisation_outcome_da if covid_hosp_outcome_date2==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure!=.&days_to_covid_admission>=2&(covid_hosp_discharge_date2-covid_hosp_outcome_date2)<=1&drug==1
drop if covid_hosp_outcome_date2==covid_hosp_discharge_date2&covid_hosp_outcome_date2!=.&days_to_covid_admission>=2
*replace covid_hospitalisation_outcome_da=. if covid_hosp_outcome_date2==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure!=.&days_to_covid_admission>=2&(covid_hosp_discharge_date2-covid_hosp_outcome_date2)<=1&drug==1
*check hosp_admission_method*
*by drug days_to_covid_admission, sort: count if covid_hospitalisation_outcome_da!=covid_hosp_date_emergency&covid_hospitalisation_outcome_da!=.

*censor*
*capture and exclude COVID-hospital admission/death on the start date
drop if covid_test_positive_date>=covid_hospitalisation_outcome_da| covid_test_positive_date>=death_with_covid_date|covid_test_positive_date>=death_date|covid_test_positive_date>=dereg_date
drop if dereg_date>=covid_test_positive_date&dereg_date<=(covid_test_positive_date+30)



*covariates* 
*10 high risk groups: downs_syndrome, solid_cancer, haematological_disease, renal_disease, liver_disease, imid, 
*immunosupression, hiv_aids, solid_organ_transplant, rare_neurological_conditions, high_risk_group_combined	

replace oral_steroid_drugs_nhsd=. if oral_steroid_drug_nhsd_3m_count < 2 & oral_steroid_drug_nhsd_12m_count < 4
gen imid_nhsd=min(oral_steroid_drugs_nhsd, immunosuppresant_drugs_nhsd)
gen rare_neuro_nhsd = min(multiple_sclerosis_nhsd, motor_neurone_disease_nhsd, myasthenia_gravis_nhsd, huntingtons_disease_nhsd)
*high risk group only based on codelists*
gen downs_syndrome=(downs_syndrome_nhsd<=start_date)
gen solid_cancer_new=(cancer_opensafely_snomed_new<=start_date)
gen haema_disease=( haematological_disease_nhsd <=start_date)
gen renal_disease=( ckd_stage_5_nhsd <=start_date)
gen liver_disease=( liver_disease_nhsd <=start_date)
gen imid=( imid_nhsd <=start_date)
gen immunosupression_new=( immunosupression_nhsd_new <=start_date)
gen hiv_aids=( hiv_aids_nhsd <=start_date)
gen solid_organ_new=( solid_organ_transplant_nhsd_new<=start_date)
gen rare_neuro=( rare_neuro_nhsd <=start_date)
gen high_risk_group_new=(( downs_syndrome + solid_cancer_new + haema_disease + renal_disease + liver_disease + imid + immunosupression_new + hiv_aids + solid_organ_new + rare_neuro )>0)
tab high_risk_group_new,m
keep if high_risk_group_new==1
count 
local N_untreated=r(N) 


*demo*
sum age,de
local age_untreated=r(mean)
gen age_group3=(age>=40)+(age>=60)
label define age_group3_Paxlovid 0 "18-39" 1 "40-59" 2 ">=60" 
label values age_group3 age_group3_Paxlovid
tab age_group3,m
egen age_5y_band=cut(age), at(18,25,30,35,40,45,50,55,60,65,70,75,80,85,110) label
tab age_5y_band,m

tab sex,m
rename sex sex_str
gen sex=0 if sex_str=="M"
replace sex=1 if sex_str=="F"
label define sex_Paxlovid 0 "Male" 1 "Female"
label values sex sex_Paxlovid
sum sex
local female_untreated=r(mean)


*vac and variant*
tab vaccination_status,m
rename vaccination_status vaccination_status_g6
gen vaccination_status=0 if vaccination_status_g6=="Un-vaccinated"|vaccination_status_g6=="Un-vaccinated (declined)"
replace vaccination_status=1 if vaccination_status_g6=="One vaccination"
replace vaccination_status=2 if vaccination_status_g6=="Two vaccinations"
replace vaccination_status=3 if vaccination_status_g6=="Three vaccinations"
replace vaccination_status=4 if vaccination_status_g6=="Four or more vaccinations"
label define vac_Paxlovid 0 "Un-vaccinated" 1 "One vaccination" 2 "Two vaccinations" 3 "Three vaccinations" 4 "Four or more vaccinations"
label values vaccination_status vac_Paxlovid
gen vaccination_3=1 if vaccination_status==3|vaccination_status==4
replace vaccination_3=0 if vaccination_status<3
tab vaccination_3
sum vaccination_3
local vac_3_untreated=r(mean)

tab variant_recorded ,m
tab symptomatic_covid_test ,m



*define outcome and follow-up time*
gen start_date_30=covid_test_positive_date+30
*30-day COVID hosp*
gen covid_hospitalisation_30day=(covid_hospitalisation_outcome_da!=.&covid_hospitalisation_outcome_da<=start_date_30)
tab covid_hospitalisation_30day,m
count if covid_hospitalisation_30day==1
local cov_hosp_untreated_n=r(N)
sum covid_hospitalisation_30day
local cov_hosp_untreated=r(mean)
*tab  covid_hospitalisation_30day if start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022), m
*tab  covid_hospitalisation_30day if start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022), m
*tab  covid_hospitalisation_30day if start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022),m
*30-day COVID death*
gen covid_death_30day=(death_with_covid_date!=.&death_with_covid_date<=start_date_30)
tab covid_death_30day,m
count if covid_death_30day==1
local cov_death_untreated_n=r(N)
sum covid_death_30day
local cov_death_untreated=r(mean)
*tab  covid_death_30day if start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022),m
*tab  covid_death_30day if start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022), m
*tab  covid_death_30day if start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022), m
*30-day all-cause death*
gen death_30day=(death_date!=.&death_date<=start_date_30)
tab death_30day,m
count if death_30day==1
local death_untreated_n=r(N)
sum death_30day
local death_untreated=r(mean)
*tab  death_30day if start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022), m
*tab  death_30day if start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022), m
*tab  death_30day if start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022), m




*exclude those with contraindications for Pax*
replace ckd_primis_stage=. if ckd_primis_stage_date>start_date
*egfr: adapted from https://github.com/opensafely/COVID-19-vaccine-breakthrough/blob/updates-feb/analysis/data_process.R*
tab creatinine_operator_ctv3,m
replace creatinine_ctv3 = . if !inrange(creatinine_ctv3, 20, 3000)| creatinine_ctv3_date>start_date
tab creatinine_operator_ctv3 if creatinine_ctv3!=.,m
replace creatinine_ctv3 = creatinine_ctv3/88.4
gen min_creatinine_ctv3=.
replace min_creatinine_ctv3 = (creatinine_ctv3/0.7)^-0.329 if sex==1
replace min_creatinine_ctv3 = (creatinine_ctv3/0.9)^-0.411 if sex==0
replace min_creatinine_ctv3 = 1 if min_creatinine_ctv3<1
gen max_creatinine_ctv3=.
replace max_creatinine_ctv3 = (creatinine_ctv3/0.7)^-1.209 if sex==1
replace max_creatinine_ctv3 = (creatinine_ctv3/0.9)^-1.209 if sex==0
replace max_creatinine_ctv3 = 1 if max_creatinine_ctv3>1
gen egfr_creatinine_ctv3 = min_creatinine_ctv3*max_creatinine_ctv3*141*(0.993^age_creatinine_ctv3) if age_creatinine_ctv3>0&age_creatinine_ctv3<=120
replace egfr_creatinine_ctv3 = egfr_creatinine_ctv3*1.018 if sex==1

tab creatinine_operator_snomed,m
tab creatinine_operator_snomed if creatinine_snomed!=.,m
replace creatinine_snomed = . if !inrange(creatinine_snomed, 20, 3000)| creatinine_snomed_date>start_date
replace creatinine_snomed_date = creatinine_short_snomed_date if missing(creatinine_snomed)
replace creatinine_operator_snomed = creatinine_operator_short_snomed if missing(creatinine_snomed)
replace age_creatinine_snomed = age_creatinine_short_snomed if missing(creatinine_snomed)
replace creatinine_snomed = creatinine_short_snomed if missing(creatinine_snomed)
replace creatinine_snomed = . if !inrange(creatinine_snomed, 20, 3000)| creatinine_snomed_date>start_date
replace creatinine_snomed = creatinine_snomed/88.4
gen min_creatinine_snomed=.
replace min_creatinine_snomed = (creatinine_snomed/0.7)^-0.329 if sex==1
replace min_creatinine_snomed = (creatinine_snomed/0.9)^-0.411 if sex==0
replace min_creatinine_snomed = 1 if min_creatinine_snomed<1
gen max_creatinine_snomed=.
replace max_creatinine_snomed = (creatinine_snomed/0.7)^-1.209 if sex==1
replace max_creatinine_snomed = (creatinine_snomed/0.9)^-1.209 if sex==0
replace max_creatinine_snomed = 1 if max_creatinine_snomed>1
gen egfr_creatinine_snomed = min_creatinine_snomed*max_creatinine_snomed*141*(0.993^age_creatinine_snomed) if age_creatinine_snomed>0&age_creatinine_snomed<=120
replace egfr_creatinine_snomed = egfr_creatinine_snomed*1.018 if sex==1

tab eGFR_operator if eGFR_record!=.,m
tab eGFR_short_operator if eGFR_short_record!=.,m
count if (egfr_creatinine_ctv3<60&creatinine_operator_ctv3!="<")|(egfr_creatinine_snomed<60&creatinine_operator_snomed!="<")|(eGFR_record<60&eGFR_record>0&eGFR_operator!=">"&eGFR_operator!=">=")|(eGFR_short_record<60&eGFR_short_record>0&eGFR_short_operator!=">"&eGFR_short_operator!=">=")

*drug interactions*
count if drugs_do_not_use<=start_date
count if drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-3*365.25)
count if drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-365.25)
count if drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-180)
count if drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-90)
count if drugs_consider_risk<=start_date
count if drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-3*365.25)
count if drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-365.25)
count if drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-180)
count if drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-90)
gen drugs_do_not_use_contra=(drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-180))
gen drugs_consider_risk_contra=(drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-180))


drop if solid_organ_new==1|solid_organ_transplant_snomed<=start_date
drop if advanced_decompensated_cirrhosis<=start_date|decompensated_cirrhosis_icd10<=start_date|ascitic_drainage_snomed<=start_date|liver_disease_nhsd_icd10<=start_date
drop if renal_disease==1|ckd_stages_3_5<=start_date|ckd_primis_stage==3|ckd_primis_stage==4|ckd_primis_stage==5|ckd3_icd10<=start_date|ckd4_icd10<=start_date|ckd5_icd10<=start_date
drop if kidney_transplant<=start_date|kidney_transplant_icd10<=start_date|kidney_transplant_procedure<=start_date
drop if dialysis<=start_date|dialysis_icd10<=start_date|dialysis_procedure<=start_date
drop if (egfr_creatinine_ctv3<60&creatinine_operator_ctv3!="<")|(egfr_creatinine_snomed<60&creatinine_operator_snomed!="<")|(eGFR_record<60&eGFR_record>0&eGFR_operator!=">"&eGFR_operator!=">=")|(eGFR_short_record<60&eGFR_short_record>0&eGFR_short_operator!=">"&eGFR_short_operator!=">=")
*drop if drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-365.25)
*drop if drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-365.25)
drop if drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-180)
*drop if drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-180)
count 
local N_untreated_no=r(N) 
sum age,de
local age_untreated_no=r(mean)
sum sex
local female_untreated_no=r(mean)
sum vaccination_3
local vac_3_untreated_no=r(mean)


*30-day COVID hosp*
tab covid_hospitalisation_30day,m
count if covid_hospitalisation_30day==1
local cov_hosp_untreated_no_n=r(N)
sum covid_hospitalisation_30day
local cov_hosp_untreated_no=r(mean)
*tab  covid_hospitalisation_30day if start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022), m
*tab  covid_hospitalisation_30day if start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022), m
*tab  covid_hospitalisation_30day if start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022),m
*30-day COVID death*
tab covid_death_30day,m
count if covid_death_30day==1
local cov_death_untreated_no_n=r(N)
sum covid_death_30day
local cov_death_untreated_no=r(mean)
*tab  covid_death_30day if start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022),m
*tab  covid_death_30day if start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022), m
*tab  covid_death_30day if start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022), m
*30-day all-cause death*
tab death_30day,m
count if death_30day==1
local death_untreated_no_n=r(N)
sum death_30day
local death_untreated_no=r(mean)
*tab  death_30day if start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022), m
*tab  death_30day if start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022), m
*tab  death_30day if start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022), m








****************************************************************************************************************************************************************
*
*	Purpose: The second section creates the variables required for the treated group.
*  
****************************************************************************************************************************************************************

* import dataset
import delimited ./output/input_raw_rate.csv, delimiter(comma) varnames(1) case(preserve) 
keep if registered_treated==1
drop if cancer_opensafely_snomed_new==""&immunosuppresant_drugs_nhsd==""&oral_steroid_drugs_nhsd==""&immunosupression_nhsd_new==""&solid_organ_transplant_nhsd_new==""&downs_syndrome_nhsd==""&haematological_disease_nhsd==""&ckd_stage_5_nhsd==""&liver_disease_nhsd==""&hiv_aids_nhsd==""&multiple_sclerosis_nhsd==""&motor_neurone_disease_nhsd==""&myasthenia_gravis_nhsd==""&huntingtons_disease_nhsd=="" 

*  Convert strings to dates  *
foreach var of varlist index_date sotrovimab_covid_therapeutics molnupiravir_covid_therapeutics paxlovid_covid_therapeutics remdesivir_covid_therapeutics	///
        casirivimab_covid_therapeutics date_treated start_date ///
        covid_test_positive_date  primary_covid_hospital_discharge primary_covid_hospital_admission ///
	   any_covid_hospital_discharge_dat any_covid_hospital_admission_dat death_date dereg_date  ///
	   cancer_opensafely_snomed_new   immunosuppresant_drugs_nhsd ///
	   oral_steroid_drugs_nhsd  immunosupression_nhsd_new   solid_organ_transplant_nhsd_new  ///
	   covid_hosp_outcome_date0 covid_hosp_outcome_date1 covid_hosp_outcome_date2 covid_hosp_discharge_date0 covid_hosp_discharge_date1 covid_hosp_discharge_date2 ///
	   death_with_covid_date death_with_covid_underly_date covid_hosp_outcome_date0T covid_hosp_outcome_date1T covid_hosp_outcome_date2T covid_hosp_discharge_date0T ///
	   covid_hosp_discharge_date1T covid_hosp_discharge_date2T death_with_covid_dateT  death_with_covid_underly_dateT death_dateT dereg_dateT  ///
	   downs_syndrome_nhsd haematological_disease_nhsd ckd_stage_5_nhsd liver_disease_nhsd hiv_aids_nhsd  ///
	   multiple_sclerosis_nhsd motor_neurone_disease_nhsd myasthenia_gravis_nhsd huntingtons_disease_nhsd advanced_decompensated_cirrhosis decompensated_cirrhosis_icd10 ///
	   ascitic_drainage_snomed  ckd_stages_3_5 ckd_primis_stage_date ckd3_icd10 ckd4_icd10 ckd5_icd10 dialysis dialysis_icd10 dialysis_procedure kidney_transplant kidney_transplant_icd10 ///
	   kidney_transplant_procedure creatinine_ctv3_date creatinine_snomed_date creatinine_short_snomed_date eGFR_record_date eGFR_short_record_date liver_disease_nhsd_icd10 ///
	   solid_organ_transplant_snomed drugs_do_not_use drugs_consider_risk  {
  capture confirm string variable `var'
  if _rc==0 {
  rename `var' a
  gen `var' = date(a, "YMD")
  drop a
  format %td `var'
  }
}


*exclusion criteria*
keep if sotrovimab_covid_therapeutics==start_date | paxlovid_covid_therapeutics==start_date
sum age,de
keep if age>=18 & age<110
tab sex,m
keep if sex=="F"|sex=="M"
keep if has_diedT==0
tab covid_test_positive covid_positive_previous_30_days,m
*restrict study period! NOTE: end of study period should be 30 days earlier than the latest hosp event date in current extraction*
keep if date_treated>=mdy(`start_MDY')&date_treated<=date("`end_DMY'", "DMY")
sum date_treated
*drop if stp==""
replace start_date=date_treated if start_date!=date_treated
*exclude those with multiple therapy*
drop if sotrovimab_covid_therapeutics!=. & ( molnupiravir_covid_therapeutics!=.| remdesivir_covid_therapeutics!=.| casirivimab_covid_therapeutics!=.|paxlovid_covid_therapeutics!=.)
drop if paxlovid_covid_therapeutics!=. & (  molnupiravir_covid_therapeutics!=.| remdesivir_covid_therapeutics!=.| casirivimab_covid_therapeutics!=.|sotrovimab_covid_therapeutics!=. )


*define exposure*
describe
gen drug=1 if sotrovimab_covid_therapeutics==start_date
replace drug=0 if paxlovid_covid_therapeutics ==start_date
label define drug_Paxlovid 1 "sotrovimab" 0 "Paxlovid" 
label values drug drug_Paxlovid
tab drug,m


*correcting COVID hosp events: further ignore any day cases or sotro initiators who had COVID hosp record with mab procedure codes on Day 0 or 1 *
*ignore day cases and mab procedures in day 0/1*
replace covid_hosp_outcome_date0T=. if covid_hosp_outcome_date0T==covid_hosp_discharge_date0T&covid_hosp_outcome_date0T!=.
replace covid_hosp_outcome_date1T=. if covid_hosp_outcome_date1T==covid_hosp_discharge_date1T&covid_hosp_outcome_date1T!=.
*replace covid_hosp_outcome_date0=. if covid_hosp_outcome_date0==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure!=.&drug==1
*replace covid_hosp_outcome_date1=. if covid_hosp_outcome_date1==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure!=.&drug==1

gen covid_hospitalisation_outcome_da=covid_hosp_outcome_date2T
replace covid_hospitalisation_outcome_da=covid_hosp_outcome_date1T if covid_hosp_outcome_date1T!=.
replace covid_hospitalisation_outcome_da=covid_hosp_outcome_date0T if covid_hosp_outcome_date0T!=.&drug==0
replace covid_hospitalisation_outcome_da=covid_hosp_outcome_date0T if covid_hosp_outcome_date0T!=.&drug==1

gen days_to_covid_admission=covid_hospitalisation_outcome_da-start_date if covid_hospitalisation_outcome_da!=.
by drug days_to_covid_admission, sort: count if covid_hospitalisation_outcome_da!=.

*ignore and censor day cases on or after day 2 from this analysis*
*ignore and censor admissions for mab procedure >= day 2 and with same-day or 1-day discharge*
*gen covid_hosp_date_day_cases_mab=covid_hospitalisation_outcome_da if covid_hosp_outcome_date2==covid_hosp_discharge_date2&covid_hosp_outcome_date2!=.&days_to_covid_admission>=2
*replace covid_hosp_date_day_cases_mab=covid_hospitalisation_outcome_da if covid_hosp_outcome_date2==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure!=.&days_to_covid_admission>=2&(covid_hosp_discharge_date2-covid_hosp_outcome_date2)<=1&drug==1
drop if covid_hosp_outcome_date2T==covid_hosp_discharge_date2T&covid_hosp_outcome_date2T!=.&days_to_covid_admission>=2
*replace covid_hospitalisation_outcome_da=. if covid_hosp_outcome_date2==covid_hosp_date_mabs_procedure&covid_hosp_date_mabs_procedure!=.&days_to_covid_admission>=2&(covid_hosp_discharge_date2-covid_hosp_outcome_date2)<=1&drug==1
*capture and exclude COVID-hospital admission/death on the start date
by drug, sort: count if start_date==covid_hospitalisation_outcome_da| start_date==death_with_covid_dateT
drop if start_date>=covid_hospitalisation_outcome_da| start_date>=death_with_covid_dateT|start_date>=death_dateT|start_date>=dereg_dateT
drop if dereg_dateT>=start_date&dereg_dateT<=(start_date+30)


*define outcome and follow-up time*
gen start_date_30=start_date+30
*30-day COVID hosp*
gen covid_hospitalisation_30day=(covid_hospitalisation_outcome_da!=.&covid_hospitalisation_outcome_da<=start_date_30)
tab covid_hospitalisation_30day,m
tab drug covid_hospitalisation_30day,row m
*tab drug covid_hospitalisation_30day if start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022),row m
*tab drug covid_hospitalisation_30day if start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022),row m
*tab drug covid_hospitalisation_30day if start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022),row m
*30-day COVID death*
gen covid_death_30day=(death_with_covid_dateT!=.&death_with_covid_dateT<=start_date_30)
tab covid_death_30day,m
tab drug covid_death_30day,row m
*tab drug covid_death_30day if start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022),row m
*tab drug covid_death_30day if start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022),row m
*tab drug covid_death_30day if start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022),row m
*30-day all-cause death*
gen death_30day=(death_dateT!=.&death_dateT<=start_date_30)
tab death_30day,m
tab drug death_30day,row m
*tab drug death_30day if start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022),row m
*tab drug death_30day if start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022),row m
*tab drug death_30day if start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022),row m


*covariates* 
*10 high risk groups: downs_syndrome, solid_cancer, haematological_disease, renal_disease, liver_disease, imid, 
*immunosupression, hiv_aids, solid_organ_transplant, rare_neurological_conditions, high_risk_group_combined	
replace oral_steroid_drugs_nhsd=. if oral_steroid_drug_nhsd_3m_count < 2 & oral_steroid_drug_nhsd_12m_count < 4
gen imid_nhsd=min(oral_steroid_drugs_nhsd, immunosuppresant_drugs_nhsd)
gen rare_neuro_nhsd = min(multiple_sclerosis_nhsd, motor_neurone_disease_nhsd, myasthenia_gravis_nhsd, huntingtons_disease_nhsd)

*high risk group only based on codelists*
gen downs_syndrome=(downs_syndrome_nhsd<=start_date)
gen solid_cancer_new=(cancer_opensafely_snomed_new<=start_date)
gen haema_disease=( haematological_disease_nhsd <=start_date)
gen renal_disease=( ckd_stage_5_nhsd <=start_date)
gen liver_disease=( liver_disease_nhsd <=start_date)
gen imid=( imid_nhsd <=start_date)
gen immunosupression_new=( immunosupression_nhsd_new <=start_date)
gen hiv_aids=( hiv_aids_nhsd <=start_date)
gen solid_organ_new=( solid_organ_transplant_nhsd_new<=start_date)
gen rare_neuro=( rare_neuro_nhsd <=start_date)
gen high_risk_group_new=(( downs_syndrome + solid_cancer_new + haema_disease + renal_disease + liver_disease + imid + immunosupression_new + hiv_aids + solid_organ_new + rare_neuro )>0)
tab high_risk_group_new,m
keep if high_risk_group_new==1
count
local N_treated=r(N) 
count if drug==0
local N_pax=r(N) 
count if drug==1
local N_sot=r(N) 


*Time between positive test and treatment*
gen d_postest_treat=start_date - covid_test_positive_date
tab d_postest_treat,m
replace d_postest_treat=. if d_postest_treat<0|d_postest_treat>7
by drug,sort: sum d_postest_treat ,de
tab drug d_postest_treat ,row


*demo*
sum age,de
local age_treated=r(mean)
sum age if drug==0,de
local age_pax=r(mean)
sum age if drug==1,de
local age_sot=r(mean)
gen age_group3=(age>=40)+(age>=60)
label define age_group3_Paxlovid 0 "18-39" 1 "40-59" 2 ">=60" 
label values age_group3 age_group3_Paxlovid
tab age_group3,m
egen age_5y_band=cut(age), at(18,25,30,35,40,45,50,55,60,65,70,75,80,85,110) label
tab age_5y_band,m
by drug,sort: sum age,de

tab sex,m
rename sex sex_str
gen sex=0 if sex_str=="M"
replace sex=1 if sex_str=="F"
label define sex_Paxlovid 0 "Male" 1 "Female"
label values sex sex_Paxlovid
tab drug sex,row chi
sum sex
local female_treated=r(mean)
sum sex if drug==0
local female_pax=r(mean)
sum sex if drug==1
local female_sot=r(mean)

*vac and variant*
tab vaccination_status,m
rename vaccination_status vaccination_status_g6
gen vaccination_status=0 if vaccination_status_g6=="Un-vaccinated"|vaccination_status_g6=="Un-vaccinated (declined)"
replace vaccination_status=1 if vaccination_status_g6=="One vaccination"
replace vaccination_status=2 if vaccination_status_g6=="Two vaccinations"
replace vaccination_status=3 if vaccination_status_g6=="Three vaccinations"
replace vaccination_status=4 if vaccination_status_g6=="Four or more vaccinations"
label define vac_Paxlovid 0 "Un-vaccinated" 1 "One vaccination" 2 "Two vaccinations" 3 "Three vaccinations" 4 "Four or more vaccinations"
label values vaccination_status vac_Paxlovid
gen vaccination_3=1 if vaccination_status==3|vaccination_status==4
replace vaccination_3=0 if vaccination_status<3
sum vaccination_3
local vac_3_treated=r(mean)
sum vaccination_3 if drug==0
local vac_3_pax=r(mean)
sum vaccination_3 if drug==1
local vac_3_sot=r(mean)

tab variant_recorded ,m
tab symptomatic_covid_test ,m
tab drug vaccination_status ,row chi


*30-day COVID hosp*
tab covid_hospitalisation_30day,m
tab drug covid_hospitalisation_30day,row m
sum covid_hospitalisation_30day
local cov_hosp_treated=r(mean)
sum covid_hospitalisation_30day if drug==0
local cov_hosp_pax=r(mean)
sum covid_hospitalisation_30day if drug==1
local cov_hosp_sot=r(mean)

count if covid_hospitalisation_30day==1
local cov_hosp_treated_n=r(N)
count if covid_hospitalisation_30day==1&drug==0
local cov_hosp_pax_n=r(N)
count if covid_hospitalisation_30day==1&drug==1
local cov_hosp_sot_n=r(N)
*30-day COVID death*
tab covid_death_30day,m
tab drug covid_death_30day,row m
sum covid_death_30day
local cov_death_treated=r(mean)
sum covid_death_30day if drug==0
local cov_death_pax=r(mean)
sum covid_death_30day if drug==1
local cov_death_sot=r(mean)

count if covid_death_30day==1
local cov_death_treated_n=r(N)
count if covid_death_30day==1&drug==0
local cov_death_pax_n=r(N)
count if covid_death_30day==1&drug==1
local cov_death_sot_n=r(N)
*30-day all-cause death*
tab death_30day,m
tab drug death_30day,row m
sum death_30day
local death_treated=r(mean)
sum death_30day if drug==0
local death_pax=r(mean)
sum death_30day if drug==1
local death_sot=r(mean)

count if death_30day==1
local death_treated_n=r(N)
count if death_30day==1&drug==0
local death_pax_n=r(N)
count if death_30day==1&drug==1
local death_sot_n=r(N)




*exclude those with contraindications for Pax*
replace ckd_primis_stage=. if ckd_primis_stage_date>start_date
*egfr: adapted from https://github.com/opensafely/COVID-19-vaccine-breakthrough/blob/updates-feb/analysis/data_process.R*
tab creatinine_operator_ctv3,m
replace creatinine_ctv3 = . if !inrange(creatinine_ctv3, 20, 3000)| creatinine_ctv3_date>start_date
tab creatinine_operator_ctv3 if creatinine_ctv3!=.,m
replace creatinine_ctv3 = creatinine_ctv3/88.4
gen min_creatinine_ctv3=.
replace min_creatinine_ctv3 = (creatinine_ctv3/0.7)^-0.329 if sex==1
replace min_creatinine_ctv3 = (creatinine_ctv3/0.9)^-0.411 if sex==0
replace min_creatinine_ctv3 = 1 if min_creatinine_ctv3<1
gen max_creatinine_ctv3=.
replace max_creatinine_ctv3 = (creatinine_ctv3/0.7)^-1.209 if sex==1
replace max_creatinine_ctv3 = (creatinine_ctv3/0.9)^-1.209 if sex==0
replace max_creatinine_ctv3 = 1 if max_creatinine_ctv3>1
gen egfr_creatinine_ctv3 = min_creatinine_ctv3*max_creatinine_ctv3*141*(0.993^age_creatinine_ctv3) if age_creatinine_ctv3>0&age_creatinine_ctv3<=120
replace egfr_creatinine_ctv3 = egfr_creatinine_ctv3*1.018 if sex==1

tab creatinine_operator_snomed,m
tab creatinine_operator_snomed if creatinine_snomed!=.,m
replace creatinine_snomed = . if !inrange(creatinine_snomed, 20, 3000)| creatinine_snomed_date>start_date
replace creatinine_snomed_date = creatinine_short_snomed_date if missing(creatinine_snomed)
replace creatinine_operator_snomed = creatinine_operator_short_snomed if missing(creatinine_snomed)
replace age_creatinine_snomed = age_creatinine_short_snomed if missing(creatinine_snomed)
replace creatinine_snomed = creatinine_short_snomed if missing(creatinine_snomed)
replace creatinine_snomed = . if !inrange(creatinine_snomed, 20, 3000)| creatinine_snomed_date>start_date
replace creatinine_snomed = creatinine_snomed/88.4
gen min_creatinine_snomed=.
replace min_creatinine_snomed = (creatinine_snomed/0.7)^-0.329 if sex==1
replace min_creatinine_snomed = (creatinine_snomed/0.9)^-0.411 if sex==0
replace min_creatinine_snomed = 1 if min_creatinine_snomed<1
gen max_creatinine_snomed=.
replace max_creatinine_snomed = (creatinine_snomed/0.7)^-1.209 if sex==1
replace max_creatinine_snomed = (creatinine_snomed/0.9)^-1.209 if sex==0
replace max_creatinine_snomed = 1 if max_creatinine_snomed>1
gen egfr_creatinine_snomed = min_creatinine_snomed*max_creatinine_snomed*141*(0.993^age_creatinine_snomed) if age_creatinine_snomed>0&age_creatinine_snomed<=120
replace egfr_creatinine_snomed = egfr_creatinine_snomed*1.018 if sex==1

*drug interactions*
tab drug if drugs_do_not_use<=start_date
tab drug if drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-3*365.25)
tab drug if drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-365.25)
tab drug if drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-180)
tab drug if drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-90)
tab drug if drugs_consider_risk<=start_date
tab drug if drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-3*365.25)
tab drug if drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-365.25)
tab drug if drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-180)
tab drug if drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-90)
gen drugs_do_not_use_contra=(drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-180))
gen drugs_consider_risk_contra=(drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-180))


drop if solid_organ_new==1|solid_organ_therapeutics==1|solid_organ_transplant_snomed<=start_date
drop if advanced_decompensated_cirrhosis<=start_date|decompensated_cirrhosis_icd10<=start_date|ascitic_drainage_snomed<=start_date|liver_disease_nhsd_icd10<=start_date
drop if renal_disease==1|renal_therapeutics==1|ckd_stages_3_5<=start_date|ckd_primis_stage==3|ckd_primis_stage==4|ckd_primis_stage==5|ckd3_icd10<=start_date|ckd4_icd10<=start_date|ckd5_icd10<=start_date
drop if kidney_transplant<=start_date|kidney_transplant_icd10<=start_date|kidney_transplant_procedure<=start_date
drop if dialysis<=start_date|dialysis_icd10<=start_date|dialysis_procedure<=start_date
drop if (egfr_creatinine_ctv3<60&creatinine_operator_ctv3!="<")|(egfr_creatinine_snomed<60&creatinine_operator_snomed!="<")|(eGFR_record<60&eGFR_record>0&eGFR_operator!=">"&eGFR_operator!=">=")|(eGFR_short_record<60&eGFR_short_record>0&eGFR_short_operator!=">"&eGFR_short_operator!=">=")
*drop if drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-365.25)
*drop if drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-365.25)
drop if drugs_do_not_use<=start_date&drugs_do_not_use>=(start_date-180)
*drop if drugs_consider_risk<=start_date&drugs_consider_risk>=(start_date-180)
count if drug==1
local N_sot_no=r(N) 
sum age if drug==1,de
local age_sot_no=r(mean)
sum sex if drug==1
local female_sot_no=r(mean)
sum vaccination_3 if drug==1
local vac_3_sot_no=r(mean)


*30-day COVID hosp*
tab covid_hospitalisation_30day,m
tab drug covid_hospitalisation_30day,row m
sum covid_hospitalisation_30day if drug==1
local cov_hosp_sot_no=r(mean)
count if covid_hospitalisation_30day==1&drug==1
local cov_hosp_sot_no_n=r(N)
*tab drug covid_hospitalisation_30day if start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022),row m
*tab drug covid_hospitalisation_30day if start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022),row m
*tab drug covid_hospitalisation_30day if start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022),row m
*30-day COVID death*
tab covid_death_30day,m
tab drug covid_death_30day,row m
sum covid_death_30day if drug==1
local cov_death_sot_no=r(mean)
count if covid_death_30day==1&drug==1
local cov_death_sot_no_n=r(N)
*tab drug covid_death_30day if start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022),row m
*tab drug covid_death_30day if start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022),row m
*tab drug covid_death_30day if start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022),row m
*30-day all-cause death*
tab death_30day,m
tab drug death_30day,row m
sum death_30day if drug==1
local death_sot_no=r(mean)
count if death_30day==1&drug==1
local death_sot_no_n=r(N)
*tab drug death_30day if start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022),row m
*tab drug death_30day if start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022),row m
*tab drug death_30day if start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022),row m



*add rows to Table*
post mytab  ("Overall treated (Sot/Pax)") ("`N_treated'") ("`cov_hosp_treated_n'/`cov_hosp_treated'") ("`cov_death_treated_n'/`cov_death_treated'") ("`death_treated_n'/`death_treated'") ("`age_treated'") ("`female_treated'") ("`vac_3_treated'") 
post mytab ("Paxlovid") ("`N_pax'") ("`cov_hosp_pax_n'/`cov_hosp_pax'") ("`cov_death_pax_n'/`cov_death_pax'") ("`death_pax_n'/`death_pax'") ("`age_pax'") ("`female_pax'") ("`vac_3_pax'")  
post mytab ("Sotrovimab") ("`N_sot'") ("`cov_hosp_sot_n'/`cov_hosp_sot'") ("`cov_death_sot_n'/`cov_death_sot'") ("`death_sot_n'/`death_sot'") ("`age_sot'") ("`female_sot'") ("`vac_3_sot'")  
post mytab ("Sotro without contra") ("`N_sot_no'") ("`cov_hosp_sot_no_n'/`cov_hosp_sot_no'") ("`cov_death_sot_no_n'/`cov_death_sot_no'") ("`death_sot_no_n'/`death_sot_no'") ("`age_sot_no'") ("`female_sot_no'") ("`vac_3_sot_no'")  
post mytab ("Untreated but eligible") ("`N_untreated'") ("`cov_hosp_untreated_n'/`cov_hosp_untreated'") ("`cov_death_untreated_n'/`cov_death_untreated'") ("`death_untreated_n'/`death_untreated'") ("`age_untreated'") ("`female_untreated'") ("`vac_3_untreated'")  
post mytab ("Untreated without contra") ("`N_untreated_no'") ("`cov_hosp_untreated_no_n'/`cov_hosp_untreated_no'") ("`cov_death_untreated_no_n'/`cov_death_untreated_no'") ("`death_untreated_no_n'/`death_untreated_no'") ("`age_untreated_no'") ("`female_untreated_no'") ("`vac_3_untreated_no'")  
post mytab ("Note:") ("") ("") ("") ("") ("") ("") ("")  
post mytab ("Start date: `start_DMY'") ("") ("") ("") ("") ("") ("") ("")  
post mytab ("End date: `end_DMY'") ("") ("") ("") ("") ("") ("") ("")   


postclose mytab
clear
use ./output/table.dta
export delimited using "./output/table.csv", replace
clear

log close




