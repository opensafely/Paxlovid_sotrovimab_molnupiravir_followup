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
log using ./logs/data_preparation_update_raw_rate_untreated, replace t
clear

* import dataset
import delimited ./output/input_raw_rate.csv, delimiter(comma) varnames(1) case(preserve) 
describe
keep if registered_eligible==1
drop if cancer_opensafely_snomed_new==""&immunosuppresant_drugs_nhsd==""&oral_steroid_drugs_nhsd==""&immunosupression_nhsd_new==""&solid_organ_transplant_nhsd_new==""&downs_syndrome_nhsd==""&haematological_disease_nhsd==""&ckd_stage_5_nhsd==""&liver_disease_nhsd==""&hiv_aids_nhsd==""&multiple_sclerosis_nhsd==""&motor_neurone_disease_nhsd==""&myasthenia_gravis_nhsd==""&huntingtons_disease_nhsd=="" 



*  Convert strings to dates  *
foreach var of varlist sotrovimab_covid_therapeutics molnupiravir_covid_therapeutics paxlovid_covid_therapeutics remdesivir_covid_therapeutics	///
        casirivimab_covid_therapeutics   ///
        covid_test_positive_date  primary_covid_hospital_discharge primary_covid_hospital_admission ///
	   any_covid_hospital_admission_dat death_date dereg_date  ///
	   cancer_opensafely_snomed_new   immunosuppresant_drugs_nhsd ///
	   oral_steroid_drugs_nhsd  immunosupression_nhsd_new   solid_organ_transplant_nhsd_new  ///
	   covid_hosp_outcome_date0 covid_hosp_outcome_date1 covid_hosp_outcome_date2 covid_hosp_discharge_date0 covid_hosp_discharge_date1 covid_hosp_discharge_date2 ///
	   death_with_covid_on_the_death_ce death_with_covid_underlying_date    ///
	   downs_syndrome_nhsd haematological_disease_nhsd ckd_stage_5_nhsd liver_disease_nhsd hiv_aids_nhsd  ///
	   multiple_sclerosis_nhsd motor_neurone_disease_nhsd myasthenia_gravis_nhsd huntingtons_disease_nhsd advanced_decompensated_cirrhosis decompensated_cirrhosis_icd10 ///
	   ascitic_drainage_snomed  ckd_stages_3_5 ckd_primis_stage_date ckd3_icd10 ckd4_icd10 ckd5_icd10 dialysis dialysis_icd10 dialysis_procedure kidney_transplant kidney_transplant_icd10 ///
	   kidney_transplant_procedure creatinine_ctv3_date creatinine_snomed_date creatinine_short_snomed_date eGFR_record_date eGFR_short_record_date ///
	   solid_organ_transplant_snomed drugs_do_not_use drugs_consider_risk  {
  capture confirm string variable `var'
  if _rc==0 {
  rename `var' a
  gen `var' = date(a, "YMD")
  drop a
  format %td `var'
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

*check hosp/death event date range*
codebook covid_hosp_outcome_date2  death_date

*exclusion criteria*
drop if sotrovimab_covid_therapeutics!=. | paxlovid_covid_therapeutics!=. | molnupiravir_covid_therapeutics!=. |remdesivir_covid_therapeutics!=. |casirivimab_covid_therapeutics!=.
sum age,de
keep if age>=18 & age<110
tab sex,m
keep if sex=="F"|sex=="M"
keep if has_died==0
tab covid_test_positive covid_positive_previous_30_days,m
keep if covid_test_positive==1 & covid_positive_previous_30_days==0
drop if primary_covid_hospital_discharge!=.|primary_covid_hospital_admission!=.
drop if any_covid_hospital_admission_dat!=.
*restrict start_date to 2022Feb10 to now*
*loose this restriction to increase N?*
keep if covid_test_positive_date>=mdy(12,16,2021)&covid_test_positive_date<=mdy(10,01,2022)
*drop if stp==""

gen start_date=covid_test_positive_date

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
*capture and exclude COVID-hospital admission/death on the start date
drop if covid_test_positive_date>=covid_hospitalisation_outcome_da| covid_test_positive_date>=death_with_covid_on_the_death_ce|covid_test_positive_date>=death_date|covid_test_positive_date>=dereg_date




*covariates* 
*10 high risk groups: downs_syndrome, solid_cancer, haematological_disease, renal_disease, liver_disease, imid, 
*immunosupression, hiv_aids, solid_organ_transplant, rare_neurological_conditions, high_risk_group_combined	

replace oral_steroid_drugs_nhsd=. if oral_steroid_drug_nhsd_3m_count < 2 & oral_steroid_drug_nhsd_12m_count < 4
gen imid_nhsd=min(oral_steroid_drugs_nhsd, immunosuppresant_drugs_nhsd)
gen rare_neuro_nhsd = min(multiple_sclerosis_nhsd, motor_neurone_disease_nhsd, myasthenia_gravis_nhsd, huntingtons_disease_nhsd)

*gen downs_syndrome=(downs_syndrome_nhsd<=start_date|downs_therapeutics==1)
*gen solid_cancer=(cancer_opensafely_snomed<=start_date|solid_cancer_therapeutics==1)
*gen solid_cancer_new=(cancer_opensafely_snomed_new<=start_date|solid_cancer_therapeutics==1)
*gen haema_disease=( haematological_disease_nhsd <=start_date|haema_disease_therapeutics==1)
*gen renal_disease=( ckd_stage_5_nhsd <=start_date|renal_therapeutics==1)
*gen liver_disease=( liver_disease_nhsd <=start_date|liver_therapeutics==1)
*gen imid=( imid_nhsd <=start_date|imid_therapeutics==1)
*gen immunosupression=( immunosupression_nhsd <=start_date|immunosup_therapeutics==1)
*gen immunosupression_new=( immunosupression_nhsd_new <=start_date|immunosup_therapeutics==1)
*gen hiv_aids=( hiv_aids_nhsd <=start_date|hiv_aids_therapeutics==1)
*gen solid_organ=( solid_organ_transplant_nhsd<=start_date|solid_organ_therapeutics==1)
*gen solid_organ_new=( solid_organ_transplant_nhsd_new<=start_date|solid_organ_therapeutics==1)
*gen rare_neuro=( rare_neuro_nhsd <=start_date|rare_neuro_therapeutics==1)
*gen high_risk_group=(( downs_syndrome + solid_cancer + haema_disease + renal_disease + liver_disease + imid + immunosupression + hiv_aids + solid_organ + rare_neuro )>0)
*tab high_risk_group,m
*gen high_risk_group_new=(( downs_syndrome + solid_cancer_new + haema_disease + renal_disease + liver_disease + imid + immunosupression_new + hiv_aids + solid_organ_new + rare_neuro )>0)
*tab high_risk_group_new,m
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

*demo*
gen age_group3=(age>=40)+(age>=60)
label define age_group3_Paxlovid 0 "18-39" 1 "40-59" 2 ">=60" 
label values age_group3 age_group3_Paxlovid
tab age_group3,m
egen age_5y_band=cut(age), at(18,25,30,35,40,45,50,55,60,65,70,75,80,85,110) label
tab age_5y_band,m
gen age_50=(age>=50)
gen age_55=(age>=55)
gen age_60=(age>=60)

tab sex,m
rename sex sex_str
gen sex=0 if sex_str=="M"
replace sex=1 if sex_str=="F"
label define sex_Paxlovid 0 "Male" 1 "Female"
label values sex sex_Paxlovid


*define outcome and follow-up time*
gen study_end_date=mdy(12,22,2022)
gen start_date_30=covid_test_positive_date+30
*30-day COVID hosp*
gen covid_hospitalisation_30day=(covid_hospitalisation_outcome_da!=.&covid_hospitalisation_outcome_da<=start_date_30)
tab covid_hospitalisation_30day,m
tab  covid_hospitalisation_30day if start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022), m
tab  covid_hospitalisation_30day if start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022), m
tab  covid_hospitalisation_30day if start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022),m
*30-day COVID death*
gen covid_death_30day=(death_with_covid_on_the_death_ce!=.&death_with_covid_on_the_death_ce<=start_date_30)
tab covid_death_30day,m
tab  covid_death_30day if start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022),m
tab  covid_death_30day if start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022), m
tab  covid_death_30day if start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022), m
*30-day all-cause death*
gen death_30day=(death_date!=.&death_date<=start_date_30)
tab death_30day,m
tab  death_30day if start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022), m
tab  death_30day if start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022), m
tab  death_30day if start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022), m




*UKRR*
tab ukrr_2021,m
*30-day COVID hosp*
tab covid_hospitalisation_30day if ukrr_2021==1,m
tab  covid_hospitalisation_30day if ukrr_2021==1&start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022), m
tab  covid_hospitalisation_30day if ukrr_2021==1&start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022), m
tab  covid_hospitalisation_30day if ukrr_2021==1&start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022),m
*30-day COVID death*
tab covid_death_30day if ukrr_2021==1,m
tab  covid_death_30day if ukrr_2021==1&start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022),m
tab  covid_death_30day if ukrr_2021==1&start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022), m
tab  covid_death_30day if ukrr_2021==1&start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022), m
*30-day all-cause death*
tab death_30day if ukrr_2021==1,m
tab  death_30day if ukrr_2021==1&start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022), m
tab  death_30day if ukrr_2021==1&start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022), m
tab  death_30day if ukrr_2021==1&start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022), m






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


*30-day COVID hosp*
tab covid_hospitalisation_30day,m
tab  covid_hospitalisation_30day if start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022), m
tab  covid_hospitalisation_30day if start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022), m
tab  covid_hospitalisation_30day if start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022),m
*30-day COVID death*
tab covid_death_30day,m
tab  covid_death_30day if start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022),m
tab  covid_death_30day if start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022), m
tab  covid_death_30day if start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022), m
*30-day all-cause death*
tab death_30day,m
tab  death_30day if start_date>=mdy(12,16,2021)&start_date<=mdy(2,10,2022), m
tab  death_30day if start_date>=mdy(2,11,2022)&start_date<=mdy(5,31,2022), m
tab  death_30day if start_date>=mdy(6,1,2022)&start_date<=mdy(10,1,2022), m







log close




