##copy_study_defintion## Adapted codes from https://github.com/opensafely/antibody-and-antiviral-deployment
## Import code building blocks from cohort extractor package 
from cohortextractor import (
  StudyDefinition,
  patients,
  codelist_from_csv,
  codelist,
  filter_codes_by_category,
  combine_codelists,
  Measure
)

## Import codelists from codelist.py (which pulls them from the codelist folder)
from codelists import *

# DEFINE STUDY POPULATION ----

## Define study time variables
from datetime import timedelta, date, datetime 

campaign_start = "2021-12-16"
end_date = date.today().isoformat()

## Define study population and variables
study = StudyDefinition(

  ## Configure the expectations framework
  default_expectations = {
    "date": {"earliest": "2021-11-01", "latest": "today"},
    "rate": "uniform",
    "incidence": 0.05,
  },
  
  ## Define index date
  index_date = campaign_start,
  
  # POPULATION ----
  population = patients.satisfying(
    """
    age >= 18 AND age < 110
    AND NOT has_died
    AND registered_treated 
    AND (sotrovimab_covid_therapeutics OR paxlovid_covid_therapeutics OR molnupiravir_covid_therapeutics)
    """,
  ),
  #require covid_test_positive_date<=date_treated (sensitivity analysis)
  #loose "AND (covid_test_positive AND NOT covid_positive_previous_30_days)"
  #AND NOT pregnancy (sensitivity analysis)
  #AND NOT (casirivimab_covid_therapeutics OR remdesivir_covid_therapeutics) (sensitivity analysis)

  # TREATMENT - NEUTRALISING MONOCLONAL ANTIBODIES OR ANTIVIRALS ----
  
  ## Sotrovimab
  sotrovimab_covid_therapeutics = patients.with_covid_therapeutics(
    #with_these_statuses = ["Approved", "Treatment Complete"],
    with_these_therapeutics = "Sotrovimab",
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.4
    },
  ),
  # restrict by status
  sotrovimab_covid_approved = patients.with_covid_therapeutics(
    with_these_statuses = ["Approved"],
    with_these_therapeutics = "Sotrovimab",
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.4
    },
  ),
  sotrovimab_covid_complete = patients.with_covid_therapeutics(
    with_these_statuses = ["Treatment Complete"],
    with_these_therapeutics = "Sotrovimab",
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.4
    },
  ),
  sotrovimab_covid_not_start = patients.with_covid_therapeutics(
    with_these_statuses = ["Treatment Not Started"],
    with_these_therapeutics = "Sotrovimab",
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.4
    },
  ),
  sotrovimab_covid_stopped = patients.with_covid_therapeutics(
    with_these_statuses = ["Treatment Stopped"],
    with_these_therapeutics = "Sotrovimab",
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.4
    },
  ),
 
# ###.-start
# ### Molnupiravir
  molnupiravir_covid_therapeutics = patients.with_covid_therapeutics(
    #with_these_statuses = ["Approved", "Treatment Complete"],
    with_these_therapeutics = "Molnupiravir",
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.4
    },
  ),
#   # restrict by status
  molnupiravir_covid_approved = patients.with_covid_therapeutics(
    with_these_statuses = ["Approved"],
    with_these_therapeutics = "Molnupiravir",
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.4
    },
  ),
  molnupiravir_covid_complete = patients.with_covid_therapeutics(
    with_these_statuses = ["Treatment Complete"],
    with_these_therapeutics = "Molnupiravir",
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.4
    },
  ),
  molnupiravir_covid_not_start = patients.with_covid_therapeutics(
    with_these_statuses = ["Treatment Not Started"],
    with_these_therapeutics = "Molnupiravir",
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.4
    },
  ),
  molnupiravir_covid_stopped = patients.with_covid_therapeutics(
    with_these_statuses = ["Treatment Stopped"],
    with_these_therapeutics = "Molnupiravir",
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.4
    },
  ),  
####.-end

#   ### Paxlovid
  paxlovid_covid_therapeutics = patients.with_covid_therapeutics(
    #with_these_statuses = ["Approved", "Treatment Complete"],
    with_these_therapeutics = "Paxlovid",
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.4
    },
  ),
#   # restrict by status
  paxlovid_covid_approved = patients.with_covid_therapeutics(
    with_these_statuses = ["Approved"],
    with_these_therapeutics = "Paxlovid",
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.4
    },
  ),
  paxlovid_covid_complete = patients.with_covid_therapeutics(
    with_these_statuses = ["Treatment Complete"],
    with_these_therapeutics = "Paxlovid",
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.4
    },
  ),
  paxlovid_covid_not_start = patients.with_covid_therapeutics(
    with_these_statuses = ["Treatment Not Started"],
    with_these_therapeutics = "Paxlovid",
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.4
    },
  ),
  paxlovid_covid_stopped = patients.with_covid_therapeutics(
    with_these_statuses = ["Treatment Stopped"],
    with_these_therapeutics = "Paxlovid",
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.4
    },
  ),

#   # ### Molnupiravir  ##.-commented out
#   # molnupiravir_covid_therapeutics = patients.with_covid_therapeutics(
#   #   #with_these_statuses = ["Approved", "Treatment Complete"],
#   #   with_these_therapeutics = "Molnupiravir",
#   #   with_these_indications = "non_hospitalised",
#   #   on_or_after = "index_date",
#   #   find_first_match_in_period = True,
#   #   returning = "date",
#   #   date_format = "YYYY-MM-DD",
#   #   return_expectations = {
#   #     "date": {"earliest": "2022-02-10"},
#   #     "incidence": 0.05
#   #   },
#   # ), 

#   # ## Remdesivir  ##.-commented out
#   # remdesivir_covid_therapeutics = patients.with_covid_therapeutics(
#   #   #with_these_statuses = ["Approved", "Treatment Complete"],
#   #   with_these_therapeutics = "Remdesivir",
#   #   with_these_indications = "non_hospitalised",
#   #   on_or_after = "index_date",
#   #   find_first_match_in_period = True,
#   #   returning = "date",
#   #   date_format = "YYYY-MM-DD",
#   #   return_expectations = {
#   #     "date": {"earliest": "2021-12-16"},
#   #     "incidence": 0.05
#   #   },
#   # ),
  
#   # ### Casirivimab and imdevimab ##.-commented out
#   # casirivimab_covid_therapeutics = patients.with_covid_therapeutics(
#   #   #with_these_statuses = ["Approved", "Treatment Complete"],
#   #   with_these_therapeutics = "Casirivimab and imdevimab",
#   #   with_these_indications = "non_hospitalised",
#   #   on_or_after = "index_date",
#   #   find_first_match_in_period = True,
#   #   returning = "date",
#   #   date_format = "YYYY-MM-DD",
#   #   return_expectations = {
#   #     "date": {"earliest": "2021-12-16"},
#   #     "incidence": 0.05
#   #   },
#   # ), 
  
  
#   ## Date treated
  date_treated = patients.minimum_of(
    "sotrovimab_covid_therapeutics",
    "paxlovid_covid_therapeutics",
    "molnupiravir_covid_therapeutics",
  ),
  
  registered_treated = patients.registered_as_of("date_treated"), 

#   # OVERALL ELIGIBILITY CRITERIA VARIABLES ----
  
#   ## Inclusion criteria variables
  
#   ### First positive SARS-CoV-2 test
#   # Note patients are eligible for treatment if diagnosed <=5d ago
  covid_test_positive = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    returning = "binary_flag",
    on_or_after = "date_treated - 30 days",
    find_first_match_in_period = True,
    restrict_to_earliest_specimen_date = False,
    return_expectations = {
      "incidence": 0.9
    },
  ),
  
  covid_test_positive_date = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    find_first_match_in_period = True,
    restrict_to_earliest_specimen_date = False,
    returning = "date",
    date_format = "YYYY-MM-DD",
    on_or_after = "date_treated - 30 days",
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date"},
      "incidence": 0.9
    },
  ),
  
#   ### Second positive SARS-CoV-2 test
  covid_test_positive_date2 = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    find_first_match_in_period = True,
    restrict_to_earliest_specimen_date = False,
    returning = "date",
    date_format = "YYYY-MM-DD",
    on_or_after = "covid_test_positive_date + 30 days",
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date"},
      "incidence": 0.1
    },
  ),
  
#   ### Covid test type - add in when avaliable 
  covid_positive_test_type = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    returning = "case_category",
    on_or_after = "date_treated - 5 days",
    restrict_to_earliest_specimen_date = True,
    return_expectations = {
      "category": {"ratios": {"LFT_Only": 0.4, "PCR_Only": 0.4, "LFT_WithPCR": 0.2}},
      "incidence": 0.2,
    },
  ),
  
#   ### Positive covid test last 30 days 
#   # (note this will only apply to patients who first tested positive towards the beginning
#   # of the study period)
  covid_positive_previous_30_days = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    returning = "binary_flag",
    between = ["covid_test_positive_date - 31 days", "covid_test_positive_date - 1 day"],
    find_last_match_in_period = True,
    restrict_to_earliest_specimen_date = False,
    return_expectations = {
      "incidence": 0.05
    },
  ),
  
#   # positive test history
  covid_test_positive_pre_date = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    find_last_match_in_period = True,
    restrict_to_earliest_specimen_date = False,
    returning = "date",
    date_format = "YYYY-MM-DD",
    on_or_before = "covid_test_positive_date - 30 days",
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date"},
      "incidence": 0.1
    },
  ),
  
#   ### Onset of symptoms of COVID-19
  symptomatic_covid_test = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "any",
    returning = "symptomatic",
    on_or_after = "date_treated - 5 days",
    find_first_match_in_period = True,
    restrict_to_earliest_specimen_date = False,
    return_expectations={
      "incidence": 0.1,
      "category": {
        "ratios": {
          "": 0.2,
          "N": 0.2,
          "Y": 0.6,
        }
      },
    },
  ),
  
  covid_symptoms_snomed = patients.with_these_clinical_events(
    covid_symptoms_snomed_codes,
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_first_match_in_period = True,
    on_or_after = "date_treated - 5 days",
  ),
  
  
#   ## Study start date for extracting variables
  start_date = patients.minimum_of(
    "sotrovimab_covid_therapeutics",
    "paxlovid_covid_therapeutics",
    "molnupiravir_covid_therapeutics" #, # .
  ),
  
#   ## Exclusion criteria variables
  
#   ### Require hospitalisation for COVID-19
#   ## NB this data lags behind the therapeutics/testing data so may be missing
  primary_covid_hospital_discharge_date = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    #with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_before = "start_date - 1 day",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date - 1 day"},
      "rate": "uniform",
      "incidence": 0.05
    },
  ),
  primary_covid_hospital_admission_date = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    #with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_before = "start_date - 1 day",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date - 1 day"},
      "rate": "uniform",
      "incidence": 0.05
    },
  ),    
  any_covid_hospital_discharge_date = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_these_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_before = "start_date - 1 day",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date - 1 day"},
      "rate": "uniform",
      "incidence": 0.05
    },
  ),
  
#   #line 490-1712 -.
#   ### New supplemental oxygen requirement specifically for the management of COVID-19 symptoms
#   #   (not currently possible to define/code)
  
#   ### Known hypersensitivity reaction to the active substances or to any of the excipients of treatments
#   #   (not currently possible to define/code)
  
  
  
#   # TREATMENT SPECIFIC ELIGIBILITY CRITERIA VARIABLES ----
  
#   ## Paxlovid - inclusion
  
#   ### Clinical judgement deems that an antiviral is the preferred option
#   #   (not currently possible to define/code)
  
#   ### Treatment is commenced within 5 days of symptom onset
#   #   (defined above)
  
#   ### The patient does NOT have a history of advanced decompensated liver cirrhosis or stage 3-5 chronic kidney disease
#   #   (OR also use renal and liver high rosk cohorts defined below)
#   #  advanced decompensated liver cirrhosis
  advanced_decompensated_cirrhosis = patients.with_these_clinical_events(
    advanced_decompensated_cirrhosis_snomed_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
    return_expectations={  ##added
            "incidence": 0.2, ##added
            "date": {"earliest": "1900-01-01", "latest": "today"},##added.
        },##added
  ),
  decompensated_cirrhosis_icd10 = patients.admitted_to_hospital(
        returning="date_admitted",
        find_last_match_in_period=True,
        date_format="YYYY-MM-DD",
        with_these_diagnoses=advanced_decompensated_cirrhosis_icd10_codes,
        on_or_before="start_date",
        return_expectations={
            "rate": "uniform",
            "date": {"earliest": "1900-01-01", "latest": "today"},
        },
  ),
#   #  regular ascitic drainage (opcs4_codes in hospital??)
  ascitic_drainage_snomed = patients.with_these_clinical_events(  
    ascitic_drainage_snomed_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ), 
  ascitic_drainage_snomed_pre = patients.with_these_clinical_events(  
    ascitic_drainage_snomed_codes,
    on_or_before = "ascitic_drainage_snomed - 1 day",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),    

#   #  hospital admission with liver disease (see below?)

#   ## CKD DEFINITIONS - adapted from https://github.com/opensafely/risk-factors-research
#   #  recorded 3-5 CKD
  ckd_stages_3_5 = patients.with_these_clinical_events(
    chronic_kidney_disease_stages_3_5_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),
  
  ckd_primis_stage=patients.with_these_clinical_events(
        codelist=primis_ckd_stage,
        on_or_before = "start_date",
        returning="category",   
        find_last_match_in_period = True,
        include_date_of_match=True,
        # date_format = "YYYY-MM-DD", #commented out-q
        return_expectations={
        #     "rate": "universal",
            "category": {"ratios": {"1": 0.5, "2": 0.5}},
            "incidence": 0.2,   ##added  -q
         },
  ),
  ckd3_icd10 = patients.admitted_to_hospital(
        returning="date_admitted",
        find_last_match_in_period=True,
        date_format="YYYY-MM-DD",
        with_these_diagnoses=codelist(["N183"], system="icd10"),
        on_or_before="start_date",
        return_expectations={
            "rate": "uniform",
            "date": {"earliest": "1900-01-01", "latest": "today"},
        },
  ),
  ckd4_icd10 = patients.admitted_to_hospital(
        returning="date_admitted",
        find_last_match_in_period=True,
        date_format="YYYY-MM-DD",
        with_these_diagnoses=codelist(["N184"], system="icd10"),
        on_or_before="start_date",
        return_expectations={
            "rate": "uniform",
            "date": {"earliest": "1900-01-01", "latest": "today"},
        },
  ),
  ckd5_icd10 = patients.admitted_to_hospital(
        returning="date_admitted",
        find_last_match_in_period=True,
        date_format="YYYY-MM-DD",
        with_these_diagnoses=codelist(["N185"], system="icd10"),
        on_or_before="start_date",
        return_expectations={
            "rate": "uniform",
            "date": {"earliest": "1900-01-01", "latest": "today"},
        },
  ),    

#   #  recorded dialysis
  dialysis = patients.with_these_clinical_events(
    dialysis_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),
  dialysis_icd10 = patients.admitted_to_hospital(
        returning="date_admitted",
        find_last_match_in_period=True,
        date_format="YYYY-MM-DD",
        with_these_diagnoses=dialysis_icd10_codelist,
        on_or_before="start_date",
        return_expectations={
            "rate": "uniform",
            "date": {"earliest": "1900-01-01", "latest": "today"},
        },
  ),
  dialysis_procedure = patients.admitted_to_hospital(
        returning="date_admitted",
        find_last_match_in_period=True,
        date_format="YYYY-MM-DD",
        with_these_procedures=dialysis_opcs4_codelist,
        on_or_before="start_date",
        return_expectations={
            "rate": "uniform",
            "date": {"earliest": "1900-01-01", "latest": "today"},
        },
  ),  

#   #  kidney transplant
  kidney_transplant = patients.with_these_clinical_events(
    kidney_transplant_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),
  kidney_transplant_icd10 = patients.admitted_to_hospital(
        returning="date_admitted",
        find_last_match_in_period=True,
        with_these_diagnoses=kidney_tx_icd10_codelist,
        date_format="YYYY-MM-DD",
        on_or_before="start_date",
        return_expectations={
            "rate": "uniform",
            "date": {"earliest": "1900-01-01", "latest": "today"},
        },
  ),
  kidney_transplant_procedure = patients.admitted_to_hospital(
        returning="date_admitted",
        find_last_match_in_period=True,
        date_format="YYYY-MM-DD",
        with_these_procedures=kidney_tx_opcs4_codelist,
        on_or_before="start_date",
        return_expectations={
            "rate": "uniform",
            "date": {"earliest": "1900-01-01", "latest": "today"},
        },
  ),

#   # RRT??
  RRT=patients.with_these_clinical_events(
    codelist=RRT_codelist,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
    return_expectations={
            "incidence": 0.2,
            "date": {"earliest": "1900-01-01", "latest": "today"},
    },
  ),
  RRT_icd10 = patients.admitted_to_hospital(
        returning="date_admitted",
        find_last_match_in_period=True,
        date_format="YYYY-MM-DD",
        with_these_diagnoses=RRT_icd10_codelist,
        on_or_before="start_date",
        return_expectations={
            "rate": "uniform",
            "date": {"earliest": "1900-01-01", "latest": "today"},
        },
  ),
  RRT_procedure = patients.admitted_to_hospital(
        returning="date_admitted",
        find_last_match_in_period=True,
        date_format="YYYY-MM-DD",
        with_these_procedures=RRT_opcs4_codelist,
        on_or_before="start_date",
        return_expectations={
            "rate": "uniform",
            "date": {"earliest": "1900-01-01", "latest": "today"},
        },
  ),  

#   #  3-5 CKD based on recorded creatinine value
  creatinine_ctv3 = patients.with_these_clinical_events(
    creatinine_codes_ctv3,
    find_last_match_in_period=True,
    on_or_before = "start_date",
    returning="numeric_value",
    include_date_of_match=True,
    date_format = "YYYY-MM-DD",
    return_expectations={
        "float": {"distribution": "normal", "mean": 60.0, "stddev": 15},
        "incidence": 0.95,
    },
  ),
  creatinine_operator_ctv3 = patients.comparator_from(
    "creatinine_ctv3",
    return_expectations={
       "rate": "universal",
       "category": {
         "ratios": {  # ~, =, >=, >, <, <=
            None: 0.10,
            "~": 0.05,
            "=": 0.65,
            ">=": 0.05,
            ">": 0.05,
            "<": 0.05,
            "<=": 0.05,
         },
       },
       "incidence": 0.80,
    },
  ),
  age_creatinine_ctv3 = patients.age_as_of(
    "creatinine_ctv3_date",
    return_expectations = {
      "rate": "universal",
      "int": {"distribution": "population_ages"},
    },
  ),
  creatinine_snomed=patients.with_these_clinical_events(
    codelist=creatinine_codes_snomed,
    find_last_match_in_period=True,
    on_or_before = "start_date",
    returning="numeric_value",
    include_date_of_match=True,
    date_format = "YYYY-MM-DD",
    return_expectations={
            "float": {"distribution": "normal", "mean": 45.0, "stddev": 20},
            "incidence": 0.5,
    },
  ),
  creatinine_operator_snomed=patients.comparator_from(
        "creatinine_snomed",
        return_expectations={
            "rate": "universal",
            "category": {
                "ratios": {  # ~, =, >= , > , < , <=
                    None: 0.10,
                    "~": 0.05,
                    "=": 0.65,
                    ">=": 0.05,
                    ">": 0.05,
                    "<": 0.05,
                    "<=": 0.05,
                }
            },
            "incidence": 0.80,
        },
  ),  
  age_creatinine_snomed = patients.age_as_of(
    "creatinine_snomed_date",
    return_expectations = {
      "rate": "universal",
      "int": {"distribution": "population_ages"},
    },
  ),  
  creatinine_short_snomed=patients.with_these_clinical_events(
    codelist=creatinine_codes_short_snomed,
    find_last_match_in_period=True,
    on_or_before = "start_date",
    returning="numeric_value",
    include_date_of_match=True,
    date_format = "YYYY-MM-DD",
    return_expectations={
            "float": {"distribution": "normal", "mean": 45.0, "stddev": 20},
            "incidence": 0.5,
    },
  ),
  creatinine_operator_short_snomed=patients.comparator_from(
        "creatinine_short_snomed",
        return_expectations={
            "rate": "universal",
            "category": {
                "ratios": {  # ~, =, >= , > , < , <=
                    None: 0.10,
                    "~": 0.05,
                    "=": 0.65,
                    ">=": 0.05,
                    ">": 0.05,
                    "<": 0.05,
                    "<=": 0.05,
                }
            },
            "incidence": 0.80,
        },
  ),  
  age_creatinine_short_snomed = patients.age_as_of(
    "creatinine_short_snomed_date",
    return_expectations = {
      "rate": "universal",
      "int": {"distribution": "population_ages"},
    },
  ),    
  #  3-5 CKD based on recorded eGFR value
  eGFR_record=patients.with_these_clinical_events(
    codelist=eGFR_level_codelist,
    find_last_match_in_period=True,
    on_or_before = "start_date",
    returning="numeric_value",
    include_date_of_match=True,
    date_format = "YYYY-MM-DD",
    return_expectations={
            "float": {"distribution": "normal", "mean": 70, "stddev": 30},
            "incidence": 0.2,
    },
  ),
  eGFR_operator=patients.comparator_from(
        "eGFR_record",
        return_expectations={
            "rate": "universal",
            "category": {
                "ratios": {  # ~, =, >= , > , < , <=
                    None: 0.10,
                    "~": 0.05,
                    "=": 0.65,
                    ">=": 0.05,
                    ">": 0.05,
                    "<": 0.05,
                    "<=": 0.05,
                }
            },
            "incidence": 0.80,
        },
  ),  
  eGFR_short_record=patients.with_these_clinical_events(
    codelist=eGFR_short_level_codelist,
    find_last_match_in_period=True,
    on_or_before = "start_date",
    returning="numeric_value",
    include_date_of_match=True,
    date_format = "YYYY-MM-DD",
    return_expectations={
            "float": {"distribution": "normal", "mean": 70, "stddev": 30},
            "incidence": 0.2,
    },
  ),
  eGFR_short_operator=patients.comparator_from(
        "eGFR_short_record",
        return_expectations={
            "rate": "universal",
            "category": {
                "ratios": {  # ~, =, >= , > , < , <=
                    None: 0.10,
                    "~": 0.05,
                    "=": 0.65,
                    ">=": 0.05,
                    ">": 0.05,
                    "<": 0.05,
                    "<=": 0.05,
                }
            },
            "incidence": 0.80,
        },
  ),  


#   # ## UKRR variables   #.commented out
#   # # Prevalent cohorts
  #2020
  ukrr_2020 = patients.with_record_in_ukrr(
      from_dataset="2020_prevalence",
      returning="binary_flag",
      return_expectations={
          "incidence": 0.25
      },
  ),
  ukrr_2020_mod = patients.with_record_in_ukrr(
      from_dataset="2020_prevalence",
      returning="treatment_modality_prevalence",
      return_expectations={
              "category": {"ratios": {"ICHD": 0.5, "Tx": 0.5}},
              "incidence": 0.25,
          },
  ),
  ukrr_2020_centre = patients.with_record_in_ukrr(
      from_dataset="2020_prevalence",
      returning="renal_centre",
      return_expectations={
              "category": {"ratios": {"RRK02": 0.5, "RDEE1": 0.5}},
              "incidence": 0.25,
          },
  ),
  ukrr_2020_startmod = patients.with_record_in_ukrr(
  from_dataset="2020_prevalence",
  returning="treatment_modality_start",
  return_expectations={
          "category": {"ratios": {"ICHD": 0.5, "Tx": 0.5}},
          "incidence": 0.25,
      },
  ),
  ukrr_2020_startdate = patients.with_record_in_ukrr(
      from_dataset="2020_prevalence",
      returning="rrt_start_date",
      date_format="YYYY-MM-DD",
      return_expectations={
              "date": {"earliest": "1970-01-01", "latest": "2020-12-31"},
          },
  ),
  #2021
  ukrr_2021 = patients.with_record_in_ukrr(
      from_dataset="2021_prevalence",
      returning="binary_flag",
      return_expectations={
          "incidence": 0.25
      },
  ),
  ukrr_2021_mod = patients.with_record_in_ukrr(
      from_dataset="2021_prevalence",
      returning="treatment_modality_prevalence",
      return_expectations={
              "category": {"ratios": {"ICHD": 0.5, "Tx": 0.5}},
              "incidence": 0.25,
          },
  ),
  ukrr_2021_centre = patients.with_record_in_ukrr(
      from_dataset="2021_prevalence",
      returning="renal_centre",
      return_expectations={
              "category": {"ratios": {"RRK02": 0.5, "RDEE1": 0.5}},
              "incidence": 0.25,
          },
  ),
  ukrr_2021_startmod = patients.with_record_in_ukrr(
  from_dataset="2021_prevalence",
  returning="treatment_modality_start",
  return_expectations={
          "category": {"ratios": {"ICHD": 0.5, "Tx": 0.5}},
          "incidence": 0.25,
      },
  ),
  ukrr_2021_startdate = patients.with_record_in_ukrr(
      from_dataset="2021_prevalence",
      returning="rrt_start_date",
      date_format="YYYY-MM-DD",
      return_expectations={
              "date": {"earliest": "1970-01-01", "latest": "2021-12-31"},
          },
  ),
#   # #2020 CKD
#   # #2020
  ukrr_ckd2020 = patients.with_record_in_ukrr(
      from_dataset="2020_ckd",
      returning="binary_flag",
      return_expectations={
          "incidence": 0.35
      },
  ),
  ukrr_ckd2020_creat = patients.with_record_in_ukrr(
      from_dataset="2020_ckd",
      returning="latest_creatinine",
          return_expectations={
              "int": {"distribution": "normal", "mean": 45, "stddev": 20},
              "incidence": 0.35,
          },
  ),
  ukrr_ckd2020_egfr = patients.with_record_in_ukrr(
      from_dataset="2020_ckd",
      returning="latest_egfr",
          return_expectations={
              "float": {"distribution": "normal", "mean": 20, "stddev": 10},
              "incidence": 0.2,
          },
  ),
  ukrr_ckd2020_centre = patients.with_record_in_ukrr(
      from_dataset="2020_ckd",
      returning="renal_centre",
      return_expectations={
              "category": {"ratios": {"RRK02": 0.5, "RDEE1": 0.5}},
              "incidence": 0.35,
          },
  ),
  #Incident cohort
  #2020
  ukrr_inc2020 = patients.with_record_in_ukrr(
      from_dataset="2020_incidence",
      returning="binary_flag",
      return_expectations={
          "incidence": 0.25
      },
  ),
  ukrr_inc2020_mod = patients.with_record_in_ukrr(
      from_dataset="2020_incidence",
      returning="treatment_modality_start",
      return_expectations={
              "category": {"ratios": {"ICHD": 0.5, "Tx": 0.5}},
              "incidence": 0.25,
          },
  ),
  ukrr_inc2020_centre = patients.with_record_in_ukrr(
      from_dataset="2020_incidence",
      returning="renal_centre",
      return_expectations={
              "category": {"ratios": {"RRK02": 0.5, "RDEE1": 0.5}},
              "incidence": 0.25,
          },
  ),
  ukrr_inc2020_date = patients.with_record_in_ukrr(
      from_dataset="2020_incidence",
      returning="rrt_start_date",
      date_format="YYYY-MM-DD",
      return_expectations={
              "date": {"earliest": "2020-01-01", "latest": "2020-12-31"},
          },
  ),

#   ## Paxlovid - exclusion
#   #  Solid organ transplant (plus solid_organ_transplant_nhsd_snomed defined below)
  solid_organ_transplant_snomed = patients.with_these_clinical_events(
    solid_organ_transplant_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),

  ### contraindicated medication
  # drugs do not use (drugs used in secondary/tertiary care?)
  drugs_do_not_use = patients.with_these_medications(
    codelist = drugs_do_not_use_codes,
    returning = "date",
    on_or_before = "start_date",
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),
  # drugs considering risks and benefits
  drugs_consider_risk = patients.with_these_medications(
    codelist = drugs_consider_risk_codes,
    returning = "date",
    on_or_before = "start_date",
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),

#   ### Children aged less than 18 years
#   #   (defined below)
  
#   ### Pregnancy
  
#   # pregnancy record in last 36 weeks
  preg_36wks_date = patients.with_these_clinical_events(
    pregnancy_primis_codes,
    returning = "date",
    find_last_match_in_period = True,
    between = ["start_date - 252 days", "start_date - 1 day"],
    date_format = "YYYY-MM-DD",
  ),
  
#   # pregnancy OR delivery code since latest pregnancy record:
#   # if one of these codes occurs later than the latest pregnancy code
#   #  this indicates pregnancy has ended, if they are same date assume 
#   #  pregnancy has most likely not ended yet
  pregdel = patients.with_these_clinical_events(
    pregdel_primis_codes,
    returning = "binary_flag",
    find_last_match_in_period = True,
    between = ["preg_36wks_date + 1 day", "start_date - 1 day"],
    date_format = "YYYY-MM-DD",
  ),
  
  pregnancy = patients.satisfying(
    
    """
    gender = 'F' AND preg_age <= 50
    AND (preg_36wks_date AND NOT pregdel)
    """,
    
    gender = patients.sex(
      return_expectations = {
        "rate": "universal",
        "category": {"ratios": {"M": 0.49, "F": 0.51}},
      }
    ),
    
    preg_age = patients.age_as_of(
      "preg_36wks_date",
      return_expectations = {
        "rate": "universal",
        "int": {"distribution": "population_ages"},
        "incidence" : 0.9
      },
    ),
    
  ),
  
  
  
#   # CENSORING ----
  
#   ## Death of any cause
  death_date = patients.died_from_any_cause(
    returning = "date_of_death",
    date_format = "YYYY-MM-DD",
    on_or_after = "start_date",
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date"},
      "incidence": 0.1
    },
  ),
  
  has_died = patients.died_from_any_cause(
    on_or_before = "start_date - 1 day",
    returning = "binary_flag",
  ),
  
  ## De-registration
  dereg_date = patients.date_deregistered_from_all_supported_practices(
    on_or_after = "start_date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date"},
      "incidence": 0.1
    },
  ),
  
  registered_eligible = patients.registered_as_of("covid_test_positive_date"),
  
  ## 1/2/3 months since treatment initiation
  ## AND study end date (today)
  ## end enrollment earlier to account for delay in outcome data update

  # HIGH RISK GROUPS ----
  
  ## NHSD ‘high risk’ cohort (codelist to be defined if/when data avaliable)
  # high_risk_cohort_nhsd = patients.with_these_clinical_events(
  #   high_risk_cohort_nhsd_codes,
  #   between = [campaign_start, index_date],
  #   returning = "date",
  #   date_format = "YYYY-MM-DD",
  #   find_first_match_in_period = True,
  # ),
  
  ## Blueteq ‘high risk’ cohort
  high_risk_cohort_covid_therapeutics = patients.with_covid_therapeutics(
    #with_these_statuses = ["Approved", "Treatment Complete"],
    with_these_therapeutics = ["Sotrovimab", "Paxlovid", "Molnupiravir"],
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "risk_group",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "rate": "universal",
      "incidence": 0.4,
      "category": {
        "ratios": {
          "Downs syndrome": 0.1,
          "sickle cell disease": 0.1,
          "solid cancer": 0.1,
          "haematological diseases,stem cell transplant recipients": 0.1,
          "renal disease,sickle cell disease": 0.1,
          "liver disease": 0.05,
          "IMID": 0.1,
          "IMID,solid cancer": 0.1,
          "haematological malignancies": 0.05,
          "primary immune deficiencies": 0.1,
          "HIV or AIDS": 0.05,
          "NA":0.05,},},
    },
  ),
  
  ## Down's syndrome
  downs_syndrome_nhsd_snomed = patients.with_these_clinical_events(
    downs_syndrome_nhsd_snomed_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),
  
  downs_syndrome_nhsd_icd10 = patients.admitted_to_hospital(
    returning = "date_admitted",
    on_or_before = "start_date",
    with_these_diagnoses = downs_syndrome_nhsd_icd10_codes,
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),
  
  downs_syndrome_nhsd = patients.minimum_of("downs_syndrome_nhsd_snomed", "downs_syndrome_nhsd_icd10"), 
  
  ## Solid cancer
  cancer_opensafely_snomed = patients.with_these_clinical_events(
    combine_codelists(
      non_haematological_cancer_opensafely_snomed_codes,
      lung_cancer_opensafely_snomed_codes,
      chemotherapy_radiotherapy_opensafely_snomed_codes
    ),
    between = ["start_date - 6 months", "start_date"],
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),
  ## Solid cance-updated  
  cancer_opensafely_snomed_new = patients.with_these_clinical_events(
    combine_codelists(
      non_haematological_cancer_opensafely_snomed_codes_new,
      lung_cancer_opensafely_snomed_codes,
      chemotherapy_radiotherapy_opensafely_snomed_codes
    ),
    between = ["start_date - 6 months", "start_date"],
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),    
  cancer_opensafely_snomed_ever = patients.with_these_clinical_events(
    combine_codelists(
      non_haematological_cancer_opensafely_snomed_codes_new,
      lung_cancer_opensafely_snomed_codes,
      chemotherapy_radiotherapy_opensafely_snomed_codes
    ),
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),    

  ## Haematological diseases
  haematopoietic_stem_cell_snomed = patients.with_these_clinical_events(
    haematopoietic_stem_cell_transplant_nhsd_snomed_codes,
    between = ["start_date - 12 months", "start_date"],
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),
  
  haematopoietic_stem_cell_icd10 = patients.admitted_to_hospital(
    returning = "date_admitted",
    between = ["start_date - 12 months", "start_date"],
    with_these_diagnoses = haematopoietic_stem_cell_transplant_nhsd_icd10_codes,
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),
  
  haematopoietic_stem_cell_opcs4 = patients.admitted_to_hospital(
    returning = "date_admitted",
    between = ["start_date - 12 months", "start_date"],
    with_these_procedures = haematopoietic_stem_cell_transplant_nhsd_opcs4_codes,
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2020-02-01"},
      "rate": "exponential_increase",
      "incidence": 0.01,
    },
  ),
  
  haematological_malignancies_snomed = patients.with_these_clinical_events(
    haematological_malignancies_nhsd_snomed_codes,
    between = ["start_date - 24 months", "start_date"],
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),
  
  haematological_malignancies_icd10 = patients.admitted_to_hospital(
    returning = "date_admitted",
    between = ["start_date - 24 months", "start_date"],
    with_these_diagnoses = haematological_malignancies_nhsd_icd10_codes,
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),
  
  sickle_cell_disease_nhsd_snomed = patients.with_these_clinical_events(
    sickle_cell_disease_nhsd_snomed_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),
  
  sickle_cell_disease_nhsd_icd10 = patients.admitted_to_hospital(
    returning = "date_admitted",
    on_or_before = "start_date",
    with_these_diagnoses = sickle_cell_disease_nhsd_icd10_codes,
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),
  
  haematological_disease_nhsd = patients.minimum_of("haematopoietic_stem_cell_snomed", 
                                                    "haematopoietic_stem_cell_icd10", 
                                                    "haematopoietic_stem_cell_opcs4", 
                                                    "haematological_malignancies_snomed", 
                                                    "haematological_malignancies_icd10",
                                                    "sickle_cell_disease_nhsd_snomed", 
                                                    "sickle_cell_disease_nhsd_icd10"), 
  
  haematopoietic_stem_cell_snomed_ever = patients.with_these_clinical_events(
    haematopoietic_stem_cell_transplant_nhsd_snomed_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),
  
  haematopoietic_stem_cell_icd10_ever = patients.admitted_to_hospital(
    returning = "date_admitted",
    on_or_before = "start_date",
    with_these_diagnoses = haematopoietic_stem_cell_transplant_nhsd_icd10_codes,
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),
  
  haematopoietic_stem_cell_opcs4_ever = patients.admitted_to_hospital(
    returning = "date_admitted",
    on_or_before = "start_date",
    with_these_procedures = haematopoietic_stem_cell_transplant_nhsd_opcs4_codes,
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2020-02-01"},
      "rate": "exponential_increase",
      "incidence": 0.01,
    },
  ),
  
  haematological_malignancies_snomed_ever = patients.with_these_clinical_events(
    haematological_malignancies_nhsd_snomed_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),
  
  haematological_malignancies_icd10_ever = patients.admitted_to_hospital(
    returning = "date_admitted",
    on_or_before = "start_date",
    with_these_diagnoses = haematological_malignancies_nhsd_icd10_codes,
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),

  haematological_disease_nhsd_ever = patients.minimum_of("haematopoietic_stem_cell_snomed_ever", 
                                                    "haematopoietic_stem_cell_icd10_ever", 
                                                    "haematopoietic_stem_cell_opcs4_ever", 
                                                    "haematological_malignancies_snomed_ever", 
                                                    "haematological_malignancies_icd10_ever",
                                                    "sickle_cell_disease_nhsd_snomed", 
                                                    "sickle_cell_disease_nhsd_icd10"), 
  

  ## Renal disease
  ckd_stage_5_nhsd_snomed = patients.with_these_clinical_events(
    ckd_stage_5_nhsd_snomed_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),
  
  ckd_stage_5_nhsd_icd10 = patients.admitted_to_hospital(
    returning = "date_admitted",
    on_or_before = "start_date",
    with_these_diagnoses = ckd_stage_5_nhsd_icd10_codes,
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),
  
  ckd_stage_5_nhsd = patients.minimum_of("ckd_stage_5_nhsd_snomed", "ckd_stage_5_nhsd_icd10"), 
  
  ## Liver disease
  liver_disease_nhsd_snomed = patients.with_these_clinical_events(
    liver_disease_nhsd_snomed_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),
  
  liver_disease_nhsd_icd10 = patients.admitted_to_hospital(
    returning = "date_admitted",
    on_or_before = "start_date",
    with_these_diagnoses = liver_disease_nhsd_icd10_codes,
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),
  
  liver_disease_nhsd = patients.minimum_of("liver_disease_nhsd_snomed", "liver_disease_nhsd_icd10"), 
  
  ## Immune-mediated inflammatory disorders (IMID)
  immunosuppresant_drugs_nhsd = patients.with_these_medications(
    codelist = combine_codelists(immunosuppresant_drugs_dmd_codes, immunosuppresant_drugs_snomed_codes),
    returning = "date",
    between = ["start_date - 6 months", "start_date"],
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),
  
  oral_steroid_drugs_nhsd = patients.with_these_medications(
    codelist = combine_codelists(oral_steroid_drugs_dmd_codes, oral_steroid_drugs_snomed_codes),
    returning = "date",
    between = ["start_date - 12 months", "start_date"],
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),
  
  oral_steroid_drug_nhsd_3m_count = patients.with_these_medications(
    codelist = combine_codelists(oral_steroid_drugs_dmd_codes, oral_steroid_drugs_snomed_codes),
    returning = "number_of_matches_in_period",
    between = ["start_date - 3 months", "start_date"],
    return_expectations = {"incidence": 0.1,
      "int": {"distribution": "normal", "mean": 2, "stddev": 1},
    },
  ),
  
  oral_steroid_drug_nhsd_12m_count = patients.with_these_medications(
    codelist = combine_codelists(oral_steroid_drugs_dmd_codes, oral_steroid_drugs_snomed_codes),
    returning = "number_of_matches_in_period",
    between = ["start_date - 12 months", "start_date"],
    return_expectations = {"incidence": 0.1,
      "int": {"distribution": "normal", "mean": 3, "stddev": 1},
    },
  ),
  
  # imid_nhsd = patients.minimum_of("immunosuppresant_drugs_nhsd", "oral_steroid_drugs_nhsd"), - define in processing script
  immunosuppresant_drugs_nhsd_ever = patients.with_these_medications(
    codelist = combine_codelists(immunosuppresant_drugs_dmd_codes, immunosuppresant_drugs_snomed_codes),
    returning = "date",
    on_or_before = "start_date",
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),
  
  oral_steroid_drugs_nhsd_ever = patients.with_these_medications(
    codelist = combine_codelists(oral_steroid_drugs_dmd_codes, oral_steroid_drugs_snomed_codes),
    returning = "date",
    on_or_before = "start_date",
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),  
  
  ## Primary immune deficiencies
  immunosupression_nhsd = patients.with_these_clinical_events(
    immunosupression_nhsd_codes,
    on_or_before = "start_date",
    returning = "date",
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),
  ## Primary immune deficiencies-updated
  immunosupression_nhsd_new = patients.with_these_clinical_events(
    immunosupression_nhsd_codes_new,
    on_or_before = "start_date",
    returning = "date",
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),  
  ## HIV/AIDs
  hiv_aids_nhsd_snomed = patients.with_these_clinical_events(
    hiv_aids_nhsd_snomed_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),
  
  hiv_aids_nhsd_icd10 = patients.admitted_to_hospital(
    returning = "date_admitted",
    on_or_before = "start_date",
    with_these_diagnoses = hiv_aids_nhsd_icd10_codes,
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),
  
  hiv_aids_nhsd = patients.minimum_of("hiv_aids_nhsd_snomed", "hiv_aids_nhsd_icd10"),
  
  ## Solid organ transplant
  solid_organ_transplant_nhsd_snomed = patients.with_these_clinical_events(
    solid_organ_transplant_nhsd_snomed_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),
  solid_organ_nhsd_snomed_new = patients.with_these_clinical_events(
    solid_organ_transplant_nhsd_snomed_codes_new,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),  
  solid_organ_transplant_nhsd_opcs4 = patients.admitted_to_hospital(
    returning = "date_admitted",
    on_or_before = "start_date",
    with_these_procedures = solid_organ_transplant_nhsd_opcs4_codes,
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2020-02-01"},
      "rate": "exponential_increase",
      "incidence": 0.01,
    },
  ),
  
  transplant_all_y_codes_opcs4 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_procedures = replacement_of_organ_transplant_nhsd_opcs4_codes,
    on_or_before = "start_date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2020-02-01"},
      "rate": "exponential_increase",
      "incidence": 0.01,
    },
  ),
  
  transplant_thymus_opcs4 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_procedures = thymus_gland_transplant_nhsd_opcs4_codes,
    between = ["transplant_all_y_codes_opcs4","transplant_all_y_codes_opcs4"],
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2020-02-01"},
      "rate": "exponential_increase",
      "incidence": 0.01,
    },
  ),
  
  transplant_conjunctiva_y_code_opcs4 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_procedures = conjunctiva_y_codes_transplant_nhsd_opcs4_codes,
    on_or_before = "start_date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2020-02-01"},
      "rate": "exponential_increase",
      "incidence": 0.01,
    },
  ),
  
  transplant_conjunctiva_opcs4 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_procedures = conjunctiva_transplant_nhsd_opcs4_codes,
    between = ["transplant_conjunctiva_y_code_opcs4","transplant_conjunctiva_y_code_opcs4"],
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2020-02-01"},
      "rate": "exponential_increase",
      "incidence": 0.01,
    },
  ),
  
  transplant_stomach_opcs4 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_procedures = stomach_transplant_nhsd_opcs4_codes,
    between = ["transplant_all_y_codes_opcs4","transplant_all_y_codes_opcs4"],
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2020-02-01"},
      "rate": "exponential_increase",
      "incidence": 0.01,
    },
  ),
  
  transplant_ileum_1_Y_codes_opcs4 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_procedures = ileum_1_y_codes_transplant_nhsd_opcs4_codes,
    on_or_before = "start_date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2020-02-01"},
      "rate": "exponential_increase",
      "incidence": 0.01,
    },
  ),
  
  transplant_ileum_2_Y_codes_opcs4 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_procedures = ileum_2_y_codes_transplant_nhsd_opcs4_codes,
    on_or_before = "start_date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2020-02-01"},
      "rate": "exponential_increase",
      "incidence": 0.01,
    },
  ),
  
  transplant_ileum_1_opcs4 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_procedures = ileum_1_transplant_nhsd_opcs4_codes,
    between = ["transplant_ileum_1_Y_codes_opcs4","transplant_ileum_1_Y_codes_opcs4"],
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2020-02-01"},
      "rate": "exponential_increase",
      "incidence": 0.01,
    },
  ),
  
  transplant_ileum_2_opcs4 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_procedures = ileum_2_transplant_nhsd_opcs4_codes,
    between = ["transplant_ileum_2_Y_codes_opcs4","transplant_ileum_2_Y_codes_opcs4"],
    date_format = "YYYY-MM-DD",
    find_first_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2020-02-01"},
      "rate": "exponential_increase",
      "incidence": 0.01,
    },
  ),
  
  solid_organ_transplant_nhsd = patients.minimum_of("solid_organ_transplant_nhsd_snomed", "solid_organ_transplant_nhsd_opcs4",
                                                    "transplant_thymus_opcs4", "transplant_conjunctiva_opcs4", "transplant_stomach_opcs4",
                                                    "transplant_ileum_1_opcs4","transplant_ileum_2_opcs4"), 
  solid_organ_transplant_nhsd_new = patients.minimum_of("solid_organ_nhsd_snomed_new", "solid_organ_transplant_nhsd_opcs4",
                                                    "transplant_thymus_opcs4", "transplant_conjunctiva_opcs4", "transplant_stomach_opcs4",
                                                    "transplant_ileum_1_opcs4","transplant_ileum_2_opcs4"), 
                                                      
  ## Rare neurological conditions
  
  ### Multiple sclerosis
  multiple_sclerosis_nhsd_snomed = patients.with_these_clinical_events(
    multiple_sclerosis_nhsd_snomed_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),
  
  multiple_sclerosis_nhsd_icd10 = patients.admitted_to_hospital(
    returning = "date_admitted",
    on_or_before = "start_date",
    with_these_diagnoses = multiple_sclerosis_nhsd_icd10_codes,
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),
  
  multiple_sclerosis_nhsd = patients.minimum_of("multiple_sclerosis_nhsd_snomed", "multiple_sclerosis_nhsd_icd10"), 
  
  ### Motor neurone disease
  motor_neurone_disease_nhsd_snomed = patients.with_these_clinical_events(
    motor_neurone_disease_nhsd_snomed_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),
  
  motor_neurone_disease_nhsd_icd10 = patients.admitted_to_hospital(
    returning = "date_admitted",
    on_or_before = "start_date",
    with_these_diagnoses = motor_neurone_disease_nhsd_icd10_codes,
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),
  
  motor_neurone_disease_nhsd = patients.minimum_of("motor_neurone_disease_nhsd_snomed", "motor_neurone_disease_nhsd_icd10"),
  
  ### Myasthenia gravis
  myasthenia_gravis_nhsd_snomed = patients.with_these_clinical_events(
    myasthenia_gravis_nhsd_snomed_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),
  
  myasthenia_gravis_nhsd_icd10 = patients.admitted_to_hospital(
    returning = "date_admitted",
    on_or_before = "start_date",
    with_these_diagnoses = myasthenia_gravis_nhsd_icd10_codes,
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),
  
  myasthenia_gravis_nhsd = patients.minimum_of("myasthenia_gravis_nhsd_snomed", "myasthenia_gravis_nhsd_icd10"),
  
  ### Huntington’s disease
  huntingtons_disease_nhsd_snomed = patients.with_these_clinical_events(
    huntingtons_disease_nhsd_snomed_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ),
  
  huntingtons_disease_nhsd_icd10 = patients.admitted_to_hospital(
    returning = "date_admitted",
    on_or_before = "start_date",
    with_these_diagnoses = huntingtons_disease_nhsd_icd10_codes,
    find_last_match_in_period = True,
    date_format = "YYYY-MM-DD",
  ),
  
  huntingtons_disease_nhsd = patients.minimum_of("huntingtons_disease_nhsd_snomed", "huntingtons_disease_nhsd_icd10"),
  
  
  
#   # CLINICAL/DEMOGRAPHIC COVARIATES ----
  
#   ## Age
  age = patients.age_as_of(
    "start_date - 1 day",
    return_expectations = {
      "rate": "universal",
      "int": {"distribution": "population_ages"},
      "incidence" : 0.9
    },
  ),
  
  ## Sex
  sex = patients.sex(
    return_expectations = {
      "rate": "universal",
      "category": {"ratios": {"M": 0.49, "F": 0.51}},
    }
  ),
  
  ## Ethnicity
  ethnicity = patients.categorised_as(
            {"Missing": "DEFAULT",
            "White": "eth='1' OR (NOT eth AND ethnicity_sus='1')", 
            "Mixed": "eth='2' OR (NOT eth AND ethnicity_sus='2')", 
            "South Asian": "eth='3' OR (NOT eth AND ethnicity_sus='3')", 
            "Black": "eth='4' OR (NOT eth AND ethnicity_sus='4')",  
            "Other": "eth='5' OR (NOT eth AND ethnicity_sus='5')",
            }, 
            return_expectations={
            "category": {"ratios": {"White": 0.6, "Mixed": 0.1, "South Asian": 0.1, "Black": 0.1, "Other": 0.1}},
            "incidence": 0.4,
            },

            ethnicity_sus = patients.with_ethnicity_from_sus(
                returning="group_6",  
                use_most_frequent_code=True,
                return_expectations={
                    "category": {"ratios": {"1": 0.6, "2": 0.1, "3": 0.1, "4": 0.1, "5": 0.1}},
                    "incidence": 0.4,
                    },
            ),

            eth=patients.with_these_clinical_events(
                ethnicity_primis_snomed_codes,
                returning="category",
                find_last_match_in_period=True,
                on_or_before="today",
                return_expectations={
                    "category": {"ratios": {"1": 0.6, "2": 0.1, "3": 0.1, "4":0.1,"5": 0.1}},
                    "incidence": 0.75,
                },
            ),
    ),
  
  ## Index of multiple deprivation
  imd = patients.categorised_as(
    {     "0": "DEFAULT",
          "1": "index_of_multiple_deprivation >= 0 AND index_of_multiple_deprivation < 32800*1/5",
          "2": "index_of_multiple_deprivation >= 32800*1/5 AND index_of_multiple_deprivation < 32800*2/5",
          "3": "index_of_multiple_deprivation >= 32800*2/5 AND index_of_multiple_deprivation < 32800*3/5",
          "4": "index_of_multiple_deprivation >= 32800*3/5 AND index_of_multiple_deprivation < 32800*4/5",
          "5": "index_of_multiple_deprivation >= 32800*4/5 AND index_of_multiple_deprivation <= 32800",
    },
    index_of_multiple_deprivation = patients.address_as_of(
      "start_date",
      returning = "index_of_multiple_deprivation",
      round_to_nearest = 100,
    ),
    return_expectations = {
      "rate": "universal",
      "category": {
        "ratios": {
          "0": 0.01,
          "1": 0.20,
          "2": 0.20,
          "3": 0.20,
          "4": 0.20,
          "5": 0.19,
        }},
    },
  ),
  
  ## Region - NHS England 9 regions
  region_nhs = patients.registered_practice_as_of(
    "start_date",
    returning = "nuts1_region_name",
    return_expectations = {
      "rate": "universal",
      "category": {
        "ratios": {
          "North East": 0.1,
          "North West": 0.1,
          "Yorkshire and The Humber": 0.1,
          "East Midlands": 0.1,
          "West Midlands": 0.1,
          "East": 0.1,
          "London": 0.2,
          "South West": 0.1,
          "South East": 0.1,},},
    },
  ),
  
  region_covid_therapeutics = patients.with_covid_therapeutics(
    #with_these_statuses = ["Approved", "Treatment Complete"],
    with_these_therapeutics = ["Sotrovimab", "Paxlovid", "Molnupiravir"],
    with_these_indications = "non_hospitalised",
    on_or_after = "start_date",
    find_first_match_in_period = True,
    returning = "region",
    return_expectations = {
      "rate": "universal",
      "category": {
        "ratios": {
          "North East": 0.1,
          "North West": 0.1,
          "Yorkshire and The Humber": 0.1,
          "East Midlands": 0.1,
          "West Midlands": 0.1,
          "East": 0.1,
          "London": 0.2,
          "South West": 0.1,
          "South East": 0.1,},},
    },
  ),
  
  # STP (NHS administration region based on geography, currently closest match to CMDU)
  stp = patients.registered_practice_as_of(
    "start_date",
    returning = "stp_code",
    return_expectations = {
      "rate": "universal",
      "category": {
        "ratios": {
          "STP1": 0.1,
          "STP2": 0.1,
          "STP3": 0.1,
          "STP4": 0.1,
          "STP5": 0.1,
          "STP6": 0.1,
          "STP7": 0.1,
          "STP8": 0.1,
          "STP9": 0.1,
          "STP10": 0.1,
        }
      },
    },
  ),
  
  # Rurality
  rural_urban = patients.address_as_of(
    "start_date",
    returning = "rural_urban_classification",
    return_expectations = {
      "rate": "universal",
      "category": {"ratios": {1: 0.125, 2: 0.125, 3: 0.125, 4: 0.125, 5: 0.125, 6: 0.125, 7: 0.125, 8: 0.125}},
      "incidence": 1,
    },
  ),
  
  
  
#   # CLINICAL GROUPS ----
  
  ## Autism
  autism_nhsd = patients.with_these_clinical_events(
    autism_nhsd_snomed_codes,
    on_or_before = "start_date",
    returning = "binary_flag",
    return_expectations = {"incidence": 0.3}
  ),
  
  ## Care home 
  care_home_primis = patients.with_these_clinical_events(
    care_home_primis_snomed_codes,
    returning = "binary_flag",
    on_or_before = "start_date",
    return_expectations = {"incidence": 0.15,}
  ),
  
  ## Dementia
  dementia_nhsd = patients.satisfying(
    
    """
    dementia_all
    AND
    age > 39
    """, 
    
    return_expectations = {
      "incidence": 0.01,
    },
    
    dementia_all = patients.with_these_clinical_events(
      dementia_nhsd_snomed_codes,
      on_or_before = "start_date",
      returning = "binary_flag",
      return_expectations = {"incidence": 0.05}
    ),
    
  ),
  
  ## Housebound
  housebound_opensafely = patients.satisfying(
    """housebound_date
                AND NOT no_longer_housebound
                AND NOT moved_into_care_home""",
    return_expectations={
      "incidence": 0.01,
    },
    
    housebound_date = patients.with_these_clinical_events( 
      housebound_opensafely_snomed_codes, 
      on_or_before = "start_date",
      find_last_match_in_period = True,
      returning = "date",
      date_format = "YYYY-MM-DD",
    ),   
    
    no_longer_housebound = patients.with_these_clinical_events( 
      no_longer_housebound_opensafely_snomed_codes, 
      on_or_after = "housebound_date",
    ),
    
    moved_into_care_home = patients.with_these_clinical_events(
      care_home_primis_snomed_codes,
      on_or_after = "housebound_date",
    ),
    
  ),
  
  ## Learning disability
  learning_disability_primis = patients.with_these_clinical_events(
    wider_ld_primis_snomed_codes,
    on_or_before = "start_date",
    returning = "binary_flag",
    return_expectations = {"incidence": 0.2}
  ),
  
  ### Serious Mental Illness
  serious_mental_illness_nhsd = patients.with_these_clinical_events(
    serious_mental_illness_nhsd_snomed_codes,
    on_or_before = "start_date",
    returning = "binary_flag",
    return_expectations = {"incidence": 0.1}
  ),
  
  ## Sickle cell disease
  sickle_cell_disease_nhsd = patients.minimum_of("sickle_cell_disease_nhsd_snomed", "sickle_cell_disease_nhsd_icd10"), 
  
  ## Vaccination status
  vaccination_status = patients.categorised_as(
    {
      "Un-vaccinated": "DEFAULT",
      "Un-vaccinated (declined)": """ covid_vax_declined AND NOT (covid_vax_1 OR covid_vax_2 OR covid_vax_3 OR covid_vax_4)""",
      "One vaccination": """ covid_vax_1 AND NOT covid_vax_2 """,
      "Two vaccinations": """ covid_vax_2 AND NOT covid_vax_3 """,
      "Three vaccinations": """ covid_vax_3 AND NOT covid_vax_4 """,
      "Four or more vaccinations": """ covid_vax_4 """
    },
    
    # first vaccine from during trials and up to treatment/test date
    covid_vax_1 = patients.with_tpp_vaccination_record(
      target_disease_matches = "SARS-2 CORONAVIRUS",
      between = ["2020-06-08", "start_date"],
      find_first_match_in_period = True,
      returning = "date",
      date_format = "YYYY-MM-DD"
    ),
    
    covid_vax_2 = patients.with_tpp_vaccination_record(
      target_disease_matches = "SARS-2 CORONAVIRUS",
      between = ["covid_vax_1 + 19 days", "start_date"],
      find_first_match_in_period = True,
      returning = "date",
      date_format = "YYYY-MM-DD"
    ),
    
    covid_vax_3 = patients.with_tpp_vaccination_record(
      target_disease_matches = "SARS-2 CORONAVIRUS",
      between = ["covid_vax_2 + 56 days", "start_date"],
      find_first_match_in_period = True,
      returning = "date",
      date_format = "YYYY-MM-DD"
    ),
    # add 4th
    covid_vax_4 = patients.with_tpp_vaccination_record(
      target_disease_matches = "SARS-2 CORONAVIRUS",
      between = ["covid_vax_3 + 56 days", "start_date"],
      find_first_match_in_period = True,
      returning = "date",
      date_format = "YYYY-MM-DD"
    ),

    covid_vax_declined = patients.with_these_clinical_events(
      covid_vaccine_declined_codes,
      returning="binary_flag",
      on_or_before = "start_date",
    ),

    return_expectations = {
      "rate": "universal",
      "category": {
        "ratios": {
          "Un-vaccinated": 0.1,
          "Un-vaccinated (declined)": 0.1,
          "One vaccination": 0.1,
          "Two vaccinations": 0.1,
          "Three vaccinations": 0.5,
          "Four or more vaccinations": 0.1,
        }
      },
    },
  ),
  # latest vaccination date
  last_vaccination_date = patients.with_tpp_vaccination_record(
      target_disease_matches = "SARS-2 CORONAVIRUS",
      on_or_before = "start_date",
      find_last_match_in_period = True,
      returning = "date",
      date_format = "YYYY-MM-DD",
      return_expectations={
            "date": {"earliest": "2020-06-08", "latest": "today"},
            "incidence": 0.95,
      }
  ),

  # CLINICAL CO-MORBIDITIES TBC ----
  #BMI, diabetes, hypertension, chronic heart diseases, Chronic respiratory disease, SGTF indicator
  # adapted codes from https://github.com/opensafely/bmi-short-data-report/
  bmi=patients.most_recent_bmi(
        on_or_before="start_date",
        minimum_age_at_measurement=18,
        include_measurement_date=True,
        date_format="YYYY-MM-DD",
        return_expectations={
            "date": {"earliest": "2020-01-01", "latest": "today"},
            "float": {"distribution": "normal", "mean": 28, "stddev": 8},
            "incidence": 0.95,
        }
  ),
  #also use bmi_code_snomed and height and weight?

  # Diabetes
  diabetes=patients.with_these_clinical_events(
        diabetes_codes,
        on_or_before="start_date",
        returning="binary_flag",
        return_expectations={"incidence": 0.1, },
  ),
  # Chronic cardiac disease
  chronic_cardiac_disease=patients.with_these_clinical_events(
        chronic_cardiac_dis_codes,
        on_or_before="start_date",
        returning="binary_flag",
        return_expectations={"incidence": 0.1, },
  ),
  # Hypertension
  hypertension=patients.with_these_clinical_events(
        hypertension_codes,
        on_or_before="start_date",
        returning="binary_flag",
        return_expectations={"incidence": 0.1, },
  ),
  # Chronic respiratory disease
  chronic_respiratory_disease=patients.with_these_clinical_events(
        chronic_respiratory_dis_codes,
        on_or_before="start_date",
        returning="binary_flag",
        return_expectations={"incidence": 0.1, },
  ),

# SGTF indicator and Variant
  sgtf=patients.with_test_result_in_sgss(
       pathogen="SARS-CoV-2",
       test_result="positive",
       find_first_match_in_period=True,
       between=["covid_test_positive_date","covid_test_positive_date + 30 days"],
       returning="s_gene_target_failure",
       return_expectations={
            "rate": "universal",
            "category": {"ratios": {"0": 0.7, "1": 0.1, "9": 0.1, "": 0.1}},
       },
  ), 
  variant_recorded=patients.with_test_result_in_sgss(
       pathogen="SARS-CoV-2",
       test_result="positive",
       find_first_match_in_period=True,
       restrict_to_earliest_specimen_date=False,
       between=["covid_test_positive_date","covid_test_positive_date + 30 days"],
       returning="variant",
       return_expectations={
            "rate": "universal",
            "category": {"ratios": {"B.1.617.2": 0.7, "VOC-21JAN-02": 0.2, "": 0.1}},
       },
  ),     
  # new sgtf data in "all tests dataset"
  sgtf_new=patients.with_test_result_in_sgss(
       pathogen="SARS-CoV-2",
       test_result="positive",
       find_first_match_in_period=True,
       restrict_to_earliest_specimen_date=False,
       between=["covid_test_positive_date","covid_test_positive_date + 30 days"],
       returning="s_gene_target_failure",
       return_expectations={
            "rate": "universal",
            "category": {"ratios": {"0": 0.7, "1": 0.1, "9": 0.1, "": 0.1}},
       },
  ), 


  # OUTCOMES ----
  
  ## COVID re-infection
  covid_positive_test_60_days_post_elig_or_treat = patients.with_test_result_in_sgss( #.-2024-0117-renamed as covid_positive_test_60_days_post_elig_or_treat(old: covid_positive_test_30_days_post_elig_or_treat)
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    returning = "date",
    date_format = "YYYY-MM-DD",
    #on_or_after = "start_date + 60 days", ## old code: (+ 30 days)#-between(start_date + 60 days,start_date + 60 days+6 months,)
    between=["start_date + 60 days","start_date + 240 days"] ##
    find_first_match_in_period = True,
    restrict_to_earliest_specimen_date = False,
    return_expectations = {
      "date": {"earliest": "2022-01-30"},
      "rate": "exponential_increase",
      "incidence": 0.4
    },
  ),
  
  ## COVID-related hospitalisation 
  # extract multiple COVID hosp events per patient because the first hosp may be for receiving sotro or day cases (Day 0 and 1):
  covid_hosp_outcome_date0 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    between = ["start_date", "start_date"],
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-16"},
      "rate": "uniform",
      "incidence": 0.1
    },
  ),  
  # in case one patient had admission records on both day 0 and 1
  covid_hosp_outcome_date1 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    between = ["start_date + 1 day", "start_date + 1 day"],
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-17"},
      "rate": "uniform",
      "incidence": 0.1
    },
  ),  
  covid_hosp_outcome_date2_60d = patients.admitted_to_hospital( #covid_hosp_outcome_date2
    returning = "date_admitted",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "start_date + 60 days",  # .-#old: on_or_after = "start_date + 2 days",
    
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-18"},
      "rate": "uniform",
      "incidence": 0.40
    },
  ),
  # capture and exclude COVID hospital admission/death on the start date

  # return discharge date to (make sure) identify and ignore day cases
  covid_hosp_discharge_date0 = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "covid_hosp_outcome_date0",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-16"},
      "rate": "uniform",
      "incidence": 0.1
    },
  ),  
  covid_hosp_discharge_date1 = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "covid_hosp_outcome_date1",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-17"},
      "rate": "uniform",
      "incidence": 0.1
    },
  ),  
  covid_hosp_discharge_date2_60d = patients.admitted_to_hospital( ## covid_hosp_discharge_date2
    returning = "date_discharged",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "covid_hosp_outcome_date2_60d",  ##covid_hosp_outcome_date2(covid_hosp_outcome_date2_60d)
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-18"},
      "rate": "uniform",
      "incidence": 0.40
    },
  ),  

  # return admission method to identify planned admissions (for sotro injection)
  covid_hosp_admission_method = patients.admitted_to_hospital(
    returning = "admission_method",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "start_date",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "rate": "universal",
      "category": {"ratios": {"21": 0.7, "22": 0.3}},
      "incidence": 0.4,
    },
  ),  

  # sensitivity analysis: emergency admissions only to ignore incidental COVID or patients receiving sotro in hospitals (planned admission)
  # separate day 0,1,2 to identify day case
  covid_hosp_date_emergency0 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    between = ["start_date", "start_date"],
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-16"},
      "rate": "uniform",
      "incidence": 0.10
    },
  ),
  covid_hosp_date_emergency1 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    between = ["start_date + 1 day", "start_date + 1 day"],
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-16"},
      "rate": "uniform",
      "incidence": 0.10
    },
  ),
  covid_hosp_date_emergency2 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "start_date + 2 days",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-16"},
      "rate": "uniform",
      "incidence": 0.40
    },
  ),
  covid_emerg_discharge_date0 = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "covid_hosp_date_emergency0",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-16"},
      "rate": "uniform",
      "incidence": 0.1
    },
  ),  
  covid_emerg_discharge_date1 = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "covid_hosp_date_emergency1",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-17"},
      "rate": "uniform",
      "incidence": 0.1
    },
  ),  
  covid_emerg_discharge_date2 = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "covid_hosp_date_emergency2",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-18"},
      "rate": "uniform",
      "incidence": 0.40
    },
  ),  

#   # identify and ignore COVID hospital admissions for community mAbs procedure on Day 0 or Day 1*
  covid_hosp_date_mabs_procedure = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    with_these_procedures=mabs_procedure_codes,
    on_or_after = "start_date",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-16"},
      "rate": "uniform",
      "incidence": 0.1
    },
  ),
  # add mab record with covid as any diagnosis
  covid_hosp_date_mabs_not_pri = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    with_these_procedures=mabs_procedure_codes,
    on_or_after = "start_date",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-16"},
      "rate": "uniform",
      "incidence": 0.1
    },
  ),  
  # with_these_diagnoses (sensitivity analysis)
  covid_hosp_date0_not_primary = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    between = ["start_date", "start_date"],
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-16"},
      "rate": "uniform",
      "incidence": 0.14
    },
  ),
  # in case one patient had admission records on both day 0 and 1
  covid_hosp_date1_not_primary = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    between = ["start_date + 1 day", "start_date + 1 day"],
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-17"},
      "rate": "uniform",
      "incidence": 0.14
    },
  ),
  covid_hosp_date2_not_primary = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "start_date + 2 days",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-18"},
      "rate": "uniform",
      "incidence": 0.46
    },
  ),

  covid_discharge_date0_not_pri = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_these_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "covid_hosp_date0_not_primary",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-16"},
      "rate": "uniform",
      "incidence": 0.14
    },
  ),  
  covid_discharge_date1_not_pri = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_these_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "covid_hosp_date1_not_primary",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-17"},
      "rate": "uniform",
      "incidence": 0.14
    },
  ),  
  covid_discharge_date2_not_pri = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_these_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "covid_hosp_date2_not_primary",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-18"},
      "rate": "uniform",
      "incidence": 0.46
    },
  ),  

  ## Critical care days for COVID-related hospitalisation 
  covid_hosp_critical_care = patients.admitted_to_hospital(
    returning = "days_in_critical_care",
    with_these_primary_diagnoses = covid_icd10_codes,
    between = ["start_date + 1 day", "start_date + 28 days"],
    find_first_match_in_period = True,
    return_expectations = {
      "category": {"ratios": {"20": 0.5, "40": 0.5}},
      "incidence": 0.4,
    },
  ),

  ## COVID related death
  death_with_covid_on_the_death_certificate_date = patients.with_these_codes_on_death_certificate(
    covid_icd10_codes,
    returning = "date_of_death",
    date_format = "YYYY-MM-DD",
    on_or_after = "start_date",
    return_expectations = {
      "date": {"earliest": "2021-01-01", "latest" : "today"},
      "rate": "uniform",
      "incidence": 0.6},
  ),
  ## COVID related death - COVID as underlying cause
  death_with_covid_underlying_date = patients.with_these_codes_on_death_certificate(
    covid_icd10_codes,
    returning = "date_of_death",
    date_format = "YYYY-MM-DD",
    on_or_after = "start_date",
    match_only_underlying_cause=True,
    return_expectations = {
      "date": {"earliest": "2021-01-01", "latest" : "today"},
      "rate": "uniform",
      "incidence": 0.6},
  ),  

#all-cause hosp; all-cause death already defined
  hospitalisation_outcome_date0 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    between = ["start_date", "start_date"],
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-16"},
      "rate": "uniform",
      "incidence": 0.2
    },
  ),
  # in case one patient had admission records on both day 0 and 1
  hospitalisation_outcome_date1 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    between = ["start_date + 1 day", "start_date + 1 day"],
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-17"},
      "rate": "uniform",
      "incidence": 0.2
    },
  ),
  hospitalisation_outcome_date2 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "start_date + 2 days",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-18"},
      "rate": "uniform",
      "incidence": 0.6
    },
  ),
  hosp_discharge_date0 = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "hospitalisation_outcome_date0",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-16"},
      "rate": "uniform",
      "incidence": 0.2
    },
  ),
  hosp_discharge_date1 = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "hospitalisation_outcome_date1",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-17"},
      "rate": "uniform",
      "incidence": 0.2
    },
  ),
  hosp_discharge_date2 = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "hospitalisation_outcome_date2",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-18"},
      "rate": "uniform",
      "incidence": 0.6
    },
  ),
  # add mab record with all-cause hosp
  covid_hosp_date_mabs_all_cause = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    with_these_procedures=mabs_procedure_codes,
    on_or_after = "start_date",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-16"},
      "rate": "uniform",
      "incidence": 0.1
    },
  ),  


#safety outcomes? 
)

