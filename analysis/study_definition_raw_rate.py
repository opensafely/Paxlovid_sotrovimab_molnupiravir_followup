## Adapted codes from https://github.com/opensafely/antibody-and-antiviral-deployment
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
end_date = date.today().isoformat()

## Define study population and variables
study = StudyDefinition(

  ## Configure the expectations framework
  default_expectations = {
    "date": {"earliest": "2021-11-01", "latest": "today"},
    "rate": "uniform",
    "incidence": 0.05,
  },
    
  # POPULATION ----
  population = patients.satisfying(
    """
    age >= 18 AND age < 110
    AND (registered_eligible OR registered_treated)
  """,
  ),
  #require covid_test_positive_date<=date_treated (sensitivity analysis)
  #loose "AND (covid_test_positive AND NOT covid_positive_previous_30_days)"
  #AND NOT pregnancy (sensitivity analysis)
  #AND NOT (casirivimab_covid_therapeutics OR remdesivir_covid_therapeutics) (sensitivity analysis)
  
  # index date need to be 1 month earlier than the beginning of study period!
  index_date = "2021-12-16",

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
  ### Paxlovid
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
  ### Molnupiravir
  molnupiravir_covid_therapeutics = patients.with_covid_therapeutics(
    #with_these_statuses = ["Approved", "Treatment Complete"],
    with_these_therapeutics = "Molnupiravir",
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-10"},
      "incidence": 0.05
    },
  ), 

  ## Remdesivir
  remdesivir_covid_therapeutics = patients.with_covid_therapeutics(
    #with_these_statuses = ["Approved", "Treatment Complete"],
    with_these_therapeutics = "Remdesivir",
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.05
    },
  ),
  
  ### Casirivimab and imdevimab
  casirivimab_covid_therapeutics = patients.with_covid_therapeutics(
    #with_these_statuses = ["Approved", "Treatment Complete"],
    with_these_therapeutics = "Casirivimab and imdevimab",
    with_these_indications = "non_hospitalised",
    on_or_after = "index_date",
    find_first_match_in_period = True,
    returning = "date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-16"},
      "incidence": 0.05
    },
  ), 

  ## Date treated
  date_treated = patients.minimum_of(
    "sotrovimab_covid_therapeutics",
    "paxlovid_covid_therapeutics",
    "molnupiravir_covid_therapeutics",
    "remdesivir_covid_therapeutics",
    "casirivimab_covid_therapeutics",
  ),
  
  registered_treated = patients.registered_as_of("date_treated"), 

  # OVERALL ELIGIBILITY CRITERIA VARIABLES ----
  
  ## Inclusion criteria variables
  
  ### First positive SARS-CoV-2 test
  # Note patients are eligible for treatment if diagnosed <=5d ago
  covid_test_positive = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    returning = "binary_flag",
    on_or_after = "index_date",
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
    on_or_after = "index_date",
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date"},
      "incidence": 0.9
    },
  ),
  
  ### Positive covid test last 30 days 
  # (note this will only apply to patients who first tested positive towards the beginning
  # of the study period)
  covid_positive_previous_30_days = patients.with_test_result_in_sgss(
    pathogen = "SARS-CoV-2",
    test_result = "positive",
    returning = "binary_flag",
    between = ["covid_test_positive_date - 30 days", "covid_test_positive_date - 1 day"],
    find_last_match_in_period = True,
    restrict_to_earliest_specimen_date = False,
    return_expectations = {
      "incidence": 0.05
    },
  ),
  
  registered_eligible = patients.registered_as_of("covid_test_positive_date"),

  ## Exclusion criteria variables
  
  ### Require hospitalisation for COVID-19
  ## NB this data lags behind the therapeutics/testing data so may be missing
  primary_covid_hospital_discharge_date = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    #with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    between = ["covid_test_positive_date - 30 days", "covid_test_positive_date - 1 day"],
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
    between = ["covid_test_positive_date - 30 days", "covid_test_positive_date - 1 day"],
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
    #with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    between = ["covid_test_positive_date - 30 days", "covid_test_positive_date - 1 day"],
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date - 1 day"},
      "rate": "uniform",
      "incidence": 0.05
    },
  ),
  any_covid_hospital_admission_date = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    #with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    between = ["covid_test_positive_date - 30 days", "covid_test_positive_date - 1 day"],
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date - 1 day"},
      "rate": "uniform",
      "incidence": 0.05
    },
  ),
  
  ## Start date for extracting variables
  start_date = patients.minimum_of(
    "date_treated",
    "covid_test_positive_date",
  ),

  ### New supplemental oxygen requirement specifically for the management of COVID-19 symptoms
  #   (not currently possible to define/code)
  
  ### Known hypersensitivity reaction to the active substances or to any of the excipients of treatments
  #   (not currently possible to define/code)
  
  
  
  # TREATMENT SPECIFIC ELIGIBILITY CRITERIA VARIABLES ----
  
  ## Paxlovid - inclusion
  
  ### Clinical judgement deems that an antiviral is the preferred option
  #   (not currently possible to define/code)
  
  ### Treatment is commenced within 5 days of symptom onset
  #   (defined above)
  
  ### The patient does NOT have a history of advanced decompensated liver cirrhosis or stage 3-5 chronic kidney disease
  #   (OR also use renal and liver high rosk cohorts defined below)
  #  advanced decompensated liver cirrhosis
  advanced_decompensated_cirrhosis = patients.with_these_clinical_events(
    advanced_decompensated_cirrhosis_snomed_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
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
  #  regular ascitic drainage (opcs4_codes in hospital??)
  ascitic_drainage_snomed = patients.with_these_clinical_events(
    ascitic_drainage_snomed_codes,
    on_or_before = "start_date",
    returning = "date",
    date_format = "YYYY-MM-DD",
    find_last_match_in_period = True,
  ), 

  #  hospital admission with liver disease (see below?)

  ## CKD DEFINITIONS - adapted from https://github.com/opensafely/risk-factors-research
  #  recorded 3-5 CKD
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
        date_format = "YYYY-MM-DD",
        return_expectations={
            "rate": "universal",
            "category": {"ratios": {"1": 0.5, "2": 0.5}},
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

  #  recorded dialysis
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

  #  kidney transplant
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


  #  3-5 CKD based on recorded creatinine value
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

  ## Paxlovid - exclusion
  #  Solid organ transplant (plus solid_organ_transplant_nhsd_snomed defined below)
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

  ### Children aged less than 18 years
  #   (defined below)
  
  
  

  # HIGH RISK GROUPS ----
  
  ## NHSD ‘high risk’ cohort (codelist to be defined if/when data avaliable)
  # high_risk_cohort_nhsd = patients.with_these_clinical_events(
  #   high_risk_cohort_nhsd_codes,
  #   between = [campaign_start, index_date],
  #   returning = "date",
  #   date_format = "YYYY-MM-DD",
  #   find_first_match_in_period = True,
  # ),
  
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
  
  
  
  # CLINICAL/DEMOGRAPHIC COVARIATES ----
  
  ## Age
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
  
# Variant
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
# symptom
  symptomatic_covid_test=patients.with_test_result_in_sgss(
       pathogen="SARS-CoV-2",
       test_result="positive",
       find_first_match_in_period=True,
       restrict_to_earliest_specimen_date=False,
       between=["covid_test_positive_date","covid_test_positive_date + 30 days"],
       returning="symptomatic",
       return_expectations={
            "rate": "universal",
            "category": {"ratios": {"Y": 0.7, "N": 0.2, "": 0.1}},
       },
  ),     


  # OUTCOMES - extracted at covid_test_positive_date ----
  ## COVID-related hospitalisation 
  # extract multiple COVID hosp events per patient because the first hosp may be for receiving sotro or day cases (Day 0 and 1):
  covid_hosp_outcome_date0 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    between = ["covid_test_positive_date", "covid_test_positive_date"],
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
    between = ["covid_test_positive_date + 1 day", "covid_test_positive_date + 1 day"],
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-17"},
      "rate": "uniform",
      "incidence": 0.1
    },
  ),  
  covid_hosp_outcome_date2 = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "covid_test_positive_date + 2 days",
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
  covid_hosp_discharge_date2 = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "covid_hosp_outcome_date2",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-18"},
      "rate": "uniform",
      "incidence": 0.40
    },
  ),  

  ## COVID related death
  death_with_covid_date = patients.with_these_codes_on_death_certificate(
    covid_icd10_codes,
    returning = "date_of_death",
    date_format = "YYYY-MM-DD",
    on_or_after = "covid_test_positive_date",
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
    on_or_after = "covid_test_positive_date",
    match_only_underlying_cause=True,
    return_expectations = {
      "date": {"earliest": "2021-01-01", "latest" : "today"},
      "rate": "uniform",
      "incidence": 0.6},
  ),  
  # CENSORING ----
  ## Death of any cause
  death_date = patients.died_from_any_cause(
    returning = "date_of_death",
    date_format = "YYYY-MM-DD",
    on_or_after = "covid_test_positive_date",
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date"},
      "incidence": 0.1
    },
  ),
  
  has_died = patients.died_from_any_cause(
    on_or_before = "covid_test_positive_date - 1 day",
    returning = "binary_flag",
  ),
  
  ## De-registration
  dereg_date = patients.date_deregistered_from_all_supported_practices(
    on_or_after = "covid_test_positive_date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date"},
      "incidence": 0.1
    },
  ),



  # OUTCOMES - extracted at date_treated ----
  ## COVID-related hospitalisation 
  # extract multiple COVID hosp events per patient because the first hosp may be for receiving sotro or day cases (Day 0 and 1):
  covid_hosp_outcome_date0T = patients.admitted_to_hospital(
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
  covid_hosp_outcome_date1T = patients.admitted_to_hospital(
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
  covid_hosp_outcome_date2T = patients.admitted_to_hospital(
    returning = "date_admitted",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "start_date + 2 days",
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
  covid_hosp_discharge_date0T = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "covid_hosp_outcome_date0T",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-16"},
      "rate": "uniform",
      "incidence": 0.1
    },
  ),  
  covid_hosp_discharge_date1T = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "covid_hosp_outcome_date1T",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-17"},
      "rate": "uniform",
      "incidence": 0.1
    },
  ),  
  covid_hosp_discharge_date2T = patients.admitted_to_hospital(
    returning = "date_discharged",
    with_these_primary_diagnoses = covid_icd10_codes,
    with_patient_classification = ["1"], # ordinary admissions only - exclude day cases and regular attenders
    # see https://docs.opensafely.org/study-def-variables/#sus for more info
    # with_admission_method=["21", "22", "23", "24", "25", "2A", "2B", "2C", "2D", "28"], # emergency admissions only to exclude incidental COVID
    on_or_after = "covid_hosp_outcome_date2T",
    find_first_match_in_period = True,
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2022-02-18"},
      "rate": "uniform",
      "incidence": 0.40
    },
  ),  

  ## COVID related death
  death_with_covid_dateT = patients.with_these_codes_on_death_certificate(
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
  death_with_covid_underlying_dateT = patients.with_these_codes_on_death_certificate(
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
  # CENSORING ----
  ## Death of any cause
  death_dateT = patients.died_from_any_cause(
    returning = "date_of_death",
    date_format = "YYYY-MM-DD",
    on_or_after = "start_date",
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date"},
      "incidence": 0.1
    },
  ),
  
  has_diedT = patients.died_from_any_cause(
    on_or_before = "start_date - 1 day",
    returning = "binary_flag",
  ),
  
  ## De-registration
  dereg_dateT = patients.date_deregistered_from_all_supported_practices(
    on_or_after = "start_date",
    date_format = "YYYY-MM-DD",
    return_expectations = {
      "date": {"earliest": "2021-12-20", "latest": "index_date"},
      "incidence": 0.1
    },
  ),


#safety outcomes? 
)