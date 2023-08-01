* Open a log file
cap log close
log using ./logs/plot, replace t
clear

use ./output/main.dta
*follow-up time and events*
stset end_date ,  origin(start_date) failure(failure==1)
keep if _st==1
tab _t,m
tab _t drug,m col
by drug, sort: sum _t ,de
tab _t drug if failure==1,m col
tab failure drug,m col
*K-M curve*
sort drug _t
by drug, sort: gen _dn=sum(_d)
gen _tn=1 if drug==0&_dn<=5
replace _tn=1 if drug==1&_dn<=5
replace _tn=4 if drug==0&_dn>5&_dn<=10
replace _tn=2 if drug==1&_dn>5&_dn<=10
replace _tn=8 if drug==0&_dn>10&_dn<=15
replace _tn=4 if drug==1&_dn>10&_dn<=15
replace _tn=12 if drug==0&_dn>15&_dn<=20
replace _tn=7 if drug==1&_dn>15&_dn<=20
replace _tn=28 if drug==0&_dn>20
replace _tn=14 if drug==1&_dn>20&_dn<=25
replace _tn=20 if drug==1&_dn>25&_dn<=30
replace _tn=28 if drug==1&_dn>30
replace _t=_tn
sts graph, by(drug) ylabel(.95(.01)1)
graph export ./output/kmcurve.svg, as(svg) replace

clear
use ./output/main_mol.dta
*follow-up time and events*
stset end_date ,  origin(start_date) failure(failure==1)
keep if _st==1
tab _t,m
tab _t drug,m col
by drug, sort: sum _t ,de
tab _t drug if failure==1,m col
tab failure drug,m col
*K-M curve*
sort drug _t
by drug, sort: gen _dn=sum(_d)
gen _tn=4 if drug==0&_dn<=5
replace _tn=1 if drug==1&_dn<=5
replace _tn=8 if drug==0&_dn>5&_dn<=10
replace _tn=2 if drug==1&_dn>5&_dn<=10
replace _tn=28 if drug==0&_dn>10
replace _tn=4 if drug==1&_dn>10&_dn<=15
replace _tn=7 if drug==1&_dn>15&_dn<=20
replace _tn=14 if drug==1&_dn>20&_dn<=25
replace _tn=20 if drug==1&_dn>25&_dn<=30
replace _tn=28 if drug==1&_dn>30
replace _t=_tn
sts graph, by(drug) ylabel(.95(.01)1)
graph export ./output/kmcurve_mol.svg, as(svg) replace

log close
