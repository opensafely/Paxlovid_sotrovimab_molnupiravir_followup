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
gen _tn=2 if drug==0&_dn<=6
replace _tn=1 if drug==1&_dn<=6
replace _tn=5 if drug==0&_dn>6&_dn<=12
replace _tn=3 if drug==1&_dn>6&_dn<=12
replace _tn=9 if drug==0&_dn>12&_dn<=18
replace _tn=5 if drug==1&_dn>12&_dn<=18
replace _tn=28 if drug==0&_dn>18
replace _tn=13 if drug==1&_dn>18&_dn<=24
replace _tn=20 if drug==1&_dn>24&_dn<=30
replace _tn=28 if drug==1&_dn>30
replace _t=_tn
sts graph, by(drug) ylabel(.98(.01)1) xtitle("Days since treatment initiation")
graph export ./output/kmcurve.svg, as(svg) replace
*monthly count*
gen month=month(start_date)
tab month drug, col


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
gen _tn=4 if drug==0&_dn<=6
replace _tn=1 if drug==1&_dn<=6
replace _tn=20 if drug==0&_dn>6&_dn<=12
replace _tn=3 if drug==1&_dn>6&_dn<=12
replace _tn=28 if drug==0&_dn>12
replace _tn=5 if drug==1&_dn>12&_dn<=18
replace _tn=13 if drug==1&_dn>18&_dn<=24
replace _tn=20 if drug==1&_dn>24&_dn<=30
replace _tn=28 if drug==1&_dn>30
replace _t=_tn
sts graph, by(drug) ylabel(.98(.01)1) xtitle("Days since treatment initiation")
graph export ./output/kmcurve_mol.svg, as(svg) replace
*monthly count*
gen month=month(start_date)
tab month drug, col


log close
