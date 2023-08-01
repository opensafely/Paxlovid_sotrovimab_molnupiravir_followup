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
sts graph, by(drug)
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
sts graph, by(drug)
graph export ./output/kmcurve_mol.svg, as(svg) replace

log close
