* ===========================================================================
* 02_calculate_kru.do
* Objectif : Calcul du KRU (mL/min) et normalisation pour 35L (Watson)
* ===========================================================================

local path "C:\Users\dajs\OneDrive - HOPITAUX UNIVERSITAIRES DE GENEVE\recherche\RKF\UMOD\stata\main prospective study\with Claude"

use "`path'\database_umod.dta", clear
* database_umod.dta = fichier source (output de 01_merge) — jamais modifié
* database_analysis.dta = fichier d'analyse enrichi (output de ce do-file)

* ── 1. Durée de collecte précise ────────────────────────────────
* Parser urinestart et urineend (format "HH:MM")
gen start_h   = real(substr(urinestart, 1, 2))
gen start_min = real(substr(urinestart, 4, 2))
gen start_tot = start_h * 60 + start_min

gen end_h   = real(substr(urineend, 1, 2))
gen end_min = real(substr(urineend, 4, 2))
gen end_tot = end_h * 60 + end_min

* T (minutes) = toujours ~24h car récolte dans les 24h précédant la dialyse
gen T_min = 1440 + (end_tot - start_tot)

* Vérification : T doit être strictement positif
assert T_min > 0 if !missing(T_min)

drop start_h start_min start_tot end_h end_min end_tot

* ── 2. Calcul du KRU (mL/min) ───────────────────────────────────
gen kru = (urineurea * urinevolume) / (bloodurea * T_min)

label variable kru "KRU (mL/min)"

* ── 3. Contrôle qualité ─────────────────────────────────────────
* KRU doit être positif
count if kru < 0 & !missing(kru)
if r(N) > 0 {
    display as error "ATTENTION : " r(N) " valeurs de KRU négatives — à vérifier"
}

* Distribution
summarize kru, detail
histogram kru, normal title("Distribution du KRU (mL/min)") xtitle("KRU (mL/min)")

* Patients sans récolte urinaire (diuresis == 0 ou données manquantes)
count if missing(kru)
display "Patients avec KRU manquant : " r(N)
tab diuresis if missing(kru)

* Explorer les 2 patients avec diurèse mais sans KRU
list id urinevolume urineurea bloodurea urinestart urineend T_min ///
    if diuresis == 1 & missing(kru)

* Assigner KRU = 0 aux patients anuriques
replace kru = 0 if diuresis == 0
display "KRU=0 assigné aux patients anuriques"

* ── Bilan ───────────────────────────────────────────────────────
count if kru > 0 & !missing(kru)
display "KRU calculé (récolte urinaire) : " r(N)

count if kru == 0
display "KRU = 0 (anuriques)            : " r(N)

count if missing(kru)
display "KRU manquant                   : " r(N)

count if !missing(kru)
display "Total avec KRU disponible      : " r(N) "/" _N

* ── 4. Normalisation par Watson (mL/min/35L) ────────────────────
* Hommes (sex==2) : V = 2.447 - 0.09516×age + 0.1074×height + 0.3362×posthdweight
* Femmes (sex==1) : V = -2.097 + 0.1069×height + 0.2466×posthdweight
* height en cm, posthdweight en kg → V en litres

gen V_watson = .
replace V_watson = 2.447 - 0.09516*age + 0.1074*height + 0.3362*posthdweight if sex == 2
replace V_watson = -2.097 + 0.1069*height + 0.2466*posthdweight              if sex == 1
label variable V_watson "Volume de distribution urée - Watson (L)"

* Contrôle : V doit être physiologiquement plausible (5–70L)
count if (V_watson < 5 | V_watson > 70) & !missing(V_watson)
if r(N) > 0 display as error "ATTENTION : " r(N) " valeurs de V_watson hors plage [5-70L]"

* Identifier les patients avec V_watson manquant
count if missing(V_watson)
display "V_watson manquant : " r(N)
list id age sex height posthdweight if missing(V_watson)

summarize V_watson, detail

gen kru_35 = kru * (35 / V_watson)
label variable kru_35 "KRU (mL/min/35L)"

summarize kru_35, detail

* Les données restent en mémoire pour la suite de l'analyse
