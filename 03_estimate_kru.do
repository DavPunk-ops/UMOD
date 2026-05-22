* ===========================================================================
* 03_estimate_kru.do
* Prérequis : lancer 02_calculate_kru.do AVANT ce do-file dans la même
*             session Stata.
* ===========================================================================

* --- Vérification que 02_calculate_kru.do a été exécuté ---
foreach v in kru_daugirdas_35 kru_naif_35 kru_pos umod diuresis {
    capture confirm variable `v'
    if _rc {
        display as error "ERREUR : variable '`v'' absente — lance d'abord 02_calculate_kru.do"
        exit 111
    }
}

display _newline "=== Prérequis OK — variables de 02_calculate_kru.do présentes ==="

* ===========================================================================
* SECTION 1 — DESCRIPTION DE L'UROMODULINE SÉRIQUE (UMOD)
* ===========================================================================

* --- 1a. Distribution globale ---
display _newline(2) "=== 1a. Distribution de l'UMOD (ng/mL) — population entière ==="
summarize umod, detail

swilk umod

count if umod == 0
display "  UMOD = 0 : " r(N) " patients"

histogram umod, normal ///
    title("Distribution de l'uromoduline sérique") ///
    xtitle("UMOD (ng/mL)") name(hist_umod, replace)

* --- 1b. UMOD selon statut de diurèse (KRU=0 vs KRU>0) ---
display _newline(2) "=== 1b. UMOD selon statut de diurèse ==="
display              "        (kru_pos : 0 = anurique, 1 = non-anurique)"

tabstat umod, by(kru_pos) statistics(n mean sd p25 p50 p75) format(%7.2f)

ranksum umod, by(kru_pos)

* --- 1c. Variable binaire KRU ≥ 2 mL/min/35L ---
capture drop kru_ge2
gen byte kru_ge2 = (kru_daugirdas_35 >= 2) if !missing(kru_daugirdas_35)
label variable kru_ge2 "KRU >= 2 mL/min/35L"
label define kruge2 0 "KRU < 2" 1 "KRU >= 2", replace
label values kru_ge2 kruge2

display _newline(2) "=== 1c. UMOD selon KRU < 2 vs KRU >= 2 mL/min/35L ==="
display              "        (population entière, N avec KRU disponible)"

tab kru_ge2, miss

tabstat umod, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%7.2f)

ranksum umod, by(kru_ge2)

display _newline(2) "=== 1c (suite). UMOD selon KRU < 2 vs KRU >= 2 — NON-ANURIQUES uniquement ==="

tabstat umod if kru_pos == 1, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%7.2f)

ranksum umod if kru_pos == 1, by(kru_ge2)
