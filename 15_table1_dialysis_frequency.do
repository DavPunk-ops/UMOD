* ===========================================================================
* 15_table1_dialysis_frequency.do
* Prérequis : lancer 02_calculate_kru.do AVANT ce do-file.
*
* Objectif : produire la LIGNE « fréquence de dialyse » à ajouter à la Table 1
*   (section B — Paramètres de dialyse), demandée par R3 (minor m3).
*   N'ajoute rien à do-03 (originale, non modifiable) : ce do-file est la
*   source reproductible autonome de la nouvelle ligne.
*
* Variable : regimen (nombre de séances / semaine).
* Population : N=151 (kru_daugirdas_35 non manquant), comparaison KRU<2 vs ≥2.
*   → 3 valeurs de regimen manquantes attendues (ligne effective N=148).
*
* Deux présentations :
*   1. Catégoriel complet  (2× / 3× / 4× / 6×) → N(%) par groupe + Fisher exact
*   2. Collapsé  2×/sem vs ≥3×/sem            → N(%) par groupe + chi2 + Fisher
*      (plus lisible pour Table 1 ; robuste malgré les petites cellules)
* ===========================================================================

* --- Prérequis ---
foreach v in kru_daugirdas_35 regimen {
    capture confirm variable `v'
    if _rc {
        display as error "ERREUR : '`v'' absente — lance d'abord 02_calculate_kru.do"
        exit 111
    }
}
display _newline "=== Prérequis OK ==="

* --- Groupe KRU≥2 ---
capture confirm variable kru_ge2
if _rc {
    gen byte kru_ge2 = (kru_daugirdas_35 >= 2) if !missing(kru_daugirdas_35)
    label variable kru_ge2 "KRU >= 2 mL/min/35L"
}

* --- Restriction aux 151 patients analysés ---
keep if !missing(kru_ge2)
display _newline "  Restriction aux patients avec KRU calculé : N = " _N

* ===========================================================================
display _newline(2) "========================================================================"
display              "  TABLE 1 — nouvelle ligne : FRÉQUENCE DE DIALYSE (regimen)"
display              "  N=151 ; comparaison KRU < 2  vs  KRU ≥ 2 mL/min/35L"
display              "========================================================================"

* --- Codage de la variable (pour vérifier les libellés) ---
display _newline "  --- Codage de 'regimen' ---"
codebook regimen
tab regimen, missing

* --- Complétude ---
quietly count if !missing(regimen)
display _newline "  regimen non manquant : " r(N) " / " _N
quietly count if missing(regimen)
display "  regimen manquant     : " r(N)

* ===========================================================================
* 1. CATÉGORIEL COMPLET — N(%) par groupe + Fisher exact
* ===========================================================================
display _newline(2) "────────────────────────────────────────────────────────────────────────"
display              "  1. FRÉQUENCE (catégoriel complet) — colonnes = groupe KRU"
display              "────────────────────────────────────────────────────────────────────────"
tab regimen kru_ge2, col
display _newline "  --- Fisher exact (cellules < 5) ---"
tab regimen kru_ge2, exact

* Effectifs par groupe (dénominateurs)
display _newline "  --- Dénominateurs (regimen non manquant) par groupe ---"
tab kru_ge2 if !missing(regimen)

* ===========================================================================
* 2. COLLAPSÉ — 2×/sem  vs  ≥3×/sem
* ===========================================================================
display _newline(2) "────────────────────────────────────────────────────────────────────────"
display              "  2. FRÉQUENCE collapsée : 2×/sem  vs  ≥3×/sem"
display              "────────────────────────────────────────────────────────────────────────"

* ATTENTION : regimen est codé 1=3x/sem, 2=2x/sem, 3=4x/sem, 5=6x/sem
*   → le code numérique NE suit PAS la fréquence. On recode par label.
capture drop freq_ge3
gen byte freq_ge3 = .
replace freq_ge3 = 0 if regimen == 2                 // 2x/sem uniquement
replace freq_ge3 = 1 if inlist(regimen, 1, 3, 5)     // 3x, 4x, 6x /sem
label define freqlbl 0 "2x/sem" 1 ">=3x/sem", replace
label values freq_ge3 freqlbl
label variable freq_ge3 "Fréquence >=3x/sem"

display _newline "  --- Vérification du recodage (regimen -> freq_ge3) ---"
tab regimen freq_ge3, missing

display _newline "  --- N(%) par groupe KRU ---"
tab freq_ge3 kru_ge2, col chi2
display _newline "  --- Fisher exact ---"
tab freq_ge3 kru_ge2, exact

display _newline(2) "=== FIN do-file 15 ==="
