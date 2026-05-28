* ===========================================================================
* 06_table1.do
* Prérequis : lancer 02_calculate_kru.do AVANT ce do-file.
*
* Objectif : Générer les statistiques pour la Table 1 du manuscrit.
*            Comparaison KRU < 2 vs KRU ≥ 2 mL/min/35L (population entière,
*            anuriques inclus avec KRU=0).
*
* Pour chaque variable continue :
*   - Shapiro-Wilk → si p<0.05 : médiane (IQR), Mann-Whitney
*                  → si p≥0.05 : moyenne ± SD, t-test
* Pour variables catégorielles : N (%), chi2 ou Fisher exact
* ===========================================================================

* --- Prérequis ---
foreach v in kru_daugirdas_35 umod labb2mprehd labcreatprehd age sex {
    capture confirm variable `v'
    if _rc {
        display as error "ERREUR : '`v'' absente — lance d'abord 02_calculate_kru.do"
        exit 111
    }
}

display _newline "=== Prérequis OK ==="

* --- Variables binaires ---
capture confirm variable kru_pos
if _rc {
    gen byte kru_pos = (kru_daugirdas_35 > 0) if !missing(kru_daugirdas_35)
    label variable kru_pos "KRU > 0"
}

capture confirm variable kru_ge2
if _rc {
    gen byte kru_ge2 = (kru_daugirdas_35 >= 2) if !missing(kru_daugirdas_35)
    label variable kru_ge2 "KRU >= 2 mL/min/35L"
}

* --- Variable female ---
capture confirm variable female
if _rc {
    gen byte female = (sex == 1) if !missing(sex)
    label variable female "Sexe féminin"
}

* ===========================================================================
* HEADER
* ===========================================================================
display _newline(2) "========================================================================"
display              "  TABLE 1 — Caractéristiques de la population"
display              "  Comparaison KRU < 2  vs  KRU ≥ 2  mL/min/35L"
display              "========================================================================"

* --- N par groupe ---
display _newline "  --- Effectifs ---"
tab kru_ge2, miss
count if missing(kru_ge2)
display "  KRU manquant : " r(N)

* ===========================================================================
* 1. DONNÉES DÉMOGRAPHIQUES
* ===========================================================================
display _newline(2) "────────────────────────────────────────────────────────────────────────"
display              "  1. DÉMOGRAPHIE"
display              "────────────────────────────────────────────────────────────────────────"

* --- Âge ---
display _newline "  AGE (ans)"
swilk age
tabstat age, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum age, by(kru_ge2)
display "  [Si normal : ttest age, by(kru_ge2)]"
ttest age, by(kru_ge2)

* --- Sexe ---
display _newline "  SEXE (% femmes)"
tab female kru_ge2, col chi2

* ===========================================================================
* 2. PARAMÈTRES DE DIALYSE
* ===========================================================================
display _newline(2) "────────────────────────────────────────────────────────────────────────"
display              "  2. PARAMÈTRES DE DIALYSE"
display              "────────────────────────────────────────────────────────────────────────"

* --- URR ---
display _newline "  URR (%)"
capture confirm variable URR
if _rc {
    display "  [URR non disponible]"
}
else {
    swilk URR
    tabstat URR, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
    ranksum URR, by(kru_ge2)
}

* --- Intervalle inter-dialytique ---
display _newline "  INTERVALLE INTER-DIALYTIQUE (jours)"
capture confirm variable interdialdays
if _rc {
    display "  [interdialdays non disponible]"
}
else {
    swilk interdialdays
    tabstat interdialdays, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
    ranksum interdialdays, by(kru_ge2)
}

* --- Poids post-HD ---
display _newline "  POIDS POST-HD (kg)"
capture confirm variable posthdweight
if _rc {
    display "  [posthdweight non disponible]"
}
else {
    swilk posthdweight
    tabstat posthdweight, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
    ranksum posthdweight, by(kru_ge2)
}

* --- Volume de distribution Watson ---
display _newline "  VOLUME DISTRIBUTION WATSON (L)"
capture confirm variable V_watson
if _rc {
    display "  [V_watson non disponible]"
}
else {
    swilk V_watson
    tabstat V_watson, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
    ranksum V_watson, by(kru_ge2)
}

* ===========================================================================
* 3. DIURÈSE ET KRU
* ===========================================================================
display _newline(2) "────────────────────────────────────────────────────────────────────────"
display              "  3. DIURÈSE ET KRU"
display              "────────────────────────────────────────────────────────────────────────"

* --- Statut anurique ---
display _newline "  STATUT ANURIQUE (KRU=0)"
tab kru_pos kru_ge2, col chi2

* --- Volume urinaire (si disponible) ---
display _newline "  VOLUME URINAIRE (mL/24h)"
capture confirm variable urinevolume
if _rc {
    display "  [urinevolume non disponible]"
}
else {
    swilk urinevolume if kru_pos == 1
    tabstat urinevolume if kru_pos == 1, by(kru_ge2) ///
        statistics(n mean sd p25 p50 p75) format(%7.0f)
    ranksum urinevolume if kru_pos == 1, by(kru_ge2)
    display "  (non-anuriques uniquement)"
}

* --- KRU Daugirdas/35L ---
display _newline "  KRU DAUGIRDAS/35L (mL/min)"
swilk kru_daugirdas_35
tabstat kru_daugirdas_35, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.2f)
display "  [Global :]"
summarize kru_daugirdas_35, detail

* ===========================================================================
* 4. BIOLOGIE
* ===========================================================================
display _newline(2) "────────────────────────────────────────────────────────────────────────"
display              "  4. BIOLOGIE"
display              "────────────────────────────────────────────────────────────────────────"

* --- Uromoduline sérique ---
display _newline "  UROMODULINE SÉRIQUE (ng/mL)"
swilk umod
summarize umod, detail
tabstat umod, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum umod, by(kru_ge2)
count if umod == 0
display "  UMOD = 0 : " r(N) " patients"

* --- Créatinine pré-HD ---
display _newline "  CRÉATININE PRÉ-HD (umol/L)"
swilk labcreatprehd
tabstat labcreatprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%7.0f)
ranksum labcreatprehd, by(kru_ge2)

* --- Urée pré-HD ---
display _newline "  URÉE PRÉ-HD (mmol/L)"
capture confirm variable labureaprehd
if _rc {
    display "  [labureaprehd non disponible]"
}
else {
    swilk labureaprehd
    tabstat labureaprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
    ranksum labureaprehd, by(kru_ge2)
}

* --- Beta-2-microglobuline ---
display _newline "  BETA-2-MICROGLOBULINE (mg/L)"
swilk labb2mprehd
summarize labb2mprehd, detail
tabstat labb2mprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum labb2mprehd, by(kru_ge2)

* ===========================================================================
* 5. VARIABLES SUPPLÉMENTAIRES À COMPLÉTER
*    Décommentez et adaptez selon les variables disponibles dans votre base.
* ===========================================================================
display _newline(2) "────────────────────────────────────────────────────────────────────────"
display              "  5. VARIABLES SUPPLÉMENTAIRES (adapter selon disponibilité)"
display              "────────────────────────────────────────────────────────────────────────"

* --- Durée de dialyse (vintage) ---
* capture confirm variable vintage
* if !_rc {
*     display _newline "  DURÉE DIALYSE (mois)"
*     swilk vintage
*     tabstat vintage, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
*     ranksum vintage, by(kru_ge2)
* }

* --- Albumine ---
* capture confirm variable labalbumin
* if !_rc {
*     display _newline "  ALBUMINE (g/L)"
*     swilk labalbumin
*     tabstat labalbumin, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%5.1f)
*     ranksum labalbumin, by(kru_ge2)
* }

* --- Hémoglobine ---
* capture confirm variable labhemoglobin
* if !_rc {
*     display _newline "  HÉMOGLOBINE (g/dL)"
*     swilk labhemoglobin
*     tabstat labhemoglobin, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%5.1f)
*     ranksum labhemoglobin, by(kru_ge2)
* }

* --- Type d'accès vasculaire ---
* capture confirm variable access_type
* if !_rc {
*     display _newline "  TYPE ACCÈS VASCULAIRE"
*     tab access_type kru_ge2, col chi2
* }

* --- Étiologie IRC ---
* capture confirm variable ckd_etiology
* if !_rc {
*     display _newline "  ÉTIOLOGIE IRC"
*     tab ckd_etiology kru_ge2, col chi2
* }

* ===========================================================================
* SYNTHÈSE FINALE
* ===========================================================================
display _newline(2) "========================================================================"
display              "  SYNTHÈSE — À vérifier avant de remplir la Table 1"
display              "========================================================================"
display "  1. Pour chaque variable continue : reporter médiane (IQR) si"
display "     Shapiro-Wilk p < 0.05, sinon moyenne ± SD."
display "  2. p-value : Mann-Whitney si non-normal, t-test si normal,"
display "     chi2 (ou Fisher si effectif < 5) pour catégorielles."
display "  3. Ajouter les variables de la section 5 selon disponibilité."
display "========================================================================"
