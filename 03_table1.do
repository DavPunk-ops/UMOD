* ===========================================================================
* 03_table1.do
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

* --- Dummy female ---
capture confirm variable female
if _rc {
    gen byte female = (sex == 1) if !missing(sex)
    label variable female "Sexe féminin"
}

* ===========================================================================
* HEADER
* ===========================================================================
display _newline(2) "========================================================================"
display              "  TABLE 1 — Caractéristiques de la population (N=151)"
display              "  Comparaison  KRU < 2  vs  KRU ≥ 2  mL/min/35L"
display              "========================================================================"

display _newline "  --- Effectifs par groupe ---"
tab kru_ge2, miss

* ===========================================================================
* 1. DÉMOGRAPHIE
* ===========================================================================
display _newline(2) "────────────────────────────────────────────────────────────────────────"
display              "  1. DÉMOGRAPHIE"
display              "────────────────────────────────────────────────────────────────────────"

* Age
display _newline "  AGE (ans)"
swilk age
tabstat age, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum age, by(kru_ge2)
ttest age, by(kru_ge2)

* Sexe
display _newline "  SEXE (% femmes)"
tab female kru_ge2, col chi2

* Race
display _newline "  RACE / ETHNICITÉ"
tab race kru_ge2, col chi2

* ===========================================================================
* 2. STATUT DIALYSE
* ===========================================================================
display _newline(2) "────────────────────────────────────────────────────────────────────────"
display              "  2. STATUT DIALYSE"
display              "────────────────────────────────────────────────────────────────────────"

* Incident vs prevalent
display _newline "  STATUT INCIDENT/PREVALENT (<3 mois vs ≥3 mois en HD)"
tab incident kru_ge2, col chi2

* Vintage (très skewed → médiane)
display _newline "  VINTAGE (mois en dialyse)"
swilk vintage
tabstat vintage, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.0f)
ranksum vintage, by(kru_ge2)

* Regimen
display _newline "  RÉGIME HD (séances/semaine)"
tab regimen kru_ge2, col chi2

* Mode HD (HDF vs HD)
display _newline "  MODE HD (HDF vs HD)"
tab mode kru_ge2, col chi2

* Accès vasculaire
display _newline "  ACCÈS VASCULAIRE (AVF vs cathéter)"
tab access kru_ge2, col chi2

* Durée de séance
display _newline "  DURÉE DE SÉANCE (min)"
swilk sessiontime
tabstat sessiontime, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.0f)
ranksum sessiontime, by(kru_ge2)
ttest sessiontime, by(kru_ge2)

* URR
display _newline "  UREA REDUCTION RATIO (%)"
swilk URR
tabstat URR, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum URR, by(kru_ge2)
ttest URR, by(kru_ge2)

* ===========================================================================
* 3. COMORBIDITÉS ET ÉTIOLOGIE
* ===========================================================================
display _newline(2) "────────────────────────────────────────────────────────────────────────"
display              "  3. COMORBIDITÉS ET ÉTIOLOGIE"
display              "────────────────────────────────────────────────────────────────────────"

* Étiologie IRC
display _newline "  ÉTIOLOGIE IRC"
tab kidneydisease kru_ge2, col chi2

* Diabète
display _newline "  DIABÈTE"
tab dm kru_ge2, col chi2

* Hypertension
display _newline "  HYPERTENSION"
tab ht kru_ge2, col chi2

* Transplantation rénale antérieure
display _newline "  TRANSPLANTATION RÉNALE ANTÉRIEURE"
display "  (attention: 48 missing)"
tab kidneytransplant kru_ge2, col chi2

* Comorbidités cardiovasculaires
display _newline "  IDM"
tab mi kru_ge2, col chi2
display _newline "  INSUFFISANCE CARDIAQUE"
tab chf kru_ge2, col chi2
display _newline "  ARTÉRIOPATHIE PÉRIPHÉRIQUE"
tab pvd kru_ge2, col chi2
display _newline "  MALADIE CÉRÉBROVASCULAIRE"
tab cva kru_ge2, col chi2

* ===========================================================================
* 4. PARAMÈTRES CLINIQUES
* ===========================================================================
display _newline(2) "────────────────────────────────────────────────────────────────────────"
display              "  4. PARAMÈTRES CLINIQUES"
display              "────────────────────────────────────────────────────────────────────────"

* Poids post-HD
display _newline "  POIDS POST-HD (kg)"
swilk posthdweight
tabstat posthdweight, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum posthdweight, by(kru_ge2)
ttest posthdweight, by(kru_ge2)

* TA pré-HD
display _newline "  PA SYSTOLIQUE PRÉ-HD (mmHg)"
swilk prehdsbp
tabstat prehdsbp, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.0f)
ranksum prehdsbp, by(kru_ge2)
ttest prehdsbp, by(kru_ge2)

display _newline "  PA DIASTOLIQUE PRÉ-HD (mmHg)"
swilk prehddbp
tabstat prehddbp, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.0f)
ranksum prehddbp, by(kru_ge2)
ttest prehddbp, by(kru_ge2)

* ===========================================================================
* 5. BIOLOGIE
* ===========================================================================
display _newline(2) "────────────────────────────────────────────────────────────────────────"
display              "  5. BIOLOGIE"
display              "────────────────────────────────────────────────────────────────────────"

* Hémoglobine
display _newline "  HÉMOGLOBINE (g/L)"
swilk labhbprehd
tabstat labhbprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.0f)
ranksum labhbprehd, by(kru_ge2)
ttest labhbprehd, by(kru_ge2)

* Sodium
display _newline "  SODIUM (mmol/L)"
swilk labnaprehd
tabstat labnaprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.0f)
ranksum labnaprehd, by(kru_ge2)

* Potassium
display _newline "  POTASSIUM (mmol/L)"
swilk labkprehd
tabstat labkprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum labkprehd, by(kru_ge2)

* Urée pré-HD
display _newline "  URÉE PRÉ-HD (mmol/L)"
swilk labureaprehd
tabstat labureaprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum labureaprehd, by(kru_ge2)
ttest labureaprehd, by(kru_ge2)

* Créatinine pré-HD
display _newline "  CRÉATININE PRÉ-HD (umol/L)"
swilk labcreatprehd
tabstat labcreatprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%7.0f)
ranksum labcreatprehd, by(kru_ge2)
ttest labcreatprehd, by(kru_ge2)

* Cystatine C
display _newline "  CYSTATINE C (mg/L)"
swilk labcyscprehd
tabstat labcyscprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.2f)
ranksum labcyscprehd, by(kru_ge2)

* Beta-2-microglobuline
display _newline "  BETA-2-MICROGLOBULINE (mg/L)"
swilk labb2mprehd
tabstat labb2mprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum labb2mprehd, by(kru_ge2)
ttest labb2mprehd, by(kru_ge2)

* Bicarbonate
display _newline "  BICARBONATE / CO2 TOTAL (mmol/L)"
swilk labco2prehd
tabstat labco2prehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum labco2prehd, by(kru_ge2)

* Calcium
display _newline "  CALCIUM TOTAL (mmol/L)"
swilk labcaprehd
tabstat labcaprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.2f)
ranksum labcaprehd, by(kru_ge2)

* Phosphate
display _newline "  PHOSPHATE (mmol/L)"
swilk labpoprehd
tabstat labpoprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.2f)
ranksum labpoprehd, by(kru_ge2)

* Albumine (49 missing — à mentionner)
display _newline "  ALBUMINE (g/L) — [49 valeurs manquantes]"
swilk labalbprehd
tabstat labalbprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum labalbprehd, by(kru_ge2)

* CRP (très skewed)
display _newline "  CRP (mg/L) — [distribution très asymétrique]"
swilk labcrpprehd
tabstat labcrpprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.0f)
ranksum labcrpprehd, by(kru_ge2)

* PTH (37 missing)
display _newline "  PTH (pg/mL) — [37 valeurs manquantes]"
swilk labpthprehd
tabstat labpthprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.0f)
ranksum labpthprehd, by(kru_ge2)

* ===========================================================================
* 6. DIURÈSE ET KRU
* ===========================================================================
display _newline(2) "────────────────────────────────────────────────────────────────────────"
display              "  6. DIURÈSE ET KRU"
display              "────────────────────────────────────────────────────────────────────────"

* Diurèse (oui/non)
display _newline "  DIURÈSE (KRU > 0 = non-anurique)"
tab kru_pos kru_ge2, col chi2

* Volume urinaire (non-anuriques seulement)
display _newline "  VOLUME URINAIRE (mL/24h — non-anuriques uniquement)"
swilk urinevolume if kru_pos == 1
tabstat urinevolume if kru_pos == 1, by(kru_ge2) ///
    statistics(n mean sd p25 p50 p75) format(%7.0f)
ranksum urinevolume if kru_pos == 1, by(kru_ge2)

* KRU Daugirdas/35L global
display _newline "  KRU DAUGIRDAS/35L (mL/min) — population entière"
summarize kru_daugirdas_35, detail
swilk kru_daugirdas_35

* ===========================================================================
* 7. UROMODULINE SÉRIQUE (biomarqueur principal)
* ===========================================================================
display _newline(2) "────────────────────────────────────────────────────────────────────────"
display              "  7. UROMODULINE SÉRIQUE"
display              "────────────────────────────────────────────────────────────────────────"

display _newline "  UMOD GLOBAL (ng/mL)"
summarize umod, detail
swilk umod

display _newline "  UMOD PAR GROUPE KRU"
tabstat umod, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.2f)
ranksum umod, by(kru_ge2)

count if umod == 0
display "  UMOD = 0 : " r(N) " patients (tous anuriques)"

* ===========================================================================
* SYNTHÈSE
* ===========================================================================
display _newline(2) "========================================================================"
display              "  SYNTHÈSE — Règles pour remplir la Table 1"
display              "========================================================================"
display "  Continue normale  → moyenne ± SD  + t-test"
display "  Continue skewed   → médiane [IQR] + Mann-Whitney"
display "  Catégorielle      → N (%)         + chi2 (ou Fisher si effectif<5)"
display ""
display "  Variables à EXCLURE ou MENTIONNER dans les limites :"
display "  - Albumine  : 49/153 missing (31.9%)"
display "  - PTH       : 37/153 missing (24.2%)"
display "  - Tabagisme : 37/153 missing"
display "  - Transplantation antérieure : 48/153 missing"
display "========================================================================"
