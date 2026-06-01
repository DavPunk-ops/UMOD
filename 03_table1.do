* ===========================================================================
* 03_table1.do
* Prérequis : lancer 02_calculate_kru.do AVANT ce do-file.
*
* Objectif : statistiques descriptives et comparaisons pour la Table 1.
*
* → manuscrit : "Description of the study cohort"
* → output    : Table 1
*
* Population : N=151 (restriction à kru_daugirdas_35 non manquant).
* Comparaison : KRU <2 vs KRU ≥2 mL/min/35L.
*
* Structure :
*   A. Paramètres cliniques (âge, sexe, IMC, BSA, Charlson, diabète, PA)
*   B. Paramètres de dialyse (HDF, UF, spKt/V, vintage, diurétiques)
*   C. Biomarqueurs (UMOD, β2M)
*
* Méthodes :
*   - Variables continues : Shapiro-Wilk → médiane [IQR] + Mann-Whitney
*                           si normal → moyenne ± SD + t-test
*   - Variables catégorielles : N (%) + chi2 ou Fisher exact
* ===========================================================================

* --- Prérequis ---
foreach v in kru_daugirdas_35 kru_daugirdas umod labb2mprehd age sex ///
             bmi bsa V_watson charlson spktv {
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

* --- Restriction aux 151 patients analysés (KRU calculable) ---
keep if !missing(kru_ge2)
display _newline "  Restriction aux patients avec KRU calculé : N = " _N

* ===========================================================================
* HEADER
* ===========================================================================
display _newline(2) "========================================================================"
display              "  TABLE 1 — Caractéristiques de la population (N=151)"
display              "  Comparaison  KRU < 2  vs  KRU ≥ 2  mL/min/35L"
display              "========================================================================"

display _newline "  --- Effectifs par groupe ---"
tab kru_ge2

* ===========================================================================
* A. PARAMÈTRES CLINIQUES
* ===========================================================================
display _newline(2) "────────────────────────────────────────────────────────────────────────"
display              "  A. PARAMÈTRES CLINIQUES"
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

* Poids
display _newline "  POIDS POST-HD (kg)"
swilk posthdweight
tabstat posthdweight, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum posthdweight, by(kru_ge2)
ttest posthdweight, by(kru_ge2)

* IMC
display _newline "  IMC (kg/m²)"
swilk bmi
tabstat bmi, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum bmi, by(kru_ge2)
ttest bmi, by(kru_ge2)

* Surface corporelle
display _newline "  SURFACE CORPORELLE Mosteller (m²)"
swilk bsa
tabstat bsa, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.2f)
ranksum bsa, by(kru_ge2)
ttest bsa, by(kru_ge2)

* Volume Watson
display _newline "  VOLUME DE DISTRIBUTION WATSON (L)"
swilk V_watson
tabstat V_watson, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum V_watson, by(kru_ge2)
ttest V_watson, by(kru_ge2)

* Ethnicité
display _newline "  RACE / ETHNICITÉ"
tab race kru_ge2, col chi2

* Étiologie IRC
display _newline "  ÉTIOLOGIE IRC"
tab kidneydisease kru_ge2, col chi2

* Score de Charlson
display _newline "  SCORE DE CHARLSON MODIFIÉ (avec ajustement âge)"
display "  (IDM, IC, AOMI, AVC/AIT, démence, BPCO, rhumato, UGD, DM±complic, ESRD, hémiplégie,"
display "   leucémie, lymphome, hépatopathie légère/sévère, cancer, métastases, AIDS)"
swilk charlson
tabstat charlson, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum charlson, by(kru_ge2)
ttest charlson, by(kru_ge2)

* Diabète
display _newline "  DIABÈTE"
tab dm kru_ge2, col chi2

* Pression artérielle
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

* KRU mL/min (distribution totale — variable de définition des groupes, pas de p)
display _newline "  KRU DAUGIRDAS (mL/min) — distribution globale (N=151)"
summarize kru_daugirdas, detail
swilk kru_daugirdas
tabstat kru_daugirdas, by(kru_ge2) statistics(n p25 p50 p75) format(%6.2f)

* KRU mL/min/35L (idem)
display _newline "  KRU DAUGIRDAS/35L (mL/min/35L) — distribution globale (N=151)"
summarize kru_daugirdas_35, detail
swilk kru_daugirdas_35
tabstat kru_daugirdas_35, by(kru_ge2) statistics(n p25 p50 p75) format(%6.2f)

* Volume urinaire (non-anuriques seulement)
display _newline "  VOLUME URINAIRE (mL/24h — non-anuriques uniquement)"
swilk urinevolume if kru_pos == 1
tabstat urinevolume if kru_pos == 1, by(kru_ge2) ///
    statistics(n mean sd p25 p50 p75) format(%7.0f)
ranksum urinevolume if kru_pos == 1, by(kru_ge2)

* % anuriques
display _newline "  STATUT ANURIQUE (KRU = 0)"
tab kru_pos kru_ge2, col chi2

* Diurétiques
display _newline "  UTILISATION DE DIURÉTIQUES"
capture confirm variable diuretic
if _rc {
    display as error "  ATTENTION : variable 'diuretic' absente — vérifier le merge 01"
}
else {
    tab diuretic kru_ge2, col chi2
}

* ===========================================================================
* B. PARAMÈTRES DE DIALYSE
* ===========================================================================
display _newline(2) "────────────────────────────────────────────────────────────────────────"
display              "  B. PARAMÈTRES DE DIALYSE"
display              "────────────────────────────────────────────────────────────────────────"

* Mode HD (HDF vs HD)
display _newline "  MODE HD (HDF vs HD)"
tab mode kru_ge2, col chi2

* Volume d'UF
display _newline "  VOLUME D'ULTRAFILTRATION (mL)"
capture confirm variable uf
if _rc {
    display as error "  ATTENTION : variable 'uf' absente"
}
else {
    swilk uf
    tabstat uf, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%7.0f)
    ranksum uf, by(kru_ge2)
    ttest uf, by(kru_ge2)
}

* spKt/V
display _newline "  spKt/V (Daugirdas)"
swilk spktv
tabstat spktv, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.2f)
ranksum spktv, by(kru_ge2)
ttest spktv, by(kru_ge2)

* Vintage
display _newline "  ANCIENNETÉ DE DIALYSE — vintage (mois)"
swilk vintage
tabstat vintage, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.0f)
ranksum vintage, by(kru_ge2)

* ===========================================================================
* C. BIOMARQUEURS
* ===========================================================================
display _newline(2) "────────────────────────────────────────────────────────────────────────"
display              "  C. BIOMARQUEURS"
display              "────────────────────────────────────────────────────────────────────────"

* UMOD
display _newline "  UROMODULINE SÉRIQUE (ng/mL)"
summarize umod, detail
swilk umod
tabstat umod, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.2f)
ranksum umod, by(kru_ge2)
count if umod == 0
display "  UMOD = 0 : " r(N) " patients (tous anuriques)"

* B2M
display _newline "  BETA-2-MICROGLOBULINE (mg/L)"
swilk labb2mprehd
tabstat labb2mprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum labb2mprehd, by(kru_ge2)
ttest labb2mprehd, by(kru_ge2)

* ===========================================================================
* SYNTHÈSE
* ===========================================================================
display _newline(2) "========================================================================"
display              "  SYNTHÈSE"
display              "========================================================================"
display "  Continue normale  (SW p≥0.05) → moyenne ± SD  + t-test"
display "  Continue skewed   (SW p<0.05) → médiane [IQR] + Mann-Whitney"
display "  Catégorielle                  → N (%)         + chi2 / Fisher"
display ""
display "  Note Charlson : ESRD (+2) constant pour tous. Toutes composantes REDCap incluses."
display "  Hépatopathie : légère +1, modérée/sévère +3 (non cumulatif). Cancer +2, métastases +6."
display "  Note KRU : variable de définition des groupes — pas de p-value reporté."
display "  Note Volume urinaire : parmi les 89 patients non-anuriques uniquement."
display "========================================================================"
