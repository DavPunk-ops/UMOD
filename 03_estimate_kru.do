* ===========================================================================
* 03_estimate_kru.do
* Objectif : Modélisation de kru_daugirdas_35 (mL/min/35L) à partir
*            de la uromoduline sérique (umod, ng/mL)
*
* Prérequis : lancer 02_calculate_kru.do AVANT ce do-file dans la même
*             session Stata. Les variables suivantes doivent être présentes
*             en mémoire : umod, kru_daugirdas_35, kru_naif_35, diuresis.
*
* PLAN :
*   PARTIE A — EXPLORATION & PRÉPARATION              [sections 1-2]
*   PARTIE B — MODÉLISATION DE KRU CONTINU            [sections 3-5]
*              two-part model (logit + OLS)
*   PARTIE C — DÉCISION CLINIQUE AU SEUIL KRU ≥ 2     [sections 6-9]
*              logistique sur tous (6), puis sur non-anuriques (7),
*              ajout des biomarqueurs prédialyse (8) et modèle
*              parcimonieux UMOD + B2M (9)
* ===========================================================================

* Vérification que 02 a bien été exécuté
capture confirm variable kru_daugirdas_35
if _rc {
    display as error "ERREUR : kru_daugirdas_35 absent — lance d'abord 02_calculate_kru.do"
    exit 111
}
capture confirm variable umod
if _rc {
    display as error "ERREUR : umod absent — vérifie le merge dans 01_merge"
    exit 111
}

* ###########################################################################
* #                                                                         #
* #   PARTIE A — EXPLORATION & PRÉPARATION                                  #
* #                                                                         #
* ###########################################################################

* ===========================================================================
*  SECTION 1 — EXPLORATION DESCRIPTIVE : UMOD ET KRU
* ===========================================================================

* ─── 1a. Distribution de umod ───────────────────────────────────
display _newline "=== Distribution umod (ng/mL) ==="
summarize umod, detail

histogram umod, normal ///
    title("Distribution de l'uromoduline sérique") ///
    xtitle("UMOD (ng/mL)") name(hist_umod, replace)

* Test de normalité (Shapiro-Wilk, valide pour N < 2000)
swilk umod

* ─── 1b. Distribution de log(umod) ──────────────────────────────
gen log_umod = log(umod)
label variable log_umod "log(UMOD) [ln, ng/mL]"

display _newline "=== Distribution log(umod) ==="
summarize log_umod, detail

histogram log_umod, normal ///
    title("Distribution log(uromoduline)") ///
    xtitle("ln(UMOD)") name(hist_log_umod, replace)

swilk log_umod

* ─── 1c. Valeurs manquantes et extrêmes ─────────────────────────
count if missing(umod)
display "UMOD manquant : " r(N)

count if umod == 0
display "UMOD = 0     : " r(N)

* Valeurs potentiellement aberrantes (< 1 ou > 2000 ng/mL)
count if umod < 1 & !missing(umod)
display "UMOD < 1 ng/mL : " r(N)
count if umod > 2000 & !missing(umod)
display "UMOD > 2000 ng/mL : " r(N)

* ─── 1d. Relation umod ~ kru_daugirdas_35 ───────────────────────
display _newline "=== Corrélations umod / kru_daugirdas_35 ==="

* Sur tous les patients avec les deux valeurs disponibles
corr umod kru_daugirdas_35
corr log_umod kru_daugirdas_35

* Scatter : tous patients
twoway (scatter kru_daugirdas_35 umod, msize(small) jitter(2)) ///
       (lfit kru_daugirdas_35 umod), ///
    title("UMOD vs KRU Daugirdas/35L — tous patients") ///
    xtitle("UMOD (ng/mL)") ytitle("KRU Daugirdas (mL/min/35L)") ///
    name(scatter_all, replace)

* Scatter log(umod) — linéarise la relation si UMOD log-normal
twoway (scatter kru_daugirdas_35 log_umod, msize(small) jitter(2)) ///
       (lfit kru_daugirdas_35 log_umod), ///
    title("log(UMOD) vs KRU Daugirdas/35L — tous patients") ///
    xtitle("ln(UMOD)") ytitle("KRU Daugirdas (mL/min/35L)") ///
    name(scatter_log, replace)

* ─── 1e. Stratification anuriques / non-anuriques ───────────────
display _newline "=== UMOD selon statut de diurèse ==="
tabstat umod, by(diuresis) statistics(n mean sd p25 p50 p75) format(%7.1f)

* Scatter séparé non-anuriques (KRU > 0)
twoway (scatter kru_daugirdas_35 log_umod if kru_daugirdas_35 > 0, ///
        msize(small)) ///
       (lfit kru_daugirdas_35 log_umod if kru_daugirdas_35 > 0), ///
    title("log(UMOD) vs KRU Daugirdas/35L — non-anuriques") ///
    xtitle("ln(UMOD)") ytitle("KRU Daugirdas (mL/min/35L)") ///
    name(scatter_nonanuric, replace)

* ─── 1f. Résumé pour orienter la modélisation ───────────────────
display _newline "========================================"
display         "  BILAN EXPLORATOIRE"
display         "========================================"
quietly summarize umod
display "  UMOD    : N=" r(N) "  médiane=" %5.1f r(mean) " ng/mL (voir detail)"
quietly count if kru_daugirdas_35 == 0 & !missing(kru_daugirdas_35)
display "  Anuriques (KRU=0) : " r(N)
quietly count if kru_daugirdas_35 > 0 & !missing(kru_daugirdas_35)
display "  Non-anuriques     : " r(N)
display "========================================"

* ===========================================================================
*  SECTION 2 — PRÉPARATION POUR MODÉLISATION
*    - UMOD utilisé brut (pas d'imputation des zéros)
*      Rationale : le dataset contient des valeurs <LLOQ déclarée (jusqu'à
*      0.38 ng/mL), donc le labo rapporte sous LLOQ quand détectable. Les
*      UMOD=0 sont donc de VRAIS zéros (indétectables au seuil absolu de
*      l'assay), pas des valeurs censurées. Les conserver à 0 reflète
*      l'extinction biologique de la masse néphronique fonctionnelle.
*    - Création de la variable binaire kru_pos = 1 si KRU > 0
* ===========================================================================

* ─── 2a. Vérification UMOD brut ─────────────────────────────────
display _newline "=== Distribution UMOD (brut, sans imputation) ==="
tabstat umod, statistics(n min p1 p5 p25 p50 mean) format(%6.2f)
count if umod == 0
display "  → " r(N) " patients avec UMOD=0 (vrais zéros, conservés)"

* ─── 2b. Variable binaire kru_pos (anurique vs non-anurique) ────
gen byte kru_pos = (kru_daugirdas_35 > 0) if !missing(kru_daugirdas_35)
label variable kru_pos "KRU Daugirdas > 0 (1=non-anurique)"
label define krupos 0 "Anurique (KRU=0)" 1 "Non-anurique (KRU>0)"
label values kru_pos krupos
display _newline "=== Distribution kru_pos ==="
tab kru_pos

* ###########################################################################
* #                                                                         #
* #   PARTIE B — MODÉLISATION DE KRU CONTINU                                #
* #   Two-part model (sections 3-5) : logit P(KRU>0) + OLS KRU|KRU>0        #
* #                                                                         #
* ###########################################################################

* ===========================================================================
*  SECTION 3 — TWO-PART [PARTIE 1] : Logit P(KRU>0) ~ UMOD
*    Modélise la probabilité d'avoir une fonction rénale résiduelle.
* ===========================================================================
display _newline(2) "=============================================="
display              "  PARTIE 1 — Logit  P(KRU>0)  ~  UMOD"
display              "=============================================="

logit kru_pos umod
estimates store logit_umod

* Odds ratios + IC95
display _newline "=== Odds ratios ==="
logit kru_pos umod, or

* Probabilité prédite et AUC
predict p_pos, pr
label variable p_pos "P(KRU>0 | UMOD) — logit"

display _newline "=== Courbe ROC ==="
lroc, name(roc_part1, replace) title("ROC — P(KRU>0) ~ UMOD")

* Calibration : décile de risque vs proportion observée
display _newline "=== Calibration Hosmer-Lemeshow ==="
estat gof, group(10) table

* ===========================================================================
*  SECTION 4 — TWO-PART [PARTIE 2] : OLS KRU ~ UMOD | KRU > 0
*    Modélise l'intensité de KRU parmi les non-anuriques uniquement.
* ===========================================================================
display _newline(2) "=============================================="
display              "  PARTIE 2 — OLS  KRU  ~  UMOD  | KRU > 0"
display              "=============================================="

regress kru_daugirdas_35 umod if kru_pos == 1
estimates store ols_umod_nonanuric

* R², coefficients
display _newline "=== R² = " %5.3f e(r2) "  Adj R² = " %5.3f e(r2_a)

* Diagnostic résidus
predict yhat_p2 if e(sample), xb
predict resid_p2 if e(sample), resid

display _newline "=== Normalité des résidus (Shapiro-Wilk) ==="
swilk resid_p2

twoway (scatter resid_p2 yhat_p2, msize(small)) ///
       (lowess resid_p2 yhat_p2, lcolor(red)), ///
    yline(0, lpattern(dash)) ///
    title("Résidus vs prédits — non-anuriques") ///
    xtitle("KRU prédit (mL/min/35L)") ytitle("Résidus") ///
    legend(off) name(resid_fit_p2, replace)

* QQ-plot des résidus
qnorm resid_p2, title("QQ-plot résidus partie 2") name(qq_p2, replace)

* ===========================================================================
*  SECTION 5 — TWO-PART [COMBINAISON] : prédiction & performance
*    E[KRU | UMOD] = P(KRU>0 | UMOD) × E[KRU | KRU>0, UMOD]
*    Évaluation RMSE/MAE/corrélation + comparaison vs OLS naïf.
* ===========================================================================
display _newline(2) "=============================================="
display              "  PRÉDICTION TWO-PART  &  PERFORMANCE"
display              "=============================================="

* Prédiction OLS étendue à tous les patients
estimates restore ols_umod_nonanuric
predict kru_cond_all, xb
label variable kru_cond_all "E[KRU | KRU>0, UMOD] (OLS, tous patients)"

* Contraindre à >= 0 (pas de KRU négatif possible)
replace kru_cond_all = 0 if kru_cond_all < 0

* Prédiction finale two-part
gen kru_pred_2p = p_pos * kru_cond_all if !missing(p_pos, kru_cond_all)
label variable kru_pred_2p "KRU prédit (two-part)"

* Performance globale
display _newline "=== Performance two-part vs KRU observé ==="
corr kru_daugirdas_35 kru_pred_2p
spearman kru_daugirdas_35 kru_pred_2p

* RMSE et MAE
gen resid_2p = kru_daugirdas_35 - kru_pred_2p
gen resid_2p_sq = resid_2p^2
gen resid_2p_abs = abs(resid_2p)
quietly summarize resid_2p_sq
local rmse_2p = sqrt(r(mean))
quietly summarize resid_2p_abs
local mae_2p = r(mean)
display _newline "  RMSE two-part = " %5.2f `rmse_2p' " mL/min/35L"
display         "  MAE  two-part = " %5.2f `mae_2p' " mL/min/35L"

* Comparaison : OLS simple sur tous (KRU=0 inclus)
display _newline "=== Comparaison : OLS naïf sur tous (référence) ==="
regress kru_daugirdas_35 umod
estimates store ols_naif_all
predict kru_pred_naif, xb
gen resid_naif_sq = (kru_daugirdas_35 - kru_pred_naif)^2
quietly summarize resid_naif_sq
display "  RMSE OLS naïf = " %5.2f sqrt(r(mean)) " mL/min/35L"

* Table comparative
display _newline "=== Comparaison des modèles ==="
estimates table logit_umod ols_umod_nonanuric ols_naif_all, ///
    b(%6.3f) se(%6.3f) stats(N r2 r2_a aic)

* Graphiques observé vs prédit
twoway (scatter kru_daugirdas_35 kru_pred_2p, msize(small)) ///
       (function y=x, range(0 12) lcolor(red) lpattern(dash)), ///
    title("KRU observé vs prédit — two-part") ///
    xtitle("KRU prédit (mL/min/35L)") ///
    ytitle("KRU observé (mL/min/35L)") ///
    legend(off) name(obs_pred_2p, replace)

* Bilan final
display _newline "========================================"
display         "  BILAN MODÈLE TWO-PART"
display         "========================================"
display "  N total          : " _N
quietly count if !missing(kru_pred_2p)
display "  N avec prédiction: " r(N)
quietly corr kru_daugirdas_35 kru_pred_2p
display "  Corrélation Pearson : " %5.3f r(rho)
display "  RMSE two-part     : " %5.2f `rmse_2p' " mL/min/35L"
display "  MAE  two-part     : " %5.2f `mae_2p' " mL/min/35L"
display "========================================"

* ###########################################################################
* #                                                                         #
* #   PARTIE C — DÉCISION CLINIQUE AU SEUIL KRU ≥ 2 mL/min/35L              #
* #   Section 6 : logistique sur TOUS les patients (n=151)                  #
* #   Section 7 : logistique sur les NON-ANURIQUES uniquement (n=89)        #
* #               → analyse principale (vraie question clinique)            #
* #   Section 8 : ajout des biomarqueurs prédialyse (créat, urée, B2M)      #
* #   Section 9 : modèle parcimonieux UMOD + B2M                            #
* #                                                                         #
* ###########################################################################

* ===========================================================================
*  SECTION 6 — Logit P(KRU≥2) ~ UMOD sur la population entière
*    Justification clinique : à KRU ≥ 2 mL/min/35L, la dose de dialyse
*    doit être adaptée (réduction Kt/V cible). C'est un seuil de décision
*    thérapeutique. Cette section inclut les anuriques (KRU=0 par définition,
*    donc < 2), ce qui amplifie artificiellement l'AUC. Voir section 7 pour
*    l'analyse clinique pertinente.
* ===========================================================================
display _newline(2) "=============================================="
display              "  SECTION 6 — SEUIL CLINIQUE KRU ≥ 2"
display              "=============================================="

* ─── 6a. Variable binaire kru_ge2 ─────────────────
gen byte kru_ge2 = (kru_daugirdas_35 >= 2) if !missing(kru_daugirdas_35)
label variable kru_ge2 "KRU >= 2 mL/min/35L (seuil clinique)"
label define kruge2 0 "< 2 (dose standard)" 1 ">= 2 (adapter dose)"
label values kru_ge2 kruge2

display _newline "=== Distribution kru_ge2 ==="
tab kru_ge2
* Détail par statut anurique
display _newline "=== Croisement kru_pos × kru_ge2 ==="
tab kru_pos kru_ge2, row

* ─── 6b. UMOD selon le seuil ─────────────────
display _newline "=== UMOD par statut KRU≥2 ==="
tabstat umod, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.2f)

* ─── 6c. Régression logistique P(KRU≥2) ~ UMOD ─────────────────
display _newline "=============================================="
display         "  6c. Logit : P(KRU≥2) ~ UMOD"
display         "=============================================="

logit kru_ge2 umod
estimates store logit_ge2

* OR + IC95
display _newline "=== Odds ratios ==="
logit kru_ge2 umod, or

* Probabilité prédite
predict p_ge2, pr
label variable p_ge2 "P(KRU≥2 | UMOD) — logit"

* ROC / AUC
display _newline "=== Courbe ROC ==="
lroc, name(roc_ge2, replace) title("ROC — P(KRU≥2) ~ UMOD")

* Calibration
display _newline "=== Calibration Hosmer-Lemeshow ==="
estat gof, group(10) table

* ─── 6d. Recherche du seuil optimal d'UMOD ─────────────────
display _newline "=============================================="
display         "  6d. Seuil optimal d'UMOD (indice de Youden)"
display         "=============================================="

* Table sensibilité / spécificité pour chaque seuil de UMOD
display _newline "=== Table Se/Sp selon le seuil UMOD ==="
roctab kru_ge2 umod, detail

* Indice de Youden : maximise Se + Sp − 1
* On utilise senspec (Stata 16+) ou calcul manuel
capture which cutpt
if _rc {
    display "  (Pour seuil optimal automatique : ssc install cutpt)"
    display "  Calcul manuel ci-dessous"
}

* Calcul manuel du Youden sur grille UMOD
quietly summarize umod
local umod_max = r(max)

tempname Youden
matrix `Youden' = J(50, 4, .)
local i = 1
forvalues u = 1(1)50 {
    if `u' <= `umod_max' {
        quietly count if umod >= `u' & kru_ge2 == 1 & !missing(kru_ge2)
        local TP = r(N)
        quietly count if umod < `u' & kru_ge2 == 0 & !missing(kru_ge2)
        local TN = r(N)
        quietly count if umod >= `u' & kru_ge2 == 0 & !missing(kru_ge2)
        local FP = r(N)
        quietly count if umod < `u' & kru_ge2 == 1 & !missing(kru_ge2)
        local FN = r(N)
        local Se = `TP' / (`TP' + `FN')
        local Sp = `TN' / (`TN' + `FP')
        local J = `Se' + `Sp' - 1
        matrix `Youden'[`i', 1] = `u'
        matrix `Youden'[`i', 2] = `Se'
        matrix `Youden'[`i', 3] = `Sp'
        matrix `Youden'[`i', 4] = `J'
        local i = `i' + 1
    }
}
matrix colnames `Youden' = "Cutoff" "Sens" "Spec" "Youden"
matlist `Youden', format(%6.3f) title("Se/Sp/Youden par seuil UMOD (ng/mL)")

* ─── 6e. Performance à des seuils cliniquement utiles ──────────
display _newline "=============================================="
display         "  6e. Performance aux seuils UMOD candidats"
display         "=============================================="

* Tester 3 seuils : Youden optimal sera identifié visuellement,
* + 2 seuils orientés Se 90% et Sp 90% pour usage clinique
foreach cut in 5 8 12 15 {
    display _newline "  ─── Seuil UMOD ≥ `cut' ng/mL ───"
    quietly count if umod >= `cut' & kru_ge2 == 1 & !missing(kru_ge2)
    local TP = r(N)
    quietly count if umod < `cut' & kru_ge2 == 0 & !missing(kru_ge2)
    local TN = r(N)
    quietly count if umod >= `cut' & kru_ge2 == 0 & !missing(kru_ge2)
    local FP = r(N)
    quietly count if umod < `cut' & kru_ge2 == 1 & !missing(kru_ge2)
    local FN = r(N)
    local Se = `TP' / (`TP' + `FN')
    local Sp = `TN' / (`TN' + `FP')
    local VPP = `TP' / (`TP' + `FP')
    local VPN = `TN' / (`TN' + `FN')
    local LRp = `Se' / (1 - `Sp')
    local LRn = (1 - `Se') / `Sp'
    display "    Se  = " %5.1f 100*`Se'  " %    Sp  = " %5.1f 100*`Sp' " %"
    display "    VPP = " %5.1f 100*`VPP' " %    VPN = " %5.1f 100*`VPN' " %"
    display "    LR+ = " %5.2f `LRp'    "      LR- = " %5.2f `LRn'
    display "    TP=`TP'  FP=`FP'  TN=`TN'  FN=`FN'"
}

* ─── 6f. Graphique : P(KRU≥2) selon UMOD ──────────────────────
estimates restore logit_ge2
quietly summarize umod
twoway (line p_ge2 umod, sort lcolor(red) lwidth(medium)) ///
       (scatter kru_ge2 umod, msize(small) mcolor(%40) jitter(2)), ///
    yline(0.5, lpattern(dot) lcolor(gray)) ///
    title("Probabilité prédite d'avoir KRU ≥ 2 mL/min/35L") ///
    xtitle("UMOD (ng/mL)") ytitle("P(KRU≥2 | UMOD)") ///
    legend(order(1 "Logit" 2 "Observé") position(11) ring(0)) ///
    name(logit_ge2_curve, replace)

* ─── 6g. Comparaison avec la prédiction two-part seuillée ──────
display _newline "=============================================="
display         "  6g. Comparaison : logit direct vs two-part seuillé"
display         "=============================================="

* Prédiction binaire à partir du two-part (KRU prédit ≥ 2 ?)
gen byte pred_2p_ge2 = (kru_pred_2p >= 2) if !missing(kru_pred_2p)
label variable pred_2p_ge2 "Two-part prédiction ≥ 2"

* Prédiction binaire à partir du logit direct (P ≥ 0.5)
gen byte pred_logit_ge2 = (p_ge2 >= 0.5) if !missing(p_ge2)
label variable pred_logit_ge2 "Logit direct prédiction ≥ 2"

display _newline "=== Two-part seuillé à 2 ==="
tab kru_ge2 pred_2p_ge2, row
display _newline "=== Logit direct (P≥0.5) ==="
tab kru_ge2 pred_logit_ge2, row

* ─── 6h. Bilan ─────────────────
display _newline(2) "========================================"
display              "  BILAN SEUIL CLINIQUE KRU ≥ 2"
display              "========================================"
quietly count if kru_ge2 == 1 & !missing(kru_ge2)
local n_ge2 = r(N)
quietly count if kru_ge2 == 0 & !missing(kru_ge2)
local n_lt2 = r(N)
display "  N total avec KRU mesuré : " `n_ge2' + `n_lt2'
display "  KRU ≥ 2 (à adapter)     : " `n_ge2'
display "  KRU < 2 (dose standard) : " `n_lt2'
display "  → voir AUC et seuil UMOD optimal ci-dessus"
display "========================================"

* ===========================================================================
*  SECTION 7 — Logit P(KRU≥2) ~ UMOD chez les NON-ANURIQUES  [PRINCIPAL]
*    Vraie question clinique : chez les patients qui urinent (kru_pos=1),
*    UMOD permet-il de distinguer KRU<2 (dose standard) vs KRU≥2 (adapter) ?
*    Les anuriques sont exclus car par définition KRU=0 → pas besoin de
*    doser UMOD pour la décision clinique chez eux.
* ===========================================================================
display _newline(2) "=============================================="
display              "  SECTION 7 — KRU ≥ 2 CHEZ LES NON-ANURIQUES"
display              "=============================================="

* ─── 7a. Description de la population concernée ─────────────
display _newline "=== Distribution kru_ge2 chez les non-anuriques ==="
tab kru_ge2 if kru_pos == 1
display _newline "=== UMOD selon KRU<2 vs KRU>=2 chez non-anuriques ==="
tabstat umod if kru_pos == 1, by(kru_ge2) ///
    statistics(n mean sd p25 p50 p75) format(%6.2f)

* Test non paramétrique (Mann-Whitney)
display _newline "=== Test Mann-Whitney sur UMOD ==="
ranksum umod if kru_pos == 1, by(kru_ge2)

* ─── 7b. Régression logistique restreinte aux non-anuriques ──
display _newline "=============================================="
display         "  7b. Logit : P(KRU≥2) ~ UMOD | non-anurique"
display         "=============================================="

logit kru_ge2 umod if kru_pos == 1
estimates store logit_ge2_na

* OR + IC95
display _newline "=== Odds ratios ==="
logit kru_ge2 umod if kru_pos == 1, or

* Probabilité prédite (uniquement sur les non-anuriques)
capture drop p_ge2_na
predict p_ge2_na if kru_pos == 1, pr
label variable p_ge2_na "P(KRU≥2 | UMOD, non-anurique)"

* ROC / AUC
display _newline "=== Courbe ROC ==="
lroc, name(roc_ge2_na, replace) ///
    title("ROC — P(KRU≥2) ~ UMOD chez non-anuriques")

* Calibration
display _newline "=== Calibration Hosmer-Lemeshow ==="
estat gof, group(10) table

* ─── 7c. Recherche du seuil UMOD optimal (Youden) ─────────────
display _newline "=============================================="
display         "  7c. Seuil optimal UMOD (Youden, non-anuriques)"
display         "=============================================="

display _newline "=== Table Se/Sp détaillée (roctab) ==="
roctab kru_ge2 umod if kru_pos == 1, detail

* Calcul manuel du Youden par grille (entiers 1-45)
quietly summarize umod if kru_pos == 1
local umod_max_na = r(max)

tempname YoudenNA
matrix `YoudenNA' = J(45, 4, .)
local i = 1
forvalues u = 1(1)45 {
    if `u' <= `umod_max_na' {
        quietly count if umod >= `u' & kru_ge2 == 1 & kru_pos == 1
        local TP = r(N)
        quietly count if umod < `u' & kru_ge2 == 0 & kru_pos == 1
        local TN = r(N)
        quietly count if umod >= `u' & kru_ge2 == 0 & kru_pos == 1
        local FP = r(N)
        quietly count if umod < `u' & kru_ge2 == 1 & kru_pos == 1
        local FN = r(N)
        if (`TP' + `FN') > 0 & (`TN' + `FP') > 0 {
            local Se = `TP' / (`TP' + `FN')
            local Sp = `TN' / (`TN' + `FP')
            local J = `Se' + `Sp' - 1
            matrix `YoudenNA'[`i', 1] = `u'
            matrix `YoudenNA'[`i', 2] = `Se'
            matrix `YoudenNA'[`i', 3] = `Sp'
            matrix `YoudenNA'[`i', 4] = `J'
        }
        local i = `i' + 1
    }
}
matrix colnames `YoudenNA' = "Cutoff" "Sens" "Spec" "Youden"
matlist `YoudenNA', format(%6.3f) ///
    title("Se/Sp/Youden par seuil UMOD chez non-anuriques")

* ─── 7d. Performance aux seuils candidats ─────────────
display _newline "=============================================="
display         "  7d. Performance aux seuils UMOD candidats"
display         "  (population : non-anuriques uniquement)"
display         "=============================================="

foreach cut in 5 8 10 12 15 {
    display _newline "  ─── Seuil UMOD ≥ `cut' ng/mL ───"
    quietly count if umod >= `cut' & kru_ge2 == 1 & kru_pos == 1
    local TP = r(N)
    quietly count if umod < `cut' & kru_ge2 == 0 & kru_pos == 1
    local TN = r(N)
    quietly count if umod >= `cut' & kru_ge2 == 0 & kru_pos == 1
    local FP = r(N)
    quietly count if umod < `cut' & kru_ge2 == 1 & kru_pos == 1
    local FN = r(N)
    if (`TP' + `FN') > 0 & (`TN' + `FP') > 0 {
        local Se = `TP' / (`TP' + `FN')
        local Sp = `TN' / (`TN' + `FP')
        local VPP = `TP' / (`TP' + `FP')
        local VPN = `TN' / (`TN' + `FN')
        local LRp = `Se' / (1 - `Sp')
        local LRn = (1 - `Se') / `Sp'
        display "    Se  = " %5.1f 100*`Se'  " %    Sp  = " %5.1f 100*`Sp' " %"
        display "    VPP = " %5.1f 100*`VPP' " %    VPN = " %5.1f 100*`VPN' " %"
        display "    LR+ = " %5.2f `LRp'    "      LR- = " %5.2f `LRn'
        display "    TP=`TP'  FP=`FP'  TN=`TN'  FN=`FN'"
    }
}

* ─── 7e. Graphique : P(KRU≥2 | non-anurique) selon UMOD ────────
estimates restore logit_ge2_na
twoway (line p_ge2_na umod if kru_pos == 1, sort lcolor(red) lwidth(medium)) ///
       (scatter kru_ge2 umod if kru_pos == 1, msize(small) ///
        mcolor(%40) jitter(2)), ///
    yline(0.5, lpattern(dot) lcolor(gray)) ///
    title("P(KRU≥2) chez les non-anuriques") ///
    xtitle("UMOD (ng/mL)") ytitle("P(KRU≥2 | UMOD, non-anurique)") ///
    legend(order(1 "Logit" 2 "Observé") position(11) ring(0)) ///
    name(logit_ge2_na_curve, replace)

* ─── 7f. Comparaison avec la section 6 (population entière) ────
display _newline "=============================================="
display         "  7f. Comparaison restreint vs population totale"
display         "=============================================="
display "  Population entière (section 6) :"
quietly logit kru_ge2 umod
quietly lroc, nograph
display "     AUC = " %5.3f r(area)
display "  Non-anuriques uniquement (section 7) :"
quietly logit kru_ge2 umod if kru_pos == 1
quietly lroc, nograph
display "     AUC = " %5.3f r(area)
display _newline "  → Cette section répond à la VRAIE question clinique :"
display "    chez un patient qui urine, UMOD prédit-il KRU≥2 ?"

* ─── 7g. Bilan ─────────────────
display _newline(2) "========================================"
display              "  BILAN — DÉCISION CLINIQUE CHEZ NON-ANURIQUES"
display              "========================================"
quietly count if kru_pos == 1 & !missing(kru_ge2)
display "  N non-anuriques avec KRU mesuré : " r(N)
quietly count if kru_pos == 1 & kru_ge2 == 1
display "  KRU ≥ 2 (à adapter)             : " r(N)
quietly count if kru_pos == 1 & kru_ge2 == 0
display "  0 < KRU < 2 (dose standard)     : " r(N)
display "========================================"

* ===========================================================================
*  SECTION 8 — UMOD + BIOMARQUEURS PRÉDIALYSE
*    Biomarqueurs explorés : créatinine, urée, β2-microglobuline prédialyse.
*    Rationale physiologique : tous ces marqueurs s'accumulent quand la
*    clairance rénale résiduelle baisse (creat & urée filtration ; β2M
*    extrêmement peu clairée par l'HD, donc reflet quasi pur de la RKF).
*    Hypothèse directionnelle : OR < 1 attendu (marqueur ↑ → KRU<2).
*
*    Stratégie :
*      8a. Description des 3 biomarqueurs (valeurs manquantes, par groupe)
*      8b. Modèle principal : UMOD + créat + urée + B2M | non-anuriques
*      8c. LR tests d'apport (vs UMOD seul, et par biomarqueur)
*      8d. Sensibilité : entraîné sur tous, testé sur NA
*      8e. Comparaison AUC tous modèles
*      8f. Seuils Youden sur le meilleur score
*      8g. Bilan
* ===========================================================================
display _newline(2) "=============================================="
display              "  SECTION 8 — UMOD + BIOMARQUEURS PRÉDIALYSE"
display              "=============================================="

* Vérification présence des variables
foreach v in labcreatprehd labureaprehd labb2mprehd {
    capture confirm variable `v'
    if _rc {
        display as error "ERREUR : variable '`v'' absente"
        exit 111
    }
}

* ─── 8a. Description des biomarqueurs ──────────────────────────
display _newline "=== Distribution des biomarqueurs (échantillon entier) ==="
tabstat labcreatprehd labureaprehd labb2mprehd, ///
    statistics(n mean sd min p25 p50 p75 max) format(%7.1f)

display _newline "=== Manquants chez les non-anuriques ==="
foreach v in labcreatprehd labureaprehd labb2mprehd {
    quietly count if missing(`v') & kru_pos == 1
    display "  `v' manquant chez NA : " r(N) "/89"
}

display _newline "=== Biomarqueurs selon KRU<2 vs KRU≥2 (non-anuriques) ==="
foreach v in labcreatprehd labureaprehd labb2mprehd {
    display _newline "--- `v' ---"
    tabstat `v' if kru_pos == 1, by(kru_ge2) ///
        statistics(n mean sd p25 p50 p75) format(%6.1f)
    ranksum `v' if kru_pos == 1, by(kru_ge2)
}

* Corrélations avec KRU continu
display _newline "=== Corrélations biomarqueurs / KRU continu (non-anuriques) ==="
foreach v in labcreatprehd labureaprehd labb2mprehd {
    quietly corr kru_daugirdas_35 `v' if kru_pos == 1
    local r_pearson = r(rho)
    quietly spearman kru_daugirdas_35 `v' if kru_pos == 1
    local r_spear = r(rho)
    display "  `v' : Pearson r = " %6.3f `r_pearson' "    Spearman = " %6.3f `r_spear'
}

* Corrélations entre biomarqueurs (collinéarité)
display _newline "=== Corrélations entre biomarqueurs ==="
corr labcreatprehd labureaprehd labb2mprehd umod if kru_pos == 1

* ─── 8b. Modèle principal : UMOD + 3 biomarqueurs | NA ────────
display _newline "=============================================="
display         "  8b. Logit : UMOD + créat + urée + B2M | NA"
display         "        (PRINCIPAL)"
display         "=============================================="

logit kru_ge2 umod labcreatprehd labureaprehd labb2mprehd if kru_pos == 1
estimates store logit_bio_na

display _newline "=== Odds ratios ==="
logit kru_ge2 umod labcreatprehd labureaprehd labb2mprehd if kru_pos == 1, or

capture drop p_bio_na
predict p_bio_na if kru_pos == 1, pr
label variable p_bio_na "P(KRU≥2) UMOD + biomarqueurs | NA"

display _newline "=== AUC ==="
lroc, name(roc_bio_na, replace) ///
    title("ROC UMOD + biomarqueurs — non-anuriques")

display _newline "=== Calibration Hosmer-Lemeshow ==="
estat gof, group(10) table

* ─── 8c. LR tests d'apport ──────────────────────────────────
display _newline "=============================================="
display         "  8c. Apport des biomarqueurs (LR tests)"
display         "=============================================="

estimates restore logit_bio_na
gen byte _sbio = e(sample)

quietly logit kru_ge2 umod if _sbio == 1
estimates store logit_u_only_sb

display _newline "=== LR test : UMOD seul vs UMOD + 3 biomarqueurs ==="
lrtest logit_u_only_sb logit_bio_na

display _newline "=== Apport individuel : UMOD + chaque biomarqueur seul ==="
foreach v in labcreatprehd labureaprehd labb2mprehd {
    quietly logit kru_ge2 umod `v' if _sbio == 1
    estimates store logit_u_`v'_sb
    display _newline "--- UMOD + `v' ---"
    quietly lroc, nograph
    display "    AUC = " %5.3f r(area)
    lrtest logit_u_only_sb logit_u_`v'_sb
}

display _newline "=== AIC/BIC : tous les modèles emboîtés ==="
estimates stats logit_u_only_sb logit_u_labcreatprehd_sb ///
    logit_u_labureaprehd_sb logit_u_labb2mprehd_sb logit_bio_na

drop _sbio

* ─── 8d. SENSIBILITÉ : entraîné sur tous, testé sur NA ───────
display _newline "=============================================="
display         "  8d. Sensibilité : entraîné sur tous"
display         "=============================================="

logit kru_ge2 umod labcreatprehd labureaprehd labb2mprehd
estimates store logit_bio_all

display _newline "=== Odds ratios (modèle sur tous) ==="
logit kru_ge2 umod labcreatprehd labureaprehd labb2mprehd, or

capture drop p_bio_all
predict p_bio_all, pr
label variable p_bio_all "P(KRU≥2) UMOD+biomark | tous"

display _newline "=== AUC sur tous (entraînement) ==="
lroc, nograph
display "    AUC sur tous : " %5.3f r(area)

display _newline "=== AUC restreint aux non-anuriques (test) ==="
roctab kru_ge2 p_bio_all if kru_pos == 1

* ─── 8e. Comparaison AUC sur les NA ──────────────────────────
display _newline "=============================================="
display         "  8e. Comparaison AUC | non-anuriques"
display         "=============================================="

display _newline "  ① UMOD seul (section 7) :"
quietly roctab kru_ge2 umod if kru_pos == 1
display "      AUC = " %5.3f r(area)

display "  ② UMOD + 3 biomarqueurs (8b) :"
quietly roctab kru_ge2 p_bio_na if kru_pos == 1
display "      AUC = " %5.3f r(area)

display "  ③ UMOD + 3 biomarqueurs entraîné sur tous (8d) :"
quietly roctab kru_ge2 p_bio_all if kru_pos == 1
display "      AUC = " %5.3f r(area)

display _newline "=== Test DeLong : UMOD seul vs UMOD+biomark ==="
roccomp kru_ge2 umod p_bio_na if kru_pos == 1, graph summary ///
    name(roccomp_bio, replace)

* ─── 8f. Seuils Youden sur le modèle UMOD + biomarqueurs ─────
display _newline "=============================================="
display         "  8f. Seuil Youden sur P(KRU≥2) - biomarqueurs"
display         "=============================================="

estimates restore logit_bio_na

foreach cut in 0.3 0.4 0.5 0.6 0.7 0.8 {
    display _newline "  ─── Seuil P(KRU≥2) ≥ `cut' ───"
    quietly count if p_bio_na >= `cut' & kru_ge2 == 1 & kru_pos == 1
    local TP = r(N)
    quietly count if p_bio_na < `cut' & kru_ge2 == 0 & kru_pos == 1
    local TN = r(N)
    quietly count if p_bio_na >= `cut' & kru_ge2 == 0 & kru_pos == 1
    local FP = r(N)
    quietly count if p_bio_na < `cut' & kru_ge2 == 1 & kru_pos == 1
    local FN = r(N)
    if (`TP' + `FN') > 0 & (`TN' + `FP') > 0 {
        local Se = `TP' / (`TP' + `FN')
        local Sp = `TN' / (`TN' + `FP')
        local VPP = `TP' / (`TP' + `FP')
        local VPN = `TN' / (`TN' + `FN')
        local J = `Se' + `Sp' - 1
        display "    Se = " %5.1f 100*`Se' " %    Sp = " %5.1f 100*`Sp' " %    J = " %5.3f `J'
        display "    VPP = " %5.1f 100*`VPP' " %    VPN = " %5.1f 100*`VPN' " %"
    }
}

* ─── 8g. Bilan ──────────────────────────────────────────────
display _newline(2) "========================================"
display              "  BILAN — UMOD + BIOMARQUEURS"
display              "========================================"
display "  Hypothèse : créat, urée, B2M ↑ → KRU<2 plus probable (OR<1)"
display "  Voir AUC, LR tests, OR et seuils Youden ci-dessus"
display "========================================"

* ===========================================================================
*  SECTION 9 — MODÈLE PARCIMONIEUX : UMOD + β2-MICROGLOBULINE
*
*    Rationale :
*      - La section 8 montre que B2M est le seul biomarqueur réellement
*        utile (LR p=0.0008, AUC +0.062 vs UMOD seul).
*      - La créatinine n'apporte rien (LR p=0.49).
*      - L'urée présente une circularité mathématique : labureaprehd ≈
*        bloodurea, variable au dénominateur du KRU Daugirdas (TAC).
*        Son coefficient inversé (OR=1.13, KRU≥2 → urée ↑) confirme
*        la suppression par collinéarité.
*      - Le modèle UMOD + B2M est biologiquement défendable (B2M =
*        reflet quasi-pur de la RKF, non clairé par HD), parcimonieux,
*        et exempt de circularité.
*
*    Stratégie :
*      9a. Modèle principal : UMOD + B2M | non-anuriques (N≈87)
*      9b. LR test d'apport de B2M vs UMOD seul
*      9c. Sensibilité : entraîné sur tous, testé sur NA
*      9d. Comparaison AUC — tableau final complet
*      9e. Seuils Youden sur le score UMOD + B2M
*      9f. Bilan
* ===========================================================================
display _newline(2) "=============================================="
display              "  SECTION 9 — MODÈLE PARCIMONIEUX : UMOD + B2M"
display              "=============================================="

* Vérification présence des variables
foreach v in umod labb2mprehd {
    capture confirm variable `v'
    if _rc {
        display as error "ERREUR : variable '`v'' absente"
        exit 111
    }
}

* ─── 9a. Modèle principal UMOD + B2M | non-anuriques ──────────
display _newline "=============================================="
display         "  9a. Modèle UMOD + B2M (non-anuriques)"
display         "=============================================="

display _newline "=== Manquants pour B2M chez les non-anuriques ==="
quietly count if missing(labb2mprehd) & kru_pos == 1
display "  labb2mprehd manquant chez NA : " r(N) "/89"
quietly count if !missing(labb2mprehd) & kru_pos == 1
display "  N analysable (UMOD + B2M, non-anuriques) : " r(N)

logit kru_ge2 umod labb2mprehd if kru_pos == 1
estimates store logit_b2m_na

display _newline "=== Odds ratios ==="
logit kru_ge2 umod labb2mprehd if kru_pos == 1, or

capture drop p_b2m_na
predict p_b2m_na if kru_pos == 1 & e(sample), pr
label variable p_b2m_na "P(KRU≥2) UMOD+B2M | NA"

display _newline "=== AUC (UMOD + B2M | non-anuriques) ==="
lroc, nograph
display "    AUC = " %5.3f r(area)

display _newline "=== Calibration Hosmer-Lemeshow (UMOD + B2M) ==="
estat gof, group(10) table

* ─── 9b. LR test : apport de B2M par-dessus UMOD seul ─────────
display _newline "=============================================="
display         "  9b. LR test : B2M au-delà de UMOD seul"
display         "=============================================="

display _newline "=== Modèle nul : UMOD seul (sur même N que 9a) ==="
quietly logit kru_ge2 umod if kru_pos == 1 & !missing(labb2mprehd)
estimates store logit_umod_b2m_n
display "    N = " e(N) "  (même sous-échantillon avec B2M non-manquant)"

display _newline "=== LR test : UMOD seul vs UMOD + B2M ==="
lrtest logit_umod_b2m_n logit_b2m_na

display _newline "=== Comparaison AIC/BIC ==="
estimates stats logit_umod_b2m_n logit_b2m_na

* ─── 9c. Sensibilité : entraîné sur tous, testé sur NA ─────────
display _newline "=============================================="
display         "  9c. Sensibilité : entraîné sur tous"
display         "=============================================="

logit kru_ge2 umod labb2mprehd
estimates store logit_b2m_all

display _newline "=== Odds ratios (modèle entraîné sur tous) ==="
logit kru_ge2 umod labb2mprehd, or

capture drop p_b2m_all
predict p_b2m_all, pr
label variable p_b2m_all "P(KRU≥2) UMOD+B2M | tous"

display _newline "=== AUC sur tous (entraînement) ==="
lroc, nograph
display "    AUC (entraînement) = " %5.3f r(area)

display _newline "=== AUC restreint aux non-anuriques (test) ==="
roctab kru_ge2 p_b2m_all if kru_pos == 1

* ─── 9d. Comparaison AUC — tableau final ───────────────────────
display _newline "=============================================="
display         "  9d. Tableau AUC — tous modèles | NA"
display         "=============================================="

display _newline "  ① UMOD seul (section 7) :"
quietly roctab kru_ge2 umod if kru_pos == 1
display "      AUC = " %5.3f r(area)

display "  ② UMOD + 3 biomarqueurs (8b) :"
quietly roctab kru_ge2 p_bio_na if kru_pos == 1
display "      AUC = " %5.3f r(area)

display "  ③ UMOD + 3 biomarqueurs entraîné sur tous (8d) :"
quietly roctab kru_ge2 p_bio_all if kru_pos == 1
display "      AUC = " %5.3f r(area)

display "  ④ UMOD + B2M | NA (9a) :"
quietly roctab kru_ge2 p_b2m_na if kru_pos == 1
display "      AUC = " %5.3f r(area)

display "  ⑤ UMOD + B2M entraîné sur tous (9c) :"
quietly roctab kru_ge2 p_b2m_all if kru_pos == 1
display "      AUC = " %5.3f r(area)

display _newline "=== Test DeLong : UMOD seul vs UMOD+B2M (modèle NA) ==="
roccomp kru_ge2 umod p_b2m_na if kru_pos == 1, graph summary ///
    name(roccomp_b2m, replace)

display _newline "=== Test DeLong : UMOD+B2M vs UMOD+3 biomarqueurs ==="
capture roccomp kru_ge2 p_b2m_na p_bio_na if kru_pos == 1, summary

* ─── 9e. Seuils Youden sur le score UMOD + B2M ─────────────────
display _newline "=============================================="
display         "  9e. Seuils Youden sur P(KRU≥2) — UMOD + B2M"
display         "=============================================="

display "  (Seuils de probabilité prédite ; non-anuriques uniquement)"
display _newline "  Cut  |  Se(%)  Sp(%)  VPP(%)  VPN(%)  J"
display          "  -----|--------------------------------------------"

foreach cut in 0.3 0.4 0.5 0.6 0.7 0.8 {
    quietly count if p_b2m_na >= `cut' & kru_ge2 == 1 & kru_pos == 1 & !missing(p_b2m_na)
    local TP = r(N)
    quietly count if p_b2m_na < `cut'  & kru_ge2 == 0 & kru_pos == 1 & !missing(p_b2m_na)
    local TN = r(N)
    quietly count if p_b2m_na >= `cut' & kru_ge2 == 0 & kru_pos == 1 & !missing(p_b2m_na)
    local FP = r(N)
    quietly count if p_b2m_na < `cut'  & kru_ge2 == 1 & kru_pos == 1 & !missing(p_b2m_na)
    local FN = r(N)
    if (`TP' + `FN') > 0 & (`TN' + `FP') > 0 {
        local Se  = `TP' / (`TP' + `FN')
        local Sp  = `TN' / (`TN' + `FP')
        local VPP = cond((`TP' + `FP') > 0, `TP' / (`TP' + `FP'), .)
        local VPN = cond((`TN' + `FN') > 0, `TN' / (`TN' + `FN'), .)
        local J   = `Se' + `Sp' - 1
        display "  ≥`cut' |  " %5.1f 100*`Se' "   " %5.1f 100*`Sp' ///
                "   " %5.1f 100*`VPP' "   " %5.1f 100*`VPN' ///
                "   " %5.3f `J'
    }
}

* ─── 9f. Bilan ──────────────────────────────────────────────────
display _newline(2) "========================================"
display              "  BILAN — SECTION 9"
display              "========================================"
display "  Modèle parcimonieux : logit P(KRU≥2) ~ UMOD + B2M"
display "  → B2M apporte un gain significatif au-delà de UMOD seul"
display "    (LR test, DeLong) sans circularité ni variable redondante."
display "  → AUC UMOD+B2M vs UMOD seul : voir 9d."
display "  → Seuil P≥0.5 recommandé si Se/Sp équilibrés requis."
display "========================================"
