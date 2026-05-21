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
*   PARTIE B — MODÉLISATION DE KRU CONTINU            [sections 3-6]
*              two-part model + test de non-linéarité
*   PARTIE C — DÉCISION CLINIQUE AU SEUIL KRU ≥ 2     [sections 7-8]
*              logistique sur tous (7), puis sur non-anuriques (8)
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
* #   Two-part model (sections 3-5) + test de non-linéarité (section 6)     #
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

* ===========================================================================
*  SECTION 6 — TEST DE NON-LINÉARITÉ (Restricted Cubic Splines)
*    Méthode Harrell : 4 nœuds aux percentiles 5/35/65/95.
*    Test LR du modèle RCS vs modèle linéaire pour les deux parties.
* ===========================================================================
display _newline(2) "=============================================="
display              "  SECTION 6 — NON-LINÉARITÉ (RCS 4 nœuds)"
display              "=============================================="

* Création des splines (3 termes pour 4 nœuds : umod_sp1 = linéaire,
* umod_sp2 et umod_sp3 = composantes non-linéaires)
capture drop umod_sp*
mkspline umod_sp = umod, cubic nknots(4) displayknots

* ─── 6a. PARTIE 1 — Logit RCS ─────────────────
display _newline "=============================================="
display         "  6a. Logit RCS : P(KRU>0) ~ rcs(UMOD)"
display         "=============================================="

logit kru_pos umod_sp1 umod_sp2 umod_sp3
estimates store logit_rcs

* Test de non-linéarité : H0 = composantes non-linéaires nulles
display _newline "=== Test de non-linéarité (Wald) ==="
test umod_sp2 umod_sp3

* LR test vs modèle linéaire
display _newline "=== LR test : RCS vs linéaire ==="
lrtest logit_umod logit_rcs

* Comparaison AIC/BIC
display _newline "=== AIC/BIC comparatifs ==="
estimates stats logit_umod logit_rcs

* Effet prédit : prédiction directe sur l'échantillon observé
* (margins ne peut pas extrapoler sur umod car le modèle contient les splines)
capture drop p_rcs_p1
predict p_rcs_p1, pr

twoway (line p_rcs_p1 umod, sort lcolor(red) lwidth(medium)) ///
       (scatter kru_pos umod, msize(small) mcolor(%40) jitter(2)), ///
    title("P(KRU>0) selon UMOD — spline cubique 4 nœuds") ///
    xtitle("UMOD (ng/mL)") ytitle("Probabilité prédite / observé") ///
    legend(order(1 "Spline RCS" 2 "Observé (jitter)") position(11) ring(0)) ///
    name(margins_logit_rcs, replace)

* ─── 6b. PARTIE 2 — OLS RCS sur non-anuriques ─────────────────
display _newline(2) "=============================================="
display              "  6b. OLS RCS : KRU ~ rcs(UMOD) | KRU > 0"
display              "=============================================="

regress kru_daugirdas_35 umod_sp1 umod_sp2 umod_sp3 if kru_pos == 1
estimates store ols_rcs_nonanuric

display _newline "=== R² RCS = " %5.3f e(r2) "  vs R² linéaire = 0.164"

* Test de non-linéarité
display _newline "=== Test de non-linéarité (F) ==="
test umod_sp2 umod_sp3

* LR test (via ftest puisque non-emboîtés en termes ML ; on utilise nestreg)
display _newline "=== Comparaison nested : linéaire vs RCS ==="
nestreg: regress kru_daugirdas_35 (umod) (umod_sp2 umod_sp3) if kru_pos == 1

* AIC/BIC
display _newline "=== AIC/BIC comparatifs ==="
estimates stats ols_umod_nonanuric ols_rcs_nonanuric

* Effet prédit : prédiction directe (idem partie 1)
estimates restore ols_rcs_nonanuric
capture drop yhat_rcs_p2_only
predict yhat_rcs_p2_only if kru_pos == 1, xb

twoway (scatter kru_daugirdas_35 umod if kru_pos == 1, ///
        msize(small) mcolor(%40)) ///
       (line yhat_rcs_p2_only umod if kru_pos == 1, ///
        sort lcolor(red) lwidth(medium)), ///
    title("KRU prédit selon UMOD (non-anuriques) — spline cubique") ///
    xtitle("UMOD (ng/mL)") ytitle("KRU (mL/min/35L)") ///
    legend(order(2 "Spline RCS" 1 "Observé") position(11) ring(0)) ///
    name(margins_ols_rcs, replace)

* ─── 6c. Comparaison visuelle linéaire vs RCS (non-anuriques) ──
estimates restore ols_umod_nonanuric
predict yhat_lin_p2 if kru_pos == 1, xb
estimates restore ols_rcs_nonanuric
predict yhat_rcs_p2 if kru_pos == 1, xb

twoway (scatter kru_daugirdas_35 umod if kru_pos == 1, msize(small) mcolor(%40)) ///
       (line yhat_lin_p2 umod if kru_pos == 1, sort lcolor(blue) lwidth(medium)) ///
       (line yhat_rcs_p2 umod if kru_pos == 1, sort lcolor(red) lwidth(medium)), ///
    title("KRU ~ UMOD chez non-anuriques : linéaire vs RCS") ///
    xtitle("UMOD (ng/mL)") ytitle("KRU (mL/min/35L)") ///
    legend(order(2 "Linéaire" 3 "Spline cubique") position(11) ring(0)) ///
    name(compare_lin_rcs, replace)

* ─── 6d. Bilan non-linéarité ─────────────────
display _newline(2) "========================================"
display              "  BILAN NON-LINÉARITÉ"
display              "========================================"
display "  → Si test umod_sp2 umod_sp3 p<0.05 : non-linéarité significative"
display "  → Si LR test / ΔAIC favorise RCS : préférer RCS au linéaire"
display "  → Sinon : conserver le modèle linéaire (parcimonie)"
display "========================================"

* ###########################################################################
* #                                                                         #
* #   PARTIE C — DÉCISION CLINIQUE AU SEUIL KRU ≥ 2 mL/min/35L              #
* #   Section 7 : logistique sur TOUS les patients (n=151)                  #
* #   Section 8 : logistique sur les NON-ANURIQUES uniquement (n=89)        #
* #               → analyse principale (vraie question clinique)            #
* #                                                                         #
* ###########################################################################

* ===========================================================================
*  SECTION 7 — Logit P(KRU≥2) ~ UMOD sur la population entière
*    Justification clinique : à KRU ≥ 2 mL/min/35L, la dose de dialyse
*    doit être adaptée (réduction Kt/V cible). C'est un seuil de décision
*    thérapeutique. Cette section inclut les anuriques (KRU=0 par définition,
*    donc < 2), ce qui amplifie artificiellement l'AUC. Voir section 8 pour
*    l'analyse clinique pertinente.
* ===========================================================================
display _newline(2) "=============================================="
display              "  SECTION 7 — SEUIL CLINIQUE KRU ≥ 2"
display              "=============================================="

* ─── 7a. Variable binaire kru_ge2 ─────────────────
gen byte kru_ge2 = (kru_daugirdas_35 >= 2) if !missing(kru_daugirdas_35)
label variable kru_ge2 "KRU >= 2 mL/min/35L (seuil clinique)"
label define kruge2 0 "< 2 (dose standard)" 1 ">= 2 (adapter dose)"
label values kru_ge2 kruge2

display _newline "=== Distribution kru_ge2 ==="
tab kru_ge2
* Détail par statut anurique
display _newline "=== Croisement kru_pos × kru_ge2 ==="
tab kru_pos kru_ge2, row

* ─── 7b. UMOD selon le seuil ─────────────────
display _newline "=== UMOD par statut KRU≥2 ==="
tabstat umod, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.2f)

* ─── 7c. Régression logistique P(KRU≥2) ~ UMOD ─────────────────
display _newline "=============================================="
display         "  7c. Logit : P(KRU≥2) ~ UMOD"
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

* ─── 7d. Recherche du seuil optimal d'UMOD ─────────────────
display _newline "=============================================="
display         "  7d. Seuil optimal d'UMOD (indice de Youden)"
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

* ─── 7e. Performance à des seuils cliniquement utiles ──────────
display _newline "=============================================="
display         "  7e. Performance aux seuils UMOD candidats"
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

* ─── 7f. Graphique : P(KRU≥2) selon UMOD ──────────────────────
estimates restore logit_ge2
quietly summarize umod
twoway (line p_ge2 umod, sort lcolor(red) lwidth(medium)) ///
       (scatter kru_ge2 umod, msize(small) mcolor(%40) jitter(2)), ///
    yline(0.5, lpattern(dot) lcolor(gray)) ///
    title("Probabilité prédite d'avoir KRU ≥ 2 mL/min/35L") ///
    xtitle("UMOD (ng/mL)") ytitle("P(KRU≥2 | UMOD)") ///
    legend(order(1 "Logit" 2 "Observé") position(11) ring(0)) ///
    name(logit_ge2_curve, replace)

* ─── 7g. Comparaison avec la prédiction two-part seuillée ──────
display _newline "=============================================="
display         "  7g. Comparaison : logit direct vs two-part seuillé"
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

* ─── 7h. Bilan ─────────────────
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
*  SECTION 8 — Logit P(KRU≥2) ~ UMOD chez les NON-ANURIQUES  [PRINCIPAL]
*    Vraie question clinique : chez les patients qui urinent (kru_pos=1),
*    UMOD permet-il de distinguer KRU<2 (dose standard) vs KRU≥2 (adapter) ?
*    Les anuriques sont exclus car par définition KRU=0 → pas besoin de
*    doser UMOD pour la décision clinique chez eux.
* ===========================================================================
display _newline(2) "=============================================="
display              "  SECTION 8 — KRU ≥ 2 CHEZ LES NON-ANURIQUES"
display              "=============================================="

* ─── 8a. Description de la population concernée ─────────────
display _newline "=== Distribution kru_ge2 chez les non-anuriques ==="
tab kru_ge2 if kru_pos == 1
display _newline "=== UMOD selon KRU<2 vs KRU>=2 chez non-anuriques ==="
tabstat umod if kru_pos == 1, by(kru_ge2) ///
    statistics(n mean sd p25 p50 p75) format(%6.2f)

* Test non paramétrique (Mann-Whitney)
display _newline "=== Test Mann-Whitney sur UMOD ==="
ranksum umod if kru_pos == 1, by(kru_ge2)

* ─── 8b. Régression logistique restreinte aux non-anuriques ──
display _newline "=============================================="
display         "  8b. Logit : P(KRU≥2) ~ UMOD | non-anurique"
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

* ─── 8c. Recherche du seuil UMOD optimal (Youden) ─────────────
display _newline "=============================================="
display         "  8c. Seuil optimal UMOD (Youden, non-anuriques)"
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

* ─── 8d. Performance aux seuils candidats ─────────────
display _newline "=============================================="
display         "  8d. Performance aux seuils UMOD candidats"
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

* ─── 8e. Graphique : P(KRU≥2 | non-anurique) selon UMOD ────────
estimates restore logit_ge2_na
twoway (line p_ge2_na umod if kru_pos == 1, sort lcolor(red) lwidth(medium)) ///
       (scatter kru_ge2 umod if kru_pos == 1, msize(small) ///
        mcolor(%40) jitter(2)), ///
    yline(0.5, lpattern(dot) lcolor(gray)) ///
    title("P(KRU≥2) chez les non-anuriques") ///
    xtitle("UMOD (ng/mL)") ytitle("P(KRU≥2 | UMOD, non-anurique)") ///
    legend(order(1 "Logit" 2 "Observé") position(11) ring(0)) ///
    name(logit_ge2_na_curve, replace)

* ─── 8f. Comparaison avec la section 7 (population entière) ────
display _newline "=============================================="
display         "  8f. Comparaison restreint vs population totale"
display         "=============================================="
display "  Population entière (section 7) :"
quietly logit kru_ge2 umod
quietly lroc, nograph
display "     AUC = " %5.3f r(area)
display "  Non-anuriques uniquement (section 8) :"
quietly logit kru_ge2 umod if kru_pos == 1
quietly lroc, nograph
display "     AUC = " %5.3f r(area)
display _newline "  → Cette section répond à la VRAIE question clinique :"
display "    chez un patient qui urine, UMOD prédit-il KRU≥2 ?"

* ─── 8g. Bilan ─────────────────
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
*  SECTION 9 — MULTIVARIÉ : UMOD + âge + sexe  [PRINCIPAL = sur non-anuriques]
*    Objectif : tester si âge et sexe apportent de l'information indépendante
*    à UMOD pour prédire KRU≥2 chez les patients qui urinent.
*
*    Stratégie :
*      9a. Description des covariables
*      9b. Modèle principal : entraîné ET testé sur non-anuriques (N=89)
*      9c. LR test : âge + sexe ajoutent-ils à UMOD seul ?
*      9d. Sensibilité : modèle entraîné sur N=151, testé sur N=89
*      9e. Comparaison AUC des trois modèles (univarié, multivar NA, multivar all)
* ===========================================================================
display _newline(2) "=============================================="
display              "  SECTION 9 — MULTIVARIÉ (UMOD + âge + sexe)"
display              "=============================================="

* ─── 9a. Description des covariables ──────────────────────────
display _newline "=== Âge selon kru_ge2 chez les non-anuriques ==="
tabstat age if kru_pos == 1, by(kru_ge2) ///
    statistics(n mean sd p25 p50 p75) format(%6.1f)

display _newline "=== Sexe selon kru_ge2 chez les non-anuriques ==="
tab sex kru_ge2 if kru_pos == 1, row

* Tests bivariés
display _newline "=== Tests bivariés (non-anuriques) ==="
ranksum age if kru_pos == 1, by(kru_ge2)
tab sex kru_ge2 if kru_pos == 1, chi2 exact

* ─── 9b. Modèle multivarié — entraîné sur non-anuriques [PRINCIPAL] ──
display _newline "=============================================="
display         "  9b. Logit multivarié | non-anuriques (PRINCIPAL)"
display         "=============================================="

logit kru_ge2 umod age i.sex if kru_pos == 1
estimates store logit_mvar_na

display _newline "=== Odds ratios ==="
logit kru_ge2 umod age i.sex if kru_pos == 1, or

* Probabilité prédite
capture drop p_mvar_na
predict p_mvar_na if kru_pos == 1, pr
label variable p_mvar_na "P(KRU≥2) multivar | non-anuriques"

* ROC / AUC
display _newline "=== AUC ==="
lroc, name(roc_mvar_na, replace) ///
    title("ROC multivarié — non-anuriques")

* Calibration
display _newline "=== Calibration Hosmer-Lemeshow ==="
estat gof, group(10) table

* ─── 9c. LR test vs UMOD seul (apport d'âge + sexe) ──────────
display _newline "=============================================="
display         "  9c. Apport d'âge + sexe vs UMOD seul"
display         "=============================================="

* Refit UMOD seul sur exactement le même échantillon (e(sample) du multivar)
estimates restore logit_mvar_na
gen byte _smvar = e(sample)
quietly logit kru_ge2 umod if _smvar == 1
estimates store logit_umod_only_na

display _newline "=== LR test : UMOD seul vs UMOD + âge + sexe ==="
lrtest logit_umod_only_na logit_mvar_na

display _newline "=== AIC/BIC comparatifs ==="
estimates stats logit_umod_only_na logit_mvar_na

drop _smvar

* ─── 9d. SENSIBILITÉ : entraîné sur TOUS, testé sur non-anuriques ──
display _newline "=============================================="
display         "  9d. Sensibilité : modèle entraîné sur tous"
display         "       et appliqué aux non-anuriques"
display         "=============================================="

logit kru_ge2 umod age i.sex
estimates store logit_mvar_all

display _newline "=== Odds ratios (modèle sur tous) ==="
logit kru_ge2 umod age i.sex, or

capture drop p_mvar_all
predict p_mvar_all, pr
label variable p_mvar_all "P(KRU≥2) multivar | échantillon entier"

* AUC sur tous (entraînement)
display _newline "=== AUC sur tous (entraînement) ==="
lroc, nograph
display "    AUC sur N=151 : " %5.3f r(area)

* AUC restreint aux non-anuriques (test pertinent)
display _newline "=== AUC restreint aux non-anuriques (test) ==="
roctab kru_ge2 p_mvar_all if kru_pos == 1

* ─── 9e. Comparaison des trois modèles ────────────────────────
display _newline "=============================================="
display         "  9e. Comparaison AUC sur les non-anuriques"
display         "=============================================="

display _newline "  Tous les AUC sont calculés sur le MÊME échantillon test"
display "  (les 89 non-anuriques) :"
display _newline "  ① UMOD seul (section 8) :"
quietly roctab kru_ge2 umod if kru_pos == 1
display "      AUC = " %5.3f r(area)

display "  ② Multivar entraîné sur non-anuriques (9b) :"
quietly roctab kru_ge2 p_mvar_na if kru_pos == 1
display "      AUC = " %5.3f r(area)

display "  ③ Multivar entraîné sur tous, testé NA (9d) :"
quietly roctab kru_ge2 p_mvar_all if kru_pos == 1
display "      AUC = " %5.3f r(area)

* Test statistique de différence (DeLong) entre ② et ③
display _newline "=== Test DeLong : multivar NA vs multivar all ==="
roccomp kru_ge2 p_mvar_na p_mvar_all if kru_pos == 1, graph summary ///
    name(roccomp_mvar, replace)

* ─── 9f. Seuil optimal (Youden) sur le score multivarié ──────
display _newline "=============================================="
display         "  9f. Seuil optimal sur le score multivarié"
display         "=============================================="

estimates restore logit_mvar_na

* Le seuil opère sur la probabilité prédite (et non plus sur UMOD direct)
* car le score multivarié intègre âge et sexe
foreach cut in 0.3 0.4 0.5 0.6 0.7 {
    display _newline "  ─── Seuil P(KRU≥2) ≥ `cut' ───"
    quietly count if p_mvar_na >= `cut' & kru_ge2 == 1 & kru_pos == 1
    local TP = r(N)
    quietly count if p_mvar_na < `cut' & kru_ge2 == 0 & kru_pos == 1
    local TN = r(N)
    quietly count if p_mvar_na >= `cut' & kru_ge2 == 0 & kru_pos == 1
    local FP = r(N)
    quietly count if p_mvar_na < `cut' & kru_ge2 == 1 & kru_pos == 1
    local FN = r(N)
    if (`TP' + `FN') > 0 & (`TN' + `FP') > 0 {
        local Se = `TP' / (`TP' + `FN')
        local Sp = `TN' / (`TN' + `FP')
        local J = `Se' + `Sp' - 1
        display "    Se = " %5.1f 100*`Se' " %    Sp = " %5.1f 100*`Sp' " %    J = " %5.3f `J'
    }
}

* ─── 9g. Bilan ──────────────────────────────────────────────
display _newline(2) "========================================"
display              "  BILAN — MULTIVARIÉ CHEZ NON-ANURIQUES"
display              "========================================"
display "  Modèle principal : logit P(KRU≥2) ~ UMOD + âge + i.sex"
display "                     entraîné sur les 89 non-anuriques"
display "  → voir AUC, LR test et seuil Youden ci-dessus"
display "========================================"

* ===========================================================================
*  SECTION 10 — MULTIVARIÉ + VINTAGE  [extension de la section 9]
*    Ajoute la vintage en dialyse au modèle UMOD + âge + sexe.
*    Hypothèse : vintage plus longue → moins de RKF résiduelle → KRU<2 plus
*    probable. C'est typiquement la covariable la plus prédictive du déclin
*    de la fonction rénale résiduelle.
*
*    Stratégie identique à la section 9 :
*      10a. Description de vintage
*      10b. Modèle principal : UMOD + âge + sexe + vintage sur non-anuriques
*      10c. LR tests d'apport (vs UMOD seul, vs UMOD+âge+sexe)
*      10d. Sensibilité : entraîné sur tous, testé sur non-anuriques
*      10e. Comparaison AUC : 4 modèles côte-à-côte
*      10f. Seuils Youden sur le score complet
* ===========================================================================
display _newline(2) "=============================================="
display              "  SECTION 10 — MULTIVARIÉ + VINTAGE"
display              "=============================================="

* Vérification présence et type de vintage
capture confirm variable vintage
if _rc {
    display as error "ERREUR : variable 'vintage' absente — vérifie le merge"
    exit 111
}

* ─── 10a. Description de vintage ──────────────────────────────
display _newline "=== Distribution de vintage (échantillon entier) ==="
summarize vintage, detail

display _newline "=== Vintage selon kru_ge2 chez les non-anuriques ==="
tabstat vintage if kru_pos == 1, by(kru_ge2) ///
    statistics(n mean sd p25 p50 p75) format(%6.1f)

display _newline "=== Test Mann-Whitney sur vintage ==="
ranksum vintage if kru_pos == 1, by(kru_ge2)

* Corrélation vintage / KRU continu pour info
display _newline "=== Corrélation vintage / KRU continu (non-anuriques) ==="
corr kru_daugirdas_35 vintage if kru_pos == 1
spearman kru_daugirdas_35 vintage if kru_pos == 1

* ─── 10b. Modèle principal : UMOD + âge + sexe + vintage | NA ──
display _newline "=============================================="
display         "  10b. Logit multivar+vintage | non-anuriques"
display         "       (PRINCIPAL)"
display         "=============================================="

logit kru_ge2 umod age i.sex vintage if kru_pos == 1
estimates store logit_mvar_vin_na

display _newline "=== Odds ratios ==="
logit kru_ge2 umod age i.sex vintage if kru_pos == 1, or

capture drop p_mvar_vin_na
predict p_mvar_vin_na if kru_pos == 1, pr
label variable p_mvar_vin_na "P(KRU≥2) UMOD+age+sex+vintage | NA"

display _newline "=== AUC ==="
lroc, name(roc_mvar_vin_na, replace) ///
    title("ROC multivar + vintage — non-anuriques")

display _newline "=== Calibration Hosmer-Lemeshow ==="
estat gof, group(10) table

* ─── 10c. LR tests d'apport ──────────────────────────────────
display _newline "=============================================="
display         "  10c. Apport de vintage (LR tests)"
display         "=============================================="

* Fixer l'échantillon du multivar+vintage
estimates restore logit_mvar_vin_na
gen byte _smv2 = e(sample)

* Refit UMOD seul sur le même échantillon
quietly logit kru_ge2 umod if _smv2 == 1
estimates store logit_u_only_s2

* Refit UMOD + âge + sexe sur le même échantillon
quietly logit kru_ge2 umod age i.sex if _smv2 == 1
estimates store logit_uas_s2

display _newline "=== LR test : UMOD seul vs UMOD+age+sex+vintage ==="
lrtest logit_u_only_s2 logit_mvar_vin_na

display _newline "=== LR test : UMOD+age+sex vs +vintage ==="
lrtest logit_uas_s2 logit_mvar_vin_na

display _newline "=== AIC/BIC : 3 modèles emboîtés ==="
estimates stats logit_u_only_s2 logit_uas_s2 logit_mvar_vin_na

drop _smv2

* ─── 10d. SENSIBILITÉ : entraîné sur TOUS, testé sur NA ───────
display _newline "=============================================="
display         "  10d. Sensibilité : entraîné sur tous"
display         "=============================================="

logit kru_ge2 umod age i.sex vintage
estimates store logit_mvar_vin_all

display _newline "=== Odds ratios (modèle sur tous) ==="
logit kru_ge2 umod age i.sex vintage, or

capture drop p_mvar_vin_all
predict p_mvar_vin_all, pr
label variable p_mvar_vin_all "P(KRU≥2) UMOD+age+sex+vintage | tous"

display _newline "=== AUC sur tous (entraînement) ==="
lroc, nograph
display "    AUC sur N=151 : " %5.3f r(area)

display _newline "=== AUC restreint aux non-anuriques (test) ==="
roctab kru_ge2 p_mvar_vin_all if kru_pos == 1

* ─── 10e. Comparaison AUC des 4 modèles sur les NA ────────────
display _newline "=============================================="
display         "  10e. Comparaison AUC sur les non-anuriques"
display         "=============================================="
display _newline "  AUC calculés sur les MÊMES 89 non-anuriques :"

display _newline "  ① UMOD seul (section 8) :"
quietly roctab kru_ge2 umod if kru_pos == 1
display "      AUC = " %5.3f r(area)

display "  ② UMOD + âge + sexe | NA (section 9b) :"
quietly roctab kru_ge2 p_mvar_na if kru_pos == 1
display "      AUC = " %5.3f r(area)

display "  ③ UMOD + âge + sexe + vintage | NA (10b) :"
quietly roctab kru_ge2 p_mvar_vin_na if kru_pos == 1
display "      AUC = " %5.3f r(area)

display "  ④ UMOD + âge + sexe + vintage | tous (10d, testé NA) :"
quietly roctab kru_ge2 p_mvar_vin_all if kru_pos == 1
display "      AUC = " %5.3f r(area)

* Test DeLong sur les paires intéressantes
display _newline "=== Test DeLong : UMOD seul vs +vintage (NA) ==="
roccomp kru_ge2 umod p_mvar_vin_na if kru_pos == 1, graph summary ///
    name(roccomp_vintage, replace)

display _newline "=== Test DeLong : +vintage NA vs +vintage all ==="
roccomp kru_ge2 p_mvar_vin_na p_mvar_vin_all if kru_pos == 1, summary

* ─── 10f. Seuils Youden sur le score complet ─────────────────
display _newline "=============================================="
display         "  10f. Seuil Youden sur P(KRU≥2) du modèle complet"
display         "=============================================="

estimates restore logit_mvar_vin_na

foreach cut in 0.3 0.4 0.5 0.6 0.7 0.8 {
    display _newline "  ─── Seuil P(KRU≥2) ≥ `cut' ───"
    quietly count if p_mvar_vin_na >= `cut' & kru_ge2 == 1 & kru_pos == 1
    local TP = r(N)
    quietly count if p_mvar_vin_na < `cut' & kru_ge2 == 0 & kru_pos == 1
    local TN = r(N)
    quietly count if p_mvar_vin_na >= `cut' & kru_ge2 == 0 & kru_pos == 1
    local FP = r(N)
    quietly count if p_mvar_vin_na < `cut' & kru_ge2 == 1 & kru_pos == 1
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

* ─── 10g. Bilan ──────────────────────────────────────────────
display _newline(2) "========================================"
display              "  BILAN — MULTIVARIÉ + VINTAGE"
display              "========================================"
display "  Modèle complet : logit P(KRU≥2) ~ UMOD + age + i.sex + vintage"
display "                   entraîné sur les non-anuriques (N=89)"
display "  → voir AUC, LR tests et seuil Youden ci-dessus"
display "========================================"

* ===========================================================================
*  SECTION 11 — UMOD + BIOMARQUEURS PRÉDIALYSE
*    Biomarqueurs explorés : créatinine, urée, β2-microglobuline prédialyse.
*    Rationale physiologique : tous ces marqueurs s'accumulent quand la
*    clairance rénale résiduelle baisse (creat & urée filtration ; β2M
*    extrêmement peu clairée par l'HD, donc reflet quasi pur de la RKF).
*    Hypothèse directionnelle : OR < 1 attendu (marqueur ↑ → KRU<2).
*
*    Stratégie :
*      11a. Description des 3 biomarqueurs (valeurs manquantes, par groupe)
*      11b. Modèle principal : UMOD + créat + urée + B2M | non-anuriques
*      11c. LR tests d'apport (vs UMOD seul, et par biomarqueur)
*      11d. Modèle étendu : + vintage (synthèse)
*      11e. Sensibilité : entraîné sur tous, testé sur NA
*      11f. Comparaison AUC tous modèles
*      11g. Seuils Youden sur le meilleur score
* ===========================================================================
display _newline(2) "=============================================="
display              "  SECTION 11 — UMOD + BIOMARQUEURS PRÉDIALYSE"
display              "=============================================="

* Vérification présence des variables
foreach v in labcreatprehd labureaprehd labb2mprehd {
    capture confirm variable `v'
    if _rc {
        display as error "ERREUR : variable '`v'' absente"
        exit 111
    }
}

* ─── 11a. Description des biomarqueurs ──────────────────────────
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

* ─── 11b. Modèle principal : UMOD + 3 biomarqueurs | NA ────────
display _newline "=============================================="
display         "  11b. Logit : UMOD + créat + urée + B2M | NA"
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

* ─── 11c. LR tests d'apport ──────────────────────────────────
display _newline "=============================================="
display         "  11c. Apport des biomarqueurs (LR tests)"
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

* ─── 11d. Modèle étendu : UMOD + biomarqueurs + vintage ──────
display _newline "=============================================="
display         "  11d. Modèle étendu : + vintage"
display         "=============================================="

logit kru_ge2 umod labcreatprehd labureaprehd labb2mprehd vintage if kru_pos == 1
estimates store logit_bio_vin_na

display _newline "=== Odds ratios ==="
logit kru_ge2 umod labcreatprehd labureaprehd labb2mprehd vintage if kru_pos == 1, or

capture drop p_bio_vin_na
predict p_bio_vin_na if kru_pos == 1, pr
label variable p_bio_vin_na "P(KRU≥2) UMOD+biomark+vintage | NA"

display _newline "=== AUC ==="
lroc, nograph
display "    AUC = " %5.3f r(area)

* LR test apport de vintage par-dessus les biomarqueurs
display _newline "=== LR test : biomarqueurs seuls vs + vintage ==="
estimates restore logit_bio_vin_na
gen byte _sbv = e(sample)
quietly logit kru_ge2 umod labcreatprehd labureaprehd labb2mprehd if _sbv == 1
estimates store logit_bio_na_sbv
lrtest logit_bio_na_sbv logit_bio_vin_na
drop _sbv

* ─── 11e. SENSIBILITÉ : entraîné sur tous, testé sur NA ───────
display _newline "=============================================="
display         "  11e. Sensibilité : entraîné sur tous"
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

* ─── 11f. Comparaison de tous les modèles sur les NA ──────────
display _newline "=============================================="
display         "  11f. Comparaison AUC tous modèles | NA"
display         "=============================================="

display _newline "  ① UMOD seul (section 8) :"
quietly roctab kru_ge2 umod if kru_pos == 1
display "      AUC = " %5.3f r(area)

display "  ② UMOD + âge + sexe (9b) :"
quietly roctab kru_ge2 p_mvar_na if kru_pos == 1
display "      AUC = " %5.3f r(area)

display "  ③ UMOD + âge + sexe + vintage (10b) :"
quietly roctab kru_ge2 p_mvar_vin_na if kru_pos == 1
display "      AUC = " %5.3f r(area)

display "  ④ UMOD + 3 biomarqueurs (11b) :"
quietly roctab kru_ge2 p_bio_na if kru_pos == 1
display "      AUC = " %5.3f r(area)

display "  ⑤ UMOD + 3 biomarqueurs + vintage (11d) :"
quietly roctab kru_ge2 p_bio_vin_na if kru_pos == 1
display "      AUC = " %5.3f r(area)

display "  ⑥ UMOD + 3 biomarqueurs entraîné sur tous (11e) :"
quietly roctab kru_ge2 p_bio_all if kru_pos == 1
display "      AUC = " %5.3f r(area)

display _newline "=== Test DeLong : UMOD seul vs UMOD+biomark ==="
roccomp kru_ge2 umod p_bio_na if kru_pos == 1, graph summary ///
    name(roccomp_bio, replace)

display _newline "=== Test DeLong : +vintage vs +biomark ==="
capture roccomp kru_ge2 p_mvar_vin_na p_bio_na if kru_pos == 1, summary

display _newline "=== Test DeLong : +biomark seul vs +biomark+vintage ==="
capture roccomp kru_ge2 p_bio_na p_bio_vin_na if kru_pos == 1, summary

* ─── 11g. Seuils Youden sur le modèle UMOD + biomarqueurs ─────
display _newline "=============================================="
display         "  11g. Seuil Youden sur P(KRU≥2) - biomarqueurs"
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

* ─── 11h. Bilan ──────────────────────────────────────────────
display _newline(2) "========================================"
display              "  BILAN — UMOD + BIOMARQUEURS"
display              "========================================"
display "  Hypothèse : créat, urée, B2M ↑ → KRU<2 plus probable (OR<1)"
display "  Voir AUC, LR tests, OR et seuils Youden ci-dessus"
display "========================================"
