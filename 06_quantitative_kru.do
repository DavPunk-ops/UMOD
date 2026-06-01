* ===========================================================================
* 06_quantitative_kru.do
* Prérequis : lancer 02_calculate_kru.do AVANT ce do-file.
*
* Objectif : relation quantitative UMOD + B2M → KRU continu (modèle two-part).
*
* → manuscrit : "Quantitative relationship between serum uromodulin,
*                β2-microglobulin, and KRU"
*
* Section 1 : corrélations Spearman UMOD / B2M / KRU (non-anuriques).
* Section 2 : modèle two-part parcimonieux UMOD + B2M.
*   2a : MFP logit P(KRU>0) → confirme la linéarité.
*   2b : logit P(KRU>0) ~ UMOD + B2M.
*   2c : MFP OLS E[KRU|KRU>0] → sélectionne (B2M/10)^-2.
*   2d : OLS E[KRU|KRU>0] ~ UMOD + (B2M/10)^-2, comparaison ΔAIC.
*   2e : combinaison two-part → KRU prédit continu.
*   2f : comparaison DeLong two-part vs logit direct (nécessite do-file 05).
* ===========================================================================

* --- Vérification des prérequis ---
foreach v in kru_daugirdas_35 kru_naif_35 umod diuresis labcreatprehd ///
             labb2mprehd age sex {
    capture confirm variable `v'
    if _rc {
        display as error "ERREUR : variable '`v'' absente — lance d'abord 02_calculate_kru.do"
        exit 111
    }
}

display _newline "=== Prérequis OK ==="

* --- Variables dérivées (recréées si absentes) ---
capture confirm variable kru_pos
if _rc {
    gen byte kru_pos = (kru_daugirdas_35 > 0) if !missing(kru_daugirdas_35)
    label variable kru_pos "KRU > 0 (1 = non-anurique)"
    label define krupos 0 "Anurique (KRU=0)" 1 "Non-anurique (KRU>0)", replace
    label values kru_pos krupos
}

capture confirm variable kru_ge2
if _rc {
    gen byte kru_ge2 = (kru_daugirdas_35 >= 2) if !missing(kru_daugirdas_35)
    label variable kru_ge2 "KRU >= 2 mL/min/35L"
    label define kruge2 0 "KRU < 2" 1 "KRU >= 2", replace
    label values kru_ge2 kruge2
}

capture confirm variable female
if _rc {
    capture drop female
    gen byte female = (sex == 1) if !missing(sex)
    label variable female "Sexe féminin (1=F, 0=H)"
}

* --- Transformation B2M^-2 (sélectionnée par MFP section 2c) ---
* Échelle X = B2M/10 (cf. output MFP : "where: X = labb2mprehd/10")
* On NE centre PAS : seul l'intercept changerait, pas les performances.
capture drop b2m_neg2
gen double b2m_neg2 = (labb2mprehd/10)^(-2) if !missing(labb2mprehd)
label variable b2m_neg2 "(B2M/10)^-2"

summarize labb2mprehd b2m_neg2

* ###########################################################################
* SECTION 1 — CORRÉLATIONS SPEARMAN UMOD / B2M / KRU (non-anuriques)
* ###########################################################################
display _newline(2) "=============================================="
display              "  1. Corrélations UMOD / B2M / KRU (non-anuriques)"
display              "=============================================="

spearman kru_daugirdas_35 umod labb2mprehd b2m_neg2 if kru_pos == 1, ///
    stats(rho p) star(0.05)

* ###########################################################################
* SECTION 2 — TWO-PART PARCIMONIEUX
*    Partie 1 : logit P(KRU>0) ~ UMOD + B2M  (linéaire)
*    Partie 2 : OLS  E[KRU|KRU>0] ~ UMOD + (B2M/10)^-2  (transformation MFP)
* ###########################################################################

* ===========================================================================
*  2a. MFP — forme fonctionnelle optimale pour le logit P(KRU>0)
*      La probabilité prédite alimentant la prédiction continue (produit
*      two-part), ses VALEURS — et non seulement leurs rangs — importent ;
*      le MFP est donc justifié comme pour la partie OLS (section 2c).
*      Résultat attendu : UMOD linéaire, B2M linéaire (puissances = 1)
* ===========================================================================
display _newline(2) "=============================================="
display              "  2a. MFP logit P(KRU>0) ~ UMOD + B2M"
display              "=============================================="

mfp logit kru_pos umod labb2mprehd

* ===========================================================================
*  2b. PARTIE 1 — Logit P(KRU>0) ~ UMOD + B2M
* ===========================================================================
display _newline(2) "=============================================="
display              "  2b. Logit P(KRU>0)  ~  UMOD + B2M"
display              "=============================================="

logit kru_pos umod labb2mprehd
estimates store ps_logit

display _newline "  --- Odds ratios ---"
logit kru_pos umod labb2mprehd, or

quietly lroc, nograph
local auc_logit_ps = r(area)
display _newline "  AUC partie 1 (parcimonieux) = " %5.3f `auc_logit_ps'

* Rappels comparatifs
quietly logit kru_pos umod
quietly lroc, nograph
local auc_logit_u = r(area)

quietly logit kru_pos umod age female labb2mprehd
quietly lroc, nograph
local auc_logit_full = r(area)

display _newline "  --- Comparaison AUC partie 1 ---"
display "    UMOD seul         : AUC = " %5.3f `auc_logit_u'
display "    UMOD + B2M        : AUC = " %5.3f `auc_logit_ps'
display "    UMOD+age+sex+B2M  : AUC = " %5.3f `auc_logit_full'

* Calibration
quietly logit kru_pos umod labb2mprehd
estat gof, group(10) table

* Prédiction
capture drop p_pos_ps
predict p_pos_ps, pr
label variable p_pos_ps "P(KRU>0 | UMOD, B2M)"

* ===========================================================================
*  2c. MFP — forme fonctionnelle optimale pour OLS E[KRU | KRU>0]
*      Résultat attendu : UMOD linéaire, B2M puissance -2 (échelle /10)
* ===========================================================================
display _newline(2) "=============================================="
display              "  2c. MFP OLS E[KRU | KRU>0] ~ UMOD + B2M"
display              "=============================================="

mfp regress kru_daugirdas_35 umod labb2mprehd if kru_pos == 1

* ===========================================================================
*  2d. OLS E[KRU | KRU>0] ~ UMOD + (B2M/10)^-2
*      (transformation sélectionnée par MFP section 2c, SE robustes)
* ===========================================================================
display _newline(2) "=============================================="
display              "  2d. OLS E[KRU | KRU>0] ~ UMOD + (B2M/10)^-2"
display              "      (transformation MFP section 2c, SE robustes)"
display              "=============================================="

regress kru_daugirdas_35 umod b2m_neg2 if kru_pos == 1, vce(robust)
estimates store ps_ols

local r2_ps    = e(r2)
local r2a_ps   = e(r2_a)
local rmse_ps  = e(rmse)
quietly estat ic
matrix _ic_ps = r(S)
local aic_ps  = _ic_ps[1,5]

display _newline "  R² = " %5.3f `r2_ps' "   R²adj = " %5.3f `r2a_ps' ///
    "   RMSE = " %5.3f `rmse_ps' "   AIC = " %6.2f `aic_ps'

* Rappels comparatifs
quietly regress kru_daugirdas_35 umod if kru_pos == 1, vce(robust)
local rmse_u = e(rmse)
local r2_u   = e(r2)
quietly estat ic
matrix _ic_u = r(S)
local aic_u  = _ic_u[1,5]

quietly regress kru_daugirdas_35 umod labb2mprehd if kru_pos == 1, vce(robust)
local rmse_lin = e(rmse)
local r2_lin   = e(r2)
quietly estat ic
matrix _ic_lin = r(S)
local aic_lin  = _ic_lin[1,5]

display _newline "  --- Comparaison OLS partie 2 ---"
display "                       R²       RMSE      AIC"
display "    UMOD seul       " %6.3f `r2_u'    "   " %6.3f `rmse_u'    "   " %6.2f `aic_u'
display "    UMOD + B2M lin. " %6.3f `r2_lin'  "   " %6.3f `rmse_lin'  "   " %6.2f `aic_lin'
display "    UMOD + B2M^-2   " %6.3f `r2_ps'   "   " %6.3f `rmse_ps'   "   " %6.2f `aic_ps'
display _newline "    ΔAIC (lin − B2M^-2) = " %5.2f (`aic_lin' - `aic_ps')
display "      > 2 : B2M^-2 préférable (confirme MFP section 2c)"

* Diagnostic résidus
capture drop resid_ps
predict resid_ps if e(sample), resid

display _newline "  Normalité des résidus :"
swilk resid_ps

* Prédiction conditionnelle
capture drop kru_cond_ps
predict kru_cond_ps, xb
label variable kru_cond_ps "E[KRU | KRU>0, UMOD, B2M^-2]"
replace kru_cond_ps = 0 if kru_cond_ps < 0

* ===========================================================================
*  2e. COMBINAISON TWO-PART
* ===========================================================================
display _newline(2) "=============================================="
display              "  2e. Prédiction two-part parcimonieuse"
display              "=============================================="

capture drop kru_pred_ps
gen double kru_pred_ps = p_pos_ps * kru_cond_ps if !missing(p_pos_ps, kru_cond_ps)
label variable kru_pred_ps "KRU prédit two-part parcimonieux"

capture drop _r_sq _r_abs
gen double _r_sq  = (kru_daugirdas_35 - kru_pred_ps)^2
gen double _r_abs = abs(kru_daugirdas_35 - kru_pred_ps)

quietly summarize _r_sq
local rmse_2p = sqrt(r(mean))
quietly summarize _r_abs
local mae_2p = r(mean)
drop _r_sq _r_abs

display _newline "  RMSE two-part parcimonieux = " %5.3f `rmse_2p' " mL/min/35L"
display         "  MAE  two-part parcimonieux = " %5.3f `mae_2p'  " mL/min/35L"

quietly corr kru_daugirdas_35 kru_pred_ps
display "  Corrélation Pearson (observé,prédit) = " %5.3f r(rho)
quietly spearman kru_daugirdas_35 kru_pred_ps
display "  Corrélation Spearman                 = " %5.3f r(rho)

* Graphique observé vs prédit
twoway (scatter kru_daugirdas_35 kru_pred_ps, msize(small)) ///
       (function y=x, range(0 10) lcolor(red) lpattern(dash)), ///
    title("KRU observé vs prédit — two-part parcimonieux") ///
    xtitle("KRU prédit (mL/min/35L)") ///
    ytitle("KRU observé (mL/min/35L)") ///
    legend(off) name(ps_obs_pred, replace)

* ===========================================================================
*  2f. COMPARAISON : two-part continu vs logit direct pour discriminer KRU≥2
*
*  ATTENTION : cette section nécessite que p_ge2_ps soit disponible en mémoire.
*              p_ge2_ps est créé par 05_multivariate_kru.do (section 2b).
*              Exécuter 05_multivariate_kru.do avant cette dofile pour activer
*              cette section. Si p_ge2_ps est absent, la section est ignorée.
* ===========================================================================
display _newline(2) "=============================================="
display              "  2f. Two-part continu vs logit direct"
display              "      pour discriminer KRU≥2 (DeLong)"
display              "=============================================="

capture confirm variable p_ge2_ps
if _rc {
    display as error "  AVERTISSEMENT : p_ge2_ps absent."
    display          "  → Exécuter d'abord 05_multivariate_kru.do, puis relancer cette section."
}
else {
    quietly roctab kru_ge2 kru_pred_ps
    local auc_cont = r(area)
    quietly roctab kru_ge2 p_ge2_ps
    local auc_dir  = r(area)

    display _newline "  AUC score continu two-part (kru_pred_ps) = " %5.3f `auc_cont'
    display         "  AUC logit direct          (p_ge2_ps)     = " %5.3f `auc_dir'
    display         "  Δ AUC (direct − continu)                 = " %6.3f (`auc_dir' - `auc_cont')

    display _newline "  --- Test de DeLong (AUC appariées) ---"
    roccomp kru_ge2 kru_pred_ps p_ge2_ps, summary
}
