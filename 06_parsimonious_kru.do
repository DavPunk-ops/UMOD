* ===========================================================================
* 06_parsimonious_kru.do
* Prérequis : lancer 02_calculate_kru.do AVANT ce do-file.
*
* Objectif : modèle parcimonieux UMOD + B2M (sans age/sex car NS dans 05)
*            avec transformation B2M^-2 dans l'OLS (suggérée par MFP dans 05).
*
* Choix sur la base des résultats de 05_multivariable_kru.do :
*   - age et female non significatifs partout (p > 0.2) → retirés
*   - MFP a sélectionné B2M^-2 dans la partie OLS (ΔAIC=12.77 vs linéaire)
*   - MFP a confirmé linéaire optimal dans le logit kru_ge2 → on garde
*     B2M linéaire dans les deux logits
* ===========================================================================

* --- Vérification des prérequis ---
foreach v in kru_daugirdas_35 kru_naif_35 umod diuresis labcreatprehd ///
             labb2mprehd {
    capture confirm variable `v'
    if _rc {
        display as error "ERREUR : variable '`v'' absente — lance d'abord 02_calculate_kru.do"
        exit 111
    }
}

display _newline "=== Prérequis OK ==="

* --- Variables binaires (recréées si absentes) ---
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

* --- Transformation B2M^-2 suggérée par MFP (do-file 04, partie OLS) ---
* Échelle X = B2M/10 (cf. MFP : "where: X = labb2mprehd/10")
* On NE centre PAS : seul l'intercept changerait, pas les performances.
capture drop b2m_neg2
gen double b2m_neg2 = (labb2mprehd/10)^(-2) if !missing(labb2mprehd)
label variable b2m_neg2 "(B2M/10)^-2"

summarize labb2mprehd b2m_neg2

* ###########################################################################
* SECTION 1 — RAPPEL DES CORRÉLATIONS UMOD / B2M / KRU
* ###########################################################################
display _newline(2) "=============================================="
display              "  1. Corrélations UMOD / B2M / KRU (non-anuriques)"
display              "=============================================="

spearman kru_daugirdas_35 umod labb2mprehd b2m_neg2 if kru_pos == 1, ///
    stats(rho p) star(0.05)

* ###########################################################################
* SECTION 2 — TWO-PART PARCIMONIEUX
*    Partie 1 : logit P(KRU>0) ~ UMOD + B2M  (linéaire)
*    Partie 2 : OLS  E[KRU|KRU>0] ~ UMOD + B2M^-2  (transformation MFP)
* ###########################################################################

* ===========================================================================
*  2a. PARTIE 1 — Logit P(KRU>0) ~ UMOD + B2M
* ===========================================================================
display _newline(2) "=============================================="
display              "  2a. Logit P(KRU>0)  ~  UMOD + B2M"
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
*  2b. PARTIE 2 — OLS E[KRU | KRU>0] ~ UMOD + B2M^-2
* ===========================================================================
display _newline(2) "=============================================="
display              "  2b. OLS E[KRU | KRU>0] ~ UMOD + (B2M/10)^-2"
display              "      (transformation MFP, SE robustes)"
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
display "      > 2 : B2M^-2 préférable (confirme MFP de 04)"

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
*  2c. COMBINAISON TWO-PART
* ===========================================================================
display _newline(2) "=============================================="
display              "  2c. Prédiction two-part parcimonieuse"
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

* ###########################################################################
* SECTION 3 — LOGIT KRU≥2 PARCIMONIEUX (UMOD + B2M linéaire)
*    MFP du do-file 04 a confirmé que le linéaire est optimal pour ce logit.
* ###########################################################################

* ===========================================================================
*  3a. LOGIT + ROC
* ===========================================================================
display _newline(2) "=============================================="
display              "  3a. Logit P(KRU≥2) ~ UMOD + B2M (N=148)"
display              "=============================================="

logit kru_ge2 umod labb2mprehd
estimates store ge2_ps

display _newline "  --- Odds ratios ---"
logit kru_ge2 umod labb2mprehd, or

quietly lroc, nograph
local auc_ge2_ps = r(area)

* Calibration
quietly logit kru_ge2 umod labb2mprehd
estat gof, group(10) table

* Prédiction
capture drop p_ge2_ps
predict p_ge2_ps, pr
label variable p_ge2_ps "P(KRU≥2 | UMOD, B2M)"

* Rappels comparatifs
quietly logit kru_ge2 umod
quietly lroc, nograph
local auc_ge2_u = r(area)

quietly logit kru_ge2 umod age female labb2mprehd
quietly lroc, nograph
local auc_ge2_full = r(area)

display _newline "  --- AUC pour KRU≥2 ---"
display "    UMOD seul         : AUC = " %5.3f `auc_ge2_u'
display "    UMOD + B2M        : AUC = " %5.3f `auc_ge2_ps'
display "    UMOD+age+sex+B2M  : AUC = " %5.3f `auc_ge2_full'

* DeLong : parcimonieux vs UMOD seul
display _newline "  --- Comparaison DeLong (parcimonieux vs UMOD seul) ---"
roccomp kru_ge2 umod p_ge2_ps, graph summary name(ps_roc_compare, replace)

* ===========================================================================
*  3a-bis. QUELLE APPROCHE CLASSIFIE LE MIEUX KRU≥2 ?
*       (1) Score continu two-part   kru_pred_ps  (sections 2a-2c)
*       (2) Logit direct             p_ge2_ps     (section 3a)
*       Les deux sont des scores continus → AUC comparées par DeLong.
*       NB : comparaison sur AUC APPARENTES. Le logit direct est entraîné
*            exactement sur kru_ge2 (avantage attendu) ; si le two-part
*            continu fait jeu égal ou mieux, c'est un argument fort.
* ===========================================================================
display _newline(2) "=============================================="
display              "  3a-bis. Two-part continu vs logit direct"
display              "          pour discriminer KRU≥2 (DeLong)"
display              "=============================================="

quietly roctab kru_ge2 kru_pred_ps
local auc_cont = r(area)
quietly roctab kru_ge2 p_ge2_ps
local auc_dir  = r(area)

display _newline "  AUC score continu two-part (kru_pred_ps) = " %5.3f `auc_cont'
display         "  AUC logit direct          (p_ge2_ps)     = " %5.3f `auc_dir'
display         "  Δ AUC (direct − continu)                 = " %6.3f (`auc_dir' - `auc_cont')

display _newline "  --- Test de DeLong (AUC appariées) ---"
roccomp kru_ge2 kru_pred_ps p_ge2_ps, summary

* ===========================================================================
*  3b. CUTOFF YOUDEN sur P(KRU≥2) — parcimonieux
* ===========================================================================
display _newline(2) "=============================================="
display              "  3b. Cutoff Youden sur P̂ (parcimonieux)"
display              "=============================================="

local best_J_ps  = -1
local best_p_ps  = .
local best_Se_ps = .
local best_Sp_ps = .

forvalues p = 0.05(0.025)0.95 {
    quietly count if p_ge2_ps >= `p' & kru_ge2 == 1 & !missing(kru_ge2, p_ge2_ps)
    local TP = r(N)
    quietly count if p_ge2_ps <  `p' & kru_ge2 == 0 & !missing(kru_ge2, p_ge2_ps)
    local TN = r(N)
    quietly count if p_ge2_ps >= `p' & kru_ge2 == 0 & !missing(kru_ge2, p_ge2_ps)
    local FP = r(N)
    quietly count if p_ge2_ps <  `p' & kru_ge2 == 1 & !missing(kru_ge2, p_ge2_ps)
    local FN = r(N)
    if (`TP' + `FN') > 0 & (`TN' + `FP') > 0 {
        local Se = `TP' / (`TP' + `FN')
        local Sp = `TN' / (`TN' + `FP')
        local J  = `Se' + `Sp' - 1
        if `J' > `best_J_ps' {
            local best_J_ps  = `J'
            local best_p_ps  = `p'
            local best_Se_ps = `Se'
            local best_Sp_ps = `Sp'
        }
    }
}

display _newline "  → Cutoff optimal Youden : P(KRU≥2) ≥ " %5.3f `best_p_ps'
display         "      Se = " %5.1f 100*`best_Se_ps' " %"
display         "      Sp = " %5.1f 100*`best_Sp_ps' " %"
display         "      J  = " %5.3f `best_J_ps'

* ###########################################################################
* SECTION 4 — BOOTSTRAP HARRELL DU MODÈLE PARCIMONIEUX FINAL
*    Modèle : logit kru_ge2 ~ UMOD + B2M (linéaire)
* ###########################################################################

* Refit pour avoir e() à jour
quietly logit kru_ge2 umod labb2mprehd

local app_AUC    = `auc_ge2_ps'
local app_p_cut  = `best_p_ps'
local app_Se     = `best_Se_ps'
local app_Sp     = `best_Sp_ps'
local app_J      = `best_J_ps'

global B = 1000
global SEED = 20260522

display _newline(2) "=============================================="
display              "  4. Bootstrap parcimonieux (B=$B, seed=$SEED)"
display              "=============================================="
display _newline "  Performance apparente :"
display "    AUC            = " %5.3f `app_AUC'
display "    Cutoff P̂ opt.  = " %5.3f `app_p_cut'
display "    Se / Sp / J    = " %5.3f `app_Se' " / " %5.3f `app_Sp' " / " %5.3f `app_J'

set seed $SEED

tempname memh
tempfile bootres
postfile `memh' double(AUC_bb AUC_bo cutoff_p Jbb Jbo Sebb Sebo Spbb Spbo) ///
    using `bootres', replace

local n_valid = 0

display _newline _continue "  Progression : "

forvalues b = 1/$B {
    if mod(`b', 100) == 0 display _continue "`b' "

    preserve
    quietly keep if !missing(kru_ge2, umod, labb2mprehd)
    quietly bsample

    * (a) Fit du modèle dans le bootstrap
    capture quietly logit kru_ge2 umod labb2mprehd
    local fit_ok = (_rc == 0)
    if !`fit_ok' {
        restore
        continue
    }

    * (b) Prédiction dans le bootstrap et AUC_bb
    capture drop p_b
    quietly predict p_b, pr
    capture quietly roctab kru_ge2 p_b, nograph
    if _rc {
        restore
        continue
    }
    local auc_bb = r(area)

    * (c) Cutoff Youden sur P_b dans le bootstrap
    local best_J_b  = -1
    local best_c_b  = .
    local best_Se_b = .
    local best_Sp_b = .

    forvalues p = 0.05(0.025)0.95 {
        quietly count if p_b >= `p' & kru_ge2 == 1
        local TP = r(N)
        quietly count if p_b <  `p' & kru_ge2 == 0
        local TN = r(N)
        quietly count if p_b >= `p' & kru_ge2 == 0
        local FP = r(N)
        quietly count if p_b <  `p' & kru_ge2 == 1
        local FN = r(N)
        if (`TP' + `FN') > 0 & (`TN' + `FP') > 0 {
            local Se_ = `TP' / (`TP' + `FN')
            local Sp_ = `TN' / (`TN' + `FP')
            local J_  = `Se_' + `Sp_' - 1
            if `J_' > `best_J_b' {
                local best_J_b  = `J_'
                local best_c_b  = `p'
                local best_Se_b = `Se_'
                local best_Sp_b = `Sp_'
            }
        }
    }

    restore

    if missing(`best_c_b') continue

    * (d) Appliquer le modèle bootstrap (encore en e()) à l'ORIGINAL
    capture drop p_bo
    quietly predict p_bo if !missing(kru_ge2, umod, labb2mprehd), pr

    capture quietly roctab kru_ge2 p_bo, nograph
    if _rc {
        drop p_bo
        continue
    }
    local auc_bo = r(area)

    quietly count if p_bo >= `best_c_b' & kru_ge2 == 1 & !missing(p_bo)
    local TPo = r(N)
    quietly count if p_bo <  `best_c_b' & kru_ge2 == 0 & !missing(p_bo)
    local TNo = r(N)
    quietly count if p_bo >= `best_c_b' & kru_ge2 == 0 & !missing(p_bo)
    local FPo = r(N)
    quietly count if p_bo <  `best_c_b' & kru_ge2 == 1 & !missing(p_bo)
    local FNo = r(N)

    drop p_bo

    if (`TPo' + `FNo') == 0 | (`TNo' + `FPo') == 0 continue

    local Se_bo = `TPo' / (`TPo' + `FNo')
    local Sp_bo = `TNo' / (`TNo' + `FPo')
    local J_bo  = `Se_bo' + `Sp_bo' - 1

    post `memh' (`auc_bb') (`auc_bo') (`best_c_b') ///
                (`best_J_b') (`J_bo') ///
                (`best_Se_b') (`Se_bo') (`best_Sp_b') (`Sp_bo')

    local n_valid = `n_valid' + 1
}

postclose `memh'

display ""
display _newline "  Itérations valides : `n_valid'/$B"

preserve
quietly use `bootres', clear

quietly summarize AUC_bb, meanonly
local m_auc_bb = r(mean)
quietly summarize AUC_bo, meanonly
local m_auc_bo = r(mean)
local opt_AUC  = `m_auc_bb' - `m_auc_bo'

quietly _pctile AUC_bo, percentiles(2.5 97.5)
local auc_lo = r(r1)
local auc_hi = r(r2)

quietly summarize Jbb, meanonly
local m_Jbb = r(mean)
quietly summarize Jbo, meanonly
local m_Jbo = r(mean)
local opt_J = `m_Jbb' - `m_Jbo'

quietly summarize Sebb, meanonly
local m_Sebb = r(mean)
quietly summarize Sebo, meanonly
local m_Sebo = r(mean)
local opt_Se = `m_Sebb' - `m_Sebo'

quietly summarize Spbb, meanonly
local m_Spbb = r(mean)
quietly summarize Spbo, meanonly
local m_Spbo = r(mean)
local opt_Sp = `m_Spbb' - `m_Spbo'

quietly summarize cutoff_p, detail
local m_cut  = r(mean)
local sd_cut = r(sd)
local cut_lo = r(p5)
local cut_hi = r(p95)

quietly count if abs(cutoff_p - `app_p_cut') < 1e-6
local pct_apparent = 100 * r(N) / `n_valid'

restore

local AUCcorr = `app_AUC' - `opt_AUC'
local Jcorr   = `app_J'   - `opt_J'
local Secorr  = `app_Se'  - `opt_Se'
local Spcorr  = `app_Sp'  - `opt_Sp'

display _newline(2) "=========================================================="
display              "  RÉSULTATS BOOTSTRAP — MODÈLE PARCIMONIEUX UMOD + B2M"
display              "=========================================================="
display _newline "  --- AUC ---"
display "    AUC apparente                       = " %5.3f `app_AUC'
display "    AUC bootstrap (mean bb)             = " %5.3f `m_auc_bb'
display "    AUC honest    (mean bo)             = " %5.3f `m_auc_bo'
display "    Optimisme (bb − bo)                 = " %6.3f `opt_AUC'
display "    AUC corrigée pour optimisme         = " %5.3f `AUCcorr'
display "    IC 95% (percentile sur AUC_bo)      = " %5.3f `auc_lo' " — " %5.3f `auc_hi'

display _newline "  --- Cutoff optimal Youden (sur P̂) ---"
display "    Cutoff apparent       = " %5.3f `app_p_cut'
display "    Cutoff bootstrap moy  = " %5.3f `m_cut' "  (SD " %5.3f `sd_cut' ")"
display "    IC bootstrap (5–95%)  = " %5.3f `cut_lo' " — " %5.3f `cut_hi'
display "    % itérations retrouvant exactement le cutoff apparent : " %4.1f `pct_apparent' " %"

display _newline "  --- Performance au cutoff (corrigée) ---"
display "                         Apparent    Optimisme    Corrigé"
display "    Sensibilité    " %6.3f `app_Se'  "       " %6.3f `opt_Se' "      " %6.3f `Secorr'
display "    Spécificité    " %6.3f `app_Sp'  "       " %6.3f `opt_Sp' "      " %6.3f `Spcorr'
display "    Youden J       " %6.3f `app_J'   "       " %6.3f `opt_J'  "      " %6.3f `Jcorr'

display _newline "  --- Comparaison des trois modèles (AUC corrigées) ---"
display "    UMOD seul              (04)  ≈ 0.89"
display "    UMOD+age+sex+B2M       (05)  ≈ 0.91"
display "    UMOD + B2M  (parcimonieux)   = " %5.3f `AUCcorr'
display "=========================================================="
