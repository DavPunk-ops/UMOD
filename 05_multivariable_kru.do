* ===========================================================================
* 05_multivariable_kru.do
* Prérequis : lancer 02_calculate_kru.do AVANT ce do-file dans la même
*             session Stata. Le do-file 04_estimate_kru.do n'est pas
*             requis mais recommandé pour comparaison.
*
* Objectif : étendre les analyses du do-file 04 (UMOD seule) à un modèle
*            multivariable incluant UMOD + age + sex + β2-microglobuline.
*            Même structure : two-part pour KRU continu, puis logit
*            multivariable + cutoff Youden + bootstrap Harrell pour KRU≥2.
* ===========================================================================

* --- Vérification des prérequis ---
foreach v in kru_daugirdas_35 kru_naif_35 umod diuresis labcreatprehd ///
             age sex labb2mprehd {
    capture confirm variable `v'
    if _rc {
        display as error "ERREUR : variable '`v'' absente — lance d'abord 02_calculate_kru.do"
        exit 111
    }
}

display _newline "=== Prérequis OK — variables nécessaires présentes ==="

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

* --- Dummy sex pour compatibilité avec mfp (n'accepte pas i.sex) ---
* sex == 1 → femme, sex == 2 → homme (cf. assert dans 02_calculate_kru.do)
capture drop female
gen byte female = (sex == 1) if !missing(sex)
label variable female "Sexe féminin (1=F, 0=H)"

* ###########################################################################
* SECTION 1 — DESCRIPTION DES COVARIABLES (age, sex, β2-microglobuline)
* ###########################################################################

* ===========================================================================
*  1a. DISTRIBUTION GLOBALE
* ===========================================================================
display _newline(2) "=============================================="
display              "  1a. Distribution age / sex / B2M — N=151"
display              "=============================================="

summarize age, detail
summarize labb2mprehd, detail
tab sex, miss

count if missing(labb2mprehd)
display _newline "  B2M manquante : " r(N) " patients"
count if missing(age)
display "  Age manquant  : " r(N) " patients"
count if missing(sex)
display "  Sex manquant  : " r(N) " patients"

* ===========================================================================
*  1b. COVARIABLES CHEZ ANURIQUES vs NON-ANURIQUES
* ===========================================================================
display _newline(2) "=============================================="
display              "  1b. Covariables : anuriques vs non-anuriques"
display              "=============================================="

display _newline "  --- Age ---"
tabstat age, by(kru_pos) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum age, by(kru_pos)

display _newline "  --- B2M ---"
tabstat labb2mprehd, by(kru_pos) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum labb2mprehd, by(kru_pos)

display _newline "  --- Sex (% femmes) ---"
tab sex kru_pos, col chi2

* ===========================================================================
*  1c. COVARIABLES CHEZ KRU<2 vs KRU>=2 — POPULATION ENTIÈRE
* ===========================================================================
display _newline(2) "=============================================="
display              "  1c. Covariables : KRU<2 vs KRU>=2 (population entière)"
display              "=============================================="

display _newline "  --- Age ---"
tabstat age, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum age, by(kru_ge2)

display _newline "  --- B2M ---"
tabstat labb2mprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%6.1f)
ranksum labb2mprehd, by(kru_ge2)

display _newline "  --- Sex ---"
tab sex kru_ge2, col chi2

* ===========================================================================
*  1d. CORRÉLATIONS UMOD / B2M / âge avec KRU (non-anuriques)
* ===========================================================================
display _newline(2) "=============================================="
display              "  1d. Corrélations avec KRU (non-anuriques, N=89)"
display              "=============================================="

spearman kru_daugirdas_35 umod labb2mprehd age if kru_pos == 1, stats(rho p) star(0.05)

* ###########################################################################
* SECTION 2 — TWO-PART MODEL MULTIVARIABLE
*              Partie 1 : logit P(KRU>0) ~ UMOD + age + sex + B2M
*              Partie 2 : OLS  E[KRU | KRU>0] ~ idem  (avec MFP)
*              Comparaison avec le modèle UMOD seul (do-file 03).
* ###########################################################################

* ===========================================================================
*  2a. PARTIE 1 — Logit multivariable P(KRU>0)
* ===========================================================================
display _newline(2) "=============================================="
display              "  2a. Logit multivariable  P(KRU>0)"
display              "      UMOD + age + female + B2M"
display              "=============================================="

* Modèle UMOD seul (référence, comme do-file 03)
quietly logit kru_pos umod
estimates store mv_logit_ref
quietly lroc, nograph
local auc_logit_ref = r(area)

* Modèle multivariable
logit kru_pos umod age female labb2mprehd
estimates store mv_logit_full

display _newline "  --- Odds ratios ---"
logit kru_pos umod age female labb2mprehd, or

* Performance
quietly lroc, nograph
local auc_logit_full = r(area)

display _newline "  --- AUC partie 1 ---"
display "    UMOD seul       : AUC = " %5.3f `auc_logit_ref'
display "    Multivariable   : AUC = " %5.3f `auc_logit_full'
display "    Δ AUC           = " %6.3f (`auc_logit_full' - `auc_logit_ref')

* Calibration
estat gof, group(10) table

* Prédiction P(KRU>0 | X)
capture drop p_pos_mv
predict p_pos_mv, pr
label variable p_pos_mv "P(KRU>0 | UMOD,age,sex,B2M)"

* ===========================================================================
*  2b. PARTIE 2 — OLS multivariable E[KRU | KRU>0]
*       Avec MFP pour vérifier la linéarité des continues
* ===========================================================================
display _newline(2) "=============================================="
display              "  2b. E[KRU | KRU>0] multivariable   (N=89)"
display              "      OLS linéaire vs MFP en parallèle"
display              "=============================================="

* --- 2b.i  OLS multivariable (SE robustes Huber-White) ---
display _newline "  --- OLS multivariable (SE robustes) ---"
regress kru_daugirdas_35 umod age female labb2mprehd if kru_pos == 1, vce(robust)
estimates store mv_ols

local r2_ols    = e(r2)
local r2_a_ols  = e(r2_a)
local rmse_ols  = e(rmse)
quietly estat ic
matrix _ic_ols = r(S)
local aic_ols  = _ic_ols[1,5]
local bic_ols  = _ic_ols[1,6]

display _newline "  R² = " %5.3f `r2_ols' "   R²adj = " %5.3f `r2_a_ols' ///
    "   RMSE = " %5.3f `rmse_ols' "   AIC = " %6.2f `aic_ols'

* Rappel : OLS UMOD seul (do-file 03)
quietly regress kru_daugirdas_35 umod if kru_pos == 1, vce(robust)
local r2_ref   = e(r2)
local rmse_ref = e(rmse)
quietly estat ic
matrix _ic_ref = r(S)
local aic_ref  = _ic_ref[1,5]

display _newline "  --- Comparaison vs UMOD seul ---"
display "                     UMOD seul    Multivariable"
display "    R²            " %7.3f `r2_ref' "       " %7.3f `r2_ols'
display "    RMSE          " %7.3f `rmse_ref' "       " %7.3f `rmse_ols'
display "    AIC           " %7.2f `aic_ref' "       " %7.2f `aic_ols'
display "    ΔAIC (réf − mv) = " %5.2f (`aic_ref' - `aic_ols')
display "      > 2 : multivariable préférable"

* Prédiction conditionnelle multivariable
capture drop kru_cond_mv
predict kru_cond_mv, xb
label variable kru_cond_mv "E[KRU | KRU>0, X] — OLS mv"
replace kru_cond_mv = 0 if kru_cond_mv < 0

* Diagnostic résidus
capture drop resid_mv
predict resid_mv if e(sample), resid
display _newline "  Normalité des résidus (Shapiro-Wilk) :"
swilk resid_mv

* --- 2b.ii  MFP multivariable ---
display _newline(2) "  --- MFP : test de non-linéarité de UMOD, age, B2M ---"
display         "  (i.sex traité comme facteur, non transformé)"

mfp: regress kru_daugirdas_35 umod age female labb2mprehd if kru_pos == 1
estimates store mv_mfp

local r2_mfp   = e(r2)
local rmse_mfp = e(rmse)
quietly estat ic
matrix _ic_mfp = r(S)
local aic_mfp  = _ic_mfp[1,5]

display _newline "  R² MFP = " %5.3f `r2_mfp' ///
    "    RMSE MFP = " %5.3f `rmse_mfp' ///
    "    AIC MFP = " %6.2f `aic_mfp'

* --- 2b.iii  Comparaison OLS vs MFP ---
display _newline(2) "  --- Comparaison OLS multivariable vs MFP ---"
estimates stats mv_ols mv_mfp

display _newline "  ───────────────────────────────────────────"
display         "                    OLS         MFP"
display         "  ───────────────────────────────────────────"
display         "  R²              " %6.3f `r2_ols' "      " %6.3f `r2_mfp'
display         "  RMSE            " %6.3f `rmse_ols' "      " %6.3f `rmse_mfp'
display         "  AIC             " %6.2f `aic_ols' "    " %6.2f `aic_mfp'
display         "  ───────────────────────────────────────────"
display         "  ΔAIC (OLS − MFP) = " %5.2f (`aic_ols' - `aic_mfp')
display         "    > 2 : MFP préférable   |   < −2 : OLS préférable"
display         "    |ΔAIC| < 2 : équivalents → préférer le plus simple (OLS)"

* ===========================================================================
*  2c. COMBINAISON TWO-PART MULTIVARIABLE
* ===========================================================================
display _newline(2) "=============================================="
display              "  2c. Prédiction two-part multivariable"
display              "=============================================="

capture drop kru_pred_mv
gen double kru_pred_mv = p_pos_mv * kru_cond_mv if !missing(p_pos_mv, kru_cond_mv)
label variable kru_pred_mv "KRU prédit two-part multivariable"

* Performance globale
capture drop _r_sq _r_abs
gen double _r_sq  = (kru_daugirdas_35 - kru_pred_mv)^2
gen double _r_abs = abs(kru_daugirdas_35 - kru_pred_mv)

quietly summarize _r_sq
local rmse_mv = sqrt(r(mean))
quietly summarize _r_abs
local mae_mv = r(mean)
drop _r_sq _r_abs

display _newline "  RMSE two-part multivariable = " %5.3f `rmse_mv' " mL/min/35L"
display         "  MAE  two-part multivariable = " %5.3f `mae_mv'  " mL/min/35L"

quietly corr kru_daugirdas_35 kru_pred_mv
display "  Corrélation Pearson (observé,prédit) = " %5.3f r(rho)
quietly spearman kru_daugirdas_35 kru_pred_mv
display "  Corrélation Spearman                 = " %5.3f r(rho)

* ===========================================================================
*  2d. GRAPHIQUES DIAGNOSTIQUES
* ===========================================================================
* Observé vs prédit two-part multivariable
twoway (scatter kru_daugirdas_35 kru_pred_mv, msize(small)) ///
       (function y=x, range(0 10) lcolor(red) lpattern(dash)), ///
    title("KRU observé vs prédit — two-part multivariable") ///
    xtitle("KRU prédit (mL/min/35L)") ///
    ytitle("KRU observé (mL/min/35L)") ///
    legend(off) name(mv_obs_pred, replace)

* ###########################################################################
* SECTION 3 — PRÉDICTION DIRECTE DE KRU ≥ 2 (logit multivariable)
*              Pendant multivariable de la section 4 du do-file 03.
*
*    Choix méthodologique : analyse sur la POPULATION ENTIÈRE (N=151),
*    anuriques inclus — cohérent avec le do-file 03 et Wong et al.
*    (Kidney International 2015).
* ###########################################################################

* ===========================================================================
*  3a. LOGIT MULTIVARIABLE : kru_ge2 ~ UMOD + age + sex + B2M
* ===========================================================================
display _newline(2) "=============================================="
display              "  3a. Logit multivariable  P(KRU≥2)"
display              "      UMOD + age + female + B2M  (N=151)"
display              "=============================================="

logit kru_ge2 umod age female labb2mprehd
estimates store ge2_full

display _newline "  --- Odds ratios ---"
logit kru_ge2 umod age female labb2mprehd, or

quietly lroc, nograph
local auc_ge2_full = r(area)
display _newline "  AUC apparente (multivariable) = " %5.3f `auc_ge2_full'

* Rappel : UMOD seul (do-file 03)
quietly logit kru_ge2 umod
quietly lroc, nograph
local auc_ge2_ref = r(area)
display "  AUC apparente (UMOD seul)     = " %5.3f `auc_ge2_ref'
display "  Δ AUC                          = " %6.3f (`auc_ge2_full' - `auc_ge2_ref')

* Calibration multivariable
quietly logit kru_ge2 umod age female labb2mprehd
estat gof, group(10) table

* Prédiction P(KRU≥2 | X)
capture drop p_ge2_mv
predict p_ge2_mv, pr
label variable p_ge2_mv "P(KRU≥2 | UMOD,age,sex,B2M)"

* ===========================================================================
*  3b. TEST DE NON-LINÉARITÉ (MFP sur le logit kru_ge2)
* ===========================================================================
display _newline(2) "=============================================="
display              "  3b. MFP sur le logit multivariable"
display              "=============================================="

mfp: logit kru_ge2 umod age female labb2mprehd
estimates store ge2_mfp

quietly lroc, nograph
local auc_ge2_mfp = r(area)
display _newline "  AUC linéaire (mv) = " %5.3f `auc_ge2_full'
display         "  AUC MFP            = " %5.3f `auc_ge2_mfp'

display _newline "  Comparaison AIC :"
estimates stats ge2_full ge2_mfp

* On revient au modèle linéaire pour les sections suivantes
quietly estimates restore ge2_full

* ===========================================================================
*  3c. ROC + CUTOFF YOUDEN sur P(KRU≥2) prédite (multivariable)
* ===========================================================================
display _newline(2) "=============================================="
display              "  3c. ROC + cutoff Youden sur P(KRU≥2) — mv"
display              "=============================================="

roctab kru_ge2 p_ge2_mv, graph summary ///
    title("ROC : modèle multivariable prédit KRU≥2 (N=151)") ///
    name(roc_mv, replace)

* Comparaison statistique des AUC : UMOD seul vs multivariable
display _newline "  --- Comparaison des AUC (test de DeLong) ---"
roccomp kru_ge2 umod p_ge2_mv, graph summary ///
    name(roc_compare, replace)

* --- Cutoff Youden optimal sur P(KRU≥2) ---
display _newline "  --- Cutoff Youden sur P(KRU≥2) prédite ---"

local best_J_mv  = -1
local best_p_mv  = .
local best_Se_mv = .
local best_Sp_mv = .

forvalues p = 0.05(0.025)0.95 {
    quietly count if p_ge2_mv >= `p' & kru_ge2 == 1 & !missing(kru_ge2, p_ge2_mv)
    local TP = r(N)
    quietly count if p_ge2_mv <  `p' & kru_ge2 == 0 & !missing(kru_ge2, p_ge2_mv)
    local TN = r(N)
    quietly count if p_ge2_mv >= `p' & kru_ge2 == 0 & !missing(kru_ge2, p_ge2_mv)
    local FP = r(N)
    quietly count if p_ge2_mv <  `p' & kru_ge2 == 1 & !missing(kru_ge2, p_ge2_mv)
    local FN = r(N)
    if (`TP' + `FN') > 0 & (`TN' + `FP') > 0 {
        local Se = `TP' / (`TP' + `FN')
        local Sp = `TN' / (`TN' + `FP')
        local J  = `Se' + `Sp' - 1
        if `J' > `best_J_mv' {
            local best_J_mv  = `J'
            local best_p_mv  = `p'
            local best_Se_mv = `Se'
            local best_Sp_mv = `Sp'
        }
    }
}

display _newline "  → Cutoff optimal Youden : P(KRU≥2) ≥ " %5.3f `best_p_mv'
display         "      Se = " %5.1f 100*`best_Se_mv' " %"
display         "      Sp = " %5.1f 100*`best_Sp_mv' " %"
display         "      J  = " %5.3f `best_J_mv'

* ===========================================================================
*  3d. SYNTHÈSE
* ===========================================================================
display _newline(2) "========================================================"
display              "  SYNTHÈSE — Prédiction multivariable de KRU ≥ 2"
display              "         (population entière, N=151)"
display              "========================================================"
display "  AUC UMOD seul             = " %5.3f `auc_ge2_ref'
display "  AUC multivariable         = " %5.3f `auc_ge2_full'
display "  Gain d'AUC                = " %6.3f (`auc_ge2_full' - `auc_ge2_ref')
display ""
display "  Cutoff optimal Youden sur P(KRU≥2) :"
display "    Seuil = " %5.3f `best_p_mv'
display "    Se    = " %4.1f 100*`best_Se_mv' " %"
display "    Sp    = " %4.1f 100*`best_Sp_mv' " %"
display "    J     = " %5.3f `best_J_mv'
display "========================================================"

* ###########################################################################
* SECTION 4 — VALIDATION BOOTSTRAP (Harrell optimism-corrected)
*              du modèle multivariable final.
*
*    Procédure :
*      1. Modèle apparent : logit kru_ge2 ~ UMOD+age+sex+B2M (N=151)
*         → AUC apparente + cutoff Youden sur P̂ apparente
*      2. Pour b=1..B :
*         a) Tirer un échantillon bootstrap (N=151 avec remise)
*         b) Refit du logit multivariable dans le bootstrap
*         c) Prédire P̂_b dans le bootstrap → AUC_bb, cutoff c_b, Se_bb, Sp_bb
*         d) Prédire P̂_b dans l'ORIGINAL → AUC_bo, Se_bo, Sp_bo
*            (en appliquant c_b)
*      3. Optimisme = mean(metric_bb − metric_bo)
*      4. Corrigé = apparent − optimisme
* ###########################################################################

* --- Performance apparente ---
local app_AUC    = `auc_ge2_full'
local app_p_cut  = `best_p_mv'
local app_Se     = `best_Se_mv'
local app_Sp     = `best_Sp_mv'
local app_J      = `best_J_mv'

* --- Configuration bootstrap ---
global B = 1000
global SEED = 20260522

display _newline(2) "=============================================="
display              "  4. Validation bootstrap multivariable"
display              "     B=$B, seed=$SEED"
display              "=============================================="
display _newline "  Performance apparente :"
display "    AUC            = " %5.3f `app_AUC'
display "    Cutoff P̂ opt.  = " %5.3f `app_p_cut'
display "    Se / Sp / J    = " %5.3f `app_Se' " / " %5.3f `app_Sp' " / " %5.3f `app_J'

set seed $SEED

* On filtre les obs avec toutes les variables nécessaires
preserve
quietly keep if !missing(kru_ge2, umod, age, sex, labb2mprehd)
quietly count
local N = r(N)
restore

display _newline "  N effectif (sans missing) = `N'"

* --- Stockage via postfile ---
tempname memh
tempfile bootres
postfile `memh' double(AUC_bb AUC_bo cutoff_p Jbb Jbo Sebb Sebo Spbb Spbo) ///
    using `bootres', replace

local n_valid = 0

display _newline _continue "  Progression : "

forvalues b = 1/$B {
    if mod(`b', 100) == 0 display _continue "`b' "

    preserve
    quietly keep if !missing(kru_ge2, umod, age, sex, labb2mprehd)
    quietly bsample

    * (a) Fit du modèle dans le bootstrap
    capture quietly logit kru_ge2 umod age female labb2mprehd
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

    * (d) Appliquer le modèle bootstrap à l'ORIGINAL
    *     Les résultats e() du dernier logit persistent à travers restore,
    *     donc predict applique les coefficients bootstrap aux données
    *     originales (predict gère correctement i.sex).
    capture drop p_bo
    quietly predict p_bo if !missing(kru_ge2, umod, age, sex, labb2mprehd), pr

    * AUC_bo : AUC du modèle bootstrap appliqué à l'original
    capture quietly roctab kru_ge2 p_bo, nograph
    if _rc {
        drop p_bo
        continue
    }
    local auc_bo = r(area)

    * Application du cutoff c_b à l'original
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

* --- Agrégation des résultats bootstrap ---
preserve
quietly use `bootres', clear

* AUC : optimisme par paires bb-bo
quietly summarize AUC_bb, meanonly
local m_auc_bb = r(mean)
quietly summarize AUC_bo, meanonly
local m_auc_bo = r(mean)
local opt_AUC  = `m_auc_bb' - `m_auc_bo'

quietly _pctile AUC_bo, percentiles(2.5 97.5)
local auc_lo = r(r1)
local auc_hi = r(r2)

* Youden
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

* --- Métriques corrigées ---
local AUCcorr = `app_AUC' - `opt_AUC'
local Jcorr   = `app_J'   - `opt_J'
local Secorr  = `app_Se'  - `opt_Se'
local Spcorr  = `app_Sp'  - `opt_Sp'

display _newline(2) "=========================================================="
display              "  RÉSULTATS BOOTSTRAP — MODÈLE MULTIVARIABLE"
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

display _newline "  Interprétation :"
display "  - AUC corrigée vs UMOD seul (do-file 04, AUC corr. ≈ 0.89) :"
display "      gain net après correction d'optimisme"
display "  - Cutoff sur P̂  → la stabilité dépend du nombre de prédicteurs"
display "  - J corrigé    → performance attendue dans une nouvelle cohorte"
display "                   (validation externe encore nécessaire)"
display "=========================================================="
