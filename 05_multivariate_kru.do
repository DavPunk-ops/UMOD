* ===========================================================================
* 05_multivariate_kru.do
* Prérequis : lancer 02_calculate_kru.do AVANT ce do-file.
*
* Objectif : modèle multivarié UMOD + B2M.
*
* Section 1 (préliminaire) teste formellement l'exclusion de l'âge et du sexe.
* Résultat : p > 0.05 dans les trois composantes → modèle parcimonieux UMOD + B2M.
* Sections 2–4 : logit direct P(KRU≥2), bootstrap Harrell, stratégie deux seuils.
*                → manuscrit : "Diagnostic performance of combined serum uromodulin
*                               and β2-microglobulin for predicting KRU ≥2 mL/min/35L"
* Sections 5–6 : corrélations et modèle two-part continu (relation quantitative).
*                → manuscrit : "Quantitative relationship between serum uromodulin,
*                               β2-microglobulin, and KRU"
* Transformation B2M^-2 dans l'OLS confirmée par ΔAIC (section 6b).
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

* --- Dummy sex (recréé si absent) — requis pour les rappels comparatifs ---
* avec le modèle complet UMOD+age+sex+B2M. sex==1 → femme (cf. 02_calculate_kru.do)
capture confirm variable female
if _rc {
    capture drop female
    gen byte female = (sex == 1) if !missing(sex)
    label variable female "Sexe féminin (1=F, 0=H)"
}

* --- Transformation B2M^-2 confirmée par ΔAIC (section 6b ci-dessous) ---
* Échelle X = B2M/10 (cf. MFP : "where: X = labb2mprehd/10")
* On NE centre PAS : seul l'intercept changerait, pas les performances.
capture drop b2m_neg2
gen double b2m_neg2 = (labb2mprehd/10)^(-2) if !missing(labb2mprehd)
label variable b2m_neg2 "(B2M/10)^-2"

summarize labb2mprehd b2m_neg2

* ###########################################################################
* SECTION 1 — JUSTIFICATION FORMELLE : EXCLUSION DE L'ÂGE ET DU SEXE
*    Trois régressions avec le modèle complet (UMOD + âge + sexe + β2M).
*    Critère : p > 0.05 pour âge ET sexe dans les trois composantes
*              → exclusion justifiée dans le modèle parcimonieux UMOD + B2M.
* ###########################################################################

display _newline(2) "=============================================="
display              "  1. Test formel âge/sexe — modèle complet"
display              "     UMOD + age + female + B2M"
display              "=============================================="

* (1) Logit P(KRU>0) ~ UMOD + age + female + B2M
display _newline "  --- (1) Logit P(KRU>0) ~ UMOD + age + female + B2M ---"
quietly logit kru_pos umod age female labb2mprehd
local p_age_1 = 2*normal(-abs(_b[age]   /_se[age]))
local p_sex_1 = 2*normal(-abs(_b[female]/_se[female]))

* (2) OLS E[KRU|KRU>0] ~ UMOD + age + female + B2M
display "  --- (2) OLS E[KRU|KRU>0] ~ UMOD + age + female + B2M ---"
quietly regress kru_daugirdas_35 umod age female labb2mprehd if kru_pos == 1, vce(robust)
local p_age_2 = 2*ttail(e(df_r), abs(_b[age]   /_se[age]))
local p_sex_2 = 2*ttail(e(df_r), abs(_b[female]/_se[female]))

* (3) Logit P(KRU≥2) ~ UMOD + age + female + B2M
display "  --- (3) Logit P(KRU≥2) ~ UMOD + age + female + B2M ---"
quietly logit kru_ge2 umod age female labb2mprehd
local p_age_3 = 2*normal(-abs(_b[age]   /_se[age]))
local p_sex_3 = 2*normal(-abs(_b[female]/_se[female]))

display _newline "  ─────────────────────────────────────────────────────"
display         "                             p (âge)     p (sexe)"
display         "  ─────────────────────────────────────────────────────"
display         "  (1) Logit P(KRU>0)       " %7.3f `p_age_1' "      " %7.3f `p_sex_1'
display         "  (2) OLS  E[KRU | KRU>0]  " %7.3f `p_age_2' "      " %7.3f `p_sex_2'
display         "  (3) Logit P(KRU≥2)       " %7.3f `p_age_3' "      " %7.3f `p_sex_3'
display         "  ─────────────────────────────────────────────────────"
display _newline "  → Âge et sexe exclus du modèle parcimonieux (p > 0.05 dans les 3 composantes)."

* ###########################################################################
* SECTION 2 — LOGIT DIRECT P(KRU≥2) ~ UMOD + B2M (N=148)
*    MFP du do-file 04_univariate_kru a confirmé que le linéaire est optimal
*    pour ce logit.
* ###########################################################################

* ===========================================================================
*  2a. LOGIT + ROC
* ===========================================================================
display _newline(2) "=============================================="
display              "  2a. Logit P(KRU≥2) ~ UMOD + B2M (N=148)"
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
*  2b. CUTOFF YOUDEN sur P(KRU≥2)
* ===========================================================================
display _newline(2) "=============================================="
display              "  2b. Cutoff Youden sur P̂"
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
* SECTION 3 — BOOTSTRAP HARRELL DU MODÈLE PARCIMONIEUX FINAL
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
display              "  3. Bootstrap parcimonieux (B=$B, seed=$SEED)"
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

* ###########################################################################
* SECTION 4 — STRATÉGIE À DEUX SEUILS (rule-out / rule-in) — UMOD + B2M
*
*    Le « score » est la PROBABILITÉ PRÉDITE p_ge2_ps du logit direct.
*
*       pc_out (rule-out) : seuil MAX de P̂ tel que Se ≥ 90%
*                           → P̂ < pc_out  →  KRU ≥2 exclu  (NPV élevée)
*       pc_in  (rule-in)  : seuil MIN de P̂ tel que Sp ≥ 90%
*                           → P̂ ≥ pc_in   →  KRU ≥2 confirmé (PPV élevée)
*
*    Trois zones :
*       [0       ; pc_out)  : rule-out  → pas de collecte urinaire
*       [pc_out  ; pc_in)   : grey zone → collecte urinaire indiquée
*       [pc_in   ; 1]       : rule-in   → pas de collecte urinaire
* ###########################################################################

* ===========================================================================
*  4a. IDENTIFICATION DES DEUX SEUILS — population entière (N=148)
* ===========================================================================
display _newline(2) "=============================================="
display              "  4a. Two-cutoff (UMOD+B2M) — cibles Se≥90% / Sp≥90%"
display              "=============================================="

* Refit pour s'assurer que p_ge2_ps reflète le modèle final
quietly logit kru_ge2 umod labb2mprehd
capture drop p_ge2_ps
quietly predict p_ge2_ps, pr

local target_se = 0.90
local target_sp = 0.90

local pc_out = .
local pc_in  = .

forvalues p = 0.01(0.01)0.99 {
    quietly count if p_ge2_ps >= `p' & kru_ge2 == 1 & !missing(kru_ge2, p_ge2_ps)
    local TP = r(N)
    quietly count if p_ge2_ps <  `p' & kru_ge2 == 0 & !missing(kru_ge2, p_ge2_ps)
    local TN = r(N)
    quietly count if p_ge2_ps >= `p' & kru_ge2 == 0 & !missing(kru_ge2, p_ge2_ps)
    local FP = r(N)
    quietly count if p_ge2_ps <  `p' & kru_ge2 == 1 & !missing(kru_ge2, p_ge2_ps)
    local FN = r(N)
    if (`TP' + `FN') > 0 & (`TN' + `FP') > 0 {
        local Se_ = `TP' / (`TP' + `FN')
        local Sp_ = `TN' / (`TN' + `FP')
        * rule-out : on garde le seuil le PLUS GRAND avec Se >= 90%
        if `Se_' >= `target_se' local pc_out = `p'
        * rule-in : on garde le seuil le PLUS PETIT avec Sp >= 90%
        if `Sp_' >= `target_sp' & missing(`pc_in') local pc_in = `p'
    }
}

display _newline "  → Rule-out cutoff (Se ≥ 90%) : P̂(KRU≥2) < " %5.3f `pc_out'
display         "  → Rule-in  cutoff (Sp ≥ 90%) : P̂(KRU≥2) ≥ " %5.3f `pc_in'

* ===========================================================================
*  4b. PERFORMANCE DES TROIS ZONES (apparent, N=148)
* ===========================================================================
display _newline(2) "=============================================="
display              "  4b. Performance des trois zones (UMOD+B2M)"
display              "=============================================="

quietly count if !missing(kru_ge2, p_ge2_ps)
local N_tot = r(N)

* --- Zone rule-out : P̂ < pc_out ---
* (!missing(p_ge2_ps) obligatoire : en Stata, . est traité comme +∞,
*  donc une valeur manquante satisferait à tort la condition >= pc_in)
quietly count if p_ge2_ps < `pc_out' & !missing(kru_ge2, p_ge2_ps)
local N_out = r(N)
quietly count if p_ge2_ps < `pc_out' & kru_ge2 == 0 & !missing(kru_ge2, p_ge2_ps)
local TN_out = r(N)
quietly count if p_ge2_ps < `pc_out' & kru_ge2 == 1 & !missing(kru_ge2, p_ge2_ps)
local FN_out = r(N)
local NPV = cond(`N_out' > 0, `TN_out' / `N_out', .)

* --- Zone rule-in : P̂ >= pc_in ---
quietly count if p_ge2_ps >= `pc_in' & !missing(kru_ge2, p_ge2_ps)
local N_in = r(N)
quietly count if p_ge2_ps >= `pc_in' & kru_ge2 == 1 & !missing(kru_ge2, p_ge2_ps)
local TP_in = r(N)
quietly count if p_ge2_ps >= `pc_in' & kru_ge2 == 0 & !missing(kru_ge2, p_ge2_ps)
local FP_in = r(N)
local PPV = cond(`N_in' > 0, `TP_in' / `N_in', .)

* --- Zone grise : pc_out <= P̂ < pc_in ---
quietly count if p_ge2_ps >= `pc_out' & p_ge2_ps < `pc_in' & !missing(kru_ge2, p_ge2_ps)
local N_grey = r(N)
quietly count if p_ge2_ps >= `pc_out' & p_ge2_ps < `pc_in' & kru_ge2 == 1 & !missing(kru_ge2, p_ge2_ps)
local KRUge2_grey = r(N)
quietly count if p_ge2_ps >= `pc_out' & p_ge2_ps < `pc_in' & kru_ge2 == 0 & !missing(kru_ge2, p_ge2_ps)
local KRUlt2_grey = r(N)

local pct_out  = 100*`N_out' /`N_tot'
local pct_grey = 100*`N_grey'/`N_tot'
local pct_in   = 100*`N_in'  /`N_tot'
local pct_class = `pct_out' + `pct_in'

display _newline "  Zone RULE-OUT  (P̂ < " %5.3f `pc_out' ")"
display         "    N = `N_out' (" %4.1f `pct_out' "%)"
display         "    TN = `TN_out' (KRU<2 correctement exclus)"
display         "    FN = `FN_out' (KRU≥2 manqués)"
display         "    NPV = " %5.1f 100*`NPV' " %"

display _newline "  Zone GREY      (" %5.3f `pc_out' " ≤ P̂ < " %5.3f `pc_in' ")"
display         "    N = `N_grey' (" %4.1f `pct_grey' "%)"
display         "    KRU<2 = `KRUlt2_grey' / KRU≥2 = `KRUge2_grey'"
display         "    → collecte urinaire indiquée"

display _newline "  Zone RULE-IN   (P̂ ≥ " %5.3f `pc_in' ")"
display         "    N = `N_in' (" %4.1f `pct_in' "%)"
display         "    TP = `TP_in' (KRU≥2 correctement confirmés)"
display         "    FP = `FP_in' (KRU<2 mal classés en KRU≥2)"
display         "    PPV = " %5.1f 100*`PPV' " %"

display _newline "  → Patients classifiés (hors zone grise) : " %4.1f `pct_class' " %"
display         "  → Collecte urinaire évitable             : " %4.1f `pct_class' " %"

local app_pc_out = `pc_out'
local app_pc_in  = `pc_in'
local app_NPV    = `NPV'
local app_PPV    = `PPV'
local app_pct_grey = `pct_grey'

* ===========================================================================
*  4c. VALIDATION BOOTSTRAP — stabilité des seuils + correction d'optimisme
*       Le logit est refité dans chaque échantillon bootstrap (Harrell).
* ===========================================================================
display _newline(2) "=============================================="
display              "  4c. Validation bootstrap two-cutoff (B=$B, seed=$SEED)"
display              "=============================================="

set seed $SEED

tempname memh2
tempfile bootres2
postfile `memh2' double(pc_out_b pc_in_b NPV_bb NPV_bo PPV_bb PPV_bo pct_grey_b) ///
    using `bootres2', replace

local n_valid2 = 0

display _newline _continue "  Progression : "

forvalues b = 1/$B {
    if mod(`b', 100) == 0 display _continue "`b' "

    preserve
    quietly keep if !missing(kru_ge2, umod, labb2mprehd)
    quietly bsample

    * (a) Refit du logit bi-marqueur dans le bootstrap
    capture quietly logit kru_ge2 umod labb2mprehd
    if _rc {
        restore
        continue
    }
    capture drop p_b
    quietly predict p_b, pr

    * (b) Recherche des deux seuils sur P̂ dans le bootstrap
    local pc_out_b = .
    local pc_in_b  = .

    forvalues p = 0.01(0.01)0.99 {
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
            if `Se_' >= 0.90 local pc_out_b = `p'
            if `Sp_' >= 0.90 & missing(`pc_in_b') local pc_in_b = `p'
        }
    }

    * (c) NPV/PPV apparents dans le bootstrap (optimistes : bb)
    if !missing(`pc_out_b') {
        quietly count if p_b < `pc_out_b'
        local N_o_b = r(N)
        quietly count if p_b < `pc_out_b' & kru_ge2 == 0
        local TN_o_b = r(N)
        local NPV_bb = cond(`N_o_b' > 0, `TN_o_b'/`N_o_b', .)
    }
    else local NPV_bb = .

    if !missing(`pc_in_b') {
        quietly count if p_b >= `pc_in_b'
        local N_i_b = r(N)
        quietly count if p_b >= `pc_in_b' & kru_ge2 == 1
        local TP_i_b = r(N)
        local PPV_bb = cond(`N_i_b' > 0, `TP_i_b'/`N_i_b', .)
    }
    else local PPV_bb = .

    if !missing(`pc_out_b', `pc_in_b') {
        quietly count if p_b >= `pc_out_b' & p_b < `pc_in_b'
        local Ng_b = r(N)
        quietly count
        local Nb = r(N)
        local pct_grey_b = 100*`Ng_b'/`Nb'
    }
    else local pct_grey_b = .

    restore

    if missing(`pc_out_b') | missing(`pc_in_b') continue

    * (d) Appliquer le modèle bootstrap (encore en e()) à l'ORIGINAL → P̂_bo
    *     puis NPV/PPV honnêtes (bo) avec les seuils bootstrap sur l'original
    capture drop p_bo
    quietly predict p_bo if !missing(kru_ge2, umod, labb2mprehd), pr

    quietly count if p_bo < `pc_out_b' & !missing(kru_ge2, p_bo)
    local N_o_orig = r(N)
    quietly count if p_bo < `pc_out_b' & kru_ge2 == 0 & !missing(kru_ge2, p_bo)
    local TN_o_orig = r(N)
    local NPV_bo = cond(`N_o_orig' > 0, `TN_o_orig'/`N_o_orig', .)

    quietly count if p_bo >= `pc_in_b' & !missing(kru_ge2, p_bo)
    local N_i_orig = r(N)
    quietly count if p_bo >= `pc_in_b' & kru_ge2 == 1 & !missing(kru_ge2, p_bo)
    local TP_i_orig = r(N)
    local PPV_bo = cond(`N_i_orig' > 0, `TP_i_orig'/`N_i_orig', .)

    drop p_bo

    post `memh2' (`pc_out_b') (`pc_in_b') (`NPV_bb') (`NPV_bo') ///
                 (`PPV_bb') (`PPV_bo') (`pct_grey_b')
    local n_valid2 = `n_valid2' + 1
}

postclose `memh2'

display ""
display _newline "  Itérations valides : `n_valid2'/$B"

* --- Agrégation ---
preserve
quietly use `bootres2', clear

quietly summarize pc_out_b, detail
local m_cout  = r(mean)
local sd_cout = r(sd)
local cout_lo = r(p5)
local cout_hi = r(p95)

quietly summarize pc_in_b, detail
local m_cin  = r(mean)
local sd_cin = r(sd)
local cin_lo = r(p5)
local cin_hi = r(p95)

quietly summarize NPV_bb, meanonly
local m_NPVbb = r(mean)
quietly summarize NPV_bo, meanonly
local m_NPVbo = r(mean)
local opt_NPV = `m_NPVbb' - `m_NPVbo'

quietly summarize PPV_bb, meanonly
local m_PPVbb = r(mean)
quietly summarize PPV_bo, meanonly
local m_PPVbo = r(mean)
local opt_PPV = `m_PPVbb' - `m_PPVbo'

quietly summarize pct_grey_b, detail
local m_grey  = r(mean)
local grey_lo = r(p5)
local grey_hi = r(p95)

quietly _pctile NPV_bo, percentiles(2.5 97.5)
local NPV_lo = r(r1)
local NPV_hi = r(r2)
quietly _pctile PPV_bo, percentiles(2.5 97.5)
local PPV_lo = r(r1)
local PPV_hi = r(r2)

restore

local NPV_corr = `app_NPV' - `opt_NPV'
local PPV_corr = `app_PPV' - `opt_PPV'

* ===========================================================================
*  4d. SYNTHÈSE + TABLE 3 + FIGURE 3
* ===========================================================================
display _newline(2) "=========================================================="
display              "  SYNTHÈSE — Two-cutoff UMOD+B2M (Se≥90% / Sp≥90%)"
display              "=========================================================="

display _newline "  --- Seuils de probabilité (apparent & bootstrap) ---"
display "    Rule-out  : apparent = " %5.3f `app_pc_out'
display "                bootstrap moy = " %5.3f `m_cout' " (SD " %5.3f `sd_cout' ")"
display "                IC bootstrap (5–95%) = " %5.3f `cout_lo' " — " %5.3f `cout_hi'

display _newline "    Rule-in   : apparent = " %5.3f `app_pc_in'
display "                bootstrap moy = " %5.3f `m_cin' " (SD " %5.3f `sd_cin' ")"
display "                IC bootstrap (5–95%) = " %5.3f `cin_lo' " — " %5.3f `cin_hi'

display _newline "  --- Performance (corrigée pour optimisme) ---"
display "                          Apparent    Optimisme    Corrigé"
display "    NPV (rule-out)   " %6.3f `app_NPV' "      " %6.3f `opt_NPV' "      " %6.3f `NPV_corr'
display "    PPV (rule-in)    " %6.3f `app_PPV' "      " %6.3f `opt_PPV' "      " %6.3f `PPV_corr'

display _newline "  --- Zone grise ---"
display "    % apparent          = " %4.1f `app_pct_grey' " %"
display "    % bootstrap moyenne = " %4.1f `m_grey' " %"
display "    IC bootstrap (5–95%) = " %4.1f `grey_lo' " — " %4.1f `grey_hi' " %"

display _newline "  --- Comparaison avec UMOD seul (do-file 04_univariate_kru, section 6) ---"
display "    UMOD seul     : grey zone ≈ 24.5%, NPV 97.1%, PPV 84.4%"
display "    UMOD + B2M    : grey zone = " %4.1f `app_pct_grey' "%, NPV " %4.1f 100*`app_NPV' "%, PPV " %4.1f 100*`app_PPV' "%"
display "=========================================================="

display _newline(2) "=================================================================="
display              "  TABLE 3 — Two-cutoff strategy: UMOD + β2M"
display              "  (corrected for optimism by Harrell bootstrap, B=$B)"
display              "=================================================================="
display "  Rule-out cutoff : P̂(KRU≥2) < " %5.3f `app_pc_out' ///
    "   (IC 5–95%: " %5.3f `cout_lo' "–" %5.3f `cout_hi' ")"
display "  Rule-in  cutoff : P̂(KRU≥2) ≥ " %5.3f `app_pc_in' ///
    "   (IC 5–95%: " %5.3f `cin_lo' "–" %5.3f `cin_hi' ")"
display _newline "  ──────────────────────────────────────────────────────────────────"
display           "  Zone             N   (%)   Metric   Apparent   Corrected   95% CI"
display           "  ──────────────────────────────────────────────────────────────────"
display "  Rule-out        " %3.0f `N_out' "  (" %4.1f `pct_out'  "%)   NPV      " ///
    %5.1f 100*`app_NPV'  "%       " %5.1f 100*`NPV_corr' "%   (" ///
    %4.1f 100*`NPV_lo' "–" %4.1f 100*`NPV_hi' "%)"
display "  Grey zone       " %3.0f `N_grey' "  (" %4.1f `pct_grey' "%)    —        —           —          —"
display "  Rule-in         " %3.0f `N_in'   "  (" %4.1f `pct_in'   "%)   PPV      " ///
    %5.1f 100*`app_PPV'  "%       " %5.1f 100*`PPV_corr' "%   (" ///
    %4.1f 100*`PPV_lo' "–" %4.1f 100*`PPV_hi' "%)"
display  "  ──────────────────────────────────────────────────────────────────"
display  "  Classified (rule-out + rule-in) : " %4.1f `pct_class' "% of patients"
display  "=================================================================="

* Figure 3 — Strip plot horizontal
preserve
keep if !missing(kru_ge2, p_ge2_ps)

set seed 20260522

gen double _y  = cond(kru_ge2==0, 1, 2)
gen double _yj = _y + (runiform()-0.5)*0.5

local cout = `app_pc_out'
local cin  = `app_pc_in'

local x_out  = `cout' / 2
local x_grey = (`cout' + `cin') / 2
local x_in   = (`cin' + 1) / 2

local lbl_npv = "(NPV " + string(round(100*`app_NPV', 0.1), "%4.1f") + "%)"
local lbl_ppv = "(PPV " + string(round(100*`app_PPV', 0.1), "%4.1f") + "%)"
local lbl_grey = string(round(`app_pct_grey', 0.1), "%3.1f") + "% of patients"

twoway ///
    (scatter _yj p_ge2_ps if kru_ge2==0, ///
        mcolor(navy%45) msize(small) msymbol(circle)) ///
    (scatter _yj p_ge2_ps if kru_ge2==1, ///
        mcolor(cranberry%45) msize(small) msymbol(circle)) ///
    , ///
    xline(`cout', lpattern(dash) lcolor(black) lwidth(medthick)) ///
    xline(`cin',  lpattern(dash) lcolor(black) lwidth(medthick)) ///
    xlabel(0(0.2)1, labsize(medium)) ///
    ylabel(1 `""KRU <2" "mL/min/35L""' 2 `""KRU ≥2" "mL/min/35L""', ///
        noticks labsize(small) angle(0)) ///
    xtitle("Predicted probability of KRU ≥2 (UMOD + {&beta}2M)", size(medlarge)) ///
    ytitle("") ///
    yscale(range(0.3 3.1)) ///
    text(2.95 `x_out'  "Rule-out",   size(small)  just(center) color(black)) ///
    text(2.80 `x_out'  "`lbl_npv'",  size(vsmall) just(center) color(black)) ///
    text(2.95 `x_grey' "Grey zone",   size(small)  just(center) color(black)) ///
    text(2.80 `x_grey' "`lbl_grey'",  size(vsmall) just(center) color(black)) ///
    text(2.95 `x_in'   "Rule-in",    size(small)  just(center) color(black)) ///
    text(2.80 `x_in'   "`lbl_ppv'",  size(vsmall) just(center) color(black)) ///
    legend(off) ///
    graphregion(color(white)) plotregion(color(white)) ///
    xsize(8) ysize(4) ///
    name(fig_twocut_mv, replace)

graph export "Figure3_twocutoff_UMOD_B2M.tif", replace width(2400)

restore

* ###########################################################################
* SECTION 5 — CORRÉLATIONS UMOD / B2M / KRU (non-anuriques)
* ###########################################################################
display _newline(2) "=============================================="
display              "  5. Corrélations UMOD / B2M / KRU (non-anuriques)"
display              "=============================================="

spearman kru_daugirdas_35 umod labb2mprehd b2m_neg2 if kru_pos == 1, ///
    stats(rho p) star(0.05)

* ###########################################################################
* SECTION 6 — TWO-PART PARCIMONIEUX
*    Partie 1 : logit P(KRU>0) ~ UMOD + B2M  (linéaire)
*    Partie 2 : OLS  E[KRU|KRU>0] ~ UMOD + B2M^-2  (transformation MFP)
* ###########################################################################

* ===========================================================================
*  6a. PARTIE 1 — Logit P(KRU>0) ~ UMOD + B2M
* ===========================================================================
display _newline(2) "=============================================="
display              "  6a. Logit P(KRU>0)  ~  UMOD + B2M"
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
*  6b. PARTIE 2 — OLS E[KRU | KRU>0] ~ UMOD + B2M^-2
* ===========================================================================
display _newline(2) "=============================================="
display              "  6b. OLS E[KRU | KRU>0] ~ UMOD + (B2M/10)^-2"
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
*  6c. COMBINAISON TWO-PART
* ===========================================================================
display _newline(2) "=============================================="
display              "  6c. Prédiction two-part parcimonieuse"
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
*  6d. COMPARAISON : two-part continu vs logit direct pour discriminer KRU≥2
*       Les deux scores (kru_pred_ps et p_ge2_ps) sont comparés par DeLong.
* ===========================================================================
display _newline(2) "=============================================="
display              "  6d. Two-part continu vs logit direct"
display              "      pour discriminer KRU≥2 (DeLong)"
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
