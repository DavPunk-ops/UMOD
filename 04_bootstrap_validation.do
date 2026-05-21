* ===========================================================================
* 04_bootstrap_validation.do
* Objectif : Validation bootstrap interne (Harrell) — correction de l'optimisme
*            de l'AUC et estimation de la pente de calibration pour les
*            modèles de prédiction de KRU≥2 chez les non-anuriques.
*
* Méthode : optimism-corrected bootstrap (Harrell, RMS book, §5.3.4)
*   1. AUC apparente : modèle ajusté ET évalué sur l'échantillon original
*   2. Pour B itérations :
*      a) Tirer un échantillon bootstrap (avec remise, taille N)
*      b) Ajuster le modèle sur le bootstrap → coefficients b*
*      c) AUC_bb = AUC(b*) appliqué au bootstrap (optimiste)
*      d) AUC_bo = AUC(b*) appliqué à l'original (honnête)
*      e) Pente_b = pente de logit(y ~ lp_b*) dans l'original
*   3. Optimism = moyenne(AUC_bb − AUC_bo)
*   4. AUC corrigée = AUC apparente − Optimism
*   5. Pente corrigée = moyenne(Pente_b) ; idéal = 1.0
*
* Population : non-anuriques (kru_pos == 1)
* Modèles évalués :
*   M1. UMOD seul                    (référence,    section 7)
*   M2. UMOD + B2M                   (parcimonieux, section 9)
*   M3. UMOD + créat + urée + B2M    (complet,      section 8)
*
* Prérequis : exécuter 03_estimate_kru.do AVANT ce do-file (variables
*             kru_ge2, kru_pos, umod, labb2mprehd, labcreatprehd,
*             labureaprehd doivent exister en mémoire).
* ===========================================================================

* --- Vérifications de prérequis ---
foreach v in kru_ge2 kru_pos umod labb2mprehd labcreatprehd labureaprehd {
    capture confirm variable `v'
    if _rc {
        display as error "ERREUR : variable '`v'' absente — lance d'abord 03_estimate_kru.do"
        exit 111
    }
}

* --- Configuration ---
local B = 1000           // nombre d'itérations bootstrap
local SEED = 20260521    // graine pour reproductibilité

* ===========================================================================
*  PROGRAMME : bootval
*    Validation bootstrap d'optimisme pour un modèle logit.
*    Syntaxe : bootval depvar indepvars [if], b(B) seed(S)
*    Retourne (r()) :
*      auc_apparent, optimism, optimism_sd, auc_corrected,
*      auc_bo_mean, auc_bo_lo, auc_bo_hi, slope_corrected,
*      n_valid, B, N
* ===========================================================================

capture program drop bootval
program define bootval, rclass
    syntax varlist(min=2) [if] [, B(integer 500) Seed(integer 42)]

    marksample touse
    gettoken depvar indepvars : varlist

    * --- 1. Performance apparente ---
    quietly logit `depvar' `indepvars' if `touse'
    if !e(converged) {
        display as error "Modèle apparent non-convergé"
        exit 198
    }

    capture drop _lp_app
    quietly predict double _lp_app if `touse', xb
    quietly roctab `depvar' _lp_app if `touse', nograph
    local auc_app = r(area)
    drop _lp_app

    * --- 2. Boucle bootstrap ---
    set seed `seed'
    quietly count if `touse'
    local N = r(N)

    display _newline "  Bootstrap : B=`B' itérations sur N=`N' patients"
    display _newline _continue "  Progression : "

    * Stockage des résultats par itération
    tempname Mres
    matrix `Mres' = J(`B', 4, .)
    matrix colnames `Mres' = AUC_bb AUC_bo Optimism Slope

    local n_valid = 0

    forvalues b = 1/`B' {
        if mod(`b', 100) == 0 display _continue "`b' "

        preserve
        quietly keep if `touse'
        quietly bsample

        capture quietly logit `depvar' `indepvars'
        local fit_ok = (_rc == 0)
        if `fit_ok' local fit_ok = e(converged)

        if !`fit_ok' {
            restore
            continue
        }

        matrix bmat = e(b)

        capture drop _lp_bb
        quietly predict double _lp_bb, xb
        quietly roctab `depvar' _lp_bb, nograph
        local auc_bb = r(area)

        restore

        capture drop _lp_bo
        matrix score double _lp_bo = bmat if `touse'
        quietly roctab `depvar' _lp_bo if `touse', nograph
        local auc_bo = r(area)

        * Pente de calibration : régresser y sur lp_bo (original)
        capture quietly logit `depvar' _lp_bo if `touse'
        local slope_ok = (_rc == 0)
        if `slope_ok' local slope_ok = e(converged)

        local slope = .
        if `slope_ok' local slope = _b[_lp_bo]

        drop _lp_bo

        local opt = `auc_bb' - `auc_bo'
        local n_valid = `n_valid' + 1

        matrix `Mres'[`b', 1] = `auc_bb'
        matrix `Mres'[`b', 2] = `auc_bo'
        matrix `Mres'[`b', 3] = `opt'
        matrix `Mres'[`b', 4] = `slope'
    }

    display ""

    if `n_valid' == 0 {
        display as error "Aucune itération bootstrap valide"
        exit 198
    }

    * --- 3. Synthèse via svmat dans un dataset temporaire ---
    preserve
    quietly drop _all
    quietly svmat double `Mres', names(__br)

    quietly summarize __br1, meanonly
    local m_auc_bb = r(mean)

    quietly summarize __br2
    local m_auc_bo = r(mean)
    local sd_auc_bo = r(sd)

    quietly summarize __br3
    local m_opt = r(mean)
    local sd_opt = r(sd)

    quietly summarize __br4, meanonly
    local m_slope = r(mean)

    quietly _pctile __br2, percentiles(2.5 97.5)
    local auc_bo_lo = r(r1)
    local auc_bo_hi = r(r2)

    restore

    local auc_corr = `auc_app' - `m_opt'

    return scalar auc_apparent = `auc_app'
    return scalar optimism = `m_opt'
    return scalar optimism_sd = `sd_opt'
    return scalar auc_corrected = `auc_corr'
    return scalar auc_bo_mean = `m_auc_bo'
    return scalar auc_bo_lo = `auc_bo_lo'
    return scalar auc_bo_hi = `auc_bo_hi'
    return scalar slope_corrected = `m_slope'
    return scalar n_valid = `n_valid'
    return scalar B = `B'
    return scalar N = `N'
end

* ###########################################################################
* #                                                                         #
* #   APPLICATION AUX 3 MODÈLES                                             #
* #                                                                         #
* ###########################################################################

* ===========================================================================
*  M1. UMOD seul (section 7)
* ===========================================================================
display _newline(2) "=============================================="
display              "  M1. Bootstrap : UMOD seul"
display              "      (non-anuriques)"
display              "=============================================="
bootval kru_ge2 umod if kru_pos == 1, b(`B') seed(`SEED')

local m1_app   = r(auc_apparent)
local m1_opt   = r(optimism)
local m1_corr  = r(auc_corrected)
local m1_slope = r(slope_corrected)
local m1_lo    = r(auc_bo_lo)
local m1_hi    = r(auc_bo_hi)
local m1_n     = r(n_valid)

display _newline "  --- Résultats M1 (UMOD seul) ---"
display "    AUC apparente   = " %5.3f `m1_app'
display "    Optimism (moy.) = " %5.3f `m1_opt'
display "    AUC corrigée    = " %5.3f `m1_corr'
display "    Pente calib.    = " %5.3f `m1_slope'
display "    AUC honnête (centiles 2.5-97.5) : " %5.3f `m1_lo' " — " %5.3f `m1_hi'
display "    Itérations valides : `m1_n'/`B'"

* ===========================================================================
*  M2. UMOD + B2M (section 9) — PARCIMONIEUX
* ===========================================================================
display _newline(2) "=============================================="
display              "  M2. Bootstrap : UMOD + B2M (parcimonieux)"
display              "      (non-anuriques)"
display              "=============================================="
bootval kru_ge2 umod labb2mprehd if kru_pos == 1, b(`B') seed(`SEED')

local m2_app   = r(auc_apparent)
local m2_opt   = r(optimism)
local m2_corr  = r(auc_corrected)
local m2_slope = r(slope_corrected)
local m2_lo    = r(auc_bo_lo)
local m2_hi    = r(auc_bo_hi)
local m2_n     = r(n_valid)

display _newline "  --- Résultats M2 (UMOD + B2M) ---"
display "    AUC apparente   = " %5.3f `m2_app'
display "    Optimism (moy.) = " %5.3f `m2_opt'
display "    AUC corrigée    = " %5.3f `m2_corr'
display "    Pente calib.    = " %5.3f `m2_slope'
display "    AUC honnête (centiles 2.5-97.5) : " %5.3f `m2_lo' " — " %5.3f `m2_hi'
display "    Itérations valides : `m2_n'/`B'"

* ===========================================================================
*  M3. UMOD + créat + urée + B2M (section 8) — COMPLET
* ===========================================================================
display _newline(2) "=============================================="
display              "  M3. Bootstrap : UMOD + créat + urée + B2M"
display              "      (non-anuriques)"
display              "=============================================="
bootval kru_ge2 umod labcreatprehd labureaprehd labb2mprehd if kru_pos == 1, ///
    b(`B') seed(`SEED')

local m3_app   = r(auc_apparent)
local m3_opt   = r(optimism)
local m3_corr  = r(auc_corrected)
local m3_slope = r(slope_corrected)
local m3_lo    = r(auc_bo_lo)
local m3_hi    = r(auc_bo_hi)
local m3_n     = r(n_valid)

display _newline "  --- Résultats M3 (UMOD + 3 biomarqueurs) ---"
display "    AUC apparente   = " %5.3f `m3_app'
display "    Optimism (moy.) = " %5.3f `m3_opt'
display "    AUC corrigée    = " %5.3f `m3_corr'
display "    Pente calib.    = " %5.3f `m3_slope'
display "    AUC honnête (centiles 2.5-97.5) : " %5.3f `m3_lo' " — " %5.3f `m3_hi'
display "    Itérations valides : `m3_n'/`B'"

* ===========================================================================
*  BILAN COMPARATIF
* ===========================================================================
display _newline(2) "=========================================================================="
display              "  BILAN — Validation bootstrap (B=`B')"
display              "=========================================================================="
display "  Modèle                      | AUC app | Optim | AUC corr | Slope  | IC 95% honnête"
display "  ----------------------------|---------|-------|----------|--------|----------------"
display "  M1. UMOD seul               |  " %5.3f `m1_app' "  | " %5.3f `m1_opt' " |  " %5.3f `m1_corr' "   | " %5.3f `m1_slope' "  | " %5.3f `m1_lo' "-" %5.3f `m1_hi'
display "  M2. UMOD + B2M              |  " %5.3f `m2_app' "  | " %5.3f `m2_opt' " |  " %5.3f `m2_corr' "   | " %5.3f `m2_slope' "  | " %5.3f `m2_lo' "-" %5.3f `m2_hi'
display "  M3. UMOD + 3 biomarqueurs   |  " %5.3f `m3_app' "  | " %5.3f `m3_opt' " |  " %5.3f `m3_corr' "   | " %5.3f `m3_slope' "  | " %5.3f `m3_lo' "-" %5.3f `m3_hi'
display "=========================================================================="
display "  Interprétation :"
display "  - Optimism élevé   → le modèle apparent surestime sa performance future"
display "  - AUC corrigée     → estimation honnête de l'AUC sur de nouveaux patients"
display "  - Slope ≈ 1        → calibration optimale des coefficients"
display "  - Slope << 1       → overfitting (coefficients trop extrêmes)"
display "  - IC honnête       → distribution des AUC du modèle bootstrap appliqué"
display "                       à l'original (centiles 2.5 / 97.5)"
display "=========================================================================="
