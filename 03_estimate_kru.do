* ===========================================================================
* 03_estimate_kru.do
* Prérequis : lancer 02_calculate_kru.do AVANT ce do-file dans la même
*             session Stata.
* ===========================================================================

* --- Vérification que 02_calculate_kru.do a été exécuté ---
foreach v in kru_daugirdas_35 kru_naif_35 umod diuresis labcreatprehd {
    capture confirm variable `v'
    if _rc {
        display as error "ERREUR : variable '`v'' absente — lance d'abord 02_calculate_kru.do"
        exit 111
    }
}

display _newline "=== Prérequis OK — variables de 02_calculate_kru.do présentes ==="

* --- Variables binaires utilisées dans tout le do-file ---
capture drop kru_pos
gen byte kru_pos = (kru_daugirdas_35 > 0) if !missing(kru_daugirdas_35)
label variable kru_pos "KRU > 0 (1 = non-anurique)"
label define krupos 0 "Anurique (KRU=0)" 1 "Non-anurique (KRU>0)", replace
label values kru_pos krupos

capture drop kru_ge2
gen byte kru_ge2 = (kru_daugirdas_35 >= 2) if !missing(kru_daugirdas_35)
label variable kru_ge2 "KRU >= 2 mL/min/35L"
label define kruge2 0 "KRU < 2" 1 "KRU >= 2", replace
label values kru_ge2 kruge2

* ###########################################################################
* SECTION 1 — DESCRIPTION DE L'UROMODULINE SÉRIQUE (UMOD)
* ###########################################################################

* ===========================================================================
*  1a. DESCRIPTION GLOBALE DE L'UMOD
* ===========================================================================
display _newline(2) "=============================================="
display              "  1a. Distribution UMOD — population entière"
display              "=============================================="

summarize umod, detail

swilk umod

count if umod == 0
display _newline "  UMOD = 0 : " r(N) " patients"

display _newline "  Répartition des UMOD=0 selon kru_pos :"
tab kru_pos if umod == 0, miss

histogram umod, normal ///
    title("Distribution de l'uromoduline sérique") ///
    xtitle("UMOD (ng/mL)") name(hist_umod, replace)

* ===========================================================================
*  1b. UMOD CHEZ ANURIQUES vs NON-ANURIQUES
* ===========================================================================
display _newline(2) "=============================================="
display              "  1b. UMOD : anuriques vs non-anuriques"
display              "=============================================="

tabstat umod, by(kru_pos) statistics(n mean sd p25 p50 p75) format(%7.2f)

ranksum umod, by(kru_pos)

* ===========================================================================
*  1c. UMOD CHEZ KRU<2 vs KRU>=2 — POPULATION ENTIÈRE
* ===========================================================================
display _newline(2) "=============================================="
display              "  1c. UMOD : KRU<2 vs KRU>=2 (population entière)"
display              "=============================================="

tab kru_ge2, miss

tabstat umod, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%7.2f)

ranksum umod, by(kru_ge2)

* ===========================================================================
*  1d. UMOD CHEZ KRU<2 vs KRU>=2 — NON-ANURIQUES UNIQUEMENT
* ===========================================================================
display _newline(2) "=============================================="
display              "  1d. UMOD : KRU<2 vs KRU>=2 (non-anuriques)"
display              "=============================================="

tabstat umod if kru_pos == 1, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%7.2f)

ranksum umod if kru_pos == 1, by(kru_ge2)

* ###########################################################################
* SECTION 2 — DESCRIPTION DE LA CRÉATININE PRÉDIALYSE (comparaison)
* ###########################################################################

* ===========================================================================
*  2a. DESCRIPTION GLOBALE DE LA CRÉATININE
* ===========================================================================
display _newline(2) "=============================================="
display              "  2a. Distribution créatinine — population entière"
display              "=============================================="

summarize labcreatprehd, detail

swilk labcreatprehd

count if missing(labcreatprehd)
display _newline "  Créatinine manquante : " r(N) " patients"

histogram labcreatprehd, normal ///
    title("Distribution de la créatinine prédialyse") ///
    xtitle("Créatinine (umol/L)") name(hist_creat, replace)

* ===========================================================================
*  2b. CRÉATININE CHEZ ANURIQUES vs NON-ANURIQUES
* ===========================================================================
display _newline(2) "=============================================="
display              "  2b. Créatinine : anuriques vs non-anuriques"
display              "=============================================="

tabstat labcreatprehd, by(kru_pos) statistics(n mean sd p25 p50 p75) format(%7.1f)

ranksum labcreatprehd, by(kru_pos)

* ===========================================================================
*  2c. CRÉATININE CHEZ KRU<2 vs KRU>=2 — POPULATION ENTIÈRE
* ===========================================================================
display _newline(2) "=============================================="
display              "  2c. Créatinine : KRU<2 vs KRU>=2 (population entière)"
display              "=============================================="

tabstat labcreatprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%7.1f)

ranksum labcreatprehd, by(kru_ge2)

* ===========================================================================
*  2d. CRÉATININE CHEZ KRU<2 vs KRU>=2 — NON-ANURIQUES UNIQUEMENT
* ===========================================================================
display _newline(2) "=============================================="
display              "  2d. Créatinine : KRU<2 vs KRU>=2 (non-anuriques)"
display              "=============================================="

tabstat labcreatprehd if kru_pos == 1, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%7.1f)

ranksum labcreatprehd if kru_pos == 1, by(kru_ge2)

* ###########################################################################
* SECTION 3 — TWO-PART MODEL : KRU ~ UMOD
*              Partie 1 : logit P(KRU>0)        — sur N=151
*              Partie 2 : OLS  E[KRU | KRU>0]   — sur N=89 (non-anuriques)
*              Combiné  : E[KRU | UMOD] = P(KRU>0|UMOD) × E[KRU|KRU>0,UMOD]
* ###########################################################################

* ===========================================================================
*  3a. PARTIE 1 — Logit P(KRU>0) ~ UMOD
* ===========================================================================
display _newline(2) "=============================================="
display              "  3a. Logit  P(KRU>0)  ~  UMOD     (N=151)"
display              "=============================================="

logit kru_pos umod
estimates store tp_logit

display _newline "  --- Odds ratios ---"
logit kru_pos umod, or

capture drop p_pos
predict p_pos, pr
label variable p_pos "P(KRU>0 | UMOD)"

display _newline "  --- Performance (ROC / calibration) ---"
lroc, nograph
display "    AUC = " %5.3f r(area)

estat gof, group(10) table

* ===========================================================================
*  3b. PARTIE 2 — E[KRU | KRU>0] ~ UMOD  (non-anuriques, N=89)
*       OLS (linéaire) et MFP (fractional polynomial) en parallèle
* ===========================================================================
display _newline(2) "=============================================="
display              "  3b. E[KRU | KRU>0]  ~  UMOD   (N=89)"
display              "      OLS linéaire vs MFP en parallèle"
display              "=============================================="

* --- 3b.i  OLS linéaire (modèle principal pour la prédiction two-part) ---
*           SE robustes (Huber-White) — coefficients inchangés, inférence
*           corrigée pour l'hétéroscédasticité résiduelle.
display _newline "  --- OLS linéaire (SE robustes Huber-White) ---"
regress kru_daugirdas_35 umod if kru_pos == 1, vce(robust)
estimates store tp_ols

local r2_ols    = e(r2)
local rmse_ols  = e(rmse)
local aic_ols   = .
quietly estat ic
matrix _ic_ols = r(S)
local aic_ols = _ic_ols[1,5]

display _newline "  R² OLS   = " %5.3f `r2_ols' ///
    "    RMSE OLS = " %5.3f `rmse_ols' ///
    "    AIC OLS  = " %6.2f `aic_ols'

* Prédiction conditionnelle (utilisée dans 3c)
capture drop kru_cond
predict kru_cond, xb
label variable kru_cond "E[KRU | KRU>0, UMOD] — OLS"
replace kru_cond = 0 if kru_cond < 0

* Diagnostic des résidus OLS
capture drop resid_p2
predict resid_p2 if e(sample), resid

display _newline "  Normalité des résidus OLS (Shapiro-Wilk) :"
swilk resid_p2

* --- 3b.ii  MFP (fractional polynomial multivariable) ---
display _newline(2) "  --- MFP (fractional polynomial) ---"
display         "  Recherche automatique de la meilleure transformation"
display         "  de UMOD parmi puissances {−2, −1, −0.5, 0=log, 0.5, 1, 2, 3}"

mfp: regress kru_daugirdas_35 umod if kru_pos == 1
estimates store tp_mfp

local r2_mfp   = e(r2)
local rmse_mfp = e(rmse)
quietly estat ic
matrix _ic_mfp = r(S)
local aic_mfp = _ic_mfp[1,5]

display _newline "  R² MFP   = " %5.3f `r2_mfp' ///
    "    RMSE MFP = " %5.3f `rmse_mfp' ///
    "    AIC MFP  = " %6.2f `aic_mfp'

* Prédiction MFP (uniquement à titre de comparaison ; non utilisée dans le two-part)
capture drop kru_cond_mfp
predict kru_cond_mfp, xb
label variable kru_cond_mfp "E[KRU | KRU>0, UMOD] — MFP"

* --- 3b.iii  Comparaison OLS vs MFP ---
display _newline(2) "  --- Comparaison OLS vs MFP (N=89 non-anuriques) ---"
estimates stats tp_ols tp_mfp

display _newline "  ───────────────────────────────────────────"
display         "                    OLS         MFP"
display         "  ───────────────────────────────────────────"
display         "  R²              " %6.3f `r2_ols'   "      " %6.3f `r2_mfp'
display         "  RMSE            " %6.3f `rmse_ols' "      " %6.3f `rmse_mfp'
display         "  AIC             " %6.2f `aic_ols'  "    " %6.2f `aic_mfp'
display         "  ───────────────────────────────────────────"
display         "  ΔAIC (OLS − MFP) = " %5.2f (`aic_ols' - `aic_mfp')
display         "    > 2 : MFP préférable"
display         "    < −2 : OLS préférable"
display         "    |ΔAIC| < 2 : équivalents → préférer le plus simple (OLS)"

* ===========================================================================
*  3c. COMBINAISON TWO-PART : E[KRU | UMOD] = P(KRU>0|UMOD) × E[KRU|KRU>0,UMOD]
* ===========================================================================
display _newline(2) "=============================================="
display              "  3c. Prédiction two-part combinée"
display              "=============================================="

capture drop kru_pred_2p
gen double kru_pred_2p = p_pos * kru_cond if !missing(p_pos, kru_cond)
label variable kru_pred_2p "KRU prédit two-part (mL/min/35L)"

* Performance globale (RMSE, MAE, corrélation observé/prédit)
capture drop _r_sq _r_abs
gen double _r_sq  = (kru_daugirdas_35 - kru_pred_2p)^2
gen double _r_abs = abs(kru_daugirdas_35 - kru_pred_2p)

quietly summarize _r_sq
local rmse_2p = sqrt(r(mean))
quietly summarize _r_abs
local mae_2p = r(mean)
drop _r_sq _r_abs

display _newline "  RMSE two-part = " %5.3f `rmse_2p' " mL/min/35L"
display         "  MAE  two-part = " %5.3f `mae_2p'  " mL/min/35L"

quietly corr kru_daugirdas_35 kru_pred_2p
display "  Corrélation Pearson (observé,prédit) = " %5.3f r(rho)

quietly spearman kru_daugirdas_35 kru_pred_2p
display "  Corrélation Spearman                 = " %5.3f r(rho)

* ===========================================================================
*  3d. GRAPHIQUES DIAGNOSTIQUES
* ===========================================================================
* P(KRU>0) prédite vs UMOD
twoway (line p_pos umod, sort lcolor(red) lwidth(medium)) ///
       (scatter kru_pos umod, msize(small) mcolor(%40) jitter(2)), ///
    yline(0.5, lpattern(dot) lcolor(gray)) ///
    title("Partie 1 : P(KRU>0) selon UMOD") ///
    xtitle("UMOD (ng/mL)") ytitle("P(KRU>0)") ///
    legend(order(1 "Logit" 2 "Observé") position(11) ring(0)) ///
    name(tp_part1, replace)

* OLS sur non-anuriques (observé vs UMOD)
twoway (scatter kru_daugirdas_35 umod if kru_pos == 1, msize(small)) ///
       (lfit kru_daugirdas_35 umod if kru_pos == 1, lcolor(red)), ///
    title("Partie 2 : KRU ~ UMOD (non-anuriques)") ///
    xtitle("UMOD (ng/mL)") ytitle("KRU Daugirdas (mL/min/35L)") ///
    legend(off) name(tp_part2, replace)

* Observé vs prédit two-part (tous patients)
twoway (scatter kru_daugirdas_35 kru_pred_2p, msize(small)) ///
       (function y=x, range(0 10) lcolor(red) lpattern(dash)), ///
    title("KRU observé vs prédit — two-part") ///
    xtitle("KRU prédit (mL/min/35L)") ///
    ytitle("KRU observé (mL/min/35L)") ///
    legend(off) name(tp_obs_pred, replace)

* Prédit two-part en fonction d'UMOD (courbe de la prédiction combinée)
twoway (line kru_pred_2p umod, sort lcolor(red) lwidth(medium)) ///
       (scatter kru_daugirdas_35 umod, msize(small) mcolor(%40) jitter(1)), ///
    title("Two-part : KRU prédit selon UMOD") ///
    xtitle("UMOD (ng/mL)") ytitle("KRU (mL/min/35L)") ///
    legend(order(1 "Prédit two-part" 2 "Observé") position(11) ring(0)) ///
    name(tp_curve, replace)

