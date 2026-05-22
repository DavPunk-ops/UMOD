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
*  3b. PARTIE 2 — Modèles E[KRU | KRU>0] ~ UMOD  (OLS et GLM Gamma en parallèle)
* ===========================================================================
display _newline(2) "=============================================="
display              "  3b. E[KRU | KRU>0]  ~  UMOD   (N=89)"
display              "      OLS et GLM Gamma (log link) en parallèle"
display              "=============================================="

* --- 3b.i  OLS ---
display _newline "  --- OLS ---"
regress kru_daugirdas_35 umod if kru_pos == 1
estimates store tp_ols

local rmse_ols_p2 = e(rmse)
local r2_ols_p2   = e(r2)

display _newline "  R² OLS = " %5.3f `r2_ols_p2' "    RMSE OLS = " %5.3f `rmse_ols_p2'

capture drop kru_cond_ols
predict kru_cond_ols, xb
label variable kru_cond_ols "E[KRU | KRU>0, UMOD] — OLS"
replace kru_cond_ols = 0 if kru_cond_ols < 0

capture drop resid_ols_p2
predict resid_ols_p2 if e(sample), resid

display _newline "  Normalité résidus OLS (Shapiro-Wilk) :"
swilk resid_ols_p2

* --- 3b.ii  GLM Gamma, log link ---
display _newline(2) "  --- GLM Gamma (log link) ---"
glm kru_daugirdas_35 umod if kru_pos == 1, family(gamma) link(log)
estimates store tp_glm

* Coefficients sur l'échelle multiplicative (exp(beta))
display _newline "  --- Exp(coef) : facteur multiplicatif de KRU par +1 ng/mL UMOD ---"
glm kru_daugirdas_35 umod if kru_pos == 1, family(gamma) link(log) eform

capture drop kru_cond_glm
predict kru_cond_glm, mu
label variable kru_cond_glm "E[KRU | KRU>0, UMOD] — GLM Gamma"
* Prédictions GLM avec log link toujours >0, pas besoin de borner

* RMSE GLM sur non-anuriques (échelle originale)
capture drop _r_glm_sq
gen double _r_glm_sq = (kru_daugirdas_35 - kru_cond_glm)^2 if kru_pos == 1
quietly summarize _r_glm_sq
local rmse_glm_p2 = sqrt(r(mean))
drop _r_glm_sq

display _newline "  RMSE GLM (échelle originale, non-anuriques) = " %5.3f `rmse_glm_p2'

* --- 3b.iii  Comparaison AIC/BIC/RMSE OLS vs GLM ---
display _newline(2) "  --- Comparaison OLS vs GLM Gamma (sur N=89 non-anuriques) ---"
estimates stats tp_ols tp_glm

display _newline "  RMSE OLS     = " %5.3f `rmse_ols_p2'
display         "  RMSE GLM     = " %5.3f `rmse_glm_p2'

* Variable kru_cond qui sera utilisée en 3c (par défaut : OLS, modifiable plus loin)
capture drop kru_cond
gen double kru_cond = kru_cond_ols
label variable kru_cond "E[KRU | KRU>0, UMOD] (OLS par défaut)"

* ===========================================================================
*  3c. COMBINAISON TWO-PART : E[KRU | UMOD] = P(KRU>0|UMOD) × E[KRU|KRU>0,UMOD]
*       Calcul pour les 2 variantes (OLS et GLM) avec comparaison RMSE/MAE
* ===========================================================================
display _newline(2) "=============================================="
display              "  3c. Prédiction two-part combinée (OLS vs GLM)"
display              "=============================================="

* --- Prédiction combinée variante OLS ---
capture drop kru_pred_2p_ols
gen double kru_pred_2p_ols = p_pos * kru_cond_ols if !missing(p_pos, kru_cond_ols)
label variable kru_pred_2p_ols "KRU prédit two-part (OLS)"

* --- Prédiction combinée variante GLM ---
capture drop kru_pred_2p_glm
gen double kru_pred_2p_glm = p_pos * kru_cond_glm if !missing(p_pos, kru_cond_glm)
label variable kru_pred_2p_glm "KRU prédit two-part (GLM Gamma)"

* --- Performance variante OLS ---
capture drop _r_ols_sq _r_ols_abs
gen double _r_ols_sq  = (kru_daugirdas_35 - kru_pred_2p_ols)^2
gen double _r_ols_abs = abs(kru_daugirdas_35 - kru_pred_2p_ols)
quietly summarize _r_ols_sq
local rmse_2p_ols = sqrt(r(mean))
quietly summarize _r_ols_abs
local mae_2p_ols = r(mean)
drop _r_ols_sq _r_ols_abs

quietly corr kru_daugirdas_35 kru_pred_2p_ols
local pear_ols = r(rho)
quietly spearman kru_daugirdas_35 kru_pred_2p_ols
local spear_ols = r(rho)

* --- Performance variante GLM ---
capture drop _r_glm_sq _r_glm_abs
gen double _r_glm_sq  = (kru_daugirdas_35 - kru_pred_2p_glm)^2
gen double _r_glm_abs = abs(kru_daugirdas_35 - kru_pred_2p_glm)
quietly summarize _r_glm_sq
local rmse_2p_glm = sqrt(r(mean))
quietly summarize _r_glm_abs
local mae_2p_glm = r(mean)
drop _r_glm_sq _r_glm_abs

quietly corr kru_daugirdas_35 kru_pred_2p_glm
local pear_glm = r(rho)
quietly spearman kru_daugirdas_35 kru_pred_2p_glm
local spear_glm = r(rho)

* --- Tableau comparatif ---
display _newline "  ─────────────────────────────────────────────────────"
display         "                 Two-part OLS    Two-part GLM Gamma"
display         "  ─────────────────────────────────────────────────────"
display         "  RMSE          " %7.3f `rmse_2p_ols' "         " %7.3f `rmse_2p_glm'
display         "  MAE           " %7.3f `mae_2p_ols'  "         " %7.3f `mae_2p_glm'
display         "  Pearson r     " %7.3f `pear_ols'    "         " %7.3f `pear_glm'
display         "  Spearman ρ    " %7.3f `spear_ols'   "         " %7.3f `spear_glm'
display         "  ─────────────────────────────────────────────────────"

* Variable principale kru_pred_2p = OLS par défaut (modifiable)
capture drop kru_pred_2p
gen double kru_pred_2p = kru_pred_2p_ols
label variable kru_pred_2p "KRU prédit two-part (OLS par défaut)"

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

