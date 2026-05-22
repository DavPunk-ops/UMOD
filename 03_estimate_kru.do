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
*  3b. PARTIE 2 — OLS  E[KRU | KRU>0] ~ UMOD  (non-anuriques)
* ===========================================================================
display _newline(2) "=============================================="
display              "  3b. OLS  E[KRU | KRU>0]  ~  UMOD   (N=89)"
display              "=============================================="

regress kru_daugirdas_35 umod if kru_pos == 1
estimates store tp_ols

display _newline "  R² = " %5.3f e(r2) "    Adj R² = " %5.3f e(r2_a) "    RMSE = " %5.3f e(rmse)

* Prédiction conditionnelle E[KRU | KRU>0, UMOD] étendue à tous les patients
capture drop kru_cond
predict kru_cond, xb
label variable kru_cond "E[KRU | KRU>0, UMOD]"

* Borner à 0 (un KRU prédit négatif n'a pas de sens biologique)
replace kru_cond = 0 if kru_cond < 0

* Diagnostics résidus (chez non-anuriques)
capture drop resid_p2
predict resid_p2 if e(sample), resid

display _newline "  Normalité résidus (Shapiro-Wilk) :"
swilk resid_p2

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

