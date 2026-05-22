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
* SECTION 3 — ESTIMATION DE KRU CONTINU À PARTIR DE UMOD SEUL
*              Modèle OLS sur l'ensemble du spectre (KRU = 0 inclus)
* ###########################################################################

* ===========================================================================
*  3a. RÉGRESSION LINÉAIRE : KRU ~ UMOD (population entière, N=151)
* ===========================================================================
display _newline(2) "=============================================="
display              "  3a. OLS  KRU = b0 + b1 * UMOD"
display              "      (population entière, KRU=0 inclus)"
display              "=============================================="

regress kru_daugirdas_35 umod
estimates store ols_umod

display _newline "  R² = " %5.3f e(r2) "    Adj R² = " %5.3f e(r2_a) "    RMSE = " %5.3f e(rmse)

* ===========================================================================
*  3b. PRÉDICTIONS ET RÉSIDUS
* ===========================================================================
display _newline(2) "=============================================="
display              "  3b. Prédictions et diagnostics"
display              "=============================================="

capture drop kru_hat
predict kru_hat, xb
label variable kru_hat "KRU prédit (OLS, mL/min/35L)"

capture drop kru_resid
predict kru_resid, resid
label variable kru_resid "Résidus OLS (mL/min/35L)"

* Performance globale (RMSE, MAE, corrélation)
capture drop _r2 _ar2
gen double _r2 = kru_resid^2
gen double _ar2 = abs(kru_resid)
quietly summarize _r2
local rmse = sqrt(r(mean))
quietly summarize _ar2
local mae = r(mean)
drop _r2 _ar2

display _newline "  RMSE = " %5.3f `rmse' " mL/min/35L"
display         "  MAE  = " %5.3f `mae'  " mL/min/35L"

quietly corr kru_daugirdas_35 kru_hat
display "  Corrélation Pearson (observé,prédit) = " %5.3f r(rho)

quietly spearman kru_daugirdas_35 kru_hat
display "  Corrélation Spearman                 = " %5.3f r(rho)

* Normalité des résidus
display _newline "  Normalité des résidus (Shapiro-Wilk) :"
swilk kru_resid

* ===========================================================================
*  3c. GRAPHIQUES DIAGNOSTIQUES
* ===========================================================================
* Observé vs prédit
twoway (scatter kru_daugirdas_35 kru_hat, msize(small)) ///
       (function y=x, range(0 10) lcolor(red) lpattern(dash)), ///
    title("KRU observé vs prédit — OLS sur UMOD") ///
    xtitle("KRU prédit (mL/min/35L)") ///
    ytitle("KRU observé (mL/min/35L)") ///
    legend(off) name(obs_pred_ols, replace)

* Résidus vs prédits
twoway (scatter kru_resid kru_hat, msize(small)) ///
       (lowess kru_resid kru_hat, lcolor(red)), ///
    yline(0, lpattern(dash)) ///
    title("Résidus vs prédits — OLS sur UMOD") ///
    xtitle("KRU prédit (mL/min/35L)") ytitle("Résidus") ///
    legend(off) name(resid_fit_ols, replace)

* QQ-plot des résidus
qnorm kru_resid, title("QQ-plot résidus OLS") name(qq_ols, replace)

* Régression observé vs UMOD (visualisation directe)
twoway (scatter kru_daugirdas_35 umod, msize(small) jitter(1)) ///
       (lfit kru_daugirdas_35 umod, lcolor(red)), ///
    title("KRU observé vs UMOD") ///
    xtitle("UMOD (ng/mL)") ytitle("KRU Daugirdas (mL/min/35L)") ///
    legend(off) name(kru_vs_umod, replace)
