* ===========================================================================
* 03_estimate_kru.do
* Objectif : Modélisation de kru_daugirdas_35 (mL/min/35L) à partir
*            de la uromoduline sérique (umod, ng/mL)
*
* Prérequis : lancer 02_calculate_kru.do AVANT ce do-file dans la même
*             session Stata. Les variables suivantes doivent être présentes
*             en mémoire : umod, kru_daugirdas_35, kru_naif_35, diuresis.
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

* ── 1. Exploration descriptive : UMOD et KRU ────────────────────

* ─ 1a. Distribution de umod ────────────────────────────────────
display _newline "=== Distribution umod (ng/mL) ==="
summarize umod, detail

histogram umod, normal ///
    title("Distribution de l'uromoduline sérique") ///
    xtitle("UMOD (ng/mL)") name(hist_umod, replace)

* Test de normalité (Shapiro-Wilk, valide pour N < 2000)
swilk umod

* ─ 1b. Distribution de log(umod) ───────────────────────────────
gen log_umod = log(umod)
label variable log_umod "log(UMOD) [ln, ng/mL]"

display _newline "=== Distribution log(umod) ==="
summarize log_umod, detail

histogram log_umod, normal ///
    title("Distribution log(uromoduline)") ///
    xtitle("ln(UMOD)") name(hist_log_umod, replace)

swilk log_umod

* ─ 1c. Valeurs manquantes et extrêmes ──────────────────────────
count if missing(umod)
display "UMOD manquant : " r(N)

count if umod == 0
display "UMOD = 0     : " r(N)

* Valeurs potentiellement aberrantes (< 1 ou > 2000 ng/mL)
count if umod < 1 & !missing(umod)
display "UMOD < 1 ng/mL : " r(N)
count if umod > 2000 & !missing(umod)
display "UMOD > 2000 ng/mL : " r(N)

* ─ 1d. Relation umod ~ kru_daugirdas_35 ────────────────────────
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

* ─ 1e. Stratification anuriques / non-anuriques ─────────────────
display _newline "=== UMOD selon statut de diurèse ==="
tabstat umod, by(diuresis) statistics(n mean sd p25 p50 p75) format(%7.1f)

* Scatter séparé non-anuriques (KRU > 0)
twoway (scatter kru_daugirdas_35 log_umod if kru_daugirdas_35 > 0, ///
        msize(small)) ///
       (lfit kru_daugirdas_35 log_umod if kru_daugirdas_35 > 0), ///
    title("log(UMOD) vs KRU Daugirdas/35L — non-anuriques") ///
    xtitle("ln(UMOD)") ytitle("KRU Daugirdas (mL/min/35L)") ///
    name(scatter_nonanuric, replace)

* ─ 1f. Résumé pour orienter la modélisation ────────────────────
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
* ── 2. Préparation pour modélisation ──────────────────────────
*  - Imputation UMOD < LLOQ (LLOQ ELISA = 2.0 ng/mL → substitution LLOQ/2 = 1.0)
*  - Création de la variable binaire kru_pos = 1 si KRU > 0
* ===========================================================================

* 2a. Imputation UMOD = 0 par LLOQ/2 = 1.0 ng/mL
gen umod_imp = umod
replace umod_imp = 1.0 if umod == 0
label variable umod_imp "UMOD (ng/mL, UMOD=0 imputé à LLOQ/2=1.0)"

display _newline "=== Vérification imputation UMOD ==="
tabstat umod umod_imp, statistics(n min p1 p5 p25 p50 mean) format(%6.2f)

* 2b. Variable binaire pour la partie 1
gen byte kru_pos = (kru_daugirdas_35 > 0) if !missing(kru_daugirdas_35)
label variable kru_pos "KRU Daugirdas > 0 (1=non-anurique)"
label define krupos 0 "Anurique (KRU=0)" 1 "Non-anurique (KRU>0)"
label values kru_pos krupos
display _newline "=== Distribution kru_pos ==="
tab kru_pos

* ===========================================================================
* ── 3. PARTIE 1 — Logit : P(KRU>0) ~ UMOD ─────────────────────
* ===========================================================================
display _newline(2) "=============================================="
display              "  PARTIE 1 — Logit  P(KRU>0)  ~  UMOD"
display              "=============================================="

logit kru_pos umod_imp
estimates store logit_umod

* Odds ratios + IC95
display _newline "=== Odds ratios ==="
logit kru_pos umod_imp, or

* Probabilité prédite et AUC
predict p_pos, pr
label variable p_pos "P(KRU>0 | UMOD) — logit"

display _newline "=== Courbe ROC ==="
lroc, name(roc_part1, replace) title("ROC — P(KRU>0) ~ UMOD")

* Calibration : décile de risque vs proportion observée
display _newline "=== Calibration Hosmer-Lemeshow ==="
estat gof, group(10) table

* ===========================================================================
* ── 4. PARTIE 2 — OLS : KRU ~ UMOD parmi non-anuriques ─────────
* ===========================================================================
display _newline(2) "=============================================="
display              "  PARTIE 2 — OLS  KRU  ~  UMOD  | KRU > 0"
display              "=============================================="

regress kru_daugirdas_35 umod_imp if kru_pos == 1
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
* ── 5. PRÉDICTION TWO-PART & PERFORMANCE ───────────────────────
*  E[KRU | UMOD] = P(KRU>0 | UMOD) × E[KRU | KRU>0, UMOD]
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
regress kru_daugirdas_35 umod_imp
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
