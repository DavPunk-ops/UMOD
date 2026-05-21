* ===========================================================================
* 03_estimate_kru.do
* Objectif : Modélisation de kru_daugirdas_35 (mL/min/35L) à partir
*            de la uromoduline sérique (umod, ng/mL)
* ===========================================================================

local path "C:\Users\dajs\OneDrive - HOPITAUX UNIVERSITAIRES DE GENEVE\recherche\RKF\UMOD\stata\main prospective study\with Claude"

* Charger la base enrichie produite par 02_calculate_kru.do
* (re-exécuter 02 pour garantir l'indépendance de ce fichier)
do "`path'\02_calculate_kru.do"

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
display " → Décision de modélisation à prendre selon"
display "   distribution log(umod) et scatter plots"
display "========================================"
