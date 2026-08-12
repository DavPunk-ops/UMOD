* ===========================================================================
* 11_ae3_clinical.do
*
* Prérequis : lancer 02_calculate_kru.do AVANT (données en mémoire).
*             Autonome ; ne modifie AUCUNE do-file existante.
*
* Objectif : réponse COMPLÈTE à AE#3 (interprétation clinique prudente /
*   screening tool). Toutes les analyses utilisent le MODÈLE PUBLIÉ (logit sur
*   la cohorte complète), restreint aux N=148 avec KRU mesuré.
*
* Sections (= livrables AE#3) :
*   1. Équation de prédiction complète (coefficients + OR)          [equation]
*   2. Calibration : Hosmer-Lemeshow + pente optimism-corrected     [calibration]
*   3. Stratégie deux-seuils — COHORTE COMPLÈTE (Table 4 + mal-classés)
*   4. Stratégie deux-seuils — NON-ANURIQUES (Table 5, seuils publiés + IC Wilson)
*   5. Seuils rule-in conservateurs : spécificité a priori (Sp≥90/95/97.5%)
*      définie sur la cohorte complète, APPLIQUÉE aux non-anuriques (+ IC Wilson)
* ===========================================================================

* --- Prérequis ---
foreach v in kru_daugirdas_35 umod labb2mprehd {
    capture confirm variable `v'
    if _rc {
        display as error "ERREUR : '`v'' absente — lance d'abord 02_calculate_kru.do"
        exit 111
    }
}
display _newline "=== Prérequis OK ==="

* --- Variables dérivées ---
capture confirm variable kru_ge2
if _rc  gen byte kru_ge2 = (kru_daugirdas_35 >= 2) if !missing(kru_daugirdas_35)
capture confirm variable kru_pos
if _rc  gen byte kru_pos = (kru_daugirdas_35 > 0)  if !missing(kru_daugirdas_35)

* Modèle combiné PUBLIÉ (logit cohorte complète)
capture drop p_ge2
quietly logit kru_ge2 umod labb2mprehd
predict p_ge2, pr
label variable p_ge2 "P(KRU≥2) — modèle combiné publié"
* Restreindre aux patients avec KRU mesuré (exclut 2 récoltes incomplètes qui
* ont une proba prédite mais pas d'outcome) → N=148, cohérent avec Table 4.
replace p_ge2 = . if missing(kru_ge2)

* ###########################################################################
* SECTION 1 — ÉQUATION DE PRÉDICTION COMPLÈTE
* ###########################################################################
display _newline(2) "########################################################"
display              "  1. Équation complète : logit P(KRU≥2) ~ UMOD + β2M"
display              "########################################################"

logit kru_ge2 umod labb2mprehd
display _newline "  --- Odds ratios ---"
logit kru_ge2 umod labb2mprehd, or
display _newline "  --- Équation (à reporter) ---"
display "  logit(P) = " %6.3f _b[_cons] " + " %6.3f _b[umod] " * UMOD + (" ///
    %6.3f _b[labb2mprehd] ") * β2M"
display "  P(KRU≥2) = 1 / (1 + exp(-logit(P)))"

* ###########################################################################
* SECTION 2 — CALIBRATION
* ###########################################################################
display _newline(2) "########################################################"
display              "  2. Calibration du modèle (cohorte complète)"
display              "########################################################"

display _newline "  --- Hosmer-Lemeshow + table observé/attendu (déciles) ---"
quietly logit kru_ge2 umod labb2mprehd
estat gof, group(10) table

display _newline "  --- Pente de calibration optimism-corrected (bootstrap) ---"
preserve
keep if !missing(kru_ge2, umod, labb2mprehd)
tempfile orig
quietly save `orig'
set seed 20260522
local B = 1000
local sumslope = 0
local nok = 0
forvalues i = 1/`B' {
    quietly use `orig', clear
    bsample
    capture quietly logit kru_ge2 umod labb2mprehd
    if _rc continue
    matrix bb = e(b)
    quietly use `orig', clear
    capture drop lp_o
    matrix score lp_o = bb
    capture quietly logit kru_ge2 lp_o
    if _rc continue
    local sumslope = `sumslope' + _b[lp_o]
    local nok = `nok' + 1
}
local slope = `sumslope'/`nok'
display "  Apparente (in-sample)        = 1.000 (par construction)"
display "  Optimism-corrected (B=`nok') = " %5.3f `slope'
restore

* ###########################################################################
* PROGRAMME — stratégie deux-seuils (rule-out <0.24 / rule-in ≥0.55)
*   `cond'   : restriction d'échantillon ("1" = cohorte complète ;
*              "kru_pos==1" = non-anuriques)
*   `wilson' : 1 pour afficher les IC binomiaux de Wilson (NPV, PPV)
* ###########################################################################
capture program drop twocut
program define twocut
    args cond label wilson
    quietly count if `cond' & !missing(p_ge2)
    local N = r(N)
    * Rule-out (<0.24)
    quietly count if `cond' & p_ge2<0.24 & !missing(p_ge2)
    local ro = r(N)
    quietly count if `cond' & p_ge2<0.24 & kru_ge2==0 & !missing(p_ge2)
    local ro_ok = r(N)
    quietly count if `cond' & p_ge2<0.24 & kru_ge2==1 & !missing(p_ge2)
    local ro_bad = r(N)
    * Rule-in (≥0.55)
    quietly count if `cond' & p_ge2>=0.55 & !missing(p_ge2)
    local ri = r(N)
    quietly count if `cond' & p_ge2>=0.55 & kru_ge2==1 & !missing(p_ge2)
    local ri_ok = r(N)
    quietly count if `cond' & p_ge2>=0.55 & kru_ge2==0 & !missing(p_ge2)
    local ri_bad = r(N)
    * Grey zone
    quietly count if `cond' & p_ge2>=0.24 & p_ge2<0.55 & !missing(p_ge2)
    local gz = r(N)

    display _newline "  === `label' (N=`N') ==="
    display "  Rule-out (P<0.24): n=`ro' (" %4.1f 100*`ro'/`N' "%)   " ///
        "NPV=" %4.1f 100*`ro_ok'/`ro' "%   FALSE rule-out (KRU≥2 missed)=`ro_bad'"
    if `wilson' cii proportions `ro' `ro_ok', wilson
    display "  Grey zone (0.24–0.55): n=`gz' (" %4.1f 100*`gz'/`N' "%)"
    display "  Rule-in (P≥0.55): n=`ri' (" %4.1f 100*`ri'/`N' "%)   " ///
        "PPV=" %4.1f 100*`ri_ok'/`ri' "%   FALSE rule-in (KRU<2 misclassified)=`ri_bad'"
    if `wilson' cii proportions `ri' `ri_ok', wilson
end

* ###########################################################################
* SECTION 3 — TABLE 4 (COHORTE COMPLÈTE) + MAL-CLASSÉS
* ###########################################################################
display _newline(2) "########################################################"
display              "  3. Deux-seuils — cohorte complète (Table 4 + mal-classés)"
display              "########################################################"
twocut "1" "Full cohort" 0

* ###########################################################################
* SECTION 4 — TABLE 5 (NON-ANURIQUES) + IC WILSON
* ###########################################################################
display _newline(2) "########################################################"
display              "  4. Deux-seuils — non-anuriques (Table 5, mêmes seuils)"
display              "########################################################"
twocut "kru_pos==1" "Measured-KRU subgroup (non-anuric)" 1

* ###########################################################################
* SECTION 5 — SEUILS RULE-IN CONSERVATEURS
*   Seuils définis par spécificité a priori SUR LA COHORTE COMPLÈTE
*   (Sp≥90/95/97.5%), puis APPLIQUÉS aux non-anuriques (seuils externes →
*   IC de Wilson valides, pas d'optimisme).
* ###########################################################################
display _newline(2) "########################################################"
display              "  5. Rule-in conservateur : Sp a priori (full cohort)"
display              "     appliqué aux non-anuriques"
display              "########################################################"

capture program drop conserv_na
program define conserv_na
    args spec
    * Seuil = centile (=spec) de P chez les KRU<2 de la COHORTE COMPLÈTE
    quietly centile p_ge2 if kru_ge2==0 & !missing(p_ge2), centile(`spec')
    local thr = r(c_1)
    * Appliqué aux NON-ANURIQUES
    quietly count if kru_pos==1 & p_ge2>=`thr' & !missing(p_ge2)
    local nri = r(N)
    quietly count if kru_pos==1 & p_ge2>=`thr' & kru_ge2==1 & !missing(p_ge2)
    local ntrue = r(N)
    quietly count if kru_pos==1 & p_ge2>=`thr' & kru_ge2==0 & !missing(p_ge2)
    local nfalse = r(N)
    display _newline "  Sp≥`spec'% → full-cohort threshold P≥" %5.3f `thr' ///
        " ; applied to non-anuric :"
    display "     rule-in n=`nri'   true KRU≥2=`ntrue'   false rule-in=`nfalse'   " ///
        "PPV=" %5.1f 100*`ntrue'/`nri' "%"
    if `nri'>0  cii proportions `nri' `ntrue', wilson
end

conserv_na 90
conserv_na 95
conserv_na 97.5

display _newline(2) "=== FIN do-file 11 ==="
