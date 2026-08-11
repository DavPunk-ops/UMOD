* ===========================================================================
* 11_ae3_clinical.do
*
* Prérequis : lancer 02_calculate_kru.do AVANT (données en mémoire).
*             Autonome ; ne modifie AUCUNE do-file existante.
*
* Objectif : réponse à AE#3 (interprétation clinique prudente / screening).
*   Sur la COHORTE COMPLÈTE (modèle publié), produit :
*     1. Équation de prédiction complète (coefficients + intercept + OR)
*     2. Calibration : Hosmer-Lemeshow + pente optimism-corrected (bootstrap)
*        + table observé/attendu par décile
*     3. Faux rule-in / faux rule-out aux seuils publiés (0.24 / 0.55)
*     4. Seuils rule-in conservateurs (Sp ≥90%, ≥95%, ≥97.5%)
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

capture confirm variable kru_ge2
if _rc  gen byte kru_ge2 = (kru_daugirdas_35 >= 2) if !missing(kru_daugirdas_35)

* Modèle combiné publié (logit cohorte complète)
capture drop p_ge2
quietly logit kru_ge2 umod labb2mprehd
predict p_ge2, pr
label variable p_ge2 "P(KRU≥2) — modèle combiné publié"
* Restreindre aux patients avec KRU mesuré : exclut les 2 récoltes incomplètes
* (proba prédite disponible mais outcome absent) → N=148, cohérent avec Table 4.
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
display "  logit(P) = " %7.4f _b[_cons] " + " %7.4f _b[umod] " * UMOD" ///
    " + (" %7.4f _b[labb2mprehd] ") * β2M"
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

* --- Pente de calibration optimism-corrected (bootstrap Harrell) ---
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
display _newline "  --- Pente de calibration ---"
display "  Apparente (in-sample)        = 1.000 (par construction)"
display "  Optimism-corrected (B=`nok') = " %5.3f `slope'
display "  (< 1 = léger surajustement ; proche de 1 = bonne calibration)"
restore

* ###########################################################################
* SECTION 3 — FAUX RULE-IN / FAUX RULE-OUT (seuils publiés 0.24 / 0.55)
* ###########################################################################
display _newline(2) "########################################################"
display              "  3. Faux classifications aux seuils publiés"
display              "########################################################"

quietly count if !missing(p_ge2)
local Ntot = r(N)

* Rule-out (<0.24)
quietly count if p_ge2<0.24 & !missing(p_ge2)
local ro = r(N)
quietly count if p_ge2<0.24 & kru_ge2==0 & !missing(p_ge2)
local ro_ok = r(N)
quietly count if p_ge2<0.24 & kru_ge2==1 & !missing(p_ge2)
local ro_bad = r(N)

* Rule-in (≥0.55)
quietly count if p_ge2>=0.55 & !missing(p_ge2)
local ri = r(N)
quietly count if p_ge2>=0.55 & kru_ge2==1 & !missing(p_ge2)
local ri_ok = r(N)
quietly count if p_ge2>=0.55 & kru_ge2==0 & !missing(p_ge2)
local ri_bad = r(N)

* Grey zone
quietly count if p_ge2>=0.24 & p_ge2<0.55 & !missing(p_ge2)
local gz = r(N)

display _newline "  N (modèle disponible) = `Ntot'"
display _newline "  RULE-OUT (P<0.24) : n=`ro'"
display "     True (KRU<2)               = `ro_ok'"
display "     FALSE rule-out (KRU≥2 missed) = `ro_bad'   (" ///
    %4.1f 100*`ro_bad'/`ro' "% of ruled-out)"
if `ro'>0 display "     NPV = " %5.1f 100*`ro_ok'/`ro' "%"
display _newline "  GREY ZONE (0.24–0.55) : n=`gz'  (" %4.1f 100*`gz'/`Ntot' "%)"
display _newline "  RULE-IN (P≥0.55) : n=`ri'"
display "     True (KRU≥2)               = `ri_ok'"
display "     FALSE rule-in (KRU<2 misclassified) = `ri_bad'   (" ///
    %4.1f 100*`ri_bad'/`ri' "% of ruled-in)"
if `ri'>0 display "     PPV = " %5.1f 100*`ri_ok'/`ri' "%"

* ###########################################################################
* SECTION 4 — SEUILS RULE-IN CONSERVATEURS (Sp ≥90%, ≥95%, ≥97.5%)
* ###########################################################################
display _newline(2) "########################################################"
display              "  4. Seuils rule-in conservateurs (cohorte complète)"
display              "     seuil = centile de P chez KRU<2 pour Sp cible"
display              "########################################################"

capture program drop rulein_full
program define rulein_full
    args spec
    quietly centile p_ge2 if kru_ge2==0 & !missing(p_ge2), centile(`spec')
    local thr = r(c_1)
    quietly count if kru_ge2==1 & !missing(p_ge2)
    local npos = r(N)
    quietly count if kru_ge2==0 & !missing(p_ge2)
    local nneg = r(N)
    quietly count if p_ge2>=`thr' & !missing(p_ge2)
    local nri = r(N)
    quietly count if p_ge2>=`thr' & kru_ge2==1 & !missing(p_ge2)
    local ntrue = r(N)
    quietly count if p_ge2>=`thr' & kru_ge2==0 & !missing(p_ge2)
    local nfalse = r(N)
    local ppv = 100*`ntrue'/`nri'
    local sens = 100*`ntrue'/`npos'
    local specobs = 100*(`nneg'-`nfalse')/`nneg'
    display "  Target Sp≥`spec'% :  threshold P≥" %5.3f `thr' ///
        "  rule-in n=" %2.0f `nri' "  true=" %2.0f `ntrue' "  false=" %1.0f `nfalse' ///
        "  PPV=" %4.1f `ppv' "%  Se=" %4.1f `sens' "%  Sp=" %4.1f `specobs' "%"
end

display _newline "  --- Rendement rule-in à spécificité croissante ---"
rulein_full 90
rulein_full 95
rulein_full 97.5

display _newline(2) "=== FIN do-file 11 ==="
