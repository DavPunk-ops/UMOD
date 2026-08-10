* ===========================================================================
* 09_nonanuric_logit.do
*
* Prérequis : lancer 02_calculate_kru.do AVANT (données en mémoire).
*             Autonome ; ne modifie AUCUNE do-file existante.
*
* Objectif : réponse CIBLÉE à AE#1 (« la performance reflète-t-elle surtout la
*   séparation anurique/non-anurique ? »). Analyse développée ENTIÈREMENT dans
*   le sous-groupe non-anurique (KRU mesuré) : le modèle combiné est REFITTÉ
*   sur les non-anuriques — rien ne touche les anuriques, ni l'outcome ni le fit.
*
* Sections :
*   1. Modèle logit KRU≥2 ~ UMOD + β2M REFITTÉ sur les non-anuriques
*      (équation + OR → contribution indépendante de chaque marqueur)
*   2. AUC : UMOD seul, β2M seul, combiné (refit) — chez les non-anuriques
*   3. Comparaisons DeLong entre les trois
*   4. Validation interne du combiné : bootstrap Harrell (optimism-corrected)
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
capture confirm variable kru_pos
if _rc  gen byte kru_pos = (kru_daugirdas_35 > 0)  if !missing(kru_daugirdas_35)
capture confirm variable kru_ge2
if _rc  gen byte kru_ge2 = (kru_daugirdas_35 >= 2) if !missing(kru_daugirdas_35)

capture drop neg_b2m
gen double neg_b2m = -labb2mprehd if !missing(labb2mprehd)
label variable neg_b2m "−β2M (orienté : haut = KRU≥2)"

* ###########################################################################
* SECTION 1 — MODÈLE COMBINÉ REFITTÉ SUR LES NON-ANURIQUES
* ###########################################################################
display _newline(2) "########################################################"
display              "  1. Logit KRU≥2 ~ UMOD + β2M — REFIT sur non-anuriques"
display              "########################################################"

display _newline "  --- Coefficients ---"
logit kru_ge2 umod labb2mprehd if kru_pos==1

display _newline "  --- Odds ratios (contribution indépendante) ---"
logit kru_ge2 umod labb2mprehd if kru_pos==1, or

* Probabilité prédite du modèle refitté (pour AUC & DeLong)
capture drop p_combo_na
predict p_combo_na, pr
label variable p_combo_na "P(KRU≥2) — modèle refitté non-anuriques"

* ###########################################################################
* SECTION 2 — AUC CHEZ LES NON-ANURIQUES
* ###########################################################################
display _newline(2) "########################################################"
display              "  2. AUC (non-anuriques) — UMOD / β2M / combiné refitté"
display              "########################################################"

display _newline "  --- UMOD seul ---"
roctab kru_ge2 umod       if kru_pos==1
display _newline "  --- β2M seul ---"
roctab kru_ge2 neg_b2m    if kru_pos==1
display _newline "  --- Combiné (modèle refitté) ---"
roctab kru_ge2 p_combo_na if kru_pos==1

* ###########################################################################
* SECTION 3 — COMPARAISONS DeLong
* ###########################################################################
display _newline(2) "########################################################"
display              "  3. Comparaisons DeLong (non-anuriques)"
display              "########################################################"

display _newline "  --- UMOD vs β2M ---"
roccomp kru_ge2 umod neg_b2m ///
    if kru_pos==1 & !missing(umod, neg_b2m), summary
display _newline "  --- Combiné vs UMOD ---"
roccomp kru_ge2 p_combo_na umod ///
    if kru_pos==1 & !missing(p_combo_na, umod), summary
display _newline "  --- Combiné vs β2M ---"
roccomp kru_ge2 p_combo_na neg_b2m ///
    if kru_pos==1 & !missing(p_combo_na, neg_b2m), summary

* ###########################################################################
* SECTION 4 — VALIDATION INTERNE DU COMBINÉ (BOOTSTRAP HARRELL)
* ###########################################################################
display _newline(2) "########################################################"
display              "  4. Bootstrap Harrell du combiné refitté (optimism)"
display              "########################################################"

preserve
keep if kru_pos==1 & !missing(umod, labb2mprehd, kru_ge2)
tempfile orig
quietly save `orig'

quietly logit kru_ge2 umod labb2mprehd
quietly predict p_app, pr
quietly roctab kru_ge2 p_app
local auc_app = r(area)
drop p_app

set seed 20260522
local B = 1000
local sumopt = 0
local nok = 0
forvalues i = 1/`B' {
    quietly use `orig', clear
    bsample
    capture quietly logit kru_ge2 umod labb2mprehd
    if _rc continue
    matrix bb = e(b)
    quietly predict pb, pr
    quietly roctab kru_ge2 pb
    local auc_b = r(area)
    quietly use `orig', clear
    capture drop lp_o
    matrix score lp_o = bb
    quietly roctab kru_ge2 lp_o
    local auc_o = r(area)
    local sumopt = `sumopt' + (`auc_b' - `auc_o')
    local nok = `nok' + 1
}
local meanopt = `sumopt'/`nok'
local auc_corr = `auc_app' - `meanopt'

display _newline "  AUC apparente (combiné refit non-anuriques) = " %5.3f `auc_app'
display "  Optimisme moyen (B=`nok')                    = " %5.3f `meanopt'
display "  AUC optimism-corrected                       = " %5.3f `auc_corr'
restore

display _newline(2) "=== FIN do-file 09 ==="
