* ===========================================================================
* 09_AE1_nonanuric_logit.do
*
* Prérequis : lancer 02_calculate_kru.do AVANT (données en mémoire).
*             Autonome ; ne modifie AUCUNE do-file existante.
*
* Objectif : réponse CIBLÉE à AE#1 (« la performance reflète-t-elle surtout la
*   séparation anurique/non-anurique ? »). Discrimination dans le sous-groupe à
*   KRU mesuré (non-anuriques), présentée de façon COHÉRENTE avec la Table 5 :
*
*   (i)  ANALYSE PRIMAIRE — le modèle PUBLIÉ (logit sur la cohorte complète)
*        est ÉVALUÉ chez les non-anuriques → AUC rapportées dans le tableau
*        AE#1. C'est le MÊME modèle que la Table 5 (seuils publiés).
*   (ii) CONTRÔLE DE ROBUSTESSE — le modèle est REFITTÉ chez les non-anuriques
*        → OR spécifiques au sous-groupe (contribution indépendante de chaque
*        marqueur) + AUC (qui COÏNCIDE avec (i)) + validation bootstrap.
*
* Sections :
*   1. (i)  Modèle publié évalué chez les non-anuriques : AUC UMOD/β2M/combiné
*   2. (ii) Refit chez les non-anuriques : OR + AUC (coïncide avec (i))
*   3. Bootstrap Harrell du refit (optimism-corrected)
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

* Échantillon d'analyse : non-anuriques avec données complètes (N=87), afin que
* TOUTES les AUC (y compris UMOD seul) soient sur le même N que le tableau.
capture drop insample
gen byte insample = (kru_pos==1 & !missing(umod, labb2mprehd, kru_ge2))

* ###########################################################################
* SECTION 1 — (i) MODÈLE PUBLIÉ ÉVALUÉ CHEZ LES NON-ANURIQUES  [PRIMAIRE]
*   Modèle publié = logit ajusté sur la COHORTE COMPLÈTE ; on évalue ses
*   prédictions dans le sous-groupe non-anurique. Ce sont les AUC rapportées
*   dans le tableau de la réponse AE#1, et le même modèle que la Table 5.
* ###########################################################################
display _newline(2) "########################################################"
display              "  1. (i) Modèle PUBLIÉ évalué chez les non-anuriques"
display              "########################################################"

quietly logit kru_ge2 umod labb2mprehd           /* cohorte complète */
capture drop p_pub
predict p_pub, pr
label variable p_pub "P(KRU≥2) — modèle publié (cohorte complète)"

display _newline "  --- AUC chez les non-anuriques (N=87, données complètes) ---"
display _newline "  UMOD seul :"
roctab kru_ge2 umod    if insample
display _newline "  β2M seul :"
roctab kru_ge2 neg_b2m if insample
display _newline "  Combiné (modèle publié) :"
roctab kru_ge2 p_pub   if insample

* ###########################################################################
* SECTION 2 — (ii) REFIT CHEZ LES NON-ANURIQUES  [CONTRÔLE DE ROBUSTESSE]
*   Refit du logit UNIQUEMENT sur les non-anuriques → OR spécifiques au
*   sous-groupe (contribution indépendante). L'AUC combinée doit COÏNCIDER
*   avec la Section 1 (démonstration que (i) = (ii)).
* ###########################################################################
display _newline(2) "########################################################"
display              "  2. (ii) Refit chez les non-anuriques (robustesse)"
display              "########################################################"

display _newline "  --- Coefficients (refit) ---"
logit kru_ge2 umod labb2mprehd if kru_pos==1
display _newline "  --- Odds ratios (contribution indépendante) ---"
logit kru_ge2 umod labb2mprehd if kru_pos==1, or

capture drop p_refit
predict p_refit, pr
label variable p_refit "P(KRU≥2) — modèle refitté non-anuriques"

display _newline "  --- AUC combiné refitté (doit coïncider avec (i)) ---"
roctab kru_ge2 p_refit if insample

* ###########################################################################
* SECTION 3 — VALIDATION INTERNE DU REFIT (BOOTSTRAP HARRELL)
*   Optimism correction de l'AUC du modèle refitté sur les non-anuriques.
* ###########################################################################
display _newline(2) "########################################################"
display              "  3. Bootstrap Harrell du combiné refitté (optimism)"
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

display _newline "  AUC apparente (refit non-anuriques) = " %5.3f `auc_app'
display "  Optimisme moyen (B=`nok')            = " %5.3f `meanopt'
display "  AUC optimism-corrected               = " %5.3f `auc_corr'
restore

display _newline(2) "=== FIN do-file 09 ==="
