* ===========================================================================
* 12_nonanuric_robustness.do
*
* Prérequis : lancer 02_calculate_kru.do AVANT (données en mémoire).
*             Autonome ; ne modifie AUCUNE do-file existante.
*
* Objectif : réponse BLINDÉE au commentaire principal (KRU=0 assigné aux
*   anuriques). Cinq analyses de robustesse, toutes restreintes aux
*   NON-ANURIQUES (kru_pos==1), qui suppriment le problème à la racine.
*
* Sections :
*   A. Refit du modèle combiné SUR les non-anuriques + bootstrap Harrell
*      → le modèle tient debout SANS les anuriques (équation autonome)
*   B. Gradient dose-réponse par bandes de KRU + AUC bande borderline (1–3)
*      → UMOD/β2M discriminent 1.5 vs 2.5 (R3#1), pas juste 0 vs positif
*   C. Robustesse en excluant UMOD ≤ 2.0 (limite de détection)  → R3#6
*   D. Sensibilité au volume urinaire (exclure <400 mL)
*      → substitut de la sensibilité-seuil d'anurie (urine <200 non dispo)
*   E. Calibration du modèle chez les non-anuriques  → AE#3
* ===========================================================================

* --- Prérequis ---
foreach v in kru_daugirdas_35 umod labb2mprehd urinevolume {
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

* Modèle combiné PUBLIÉ (logit cohorte complète)
capture drop p_ge2_na
quietly logit kru_ge2 umod labb2mprehd
predict p_ge2_na, pr
label variable p_ge2_na "P(KRU≥2) — modèle combiné publié"

* ###########################################################################
* SECTION A — REFIT SUR NON-ANURIQUES + BOOTSTRAP HARRELL
* ###########################################################################
display _newline(2) "########################################################"
display              "  A. Refit du modèle combiné SUR les non-anuriques"
display              "     (équation autonome + AUC optimism-corrected)"
display              "########################################################"

* --- Équation autonome (coefficients + OR) ---
display _newline "  --- Logit KRU≥2 ~ UMOD + β2M (non-anuriques uniquement) ---"
logit kru_ge2 umod labb2mprehd if kru_pos==1
display _newline "  --- Odds ratios ---"
logit kru_ge2 umod labb2mprehd if kru_pos==1, or

* --- AUC apparente + bootstrap Harrell (optimism) ---
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

display _newline "  --- Validation interne (Harrell, B=`nok') ---"
display "  AUC apparente (refit non-anurique) = " %5.3f `auc_app'
display "  Optimisme moyen                    = " %5.3f `meanopt'
display "  AUC optimism-corrected             = " %5.3f `auc_corr'
restore

* ###########################################################################
* SECTION B — GRADIENT DOSE-RÉPONSE + AUC BANDE BORDERLINE
* ###########################################################################
display _newline(2) "########################################################"
display              "  B. Gradient par bandes de KRU + bande borderline"
display              "########################################################"

capture drop kru_band
gen byte kru_band = .
replace kru_band = 1 if kru_pos==1 & kru_daugirdas_35>0 & kru_daugirdas_35<=1
replace kru_band = 2 if kru_pos==1 & kru_daugirdas_35>1 & kru_daugirdas_35<=2
replace kru_band = 3 if kru_pos==1 & kru_daugirdas_35>2 & kru_daugirdas_35<=3
replace kru_band = 4 if kru_pos==1 & kru_daugirdas_35>3 & kru_daugirdas_35<=4
replace kru_band = 5 if kru_pos==1 & kru_daugirdas_35>4
label define bandlbl 1 "(0,1]" 2 "(1,2]" 3 "(2,3]" 4 "(3,4]" 5 ">4", replace
label values kru_band bandlbl

display _newline "  --- UMOD (ng/mL) par bande de KRU (non-anuriques) ---"
tabstat umod if kru_pos==1, by(kru_band) statistics(n p25 p50 p75) format(%7.2f)

display _newline "  --- β2M (mg/L) par bande de KRU (non-anuriques) ---"
tabstat labb2mprehd if kru_pos==1, by(kru_band) statistics(n p25 p50 p75) format(%7.2f)

display _newline "  --- Tendance monotone (Spearman, non-anuriques) ---"
spearman kru_daugirdas_35 umod        if kru_pos==1, stats(rho p)
spearman kru_daugirdas_35 labb2mprehd if kru_pos==1, stats(rho p)

display _newline "  --- AUC BANDE BORDERLINE : KRU 1–3 (discriminer <2 vs ≥2) ---"
display "      Réponse directe à R3#1 : distingue-t-on 1.5 de 2.5 ?"
quietly count if kru_pos==1 & kru_daugirdas_35>1 & kru_daugirdas_35<=3
display "      N dans la bande (1,3] = " r(N)
tab kru_ge2 if kru_pos==1 & kru_daugirdas_35>1 & kru_daugirdas_35<=3
display _newline "  UMOD seul (bande 1–3) :"
roctab kru_ge2 umod     if kru_pos==1 & kru_daugirdas_35>1 & kru_daugirdas_35<=3
display _newline "  β2M seul (bande 1–3) :"
roctab kru_ge2 neg_b2m  if kru_pos==1 & kru_daugirdas_35>1 & kru_daugirdas_35<=3
display _newline "  Combiné (bande 1–3) :"
roctab kru_ge2 p_ge2_na if kru_pos==1 & kru_daugirdas_35>1 & kru_daugirdas_35<=3

* ###########################################################################
* SECTION C — ROBUSTESSE : EXCLURE UMOD ≤ 2.0 (limite de détection)
* ###########################################################################
display _newline(2) "########################################################"
display              "  C. AUC non-anuriques en excluant UMOD ≤ 2.0 ng/mL"
display              "########################################################"

quietly count if kru_pos==1 & umod<=2.0
display "  Non-anuriques exclus (UMOD ≤ 2.0) = " r(N)
quietly count if kru_pos==1 & umod>2.0
display "  Non-anuriques analysés (UMOD > 2.0)= " r(N)

display _newline "  UMOD seul :"
roctab kru_ge2 umod     if kru_pos==1 & umod>2.0
display _newline "  β2M seul :"
roctab kru_ge2 neg_b2m  if kru_pos==1 & umod>2.0
display _newline "  Combiné :"
roctab kru_ge2 p_ge2_na if kru_pos==1 & umod>2.0

* ###########################################################################
* SECTION D — SENSIBILITÉ AU VOLUME URINAIRE (exclure <400 mL)
* ###########################################################################
display _newline(2) "########################################################"
display              "  D. AUC non-anuriques en excluant volume urinaire <400 mL"
display              "########################################################"

quietly count if kru_pos==1 & urinevolume<400 & !missing(urinevolume)
display "  Non-anuriques exclus (urine <400 mL) = " r(N)
quietly count if kru_pos==1 & urinevolume>=400 & !missing(urinevolume)
display "  Non-anuriques analysés (urine ≥400)  = " r(N)

display _newline "  UMOD seul :"
roctab kru_ge2 umod     if kru_pos==1 & urinevolume>=400
display _newline "  β2M seul :"
roctab kru_ge2 neg_b2m  if kru_pos==1 & urinevolume>=400
display _newline "  Combiné :"
roctab kru_ge2 p_ge2_na if kru_pos==1 & urinevolume>=400

* ###########################################################################
* SECTION E — CALIBRATION CHEZ LES NON-ANURIQUES
* ###########################################################################
display _newline(2) "########################################################"
display              "  E. Calibration du modèle chez les non-anuriques"
display              "########################################################"

* Prédicteur linéaire du modèle PUBLIÉ (cohorte complète)
capture drop lp_pub
quietly logit kru_ge2 umod labb2mprehd
predict lp_pub, xb

display _newline "  --- Pente de calibration (idéal = 1.0) ---"
display "      coef de lp_pub dans : logit KRU≥2 ~ lp_pub (non-anuriques)"
logit kru_ge2 lp_pub if kru_pos==1

display _newline "  --- Calibration-in-the-large (offset ; intercept idéal = 0) ---"
logit kru_ge2 if kru_pos==1, offset(lp_pub)

display _newline "  --- Hosmer-Lemeshow, modèle refitté non-anurique ---"
quietly logit kru_ge2 umod labb2mprehd if kru_pos==1
estat gof, group(10) table

display _newline(2) "=== FIN do-file 12 ==="
