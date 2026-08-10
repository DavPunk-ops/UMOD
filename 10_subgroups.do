* ===========================================================================
* 10_subgroups.do
*
* Prérequis : lancer 02_calculate_kru.do AVANT (données en mémoire).
*             Autonome ; ne modifie AUCUNE do-file existante.
*
* Objectif : réponse à AE#1 (dernière question) — « Is it different depending
*   sex? Age? comorbidities? » — et bonus R1 (sous-groupe ADPKD).
*   Performance discriminante du modèle (UMOD, β2M, combiné) pour KRU≥2,
*   évaluée dans chaque sous-groupe de la cohorte complète.
*   Le modèle combiné est le modèle PUBLIÉ (logit sur cohorte complète),
*   dont on évalue les prédictions dans chaque sous-groupe.
*
* Sous-groupes : sexe · âge (split médian) · diabète · Charlson (split médian)
*                · exclusion ADPKD (sensibilité R1)
* ===========================================================================

* --- Prérequis ---
foreach v in kru_daugirdas_35 umod labb2mprehd age {
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
capture confirm variable female
if _rc {
    capture drop female
    gen byte female = (sex == 1) if !missing(sex)
}
capture drop neg_b2m
gen double neg_b2m = -labb2mprehd if !missing(labb2mprehd)

* Modèle combiné PUBLIÉ (logit cohorte complète)
capture drop p_ge2
quietly logit kru_ge2 umod labb2mprehd
predict p_ge2, pr
label variable p_ge2 "P(KRU≥2) — modèle combiné publié"

* Médianes pour les splits
quietly summarize age, detail
local agemed = r(p50)
capture confirm variable charlson
local has_char = (_rc==0)
if `has_char' {
    quietly summarize charlson, detail
    local charmed = r(p50)
}

* Détection automatique du code ADPKD dans kidneydisease
local adpkd = .
capture confirm variable kidneydisease
if _rc==0 {
    quietly levelsof kidneydisease, local(kdcodes)
    foreach c of local kdcodes {
        local lbl : label (kidneydisease) `c'
        if strpos(lower("`lbl'"),"adpkd")>0 | strpos(lower("`lbl'"),"polycystic")>0 {
            local adpkd = `c'
        }
    }
}

* --- Programme utilitaire : AUC des 3 prédicteurs dans un sous-groupe ---
capture program drop auc_sub
program define auc_sub
    args label cond
    quietly count if `cond' & !missing(kru_ge2)
    local n = r(N)
    quietly count if `cond' & kru_ge2==1 & !missing(kru_ge2)
    local nev = r(N)
    quietly roctab kru_ge2 umod    if `cond'
    local a_u = r(area)
    quietly roctab kru_ge2 neg_b2m if `cond'
    local a_b = r(area)
    quietly roctab kru_ge2 p_ge2   if `cond'
    local a_c = r(area)
    display "  " %-24s "`label'" "  N=" %3.0f `n' " (KRU≥2=" %2.0f `nev' ")" ///
        "   UMOD=" %5.3f `a_u' "   β2M=" %5.3f `a_b' "   Combiné=" %5.3f `a_c'
end

* ###########################################################################
display _newline(2) "########################################################"
display              "  AUC pour KRU≥2 par sous-groupe (cohorte complète)"
display              "  Colonnes : UMOD seul · β2M seul · modèle combiné"
display              "########################################################"

display _newline "  --- Référence : cohorte entière ---"
auc_sub "All patients"        "1"

display _newline "  --- Sexe ---"
auc_sub "Men"                 "female==0"
auc_sub "Women"               "female==1"

display _newline "  --- Âge (split médian = `agemed' ans) ---"
auc_sub "Age < median"        "age<`agemed'"
auc_sub "Age >= median"       "age>=`agemed'"

capture confirm variable dm
if _rc==0 {
    display _newline "  --- Diabète ---"
    auc_sub "No diabetes"     "dm==0"
    auc_sub "Diabetes"        "dm==1"
}
else display _newline "  (variable 'dm' absente — sous-groupe diabète ignoré)"

if `has_char' {
    display _newline "  --- Comorbidité (Charlson, split médian = `charmed') ---"
    auc_sub "Charlson < median"  "charlson<`charmed'"
    auc_sub "Charlson >= median" "charlson>=`charmed'"
}

if `adpkd' != . {
    display _newline "  --- Sensibilité : exclusion ADPKD (code `adpkd') ---"
    auc_sub "Excluding ADPKD"  "kidneydisease!=`adpkd'"
    quietly count if kidneydisease==`adpkd'
    display "     (N ADPKD exclus = " r(N) ")"
}
else {
    display _newline "  (code ADPKD non détecté dans kidneydisease — voir tab ci-dessous)"
    capture noisily tab kidneydisease
}

display _newline(2) "=== FIN do-file 10 ==="
