* ===========================================================================
* 12_nonanuric_reanalysis.do
*
* Prérequis : lancer 02_calculate_kru.do AVANT ce do-file (données en mémoire).
*             N'utilise QUE des variables déjà créées par 02 (+ recrée localement
*             kru_pos, kru_ge2, neg_b2m, b2m_neg2, p_ge2_na si besoin).
*             NE MODIFIE AUCUNE do-file existante.
*
* Objectif : répondre au commentaire PRINCIPAL des reviewers (major revision).
*   Les 62 patients anuriques se sont vu ASSIGNER KRU=0 (non mesuré). Les
*   reviewers (AE#1 ; R3 majeurs #1-4) craignent que la discrimination reflète
*   surtout la séparation « facile » anuriques vs non-anuriques.
*   → On refait TOUTES les analyses primaires UNIQUEMENT chez les non-anuriques
*     (KRU réellement mesuré, kru_pos==1), population cliniquement pertinente.
*
* Sections :
*   A. Effectifs + description UMOD/β2M par KRU≥2 (non-anuriques)
*   B. AUC UMOD / β2M / combiné + comparaisons DeLong (non-anuriques)
*   C. Stratégie grey-zone — SEUILS PUBLIÉS appliqués aux non-anuriques :
*        NPV, PPV, faux rule-out, faux rule-in, taille des 3 zones
*        (réponse directe à R3#3 : « les 85% sont-ils surestimés ? »)
*   D. Contribution indépendante UMOD vs β2M pour le KRU CONTINU (non-anuriques)
*        corrélations partielles + OLS emboîtés
*   E. Uromoduline près de la limite de détection (2.0 ng/mL)
* ===========================================================================

* --- Vérification des prérequis ---
foreach v in kru_daugirdas_35 umod labb2mprehd {
    capture confirm variable `v'
    if _rc {
        display as error "ERREUR : '`v'' absente — lance d'abord 02_calculate_kru.do"
        exit 111
    }
}
display _newline "=== Prérequis OK ==="

* --- Variables dérivées (recréées si absentes, sans écraser l'existant) ---
capture confirm variable kru_pos
if _rc {
    gen byte kru_pos = (kru_daugirdas_35 > 0) if !missing(kru_daugirdas_35)
    label variable kru_pos "KRU > 0 (1 = non-anurique)"
    label define krupos08 0 "Anurique (KRU=0)" 1 "Non-anurique (KRU>0)", replace
    label values kru_pos krupos08
}
capture confirm variable kru_ge2
if _rc {
    gen byte kru_ge2 = (kru_daugirdas_35 >= 2) if !missing(kru_daugirdas_35)
    label variable kru_ge2 "KRU >= 2 mL/min/35L"
    label define kruge208 0 "KRU < 2" 1 "KRU >= 2", replace
    label values kru_ge2 kruge208
}

* β2M orienté pour ROC (valeur haute → KRU≥2, car β2M inversement associée)
capture drop neg_b2m
gen double neg_b2m = -labb2mprehd if !missing(labb2mprehd)
label variable neg_b2m "−β2M (orienté ROC)"

* Transformation β2M^-2 (identique au modèle continu du do-file 06)
capture drop b2m_neg2
gen double b2m_neg2 = (labb2mprehd/10)^(-2) if !missing(labb2mprehd)
label variable b2m_neg2 "(β2M/10)^-2"

* Probabilité prédite du modèle combiné PUBLIÉ (logit sur cohorte COMPLÈTE),
* évaluée ensuite dans le sous-groupe non-anurique → « le modèle publié
* fonctionne-t-il chez les patients à KRU mesuré ? »
capture drop p_ge2_na
quietly logit kru_ge2 umod labb2mprehd
predict p_ge2_na, pr
label variable p_ge2_na "P(KRU≥2) — modèle combiné publié (cohorte complète)"

* ###########################################################################
* SECTION A — DESCRIPTION CHEZ LES NON-ANURIQUES
* ###########################################################################
display _newline(2) "########################################################"
display              "  A. NON-ANURIQUES : effectifs et distributions"
display              "########################################################"

display _newline "  --- Effectifs KRU<2 vs KRU≥2 chez les non-anuriques ---"
tab kru_ge2 if kru_pos == 1, miss

display _newline "  --- UMOD (ng/mL) : KRU<2 vs KRU≥2 (non-anuriques) ---"
tabstat umod if kru_pos == 1, by(kru_ge2) ///
    statistics(n p25 p50 p75) format(%7.2f)
ranksum umod if kru_pos == 1, by(kru_ge2)

display _newline "  --- β2M (mg/L) : KRU<2 vs KRU≥2 (non-anuriques) ---"
tabstat labb2mprehd if kru_pos == 1, by(kru_ge2) ///
    statistics(n p25 p50 p75) format(%7.2f)
ranksum labb2mprehd if kru_pos == 1, by(kru_ge2)

* ###########################################################################
* SECTION B — AUC ET DeLong CHEZ LES NON-ANURIQUES
* ###########################################################################
display _newline(2) "########################################################"
display              "  B. AUC (non-anuriques) — UMOD / β2M / combiné"
display              "########################################################"

display _newline "  --- AUC UMOD seul (non-anuriques) ---"
roctab kru_ge2 umod if kru_pos == 1

display _newline "  --- AUC β2M seul (non-anuriques) ---"
roctab kru_ge2 neg_b2m if kru_pos == 1

display _newline "  --- AUC modèle combiné publié (non-anuriques) ---"
roctab kru_ge2 p_ge2_na if kru_pos == 1

display _newline "  --- DeLong : UMOD vs β2M (non-anuriques) ---"
roccomp kru_ge2 umod neg_b2m if kru_pos == 1 & !missing(umod, neg_b2m), summary

display _newline "  --- DeLong : combiné vs UMOD (non-anuriques) ---"
roccomp kru_ge2 p_ge2_na umod if kru_pos == 1 & !missing(p_ge2_na, umod), summary

display _newline "  --- DeLong : combiné vs β2M (non-anuriques) ---"
roccomp kru_ge2 p_ge2_na neg_b2m if kru_pos == 1 & !missing(p_ge2_na, neg_b2m), summary

* ###########################################################################
* SECTION C — STRATÉGIE GREY-ZONE (SEUILS PUBLIÉS) CHEZ LES NON-ANURIQUES
*   Réponse directe à R3#3 : quel est le VRAI rendement clinique de la règle
*   quand on l'applique aux seuls patients qui auraient une récolte d'urine ?
*
*   NOTE IMPORTANTE (piège Stata) : les valeurs manquantes sont traitées comme
*   +∞, donc « x >= seuil » est VRAI si x est manquant. Les 2 non-anuriques
*   sans β2M (donc sans p_ge2_na, tous deux KRU≥2) seraient comptés à tort dans
*   le rule-in. TOUS les comptages incluent donc « & !missing(...) ».
* ###########################################################################
display _newline(2) "########################################################"
display              "  C. GREY-ZONE (seuils publiés) chez les non-anuriques"
display              "########################################################"

* --- C1. UMOD seul : rule-out <7.0 ; rule-in ≥14.0 -------------------------
display _newline "  ============ UMOD seul (seuils 7.0 / 14.0) ============"
quietly count if kru_pos==1 & !missing(umod)
local N_na = r(N)

* Rule-out zone (<7.0)
quietly count if kru_pos==1 & umod<7.0 & !missing(umod)
local ro     = r(N)
quietly count if kru_pos==1 & umod<7.0 & kru_ge2==0 & !missing(umod)
local ro_ok  = r(N)
quietly count if kru_pos==1 & umod<7.0 & kru_ge2==1 & !missing(umod)
local ro_bad = r(N)

* Rule-in zone (≥14.0)
quietly count if kru_pos==1 & umod>=14.0 & !missing(umod)
local ri     = r(N)
quietly count if kru_pos==1 & umod>=14.0 & kru_ge2==1 & !missing(umod)
local ri_ok  = r(N)
quietly count if kru_pos==1 & umod>=14.0 & kru_ge2==0 & !missing(umod)
local ri_bad = r(N)

* Grey zone (7.0–<14.0)
quietly count if kru_pos==1 & umod>=7.0 & umod<14.0 & !missing(umod)
local gz = r(N)

display "  N non-anuriques (UMOD dispo) = `N_na'"
display _newline "  Rule-out (<7.0)  : n=`ro'  (" %4.1f 100*`ro'/`N_na' "%)   " ///
    "vrais KRU<2=`ro_ok'  FAUX rule-out (KRU≥2 manqués)=`ro_bad'"
if `ro'>0 display "      → NPV = " %5.1f 100*`ro_ok'/`ro' "%"
display "  Grey zone (7.0–14): n=`gz'  (" %4.1f 100*`gz'/`N_na' "%)"
display "  Rule-in (≥14.0)  : n=`ri'  (" %4.1f 100*`ri'/`N_na' "%)   " ///
    "vrais KRU≥2=`ri_ok'  FAUX rule-in (KRU<2 classés ≥2)=`ri_bad'"
if `ri'>0 display "      → PPV = " %5.1f 100*`ri_ok'/`ri' "%"
display _newline "  Classés sans récolte (rule-out+rule-in) = " ///
    %4.1f 100*(`ro'+`ri')/`N_na' "%   (vs 75.5% en cohorte complète)"

* --- C2. Modèle combiné : rule-out <0.24 ; rule-in ≥0.55 -------------------
display _newline "  ======= Combiné UMOD+β2M (seuils 0.24 / 0.55) ======="
quietly count if kru_pos==1 & !missing(p_ge2_na)
local M_na = r(N)

quietly count if kru_pos==1 & p_ge2_na<0.24 & !missing(p_ge2_na)
local cro     = r(N)
quietly count if kru_pos==1 & p_ge2_na<0.24 & kru_ge2==0 & !missing(p_ge2_na)
local cro_ok  = r(N)
quietly count if kru_pos==1 & p_ge2_na<0.24 & kru_ge2==1 & !missing(p_ge2_na)
local cro_bad = r(N)

quietly count if kru_pos==1 & p_ge2_na>=0.55 & !missing(p_ge2_na)
local cri     = r(N)
quietly count if kru_pos==1 & p_ge2_na>=0.55 & kru_ge2==1 & !missing(p_ge2_na)
local cri_ok  = r(N)
quietly count if kru_pos==1 & p_ge2_na>=0.55 & kru_ge2==0 & !missing(p_ge2_na)
local cri_bad = r(N)

quietly count if kru_pos==1 & p_ge2_na>=0.24 & p_ge2_na<0.55 & !missing(p_ge2_na)
local cgz = r(N)

display "  N non-anuriques (modèle dispo) = `M_na'"
display _newline "  Rule-out (<0.24) : n=`cro'  (" %4.1f 100*`cro'/`M_na' "%)   " ///
    "vrais KRU<2=`cro_ok'  FAUX rule-out=`cro_bad'"
if `cro'>0 display "      → NPV = " %5.1f 100*`cro_ok'/`cro' "%"
display "  Grey zone        : n=`cgz'  (" %4.1f 100*`cgz'/`M_na' "%)"
display "  Rule-in (≥0.55)  : n=`cri'  (" %4.1f 100*`cri'/`M_na' "%)   " ///
    "vrais KRU≥2=`cri_ok'  FAUX rule-in=`cri_bad'"
if `cri'>0 display "      → PPV = " %5.1f 100*`cri_ok'/`cri' "%"

* --- IC binomiaux (Wilson) pour Table 5 ------------------------------------
* Seuils EXTERNES (0.24/0.55, dérivés de la cohorte complète) appliqués au
* sous-groupe → NPV et PPV sont de simples proportions → IC analytique
* (Wilson), sans bootstrap. cii proportions #obs #succès.
if `cro'>0 {
    display _newline "  --- NPV rule-out : IC binomial de Wilson ---"
    cii proportions `cro' `cro_ok', wilson
}
if `cri'>0 {
    display "  --- PPV rule-in : IC binomial de Wilson ---"
    cii proportions `cri' `cri_ok', wilson
}

display _newline "  Classés sans récolte (rule-out+rule-in) = " ///
    %4.1f 100*(`cro'+`cri')/`M_na' "%   (vs 85.1% en cohorte complète)"

* ###########################################################################
* SECTION D — CONTRIBUTION INDÉPENDANTE POUR LE KRU CONTINU (non-anuriques)
*   UMOD garde-t-il une information indépendante de β2M une fois les
*   anuriques exclus ? (filet de sécurité du message central)
* ###########################################################################
display _newline(2) "########################################################"
display              "  D. Info indépendante UMOD vs β2M — KRU continu"
display              "     (non-anuriques uniquement)"
display              "########################################################"

display _newline "  --- Corrélations partielles (KRU continu, non-anuriques) ---"
display "      pcorr : UMOD ajusté sur β2M, et β2M ajusté sur UMOD"
pcorr kru_daugirdas_35 umod labb2mprehd if kru_pos == 1

display _newline "  --- OLS : KRU ~ UMOD seul (non-anuriques) ---"
regress kru_daugirdas_35 umod if kru_pos == 1, vce(robust)

display _newline "  --- OLS : KRU ~ (β2M/10)^-2 seul (non-anuriques) ---"
regress kru_daugirdas_35 b2m_neg2 if kru_pos == 1, vce(robust)

display _newline "  --- OLS : KRU ~ UMOD + (β2M/10)^-2 (non-anuriques) ---"
display "      → p-value de UMOD = contribution indépendante de β2M"
regress kru_daugirdas_35 umod b2m_neg2 if kru_pos == 1, vce(robust)

* ###########################################################################
* SECTION E — UROMODULINE PRÈS DE LA LIMITE DE DÉTECTION (2.0 ng/mL)
*   Réponse à R3#6 : combien de valeurs au niveau/en dessous du seuil ?
* ###########################################################################
display _newline(2) "########################################################"
display              "  E. UMOD près de la limite de détection (2.0 ng/mL)"
display              "########################################################"

quietly count if !missing(umod)
display "  N total avec UMOD dosé          = " r(N)
quietly count if umod == 0
display "  UMOD = 0                        = " r(N)
quietly count if umod <= 2.0 & !missing(umod)
display "  UMOD ≤ 2.0 (limite de détection)= " r(N)
quietly count if umod <= 2.0 & kru_pos==0
display "     dont anuriques               = " r(N)
quietly count if umod <= 2.0 & kru_pos==1
display "     dont non-anuriques           = " r(N)

display _newline "  --- Distribution UMOD par statut anurique ---"
tabstat umod, by(kru_pos) statistics(n min p25 p50 p75 max) format(%7.2f)

display _newline(2) "=== FIN do-file 12 ==="
