* ===========================================================================
* 11_anuria400_sensitivity.do
*
* Prérequis : lancer 02_calculate_kru.do AVANT (données en mémoire).
*             Autonome ; ne modifie AUCUNE do-file existante.
*
* Objectif : réponse DIRECTE à AE#1 / R3#4 — refaire l'analyse primaire en
*   redéfinissant l'anurie à un seuil PLUS HAUT (<400 mL/j au lieu de <200).
*   (Un seuil plus bas — 100 mL — est infaisable : urine non récoltée <200 mL.)
*
*   Sous la définition 400 mL, les patients à diurèse 200–399 mL sont reclassés
*   « anuriques » et se voient assigner KRU=0.
*
* Sections :
*   A. Identifier les patients reclassés (200–399 mL) et leur KRU réel
*   B. Créer les variables sous définition 400 mL + vérifier si la
*      classification binaire KRU≥2 change réellement
*   C. Refaire les AUC primaires (UMOD / β2M / combiné) sous 400 mL vs 200 mL
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

* --- Variables dérivées (définition originale, 200 mL) ---
capture confirm variable kru_pos
if _rc  gen byte kru_pos = (kru_daugirdas_35 > 0)  if !missing(kru_daugirdas_35)
capture confirm variable kru_ge2
if _rc  gen byte kru_ge2 = (kru_daugirdas_35 >= 2) if !missing(kru_daugirdas_35)

capture drop neg_b2m
gen double neg_b2m = -labb2mprehd if !missing(labb2mprehd)
label variable neg_b2m "−β2M (orienté : haut = KRU≥2)"

* ###########################################################################
* SECTION A — PATIENTS RECLASSÉS (diurèse 200–399 mL)
* ###########################################################################
display _newline(2) "########################################################"
display              "  A. Patients reclassés « anuriques » sous seuil 400 mL"
display              "     (diurèse 200–399 mL, actuellement non-anuriques)"
display              "########################################################"

quietly count if kru_pos==1 & urinevolume<400 & !missing(urinevolume)
display _newline "  N reclassés (non-anuriques avec urine <400 mL) = " r(N)

display _newline "  --- Leur KRU réel mesuré et statut KRU≥2 ---"
list id urinevolume kru_daugirdas_35 kru_ge2 ///
    if kru_pos==1 & urinevolume<400 & !missing(urinevolume), noobs sepby(kru_ge2)

display _newline "  --- Combien d'entre eux ont KRU≥2 ? ---"
tab kru_ge2 if kru_pos==1 & urinevolume<400 & !missing(urinevolume)

* ###########################################################################
* SECTION B — DÉFINITION 400 mL + VÉRIFICATION DU CHANGEMENT BINAIRE
* ###########################################################################
display _newline(2) "########################################################"
display              "  B. Définition 400 mL — la classification KRU≥2 change-t-elle ?"
display              "########################################################"

* KRU sous définition 400 mL : 0 si diurèse <400 mL
capture drop kru400
gen double kru400 = kru_daugirdas_35
replace kru400 = 0 if urinevolume<400 & !missing(urinevolume)
label variable kru400 "KRU sous définition anurie <400 mL"

capture drop kru_pos400
gen byte kru_pos400 = (kru400 > 0)  if !missing(kru400)
capture drop kru_ge2_400
gen byte kru_ge2_400 = (kru400 >= 2) if !missing(kru400)
label variable kru_ge2_400 "KRU≥2 sous définition 400 mL"

display _newline "  --- Effectifs anuriques : 200 mL vs 400 mL ---"
quietly count if kru_pos==0
local an200 = r(N)
quietly count if kru_pos400==0
local an400 = r(N)
display "  Anuriques (déf. 200 mL) = `an200'"
display "  Anuriques (déf. 400 mL) = `an400'"

display _newline "  --- Patients dont la classification KRU≥2 CHANGE ---"
quietly count if kru_ge2 != kru_ge2_400 & !missing(kru_ge2, kru_ge2_400)
display "  N avec classification KRU≥2 modifiée = " r(N)
display "  (0 = la discrimination pour KRU≥2 est inchangée par construction,"
display "   car les patients reclassés étaient déjà KRU<2)"

* ###########################################################################
* SECTION C — AUC PRIMAIRES SOUS 400 mL vs 200 mL
* ###########################################################################
display _newline(2) "########################################################"
display              "  C. AUC pour KRU≥2 : définition 200 mL vs 400 mL"
display              "########################################################"

display _newline "  ===== Définition ORIGINALE (anurie <200 mL) ====="
display "  UMOD seul :"
roctab kru_ge2 umod
display "  β2M seul :"
roctab kru_ge2 neg_b2m
display "  Combiné :"
quietly logit kru_ge2 umod labb2mprehd
capture drop p200
predict p200, pr
roctab kru_ge2 p200

display _newline "  ===== Définition ALTERNATIVE (anurie <400 mL) ====="
display "  UMOD seul :"
roctab kru_ge2_400 umod
display "  β2M seul :"
roctab kru_ge2_400 neg_b2m
display "  Combiné :"
quietly logit kru_ge2_400 umod labb2mprehd
capture drop p400
predict p400, pr
roctab kru_ge2_400 p400

display _newline(2) "=== FIN do-file 11 ==="
