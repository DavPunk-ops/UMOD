* ===========================================================================
* 16_kru_normalization.do
*
* Prérequis : lancer 02_calculate_kru.do AVANT (données en mémoire).
*             Autonome ; ne modifie AUCUNE do-file existante.
*
* Objectif : réponse à R3#5 — « l'outcome est en KRU/35L mais l'intro cite
*   /1.73 m² ; ces approches ne sont pas interchangeables. Watson est-elle
*   valide sur toute la composition corporelle ? Les résultats changent-ils
*   avec un KRU non-normalisé ou normalisé à la BSA ? »
*
* Stratégie : l'outcome primaire RESTE /35L (standard cinétique urée / HD
*   incrémentale). On montre que la discrimination des biomarqueurs est
*   ROBUSTE à la normalisation, en re-définissant KRU≥2 sous 3 échelles :
*     (a) brut (mL/min, non-normalisé)  → kru_daugirdas (existe déjà)
*     (b) /35L Watson (primaire)         → kru_daugirdas_35 (existe déjà)
*     (c) /1.73 m² BSA (DuBois)          → calculé ici
*
* Sections :
*   1. Composition corporelle (BMI, V_watson, BSA, imputations) → validité Watson
*   2. Corrélations entre les 3 échelles de KRU (attendu ρ élevé)
*   3. Reclassification au seuil ≥2 selon la normalisation
*   4. AUC (UMOD / β2M / combiné) sous chaque définition d'outcome
*   5. Tableau récapitulatif
* ===========================================================================

* --- Prérequis ---
foreach v in kru_daugirdas kru_daugirdas_35 umod labb2mprehd V_watson height weight {
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

capture drop neg_b2m
gen double neg_b2m = -labb2mprehd if !missing(labb2mprehd)
label variable neg_b2m "−β2M (orienté : haut = KRU≥2)"

* ###########################################################################
* SECTION 1 — COMPOSITION CORPORELLE (validité de Watson)
* ###########################################################################
display _newline(2) "########################################################"
display              "  1. Composition corporelle (gamme + imputations Watson)"
display              "########################################################"

capture drop bmi
gen double bmi = weight / ((height/100)^2) if !missing(weight, height)
label variable bmi "IMC (kg/m²)"

* BSA DuBois : 0.007184 × height(cm)^0.725 × weight(kg)^0.425
capture drop bsa
gen double bsa = 0.007184 * (height^0.725) * (weight^0.425) if !missing(weight, height)
label variable bsa "Surface corporelle DuBois (m²)"

display _newline "  --- Distribution IMC / V_watson / BSA ---"
summarize bmi V_watson bsa, detail

display _newline "  --- Patients avec V_watson imputé à 35L (poids/taille manquant) ---"
count if V_watson == 35
display "  N imputés (V=35L exactement) = " r(N)

display _newline "  --- Gamme de composition corporelle (pour la validité de Watson) ---"
quietly summarize bmi
display "  IMC : min=" %4.1f r(min) "  max=" %4.1f r(max)
quietly summarize V_watson
display "  V_watson (L) : min=" %4.1f r(min) "  max=" %4.1f r(max)

* ###########################################################################
* SECTION 2 — KRU SOUS 3 NORMALISATIONS + CORRÉLATIONS
* ###########################################################################
display _newline(2) "########################################################"
display              "  2. KRU brut / /35L / /1.73m²  +  corrélations"
display              "########################################################"

* (a) brut : kru_daugirdas (existe, =0 chez anuriques)
* (b) /35L : kru_daugirdas_35 (existe)
* (c) /1.73m² BSA : kru_daugirdas × 1.73/BSA   (anuriques =0 conservés)
capture drop kru_bsa
gen double kru_bsa = kru_daugirdas * (1.73 / bsa) if !missing(kru_daugirdas, bsa)
replace     kru_bsa = 0 if kru_daugirdas == 0 & !missing(bsa)
label variable kru_bsa "KRU Daugirdas (mL/min/1.73m²)"

display _newline "  --- Résumé des 3 échelles (>0 seulement) ---"
foreach v in kru_daugirdas kru_daugirdas_35 kru_bsa {
    quietly summarize `v' if `v'>0 & !missing(`v'), detail
    display "  `v' : N>0=" r(N) "  médiane=" %5.2f r(p50) "  moyenne=" %5.2f r(mean)
}

display _newline "  --- Corrélations Spearman entre les 3 échelles ---"
spearman kru_daugirdas kru_daugirdas_35 kru_bsa, stats(rho p)

* ###########################################################################
* SECTION 3 — RECLASSIFICATION AU SEUIL ≥2 SELON LA NORMALISATION
* ###########################################################################
display _newline(2) "########################################################"
display              "  3. Reclassification KRU≥2 selon la normalisation"
display              "########################################################"

capture drop ge2_raw ge2_35 ge2_bsa
gen byte ge2_raw = (kru_daugirdas    >= 2) if !missing(kru_daugirdas)
gen byte ge2_35  = (kru_daugirdas_35 >= 2) if !missing(kru_daugirdas_35)
gen byte ge2_bsa = (kru_bsa          >= 2) if !missing(kru_bsa)
label variable ge2_raw "KRU≥2 (brut)"
label variable ge2_35  "KRU≥2 (/35L, primaire)"
label variable ge2_bsa "KRU≥2 (/1.73m²)"

display _newline "  --- /35L (primaire)  vs  brut ---"
tab ge2_35 ge2_raw
display _newline "  --- /35L (primaire)  vs  /1.73m² ---"
tab ge2_35 ge2_bsa

quietly count if ge2_35 != ge2_raw & !missing(ge2_35, ge2_raw)
display _newline "  Reclassés /35L→brut     : " r(N)
quietly count if ge2_35 != ge2_bsa & !missing(ge2_35, ge2_bsa)
display "  Reclassés /35L→/1.73m²  : " r(N)

* ###########################################################################
* SECTION 4 — AUC SOUS CHAQUE DÉFINITION D'OUTCOME
* ###########################################################################
display _newline(2) "########################################################"
display              "  4. AUC (UMOD / β2M / combiné) sous chaque normalisation"
display              "########################################################"

capture program drop auc_norm
program define auc_norm
    args outc label
    display _newline "  === Outcome : `label' ==="
    quietly count if `outc'==1 & !missing(`outc')
    local nev = r(N)
    quietly count if !missing(`outc')
    local ntot = r(N)
    display "  N=`ntot'  événements(KRU≥2)=`nev'"
    display "  UMOD seul :"
    quietly roctab `outc' umod
    display "     AUC = " %5.3f r(area)
    display "  β2M seul :"
    quietly roctab `outc' neg_b2m
    display "     AUC = " %5.3f r(area)
    display "  Combiné (logit UMOD+β2M) :"
    capture drop _pnorm
    quietly logit `outc' umod labb2mprehd
    quietly predict _pnorm, pr
    quietly roctab `outc' _pnorm
    display "     AUC = " %5.3f r(area)
    capture drop _pnorm
end

auc_norm ge2_raw "KRU≥2 brut (mL/min)"
auc_norm ge2_35  "KRU≥2 /35L (primaire)"
auc_norm ge2_bsa "KRU≥2 /1.73m² (BSA)"

* ###########################################################################
* SECTION 5 — RÉCAPITULATIF
* ###########################################################################
display _newline(2) "########################################################"
display              "  5. Récapitulatif : la discrimination est-elle stable ?"
display              "########################################################"
display "  → Si les AUC combinées sont quasi identiques sous les 3 échelles"
display "    et que peu de patients sont reclassés, la performance des"
display "    biomarqueurs ne dépend PAS du choix de normalisation."

display _newline(2) "=== FIN do-file 16 ==="
