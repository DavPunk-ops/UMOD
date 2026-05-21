* ===========================================================================
* 02_calculate_kru.do
* Objectif : Calcul du KRU (mL/min) et normalisation pour 35L (Watson)
* ===========================================================================

local path "C:\Users\dajs\OneDrive - HOPITAUX UNIVERSITAIRES DE GENEVE\recherche\RKF\UMOD\stata\main prospective study\with Claude"

use "`path'\database_umod.dta", clear
* database_umod.dta = fichier source (output de 01_merge) — jamais modifié
* database_analysis.dta = fichier d'analyse enrichi (output de ce do-file)

* ── 1. Durée de collecte précise ────────────────────────────────
* Parser urinestart et urineend (format "HH:MM")
gen start_h   = real(substr(urinestart, 1, 2))
gen start_min = real(substr(urinestart, 4, 2))
gen start_tot = start_h * 60 + start_min

gen end_h   = real(substr(urineend, 1, 2))
gen end_min = real(substr(urineend, 4, 2))
gen end_tot = end_h * 60 + end_min

* T (minutes) = toujours ~24h car récolte dans les 24h précédant la dialyse
gen T_min = 1440 + (end_tot - start_tot)

* Vérification : T doit être strictement positif
assert T_min > 0 if !missing(T_min)

drop start_h start_min start_tot end_h end_min end_tot

* ── 1b. Comparaison bloodurea vs labureaprehd ───────────────────
* Objectif : vérifier que les deux mesures d'urée plasmatique ne
* diffèrent QUE lorsque les dates de prélèvement diffèrent.
*   - bloodurea       : urée au moment de la fin de récolte urinaire (urinedate)
*   - labureaprehd    : urée du laboratoire pré-dialyse        (datevisit)
*
* Si urinedate == datevisit : les deux mesures concernent le même
* prélèvement → elles doivent être (quasi) identiques.
* Si urinedate != datevisit : prélèvements à des moments différents
* → un écart est attendu et bloodurea reste la mesure cohérente
*   temporellement avec la récolte.

gen byte same_date = (urinedate == datevisit) if !missing(urinedate) & !missing(datevisit)
label variable same_date "urinedate == datevisit (1=oui, 0=non)"
label define same_date_lbl 0 "Dates différentes" 1 "Mêmes dates"
label values same_date same_date_lbl

display _newline "=== Distribution : urinedate == datevisit ? ==="
tab same_date, miss

* Différence entre les deux mesures d'urée plasmatique
gen delta_urea = labureaprehd - bloodurea
label variable delta_urea "labureaprehd − bloodurea (mmol/L)"

* ─ Vue globale (tous patients avec les deux valeurs) ──────────
display _newline "=== Comparaison globale ==="
summarize bloodurea labureaprehd delta_urea, detail
corr bloodurea labureaprehd

* ─ Stratification selon concordance des dates ────────────────
display _newline "=== Mêmes dates (urinedate == datevisit) ==="
summarize bloodurea labureaprehd delta_urea if same_date == 1, detail
corr bloodurea labureaprehd if same_date == 1

display _newline "=== Dates différentes (urinedate != datevisit) ==="
summarize bloodurea labureaprehd delta_urea if same_date == 0, detail
corr bloodurea labureaprehd if same_date == 0

* ─ Cas avec écart important (|Δ| > 1 mmol/L) ──────────────────
* Permet de vérifier que les écarts coïncident bien avec dates différentes
display _newline "=== Patients avec |labureaprehd − bloodurea| > 1 mmol/L ==="
list id urinedate datevisit same_date bloodurea labureaprehd delta_urea ///
    if abs(delta_urea) > 1 & !missing(delta_urea), noobs separator(0)

drop same_date delta_urea

* ── 2. Calcul du KRU (mL/min) ───────────────────────────────────

* --- Méthode naïve (référence de comparaison) -------------------
gen kru_naif = (urineurea * urinevolume) / (bloodurea * T_min)
label variable kru_naif "KRU naïf (mL/min)"

* --- Méthode Daugirdas (cinétique non-linéaire) -----------------
* Référence : Daugirdas, Handbook of Dialysis
* Erreur < 5% dans >98% des cas vs modélisation formelle
*
* Variables utilisées :
*   E         : urinevolume [mL], urineurea [mmol/L], T_min [min]
*   URR       : labureaprehd [mmol/L], labureaposthd [mmol/L]
*   IDI       : interdialdays × 1440 [min]
*   R_adj    : URR, T_min/IDI
*   TAC      : bloodurea × R_adj [mmol/L]
*   KRU      : E / TAC [mL/min]
*
* Note temporelle (cf. section 1b) :
*   - bloodurea     : urée au moment de la fin de récolte urinaire
*                     → utilisée pour TAC (ancrage temporel cohérent avec E)
*   - labureaprehd  : urée mesurée au démarrage de la séance d'HD
*                     → utilisée pour URR (réduction réelle PENDANT la séance)
*   Dans 98% des cas les deux valeurs coïncident (mêmes dates) ;
*   ce choix n'affecte que les 3 patients avec urinedate ≠ datevisit.

* (a) Taux d'excrétion urinaire d'urée
gen E_rate = (urineurea * urinevolume) / T_min
label variable E_rate "Taux excrétion urée (mmol·mL/L/min)"

* (b) Urea Reduction Ratio (sur la séance suivant la récolte)
* labureaprehd = urée mesurée au démarrage de la séance d'HD
* → réduction REELLE pendant la dialyse
gen URR = (labureaprehd - labureaposthd) / labureaprehd * 100
label variable URR "Urea Reduction Ratio (%)"

* (c) Intervalle inter-dialytique en minutes
gen IDI = interdialdays * 1440
label variable IDI "Intervalle inter-dialytique (min)"

* (d) Ratio d'ajustement (cinétique non-linéaire)
* 1.075 = correction eau plasmatique (1/0.93)
gen R_adj = 1.075 - (0.0038 * URR + 0.059) * (T_min / IDI)
label variable R_adj "Ratio d'ajustement Daugirdas"

* (e) Concentration plasmatique moyenne pondérée par le temps
gen TAC_urea = bloodurea * R_adj
label variable TAC_urea "TAC urée (mmol/L)"

* (f) KRU final
gen kru_daugirdas = E_rate / TAC_urea
label variable kru_daugirdas "KRU Daugirdas (mL/min)"

* ── 2b. Option A : fallback kru_naif pour les patients avec diurèse
*        mais sans labureaposthd (corrélation naïf/Daugirdas = 0.998) ──
replace kru_daugirdas = kru_naif if missing(kru_daugirdas) & diuresis == 1
label variable kru_daugirdas "KRU Daugirdas (mL/min) [naïf si labureaposthd manquant]"
count if !missing(kru_daugirdas) & diuresis == 1
display "Patients avec KRU Daugirdas récupérés via fallback naïf"

* ── 3. Comparaison naïf vs Daugirdas ────────────────────────────
summarize kru_naif kru_daugirdas, detail
corr kru_naif kru_daugirdas
gen kru_diff = kru_daugirdas - kru_naif
label variable kru_diff "Différence KRU Daugirdas - naïf (mL/min)"
summarize kru_diff, detail

* ── 4. Contrôle qualité ─────────────────────────────────────────
* KRU doit être positif
count if kru_daugirdas < 0 & !missing(kru_daugirdas)
if r(N) > 0 {
    display as error "ATTENTION : " r(N) " valeurs de KRU négatives — à vérifier"
}

* Distribution
summarize kru_daugirdas, detail
histogram kru_daugirdas, normal title("Distribution du KRU Daugirdas (mL/min)") xtitle("KRU (mL/min)")

* Patients sans récolte urinaire (diuresis == 0 ou données manquantes)
count if missing(kru_daugirdas)
display "Patients avec KRU manquant : " r(N)
tab diuresis if missing(kru_daugirdas)

* Explorer les patients avec diurèse mais sans KRU
list id urinevolume urineurea bloodurea labureaposthd urinestart urineend T_min ///
    if diuresis == 1 & missing(kru_daugirdas)

* Assigner KRU = 0 aux patients anuriques (mêmes deux versions)
replace kru_naif      = 0 if diuresis == 0
replace kru_daugirdas = 0 if diuresis == 0
display "KRU=0 assigné aux patients anuriques (naif + Daugirdas)"

* ── Bilan (sur kru_daugirdas) ───────────────────────────────────
count if kru_daugirdas > 0 & !missing(kru_daugirdas)
display "KRU calculé (récolte urinaire) : " r(N)

count if kru_daugirdas == 0
display "KRU = 0 (anuriques)            : " r(N)

count if missing(kru_daugirdas)
display "KRU manquant                   : " r(N)

count if !missing(kru_daugirdas)
display "Total avec KRU disponible      : " r(N) "/" _N

* ── 5. Normalisation par Watson (mL/min/35L) ────────────────────
* Hommes (sex==2) : V = 2.447 - 0.09516×age + 0.1074×height + 0.3362×posthdweight
* Femmes (sex==1) : V = -2.097 + 0.1069×height + 0.2466×posthdweight
* height en cm, posthdweight en kg → V en litres

* Poids : posthdweight en priorité, prehdweight en substitut si manquant
gen weight = posthdweight
replace weight = prehdweight if missing(posthdweight)
label variable weight "Poids (post-HD, ou pré-HD si manquant)"

gen V_watson = .
replace V_watson = 2.447 - 0.09516*age + 0.1074*height + 0.3362*weight if sex == 2
replace V_watson = -2.097 + 0.1069*height + 0.2466*weight              if sex == 1
label variable V_watson "Volume de distribution urée - Watson (L)"

* Contrôle : V doit être physiologiquement plausible (5–70L)
count if (V_watson < 5 | V_watson > 70) & !missing(V_watson)
if r(N) > 0 display as error "ATTENTION : " r(N) " valeurs de V_watson hors plage [5-70L]"

* Si poids totalement absent : assigner V_watson = 35L (moyenne population)
replace V_watson = 35 if missing(V_watson)

* Identifier les patients avec V_watson imputé à 35L
list id age sex height posthdweight prehdweight if V_watson == 35

summarize V_watson, detail

gen kru_naif_35      = kru_naif      * (35 / V_watson)
gen kru_daugirdas_35 = kru_daugirdas * (35 / V_watson)
label variable kru_naif_35      "KRU naïf (mL/min/35L)"
label variable kru_daugirdas_35 "KRU Daugirdas (mL/min/35L)"

summarize kru_naif_35 kru_daugirdas_35, detail

* ── 6. Résumé comparatif naïf vs Daugirdas ──────────────────────
display _newline "========================================"
display         "  RÉSUMÉ COMPARATIF KRU (N=" _N ")"
display         "========================================"

foreach suffix in "" "_35" {
    if "`suffix'" == "" local unit "mL/min"
    else                local unit "mL/min/35L"

    display _newline "--- `unit' ---"
    foreach method in naif daugirdas {
        local var "kru_`method'`suffix'"
        quietly count if !missing(`var')
        local n_tot = r(N)
        quietly count if `var' == 0 & !missing(`var')
        local n_an = r(N)
        local pct_an = string(round(`n_an' / `n_tot' * 100, 0.1))
        quietly summarize `var' if `var' > 0 & !missing(`var'), detail
        local med = string(round(r(p50), 0.01))
        local moy = string(round(r(mean), 0.01))
        display "  `method' (`unit') : N=`n_tot'  anuriques=`n_an' (`pct_an'%)  médiane=`med'  moyenne=`moy'"
    }
}

display "========================================"

* Les données restent en mémoire pour la suite de l'analyse
