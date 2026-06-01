* ===========================================================================
* 07_temporal_explore.do
* Prérequis : lancer 02_calculate_kru.do AVANT ce do-file.
*
* Objectif : EXPLORATOIRE uniquement — évaluer la faisabilité d'une
*            validation interne par temporal split (séparation chronologique
*            train/test) à partir de la date de visite (datevisit).
*
* → manuscrit : aide à la décision pour la section
*               "Diagnostic performance of combined serum uromodulin and
*                β2-microglobulin for predicting KRU ≥2 mL/min/35L"
*
* CE DO-FILE NE SAUVEGARDE RIEN ET NE MODIFIE AUCUNE DONNÉE.
* Il sert seulement à juger si N et le nombre d'events par période
* suffisent pour un split temporel exploitable.
*
* Section 1 : étendue et distribution temporelle des visites.
* Section 2 : prévalence de KRU≥2 dans le temps (année / semestre).
* Section 3 : faisabilité de différents points de coupure (médiane, 70/30).
* ===========================================================================

* --- Vérification des prérequis ---
foreach v in kru_daugirdas_35 datevisit umod labb2mprehd {
    capture confirm variable `v'
    if _rc {
        display as error "ERREUR : variable '`v'' absente — lance d'abord 02_calculate_kru.do"
        exit 111
    }
}
display _newline "=== Prérequis OK ==="

* --- Variable binaire KRU≥2 (recréée si absente) ---
capture confirm variable kru_ge2
if _rc {
    gen byte kru_ge2 = (kru_daugirdas_35 >= 2) if !missing(kru_daugirdas_35)
    label variable kru_ge2 "KRU >= 2 mL/min/35L"
}

* --- Restriction à la population analysée (KRU calculable) ---
keep if !missing(kru_ge2)
display _newline "  Population analysée (KRU calculable) : N = " _N

* --- Normalisation de datevisit en date numérique Stata (%td) ---
* Gère le cas où datevisit est déjà numérique (%td) ou une chaîne.
capture confirm numeric variable datevisit
if _rc {
    * datevisit est une chaîne → tentative de conversion (formats usuels)
    display as text "  datevisit est une chaîne — conversion en date %td"
    gen double _vdate = .
    capture replace _vdate = date(datevisit, "DMY")
    capture replace _vdate = date(datevisit, "MDY") if missing(_vdate)
    capture replace _vdate = date(datevisit, "YMD") if missing(_vdate)
    format _vdate %td
}
else {
    gen double _vdate = datevisit
    format _vdate %td
}

count if missing(_vdate)
if r(N) > 0 {
    display as error "  ATTENTION : " r(N) " date(s) de visite manquante(s)/non converties"
}

* Année et semestre de visite
gen int  _vyear = year(_vdate)
gen int  _vsem  = halfyear(_vdate)
label variable _vyear "Année de visite"
label variable _vsem  "Semestre (1 = jan-juin, 2 = juil-déc)"

* ===========================================================================
* SECTION 1 — ÉTENDUE ET DISTRIBUTION TEMPORELLE
* ===========================================================================
display _newline(2) "========================================================================"
display              "  SECTION 1 — Étendue et distribution des visites"
display              "========================================================================"

display _newline "  --- Bornes et médiane de la date de visite ---"
summarize _vdate, detail format

display _newline "  Étendue : du " %td r(min) "  au  " %td r(max)
summarize _vdate
local span_days = r(max) - r(min)
display "  Durée totale de recrutement : " %6.0f `span_days' " jours" ///
        "  (~" %4.1f `span_days'/365.25 " ans)"

display _newline "  --- Effectifs par année de visite ---"
tab _vyear

display _newline "  --- Effectifs par année × semestre ---"
tab _vyear _vsem

* ===========================================================================
* SECTION 2 — PRÉVALENCE DE KRU≥2 DANS LE TEMPS
* ===========================================================================
display _newline(2) "========================================================================"
display              "  SECTION 2 — Prévalence de KRU≥2 au fil du temps"
display              "========================================================================"

display _newline "  --- KRU≥2 par année (lignes = année, % en ligne) ---"
tab _vyear kru_ge2, row

display _newline "  --- Nombre d'events (KRU≥2) par année ---"
tabstat kru_ge2, by(_vyear) statistics(n sum mean) format(%6.3f)

display _newline "  Note : 'sum' = nb d'events KRU≥2 ; 'mean' = prévalence dans l'année."
display "  Une dérive marquée de la prévalence rend un split temporel déséquilibré."

* ===========================================================================
* SECTION 3 — FAISABILITÉ DES POINTS DE COUPURE
* ===========================================================================
display _newline(2) "========================================================================"
display              "  SECTION 3 — Faisabilité de différents splits temporels"
display              "========================================================================"

* Effectif total et nombre total d'events disponibles
count
local n_tot = r(N)
count if kru_ge2 == 1
local ev_tot = r(N)
display _newline "  Total : N = `n_tot'  |  events KRU≥2 = `ev_tot'  |  non-events = " ///
        `n_tot' - `ev_tot'

* ---- Coupure A : médiane de la date (split ~50/50) ----
display _newline "  --------------------------------------------------------------------"
display          "  COUPURE A — médiane de la date de visite (≈ 50/50)"
display          "  --------------------------------------------------------------------"
summarize _vdate, detail
local cutA = r(p50)
gen byte _splitA = (_vdate > `cutA') if !missing(_vdate)
label define splitA 0 "Train (≤ médiane)" 1 "Test (> médiane)", replace
label values _splitA splitA
display "  Date de coupure (médiane) : " %td `cutA'
tab _splitA kru_ge2, row

display _newline "  Events KRU≥2 de part et d'autre :"
tabstat kru_ge2, by(_splitA) statistics(n sum mean) format(%6.3f)

* ---- Coupure B : 70e centile de la date (split 70/30) ----
display _newline "  --------------------------------------------------------------------"
display          "  COUPURE B — 70e centile de la date (≈ 70 train / 30 test)"
display          "  --------------------------------------------------------------------"
_pctile _vdate, percentiles(70)
local cutB = r(r1)
gen byte _splitB = (_vdate > `cutB') if !missing(_vdate)
label define splitB 0 "Train (1ers 70%)" 1 "Test (derniers 30%)", replace
label values _splitB splitB
display "  Date de coupure (P70) : " %td `cutB'
tab _splitB kru_ge2, row

display _newline "  Events KRU≥2 de part et d'autre :"
tabstat kru_ge2, by(_splitB) statistics(n sum mean) format(%6.3f)

* ===========================================================================
* SYNTHÈSE — règle de décision
* ===========================================================================
display _newline(2) "========================================================================"
display              "  SYNTHÈSE — faisabilité du temporal split"
display              "========================================================================"
display "  Repère pratique pour estimer des seuils stables (Se≥90% / Sp≥90%) :"
display "    - viser ≥ ~20-25 events KRU≥2 ET ≥ ~20-25 non-events dans le TRAIN ;"
display "    - un test set < ~30 patients donne des Se/Sp à IC très larges."
display ""
display "  Si le plus petit des deux groupes (train/test) tombe sous ces repères,"
display "  le split temporel sera sous-puissant → préférer le mentionner comme"
display "  limitation (absence de validation externe) plutôt que de le réaliser."
display "========================================================================"

* --- Nettoyage des variables temporaires de travail ---
* (ce do-file n'altère pas le dataset sauvegardé ; on lance sur données en mémoire)
capture drop _vdate _vyear _vsem _splitA _splitB
display _newline "  [07] Variables exploratoires temporaires supprimées de la mémoire."
