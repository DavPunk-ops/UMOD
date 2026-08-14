* ===========================================================================
* 01_merge_baseline_umod.do
* Objectif : Extraire le baseline, merger avec les valeurs UMOD,
*            puis ajouter les données médicaments (medication_arm_1).
*
* ── Traçabilité UMOD (source labo : dosages_UMOD.xlsx) ──────────────────────
* Correspondance Excel ↔ .dta VÉRIFIÉE ligne-à-ligne : 151/151 patients
* analysés identiques (0 mismatch). La colonne UMOD de l'Excel = moyenne des
* DUPLICATS rendue par le labo (chaque dosage UMOD fait en double ; β2M non).
*   - "undetectable" / "Range?" (10 samples, tous anuriques) → codés umod = 0.
*   - "CV>20%" (flag qualité, 19 samples) : métadonnée Excel uniquement, absente
*     de la .dta ; ces samples restent dans l'analyse avec leur valeur moyenne.
*   - 2 samples Excel absents de la .dta (126-T0, 143-T0) = exclus pour collecte
*     urinaire incomplète (UMOD valide, sans lien avec l'assay).
* ===========================================================================

local path "C:\Users\dajs\OneDrive - HOPITAUX UNIVERSITAIRES DE GENEVE\recherche\RKF\UMOD\stata\main prospective study\with Claude"

local main_db  "`path'\database.dta"
local umod_db  "`path'\umod.dta"
local output   "`path'\database_umod.dta"

* ── 1. Extraire les données médicaments ─────────────────────────
use "`main_db'", clear
keep if redcap_event_name == "medication_arm_1"
* Attendu : ~152 observations

display _newline "  Event medication_arm_1 : " _N " observations"

* Vérifier que les variables médicaments sont bien présentes
foreach v in diuretic antiht ado insulin lipid epo pobinder kbinder vitdanalog bicarbonate {
    capture confirm variable `v'
    if _rc {
        display as error "  ATTENTION : variable '`v'' absente dans medication_arm_1"
    }
}

* Garder uniquement record_id + variables médicaments
keep record_id diuretic antiht ado insulin lipid epo pobinder kbinder vitdanalog bicarbonate

* Nettoyer record_id pour le merge (même format que baseline)
gen record_id_num = real(record_id)
drop record_id

* Sauvegarder temporairement
tempfile medic_data
save `medic_data'

display "  Données médicaments extraites : " _N " patients"

* ── 2. Charger le baseline et merger avec UMOD ──────────────────
use "`main_db'", clear
keep if redcap_event_name == "baseline_t0_arm_1"
* Attendu : ~154 observations
display _newline "  Event baseline_t0_arm_1 : " _N " observations"

* Construire la clé de merge avec UMOD
gen id = string(real(record_id), "%03.0f") + "-T0"

* Merger avec les valeurs UMOD
merge 1:1 id using "`umod_db'"
tab _merge
* 3 = match (attendu pour la majorité)
keep if _merge == 3
drop _merge

display _newline "  Après merge UMOD : " _N " observations"

* ── 3. Merger avec les données médicaments ──────────────────────
gen record_id_num = real(record_id)
* Option update : remplit les missing du master avec les valeurs du using
* (nécessaire car diuretic etc. existent dans baseline mais sont toutes missing)
merge m:1 record_id_num using `medic_data', ///
    keepusing(diuretic antiht ado insulin lipid epo pobinder kbinder vitdanalog bicarbonate) ///
    gen(_merge_medic) update

tab _merge_medic
* 1 = baseline sans médicaments (rare)
* 2 = médicaments sans baseline → à exclure
* 3 = match complet (attendu)

* Exclure le patient présent dans medication mais absent du baseline UMOD
keep if _merge_medic != 2

display _newline "  Après merge médicaments : " _N " observations"

drop _merge_medic record_id_num

* ── 4. Vérification rapide ──────────────────────────────────────
display _newline "  --- Médicaments (N total = " _N ") ---"
foreach v in diuretic antiht ado insulin lipid epo pobinder kbinder vitdanalog bicarbonate {
    capture confirm variable `v'
    if !_rc {
        quietly count if !missing(`v')
        display "  `v' : " r(N) " non-manquants"
    }
}

* ── 5. Sauvegarder ──────────────────────────────────────────────
save "`output'", replace
display _newline "=== Base sauvegardée : `output' ==="
display "  N final = " _N " patients"
