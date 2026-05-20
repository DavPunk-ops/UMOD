* ===========================================================================
* 01_merge_baseline_umod.do
* Objectif : Extraire le baseline et merger avec les valeurs UMOD
* ===========================================================================

local path "C:\Users\dajs\OneDrive - HOPITAUX UNIVERSITAIRES DE GENEVE\recherche\RKF\UMOD\stata\main prospective study\with Claude"

local main_db  "`path'\database.dta"
local umod_db  "`path'\umod.dta"
local output   "`path'\database_baseline.dta"

* ── 1. Charger la base principale et garder le baseline ─────────
use "`main_db'", clear

keep if redcap_event_name == "baseline_t0_arm_1"
* Attendu : ~154 observations

* ── 2. Construire la clé de merge ───────────────────────────────
* Format : "029-T0" (numéro patient sur 3 chiffres + "-T0")
gen id = string(real(record_id), "%03.0f") + "-T0"

* ── 3. Merger avec les valeurs UMOD ─────────────────────────────
merge 1:1 id using "`umod_db'"

tab _merge
* 1 = base principale seulement (pas de valeur UMOD)
* 2 = fichier UMOD seulement (patient absent du baseline)
* 3 = match des deux côtés (attendu pour la majorité)

keep if _merge == 3
drop _merge

* ── 4. Sauvegarder ──────────────────────────────────────────────
save "`output'", replace
