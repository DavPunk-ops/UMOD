* ===========================================================================
* 03_normalize_kru.do
* Objectif : Normalisation du KRU pour 35L via formule de Watson
* Formule  : KRU_norm (mL/min/35L) = KRU × (35 / V_watson)
* ===========================================================================

* ── 1. Volume de distribution individuel (Watson) ───────────────
* Hommes (sex==2) : V = 2.447 - 0.09516×age + 0.1074×height + 0.3362×posthdweight
* Femmes (sex==1) : V = -2.097 + 0.1069×height + 0.2466×posthdweight
* height en cm, posthdweight en kg → V en litres

gen V_watson = .
replace V_watson = 2.447 - 0.09516*age + 0.1074*height + 0.3362*posthdweight if sex == 2
replace V_watson = -2.097 + 0.1069*height + 0.2466*posthdweight              if sex == 1

label variable V_watson "Volume de distribution urée - Watson (L)"

* Contrôle : V doit être physiologiquement plausible (5–70L)
count if V_watson < 5 | V_watson > 70
if r(N) > 0 display as error "ATTENTION : " r(N) " valeurs de V_watson hors plage [5-70L]"

summarize V_watson, detail

* ── 2. KRU normalisé pour 35L ───────────────────────────────────
gen kru_35 = kru * (35 / V_watson)

label variable kru_35 "KRU (mL/min/35L)"

* ── 3. Distribution ─────────────────────────────────────────────
summarize kru_35, detail

* Les données restent en mémoire pour la suite de l'analyse
