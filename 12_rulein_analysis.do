* ===========================================================================
* 12_rulein_analysis.do
*
* Prérequis : lancer 02_calculate_kru.do AVANT (données en mémoire).
*             Autonome ; ne modifie AUCUNE do-file existante.
*
* QUESTION CLINIQUE (objectif principal du révisé) :
*   « Puis-je être CERTAIN qu'un patient a KRU ≥2 mL/min/35L, pour lui offrir
*     une HD incrémentale en sécurité ? »  → c'est une question de RULE-IN
*     (haute spécificité / haut PPV). Et : UMOD AJOUTE-T-IL à β2M pour cela ?
*
*   Analyse restreinte aux NON-ANURIQUES (kru_pos==1), population à KRU mesuré.
*
* Méthode :
*   Pour chaque prédicteur (UMOD seul, β2M seul, modèle combiné), on RE-DÉRIVE
*   le seuil rule-in DANS les non-anuriques comme le centile du score chez les
*   KRU<2 donnant la spécificité cible (90% puis 95%, « certitude »).
*   On rapporte : seuil, patients rule-in (rendement), vrais/faux, PPV, Se.
*   Comparaison head-to-head → β2M ajoute-t-il un rendement au-delà de UMOD ?
*
*   Orientation : on utilise des scores où HAUT = KRU≥2 probable :
*     umod (+), neg_b2m = −β2M (+), p_ge2_na = P(KRU≥2) du modèle combiné (+).
* ===========================================================================

* --- Prérequis ---
foreach v in kru_daugirdas_35 umod labb2mprehd {
    capture confirm variable `v'
    if _rc {
        display as error "ERREUR : '`v'' absente — lance d'abord 02_calculate_kru.do"
        exit 111
    }
}
display _newline "=== Prérequis OK ==="

* --- Variables dérivées ---
capture confirm variable kru_pos
if _rc {
    gen byte kru_pos = (kru_daugirdas_35 > 0) if !missing(kru_daugirdas_35)
}
capture confirm variable kru_ge2
if _rc {
    gen byte kru_ge2 = (kru_daugirdas_35 >= 2) if !missing(kru_daugirdas_35)
}
capture drop neg_b2m
gen double neg_b2m = -labb2mprehd if !missing(labb2mprehd)
label variable neg_b2m "−β2M (orienté : haut = KRU≥2)"

* Modèle combiné publié (logit cohorte complète), évalué chez non-anuriques
capture drop p_ge2_na
quietly logit kru_ge2 umod labb2mprehd
predict p_ge2_na, pr
label variable p_ge2_na "P(KRU≥2) — modèle combiné publié"

* ===========================================================================
* PROGRAMME UTILITAIRE : rule-in d'un score à une spécificité cible
*   Seuil = centile (100 - specTarget) du score chez les NON-ANURIQUES KRU<2
*   (p.ex. spéc 90% → 90e centile des KRU<2 → 10% de faux positifs)
* ===========================================================================
capture program drop rulein_report
program define rulein_report
    args score spec label
    * spec = spécificité cible en % (ex. 90) ; label = nom lisible du score
    local pctile = `spec'                 // centile du score chez KRU<2
    quietly centile `score' if kru_pos==1 & kru_ge2==0 & !missing(`score'), ///
        centile(`pctile')
    local thr = r(c_1)

    quietly count if kru_pos==1 & kru_ge2==1 & !missing(`score')
    local n_pos = r(N)                    // total KRU≥2 (non-anuriques, score dispo)
    quietly count if kru_pos==1 & kru_ge2==0 & !missing(`score')
    local n_neg = r(N)                    // total KRU<2

    quietly count if kru_pos==1 & `score'>=`thr' & !missing(`score')
    local nri = r(N)                      // rule-in total
    quietly count if kru_pos==1 & `score'>=`thr' & kru_ge2==1 & !missing(`score')
    local ntrue = r(N)                    // vrais KRU≥2
    quietly count if kru_pos==1 & `score'>=`thr' & kru_ge2==0 & !missing(`score')
    local nfalse = r(N)                   // faux rule-in (KRU<2)

    local ppv = 100*`ntrue'/`nri'
    local sens = 100*`ntrue'/`n_pos'
    local specobs = 100*(`n_neg'-`nfalse')/`n_neg'

    display "  " %-22s "`label'" "  seuil=" %7.3f `thr' ///
        "  rule-in n=" %2.0f `nri' ///
        "  vrais=" %2.0f `ntrue' "  faux=" %2.0f `nfalse' ///
        "  PPV=" %4.1f `ppv' "%" ///
        "  Se=" %4.1f `sens' "%  Sp(obs)=" %4.1f `specobs' "%"
end

* ===========================================================================
* SECTION A — RULE-IN À SPÉCIFICITÉ ≥ 90%  (règle standard)
* ===========================================================================
display _newline(2) "########################################################"
display              "  A. RULE-IN à spécificité 90% (non-anuriques)"
display              "     Seuils RE-DÉRIVÉS dans les non-anuriques"
display              "########################################################"
display _newline "  Total non-anuriques : KRU≥2 et KRU<2"
tab kru_ge2 if kru_pos==1 & !missing(p_ge2_na)

display _newline "  --- Rendement rule-in à Sp≈90% ---"
rulein_report umod     90 "UMOD seul"
rulein_report neg_b2m  90 "β2M seul"
rulein_report p_ge2_na 90 "Combiné UMOD+β2M"

* ===========================================================================
* SECTION B — RULE-IN À SPÉCIFICITÉ ≥ 95%  (« certitude » pour prescrire)
* ===========================================================================
display _newline(2) "########################################################"
display              "  B. RULE-IN à spécificité 95% (non-anuriques)"
display              "     seuil conservateur = quasi-certitude KRU≥2"
display              "########################################################"

display _newline "  --- Rendement rule-in à Sp≈95% ---"
rulein_report umod     95 "UMOD seul"
rulein_report neg_b2m  95 "β2M seul"
rulein_report p_ge2_na 95 "Combiné UMOD+β2M"

* ===========================================================================
* SECTION C — LECTURE
*   Question : à spécificité fixée, le COMBINÉ rule-in PLUS de patients
*   (rendement = n rule-in / Se) que β2M seul ? Si oui → UMOD ajoute une
*   valeur clinique pour la certitude KRU≥2. Sinon → β2M seul suffit.
* ===========================================================================
display _newline(2) "########################################################"
display              "  C. INTERPRÉTATION"
display              "########################################################"
display "  À spécificité fixée, comparer le nombre de VRAIS rule-in :"
display "    - si Combiné > β2M seul  → UMOD ajoute du rendement (certitude+)"
display "    - si Combiné ≈ β2M seul  → β2M seul suffit pour le rule-in"
display "  Le PPV à Sp=95% indique le niveau de certitude atteignable."

display _newline(2) "=== FIN do-file 12 ==="
