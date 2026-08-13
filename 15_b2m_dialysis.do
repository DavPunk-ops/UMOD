* ===========================================================================
* 15_b2m_dialysis.do
*
* Prérequis : lancer 02_calculate_kru.do AVANT (données en mémoire).
*             Autonome ; ne modifie AUCUNE do-file existante.
*
* Objectif : réponse à R3#7 — « β2M pré-dialyse est influencée par la
*   prescription de dialyse ; le modèle combiné pourrait capter la dialyse
*   plutôt que la biologie rénale. »
*
* Deux volets :
*   1. Corrélations UMOD / β2M avec les paramètres de dialyse
*      (UMOD, marqueur tubulaire non dialysé, devrait être peu lié ;
*       β2M, moyenne molécule, davantage lié — comme le note le reviewer)
*   2. Le signal survit-il à l'AJUSTEMENT sur la dialyse ?
*      logit KRU≥2 ~ UMOD + β2M   vs   + vintage + mode + UF + spKt/V
*      → UMOD et β2M restent-ils des prédicteurs indépendants ?
*      → l'ajout des paramètres de dialyse améliore-t-il l'AUC ?
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

capture confirm variable kru_ge2
if _rc  gen byte kru_ge2 = (kru_daugirdas_35 >= 2) if !missing(kru_daugirdas_35)

* Vérifier la présence des variables de dialyse
display _newline "  --- Disponibilité des paramètres de dialyse ---"
foreach v in vintage mode uf spktv sessiontime {
    capture confirm variable `v'
    if _rc  display as error "  ATTENTION : '`v'' absente — vérifier le nom"
    else {
        quietly count if !missing(`v')
        display "  `v' : N disponible = " r(N)
    }
}

* ###########################################################################
* SECTION 1 — CORRÉLATIONS BIOMARQUEURS × PARAMÈTRES DE DIALYSE
* ###########################################################################
display _newline(2) "########################################################"
display              "  1. Corrélations UMOD / β2M × paramètres de dialyse"
display              "########################################################"

display _newline "  --- UMOD vs (vintage, UF, spKt/V, durée de séance) ---"
spearman umod vintage uf spktv sessiontime, stats(rho p) star(0.05)

display _newline "  --- β2M vs (vintage, UF, spKt/V, durée de séance) ---"
spearman labb2mprehd vintage uf spktv sessiontime, stats(rho p) star(0.05)

capture confirm variable mode
if _rc==0 {
    display _newline "  --- UMOD selon le mode (HD vs HDF) ---"
    ranksum umod, by(mode)
    display _newline "  --- β2M selon le mode (HD vs HDF) ---"
    ranksum labb2mprehd, by(mode)
}

* ###########################################################################
* SECTION 2 — SURVIE DU SIGNAL APRÈS AJUSTEMENT SUR LA DIALYSE
* ###########################################################################
display _newline(2) "########################################################"
display              "  2. logit KRU≥2 ~ UMOD + β2M  (± paramètres de dialyse)"
display              "########################################################"

* --- Modèle de base (publié) ---
display _newline "  === Modèle de base : UMOD + β2M ==="
logit kru_ge2 umod labb2mprehd, or
quietly lroc, nograph
local auc_base = r(area)
display "  AUC (base) = " %5.3f `auc_base'

* --- Modèle ajusté sur la dialyse ---
display _newline "  === Modèle ajusté : + vintage + mode + UF + spKt/V + durée séance ==="
logit kru_ge2 umod labb2mprehd vintage i.mode uf spktv sessiontime, or
quietly lroc, nograph
local auc_adj = r(area)
display "  AUC (ajusté dialyse) = " %5.3f `auc_adj'
quietly count if !missing(kru_ge2, umod, labb2mprehd, vintage, mode, uf, spktv, sessiontime)
display "  N (modèle ajusté)    = " r(N)

display _newline "  --- Comparaison ---"
display "  AUC base            = " %5.3f `auc_base'
display "  AUC + dialyse       = " %5.3f `auc_adj'
display "  Δ AUC               = " %5.3f (`auc_adj' - `auc_base')
display _newline "  → Si UMOD & β2M restent significatifs et ΔAUC faible :"
display "    le signal des biomarqueurs n'est PAS expliqué par la dialyse."

display _newline(2) "=== FIN do-file 15 ==="
