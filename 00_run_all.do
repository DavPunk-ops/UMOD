* ===========================================================================
* 00_run_all.do
* Objectif : Lance la totalité de la chaîne d'analyse en une seule commande.
*
* Chaîne :
*   02_calculate_kru.do          → charge la base et calcule KRU Daugirdas/35L
*   03_estimate_kru.do           → modélisation UMOD → KRU (two-part, KRU≥2)
*   04_bootstrap_validation.do   → validation bootstrap interne des modèles
*
* Note : 01_merge_baseline_umod.do n'est PAS inclus (à lancer une seule fois
*        pour créer database_umod.dta ; les do-files 02-04 partent de ce
*        fichier déjà mergé).
*
* Usage :  do 00_run_all.do
* ===========================================================================

* Répertoire courant du do-file (à adapter si nécessaire)
local dofile_dir "C:\Users\dajs\OneDrive - HOPITAUX UNIVERSITAIRES DE GENEVE\recherche\RKF\UMOD\stata\main prospective study\with Claude"

display _newline(2) "###########################################################"
display              "#  CHAÎNE D'ANALYSE UMOD → KRU                            #"
display              "###########################################################"

display _newline "=== ÉTAPE 1/3 : 02_calculate_kru.do ==="
do "`dofile_dir'\02_calculate_kru.do"

display _newline(2) "=== ÉTAPE 2/3 : 03_estimate_kru.do ==="
do "`dofile_dir'\03_estimate_kru.do"

display _newline(2) "=== ÉTAPE 3/3 : 04_bootstrap_validation.do ==="
do "`dofile_dir'\04_bootstrap_validation.do"

display _newline(2) "###########################################################"
display              "#  CHAÎNE TERMINÉE                                        #"
display              "###########################################################"
