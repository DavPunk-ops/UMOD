* ===========================================================================
* 04_estimate_kru.do
* Prérequis : lancer 02_calculate_kru.do AVANT ce do-file dans la même
*             session Stata.
* ===========================================================================

* --- Vérification que 02_calculate_kru.do a été exécuté ---
foreach v in kru_daugirdas_35 kru_naif_35 umod diuresis labcreatprehd {
    capture confirm variable `v'
    if _rc {
        display as error "ERREUR : variable '`v'' absente — lance d'abord 02_calculate_kru.do"
        exit 111
    }
}

display _newline "=== Prérequis OK — variables de 02_calculate_kru.do présentes ==="

* --- Variables binaires utilisées dans tout le do-file ---
capture drop kru_pos
gen byte kru_pos = (kru_daugirdas_35 > 0) if !missing(kru_daugirdas_35)
label variable kru_pos "KRU > 0 (1 = non-anurique)"
label define krupos 0 "Anurique (KRU=0)" 1 "Non-anurique (KRU>0)", replace
label values kru_pos krupos

capture drop kru_ge2
gen byte kru_ge2 = (kru_daugirdas_35 >= 2) if !missing(kru_daugirdas_35)
label variable kru_ge2 "KRU >= 2 mL/min/35L"
label define kruge2 0 "KRU < 2" 1 "KRU >= 2", replace
label values kru_ge2 kruge2

* ###########################################################################
* SECTION 1 — DESCRIPTION DE L'UROMODULINE SÉRIQUE (UMOD)
* ###########################################################################

* ===========================================================================
*  1a. DESCRIPTION GLOBALE DE L'UMOD
* ===========================================================================
display _newline(2) "=============================================="
display              "  1a. Distribution UMOD — population entière"
display              "=============================================="

summarize umod, detail

swilk umod

count if umod == 0
display _newline "  UMOD = 0 : " r(N) " patients"

display _newline "  Répartition des UMOD=0 selon kru_pos :"
tab kru_pos if umod == 0, miss

histogram umod, normal ///
    title("Distribution de l'uromoduline sérique") ///
    xtitle("UMOD (ng/mL)") name(hist_umod, replace)

* ===========================================================================
*  1b. UMOD CHEZ ANURIQUES vs NON-ANURIQUES
* ===========================================================================
display _newline(2) "=============================================="
display              "  1b. UMOD : anuriques vs non-anuriques"
display              "=============================================="

tabstat umod, by(kru_pos) statistics(n mean sd p25 p50 p75) format(%7.2f)

ranksum umod, by(kru_pos)

* ===========================================================================
*  1c. UMOD CHEZ KRU<2 vs KRU>=2 — POPULATION ENTIÈRE
* ===========================================================================
display _newline(2) "=============================================="
display              "  1c. UMOD : KRU<2 vs KRU>=2 (population entière)"
display              "=============================================="

tab kru_ge2, miss

tabstat umod, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%7.2f)

ranksum umod, by(kru_ge2)

* ===========================================================================
*  1d. UMOD CHEZ KRU<2 vs KRU>=2 — NON-ANURIQUES UNIQUEMENT
* ===========================================================================
display _newline(2) "=============================================="
display              "  1d. UMOD : KRU<2 vs KRU>=2 (non-anuriques)"
display              "=============================================="

tabstat umod if kru_pos == 1, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%7.2f)

ranksum umod if kru_pos == 1, by(kru_ge2)

* ===========================================================================
*  1e. FIGURE 1 — Strip plot UMOD : anuriques vs non-anuriques (A)
*                                    KRU<2 vs KRU≥2 (B)
*      Points individuels (jitter) + barre médiane + IQR — échelle native
* ===========================================================================

preserve
keep if !missing(kru_ge2)

set seed 20260522

gen double _xA  = 1*(kru_pos==0) + 2*(kru_pos==1)
gen double _xjA = _xA + (runiform()-0.5)*0.35

gen double _xB  = 1*(kru_ge2==0) + 2*(kru_ge2==1)
gen double _xjB = _xB + (runiform()-0.5)*0.35

foreach g in 0 1 {
    quietly summarize umod if kru_pos == `g', detail
    local med_A`g' = r(p50)
    local p25_A`g' = r(p25)
    local p75_A`g' = r(p75)

    quietly summarize umod if kru_ge2 == `g', detail
    local med_B`g' = r(p50)
    local p25_B`g' = r(p25)
    local p75_B`g' = r(p75)
}

local mw = 0.20

twoway ///
    (scatter umod _xjA if kru_pos==0, ///
        mcolor(navy%35) msize(small) msymbol(circle)) ///
    (scatter umod _xjA if kru_pos==1, ///
        mcolor(cranberry%35) msize(small) msymbol(circle)) ///
    (pci `p25_A0' 1 `p75_A0' 1, lcolor(navy)     lwidth(medthick)) ///
    (pci `p25_A1' 2 `p75_A1' 2, lcolor(cranberry) lwidth(medthick)) ///
    (pci `med_A0' `=1-`mw'' `med_A0' `=1+`mw'', lcolor(navy)     lwidth(vthick)) ///
    (pci `med_A1' `=2-`mw'' `med_A1' `=2+`mw'', lcolor(cranberry) lwidth(vthick)) ///
    , ///
    xlabel(1 "Anuric" 2 "Non-anuric", noticks labsize(medlarge)) ///
    xtitle("") ///
    ytitle("Serum UMOD (ng/mL)", size(medlarge)) ///
    xscale(range(0.4 2.6)) ///
    ylabel(0(10)50, grid glcolor(gs14) labsize(medlarge)) ///
    text(48 1.5 "p<0.001", size(medlarge) color(black)) ///
    legend(off) ///
    graphregion(color(white)) plotregion(color(white)) ///
    title("A", pos(11) size(large)) ///
    name(figA, replace)

twoway ///
    (scatter umod _xjB if kru_ge2==0, ///
        mcolor(navy%35) msize(small) msymbol(circle)) ///
    (scatter umod _xjB if kru_ge2==1, ///
        mcolor(cranberry%35) msize(small) msymbol(circle)) ///
    (pci `p25_B0' 1 `p75_B0' 1, lcolor(navy)     lwidth(medthick)) ///
    (pci `p25_B1' 2 `p75_B1' 2, lcolor(cranberry) lwidth(medthick)) ///
    (pci `med_B0' `=1-`mw'' `med_B0' `=1+`mw'', lcolor(navy)     lwidth(vthick)) ///
    (pci `med_B1' `=2-`mw'' `med_B1' `=2+`mw'', lcolor(cranberry) lwidth(vthick)) ///
    , ///
    xlabel(1 "KRU <2" 2 "KRU ≥2", noticks labsize(medlarge)) ///
    xtitle("mL/min/35L", size(medsmall)) ///
    ytitle("Serum UMOD (ng/mL)", size(medlarge)) ///
    xscale(range(0.4 2.6)) ///
    ylabel(0(10)50, grid glcolor(gs14) labsize(medlarge)) ///
    text(48 1.5 "p<0.001", size(medlarge) color(black)) ///
    legend(off) ///
    graphregion(color(white)) plotregion(color(white)) ///
    title("B", pos(11) size(large)) ///
    name(figB, replace)

graph combine figA figB, ///
    cols(2) imargin(small) ///
    graphregion(color(white)) ///
    xsize(8) ysize(4.5) ///
    name(fig1_umod, replace)

graph export "Figure1_UMOD.tif", replace width(2400)

restore

* ###########################################################################
* SECTION 2 — DESCRIPTION DE LA CRÉATININE PRÉDIALYSE (comparaison)
* ###########################################################################

* ===========================================================================
*  2a. DESCRIPTION GLOBALE DE LA CRÉATININE
* ===========================================================================
display _newline(2) "=============================================="
display              "  2a. Distribution créatinine — population entière"
display              "=============================================="

summarize labcreatprehd, detail

swilk labcreatprehd

count if missing(labcreatprehd)
display _newline "  Créatinine manquante : " r(N) " patients"

histogram labcreatprehd, normal ///
    title("Distribution de la créatinine prédialyse") ///
    xtitle("Créatinine (umol/L)") name(hist_creat, replace)

* ===========================================================================
*  2b. CRÉATININE CHEZ ANURIQUES vs NON-ANURIQUES
* ===========================================================================
display _newline(2) "=============================================="
display              "  2b. Créatinine : anuriques vs non-anuriques"
display              "=============================================="

tabstat labcreatprehd, by(kru_pos) statistics(n mean sd p25 p50 p75) format(%7.1f)

ranksum labcreatprehd, by(kru_pos)

* ===========================================================================
*  2c. CRÉATININE CHEZ KRU<2 vs KRU>=2 — POPULATION ENTIÈRE
* ===========================================================================
display _newline(2) "=============================================="
display              "  2c. Créatinine : KRU<2 vs KRU>=2 (population entière)"
display              "=============================================="

tabstat labcreatprehd, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%7.1f)

ranksum labcreatprehd, by(kru_ge2)

* ===========================================================================
*  2d. CRÉATININE CHEZ KRU<2 vs KRU>=2 — NON-ANURIQUES UNIQUEMENT
* ===========================================================================
display _newline(2) "=============================================="
display              "  2d. Créatinine : KRU<2 vs KRU>=2 (non-anuriques)"
display              "=============================================="

tabstat labcreatprehd if kru_pos == 1, by(kru_ge2) statistics(n mean sd p25 p50 p75) format(%7.1f)

ranksum labcreatprehd if kru_pos == 1, by(kru_ge2)

* ###########################################################################
* SECTION 3 — TWO-PART MODEL : KRU ~ UMOD
*              Partie 1 : logit P(KRU>0)        — sur N=151
*              Partie 2 : OLS  E[KRU | KRU>0]   — sur N=89 (non-anuriques)
*              Combiné  : E[KRU | UMOD] = P(KRU>0|UMOD) × E[KRU|KRU>0,UMOD]
* ###########################################################################

* ===========================================================================
*  3a. PARTIE 1 — Logit P(KRU>0) ~ UMOD
* ===========================================================================
display _newline(2) "=============================================="
display              "  3a. Logit  P(KRU>0)  ~  UMOD     (N=151)"
display              "=============================================="

logit kru_pos umod
estimates store tp_logit

display _newline "  --- Odds ratios ---"
logit kru_pos umod, or

capture drop p_pos
predict p_pos, pr
label variable p_pos "P(KRU>0 | UMOD)"

display _newline "  --- Performance (ROC / calibration) ---"
lroc, nograph
display "    AUC = " %5.3f r(area)

estat gof, group(10) table

* ===========================================================================
*  3b. PARTIE 2 — E[KRU | KRU>0] ~ UMOD  (non-anuriques, N=89)
*       OLS (linéaire) et MFP (fractional polynomial) en parallèle
* ===========================================================================
display _newline(2) "=============================================="
display              "  3b. E[KRU | KRU>0]  ~  UMOD   (N=89)"
display              "      OLS linéaire vs MFP en parallèle"
display              "=============================================="

* --- 3b.i  OLS linéaire (modèle principal pour la prédiction two-part) ---
*           SE robustes (Huber-White) — coefficients inchangés, inférence
*           corrigée pour l'hétéroscédasticité résiduelle.
display _newline "  --- OLS linéaire (SE robustes Huber-White) ---"
regress kru_daugirdas_35 umod if kru_pos == 1, vce(robust)
estimates store tp_ols

local r2_ols    = e(r2)
local rmse_ols  = e(rmse)
local aic_ols   = .
quietly estat ic
matrix _ic_ols = r(S)
local aic_ols = _ic_ols[1,5]

display _newline "  R² OLS   = " %5.3f `r2_ols' ///
    "    RMSE OLS = " %5.3f `rmse_ols' ///
    "    AIC OLS  = " %6.2f `aic_ols'

* Prédiction conditionnelle (utilisée dans 3c)
capture drop kru_cond
predict kru_cond, xb
label variable kru_cond "E[KRU | KRU>0, UMOD] — OLS"
replace kru_cond = 0 if kru_cond < 0

* Diagnostic des résidus OLS
capture drop resid_p2
predict resid_p2 if e(sample), resid

display _newline "  Normalité des résidus OLS (Shapiro-Wilk) :"
swilk resid_p2

* --- 3b.ii  MFP (fractional polynomial multivariable) ---
display _newline(2) "  --- MFP (fractional polynomial) ---"
display         "  Recherche automatique de la meilleure transformation"
display         "  de UMOD parmi puissances {−2, −1, −0.5, 0=log, 0.5, 1, 2, 3}"

mfp: regress kru_daugirdas_35 umod if kru_pos == 1
estimates store tp_mfp

local r2_mfp   = e(r2)
local rmse_mfp = e(rmse)
quietly estat ic
matrix _ic_mfp = r(S)
local aic_mfp = _ic_mfp[1,5]

display _newline "  R² MFP   = " %5.3f `r2_mfp' ///
    "    RMSE MFP = " %5.3f `rmse_mfp' ///
    "    AIC MFP  = " %6.2f `aic_mfp'

* Prédiction MFP (uniquement à titre de comparaison ; non utilisée dans le two-part)
capture drop kru_cond_mfp
predict kru_cond_mfp, xb
label variable kru_cond_mfp "E[KRU | KRU>0, UMOD] — MFP"

* --- 3b.iii  Comparaison OLS vs MFP ---
display _newline(2) "  --- Comparaison OLS vs MFP (N=89 non-anuriques) ---"
estimates stats tp_ols tp_mfp

display _newline "  ───────────────────────────────────────────"
display         "                    OLS         MFP"
display         "  ───────────────────────────────────────────"
display         "  R²              " %6.3f `r2_ols'   "      " %6.3f `r2_mfp'
display         "  RMSE            " %6.3f `rmse_ols' "      " %6.3f `rmse_mfp'
display         "  AIC             " %6.2f `aic_ols'  "    " %6.2f `aic_mfp'
display         "  ───────────────────────────────────────────"
display         "  ΔAIC (OLS − MFP) = " %5.2f (`aic_ols' - `aic_mfp')
display         "    > 2 : MFP préférable"
display         "    < −2 : OLS préférable"
display         "    |ΔAIC| < 2 : équivalents → préférer le plus simple (OLS)"

* ===========================================================================
*  3c. COMBINAISON TWO-PART : E[KRU | UMOD] = P(KRU>0|UMOD) × E[KRU|KRU>0,UMOD]
* ===========================================================================
display _newline(2) "=============================================="
display              "  3c. Prédiction two-part combinée"
display              "=============================================="

capture drop kru_pred_2p
gen double kru_pred_2p = p_pos * kru_cond if !missing(p_pos, kru_cond)
label variable kru_pred_2p "KRU prédit two-part (mL/min/35L)"

* Performance globale (RMSE, MAE, corrélation observé/prédit)
capture drop _r_sq _r_abs
gen double _r_sq  = (kru_daugirdas_35 - kru_pred_2p)^2
gen double _r_abs = abs(kru_daugirdas_35 - kru_pred_2p)

quietly summarize _r_sq
local rmse_2p = sqrt(r(mean))
quietly summarize _r_abs
local mae_2p = r(mean)
drop _r_sq _r_abs

display _newline "  RMSE two-part = " %5.3f `rmse_2p' " mL/min/35L"
display         "  MAE  two-part = " %5.3f `mae_2p'  " mL/min/35L"

quietly corr kru_daugirdas_35 kru_pred_2p
display "  Corrélation Pearson (observé,prédit) = " %5.3f r(rho)

quietly spearman kru_daugirdas_35 kru_pred_2p
display "  Corrélation Spearman                 = " %5.3f r(rho)

* ===========================================================================
*  3d. GRAPHIQUES DIAGNOSTIQUES
* ===========================================================================
* P(KRU>0) prédite vs UMOD
twoway (line p_pos umod, sort lcolor(red) lwidth(medium)) ///
       (scatter kru_pos umod, msize(small) mcolor(%40) jitter(2)), ///
    yline(0.5, lpattern(dot) lcolor(gray)) ///
    title("Partie 1 : P(KRU>0) selon UMOD") ///
    xtitle("UMOD (ng/mL)") ytitle("P(KRU>0)") ///
    legend(order(1 "Logit" 2 "Observé") position(11) ring(0)) ///
    name(tp_part1, replace)

* OLS sur non-anuriques (observé vs UMOD)
twoway (scatter kru_daugirdas_35 umod if kru_pos == 1, msize(small)) ///
       (lfit kru_daugirdas_35 umod if kru_pos == 1, lcolor(red)), ///
    title("Partie 2 : KRU ~ UMOD (non-anuriques)") ///
    xtitle("UMOD (ng/mL)") ytitle("KRU Daugirdas (mL/min/35L)") ///
    legend(off) name(tp_part2, replace)

* Observé vs prédit two-part (tous patients)
twoway (scatter kru_daugirdas_35 kru_pred_2p, msize(small)) ///
       (function y=x, range(0 10) lcolor(red) lpattern(dash)), ///
    title("KRU observé vs prédit — two-part") ///
    xtitle("KRU prédit (mL/min/35L)") ///
    ytitle("KRU observé (mL/min/35L)") ///
    legend(off) name(tp_obs_pred, replace)

* Prédit two-part en fonction d'UMOD (courbe de la prédiction combinée)
twoway (line kru_pred_2p umod, sort lcolor(red) lwidth(medium)) ///
       (scatter kru_daugirdas_35 umod, msize(small) mcolor(%40) jitter(1)), ///
    title("Two-part : KRU prédit selon UMOD") ///
    xtitle("UMOD (ng/mL)") ytitle("KRU (mL/min/35L)") ///
    legend(order(1 "Prédit two-part" 2 "Observé") position(11) ring(0)) ///
    name(tp_curve, replace)

* ###########################################################################
* SECTION 4 — CUT-OFF UMOD POUR PRÉDIRE KRU ≥ 2 mL/min/35L
*              Analyse ROC + indice de Youden sur UMOD seul.
*
*    Choix méthodologique : analyse sur la POPULATION ENTIÈRE (N=151),
*    anuriques inclus (KRU=0 par définition). Justification :
*    de manière analogue à Wong et al. (Kidney International 2015),
*    l'objectif clinique d'un biomarqueur sérique est précisément
*    d'éviter la récolte urinaire ; stratifier a priori sur le statut
*    anurique annulerait ce bénéfice. La performance s'évalue donc
*    dans la population réelle d'application (tous les patients HD).
*    Wong et al. ont obtenu AUC = 0.903 (modeling, N=191, 34% anuriques)
*    et AUC = 0.948 (validation, N=40, 42.5% anuriques) selon cette
*    même logique avec β2M + β-trace protein.
* ###########################################################################

* ===========================================================================
*  4a. ROC : UMOD prédit KRU≥2  (population entière, N=151)
* ===========================================================================
display _newline(2) "=============================================="
display              "  4a. ROC UMOD ~ KRU≥2 — population entière"
display              "=============================================="

roctab kru_ge2 umod, graph summary ///
    title("ROC : UMOD prédit KRU≥2 (N=151)") ///
    name(roc_all, replace)
local app_AUC = r(area)

* ===========================================================================
*  4b. CUT-OFF OPTIMAL (Youden) — population entière
* ===========================================================================
display _newline(2) "=============================================="
display              "  4b. Cutoff Youden — population entière"
display              "=============================================="

* Grille fine (0.5 ng/mL) sur la plage [0 ; max]
quietly summarize umod
local umod_max = r(max)

tempname YouAll
matrix `YouAll' = J(100, 4, .)
local i = 1
local best_J_all = -1
local best_cut_all = .

forvalues u = 0.5(0.5)50 {
    if `u' <= `umod_max' {
        quietly count if umod >= `u' & kru_ge2 == 1 & !missing(kru_ge2)
        local TP = r(N)
        quietly count if umod <  `u' & kru_ge2 == 0 & !missing(kru_ge2)
        local TN = r(N)
        quietly count if umod >= `u' & kru_ge2 == 0 & !missing(kru_ge2)
        local FP = r(N)
        quietly count if umod <  `u' & kru_ge2 == 1 & !missing(kru_ge2)
        local FN = r(N)
        if (`TP' + `FN') > 0 & (`TN' + `FP') > 0 {
            local Se = `TP' / (`TP' + `FN')
            local Sp = `TN' / (`TN' + `FP')
            local J  = `Se' + `Sp' - 1
            matrix `YouAll'[`i', 1] = `u'
            matrix `YouAll'[`i', 2] = `Se'
            matrix `YouAll'[`i', 3] = `Sp'
            matrix `YouAll'[`i', 4] = `J'
            if `J' > `best_J_all' {
                local best_J_all = `J'
                local best_cut_all = `u'
                local best_Se_all = `Se'
                local best_Sp_all = `Sp'
            }
            local i = `i' + 1
        }
    }
}

display _newline "  → Cutoff optimal Youden : UMOD ≥ " %5.2f `best_cut_all' " ng/mL"
display         "      Se = " %5.1f 100*`best_Se_all' " %"
display         "      Sp = " %5.1f 100*`best_Sp_all' " %"
display         "      J  = " %5.3f `best_J_all'

display _newline "  --- Performance à des seuils cliniques candidats ---"
foreach cut in 5 8 10 12 15 {
    quietly count if umod >= `cut' & kru_ge2 == 1 & !missing(kru_ge2)
    local TP = r(N)
    quietly count if umod <  `cut' & kru_ge2 == 0 & !missing(kru_ge2)
    local TN = r(N)
    quietly count if umod >= `cut' & kru_ge2 == 0 & !missing(kru_ge2)
    local FP = r(N)
    quietly count if umod <  `cut' & kru_ge2 == 1 & !missing(kru_ge2)
    local FN = r(N)
    if (`TP' + `FN') > 0 & (`TN' + `FP') > 0 {
        local Se  = `TP' / (`TP' + `FN')
        local Sp  = `TN' / (`TN' + `FP')
        local VPP = cond((`TP' + `FP') > 0, `TP' / (`TP' + `FP'), .)
        local VPN = cond((`TN' + `FN') > 0, `TN' / (`TN' + `FN'), .)
        local J   = `Se' + `Sp' - 1
        display _newline "  Seuil UMOD ≥ `cut' ng/mL :"
        display "    Se=" %5.1f 100*`Se' "%   Sp=" %5.1f 100*`Sp' ///
                "%   VPP=" %5.1f 100*`VPP' "%   VPN=" %5.1f 100*`VPN' ///
                "%   J=" %5.3f `J'
    }
}

* ===========================================================================
*  4c. SYNTHÈSE
* ===========================================================================
display _newline(2) "========================================================"
display              "  SYNTHÈSE — Cutoff UMOD pour KRU ≥ 2 mL/min/35L"
display              "         (population entière, N=151)"
display              "========================================================"
display "  Cutoff optimal Youden = " %5.2f `best_cut_all' " ng/mL"
display "    Se = " %4.1f 100*`best_Se_all' " %"
display "    Sp = " %4.1f 100*`best_Sp_all' " %"
display "    J  = " %5.3f `best_J_all'
display "========================================================"

* ###########################################################################
* SECTION 5 — VALIDATION BOOTSTRAP (Harrell optimism-corrected)
*
*    Objectif : corriger l'optimisme de l'AUC et du cutoff Youden
*    sélectionnés et évalués sur la même population.
*
*    Méthode (Harrell, RMS §5.3) :
*      1. Performance apparente : AUC + cutoff Youden sur N=151
*      2. Pour b=1..B :
*         a) Tirer un échantillon bootstrap (N avec remise)
*         b) Trouver le cutoff Youden c_b dans le bootstrap
*         c) Calculer Se/Sp/J de c_b dans le BOOTSTRAP    (optimiste : bb)
*         d) Calculer Se/Sp/J de c_b dans l'ORIGINAL      (honnête : bo)
*         e) Calculer AUC dans bootstrap (bb) et original (bo)
*      3. Optimisme = mean(metric_bb − metric_bo)
*      4. Métrique corrigée = apparent − optimisme
*
*    Population : entière (N=151), cohérent avec section 4.
* ###########################################################################

* --- Sauvegarde des résultats apparents (section 4) ---
local app_cutoff = `best_cut_all'
local app_Se     = `best_Se_all'
local app_Sp     = `best_Sp_all'
local app_J      = `best_J_all'

* --- Configuration bootstrap ---
global B = 1000
global SEED = 20260522

display _newline(2) "=============================================="
display              "  5. Validation bootstrap (B=$B, seed=$SEED)"
display              "=============================================="
display _newline "  Performance apparente (section 4) :"
display "    AUC          = " %5.3f `app_AUC'
display "    Cutoff opt.  = " %5.2f `app_cutoff' " ng/mL"
display "    Se / Sp / J  = " %5.3f `app_Se' " / " %5.3f `app_Sp' " / " %5.3f `app_J'

set seed $SEED
quietly count if !missing(kru_ge2, umod)
local N = r(N)

* --- Stockage via postfile (compatible Stata/BE, pas de limite de taille) ---
tempname memh
tempfile bootres
postfile `memh' double(AUC_bb cutoff Jbb Jbo Sebb Sebo Spbb Spbo) using `bootres', replace

local n_valid = 0

display _newline _continue "  Progression : "

forvalues b = 1/$B {
    if mod(`b', 100) == 0 display _continue "`b' "

    preserve
    quietly keep if !missing(kru_ge2, umod)
    quietly bsample

    * (a) Cutoff Youden optimal dans le bootstrap
    local best_J_b = -1
    local best_c_b = .
    local best_Se_b = .
    local best_Sp_b = .

    quietly summarize umod
    local umax = r(max)

    forvalues u = 0.5(0.5)50 {
        if `u' <= `umax' {
            quietly count if umod >= `u' & kru_ge2 == 1
            local TP = r(N)
            quietly count if umod <  `u' & kru_ge2 == 0
            local TN = r(N)
            quietly count if umod >= `u' & kru_ge2 == 0
            local FP = r(N)
            quietly count if umod <  `u' & kru_ge2 == 1
            local FN = r(N)
            if (`TP' + `FN') > 0 & (`TN' + `FP') > 0 {
                local Se_ = `TP' / (`TP' + `FN')
                local Sp_ = `TN' / (`TN' + `FP')
                local J_  = `Se_' + `Sp_' - 1
                if `J_' > `best_J_b' {
                    local best_J_b  = `J_'
                    local best_c_b  = `u'
                    local best_Se_b = `Se_'
                    local best_Sp_b = `Sp_'
                }
            }
        }
    }

    * (b) AUC dans le bootstrap
    capture quietly roctab kru_ge2 umod, nograph
    local fit_ok = (_rc == 0)
    if `fit_ok' {
        local auc_bb = r(area)
    }
    else {
        local auc_bb = .
    }

    restore

    if missing(`best_c_b') | !`fit_ok' continue

    * (c) Application du cutoff c_b à l'ORIGINAL
    quietly count if umod >= `best_c_b' & kru_ge2 == 1 & !missing(kru_ge2, umod)
    local TPo = r(N)
    quietly count if umod <  `best_c_b' & kru_ge2 == 0 & !missing(kru_ge2, umod)
    local TNo = r(N)
    quietly count if umod >= `best_c_b' & kru_ge2 == 0 & !missing(kru_ge2, umod)
    local FPo = r(N)
    quietly count if umod <  `best_c_b' & kru_ge2 == 1 & !missing(kru_ge2, umod)
    local FNo = r(N)

    if (`TPo' + `FNo') == 0 | (`TNo' + `FPo') == 0 continue

    local Se_bo = `TPo' / (`TPo' + `FNo')
    local Sp_bo = `TNo' / (`TNo' + `FPo')
    local J_bo  = `Se_bo' + `Sp_bo' - 1

    post `memh' (`auc_bb') (`best_c_b') (`best_J_b') (`J_bo') ///
                (`best_Se_b') (`Se_bo') (`best_Sp_b') (`Sp_bo')

    local n_valid = `n_valid' + 1
}

postclose `memh'

display ""
display _newline "  Itérations valides : `n_valid'/$B"

* --- Agrégation : charger le dataset bootstrap ---
preserve
quietly use `bootres', clear

quietly summarize AUC_bb, meanonly
local m_auc_bb = r(mean)
quietly _pctile AUC_bb, percentiles(2.5 97.5)
local auc_lo = r(r1)
local auc_hi = r(r2)

quietly summarize Jbb, meanonly
local m_Jbb = r(mean)
quietly summarize Jbo, meanonly
local m_Jbo = r(mean)
local opt_J = `m_Jbb' - `m_Jbo'

quietly summarize Sebb, meanonly
local m_Sebb = r(mean)
quietly summarize Sebo, meanonly
local m_Sebo = r(mean)
local opt_Se = `m_Sebb' - `m_Sebo'

quietly summarize Spbb, meanonly
local m_Spbb = r(mean)
quietly summarize Spbo, meanonly
local m_Spbo = r(mean)
local opt_Sp = `m_Spbb' - `m_Spbo'

quietly summarize cutoff, detail
local m_cut = r(mean)
local sd_cut = r(sd)
local cut_lo = r(p5)
local cut_hi = r(p95)

quietly count if cutoff == `app_cutoff'
local pct_apparent = 100 * r(N) / `n_valid'

restore

* --- Synthèse ---
local Jcorr   = `app_J'  - `opt_J'
local Secorr  = `app_Se' - `opt_Se'
local Spcorr  = `app_Sp' - `opt_Sp'

display _newline(2) "=========================================================="
display              "  RÉSULTATS BOOTSTRAP — VALIDATION INTERNE"
display              "=========================================================="
display _newline "  --- AUC ---"
display "    AUC apparente                 = " %5.3f `app_AUC'
display "    AUC bootstrap (moyenne)       = " %5.3f `m_auc_bb'
display "    IC 95% bootstrap (percentile) = " %5.3f `auc_lo' " — " %5.3f `auc_hi'

display _newline "  --- Cutoff optimal Youden ---"
display "    Cutoff apparent      = " %5.2f `app_cutoff' " ng/mL"
display "    Cutoff bootstrap moy = " %5.2f `m_cut' " ng/mL (SD " %5.2f `sd_cut' ")"
display "    IC bootstrap (5–95%) = " %5.2f `cut_lo' " — " %5.2f `cut_hi' " ng/mL"
display "    % itérations retrouvant le cutoff apparent (" %3.1f `app_cutoff' ") : " %4.1f `pct_apparent' " %"

display _newline "  --- Performance au cutoff (corrigée pour optimisme) ---"
display "                         Apparent    Optimisme    Corrigé"
display "    Sensibilité    " %6.3f `app_Se'  "       " %6.3f `opt_Se' "      " %6.3f `Secorr'
display "    Spécificité    " %6.3f `app_Sp'  "       " %6.3f `opt_Sp' "      " %6.3f `Spcorr'
display "    Youden J       " %6.3f `app_J'   "       " %6.3f `opt_J'  "      " %6.3f `Jcorr'

display _newline "  Interprétation :"
display "  - AUC IC 95%   → précision de l'estimation"
display "  - Cutoff IC    → stabilité du seuil sélectionné"
display "  - J corrigé    → performance attendue sur de nouveaux patients"
display "  - Optimisme    → biais dû à la sélection du cutoff sur les mêmes données"
display "=========================================================="