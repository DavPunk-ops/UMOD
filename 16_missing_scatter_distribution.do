* ===========================================================================
* 16_missing_scatter_distribution.do
* Prérequis : lancer 02_calculate_kru.do AVANT ce do-file.
*
* Objectif : répondre aux minor comments R3 #4, #5, #6.
*   Autonome ; ne modifie AUCUNE do-file existante.
*
* Sections :
*   A. AUDIT DES DONNÉES MANQUANTES (R3 #4)
*      → N non manquant / manquant pour chaque variable clé (sur N=151)
*   B. SCATTERPLOTS biomarqueur vs KRU mesuré, non-anuriques (R3 #5)
*      → UMOD~KRU et β2M~KRU + lowess + Spearman ; exports PNG
*   C. DISTRIBUTION DE L'UROMODULINE (R3 #6)
*      → histogramme + nombre de valeurs ≤ 2.0 ng/mL (limite de détection)
*        et = 0 (undetectable) ; export PNG
*
* Chemin d'export des graphes :
local gpath "C:\Users\dajs\OneDrive - HOPITAUX UNIVERSITAIRES DE GENEVE\recherche\RKF\UMOD\stata\main prospective study\with Claude"
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
capture confirm variable kru_ge2
if _rc  gen byte kru_ge2 = (kru_daugirdas_35 >= 2) if !missing(kru_daugirdas_35)
capture confirm variable kru_pos
if _rc  gen byte kru_pos = (kru_daugirdas_35 > 0)  if !missing(kru_daugirdas_35)

* --- Restriction aux 151 patients analysés ---
keep if !missing(kru_daugirdas_35)
display _newline "  Population analysée : N = " _N

* ###########################################################################
* SECTION A — AUDIT DES DONNÉES MANQUANTES  (R3 #4)
* ###########################################################################
display _newline(2) "########################################################"
display              "  A. Audit des données manquantes (N=151)"
display              "########################################################"

* Liste des variables clés à auditer (adapter les noms si besoin)
local auditvars ///
    age sex race bmi bsa V_watson charlson dm kidneydisease ///
    mode uf spktv vintage regimen sessiontime ///
    prehdsbp prehddbp posthdweight ///
    kru_daugirdas_35 urinevolume ///
    umod labb2mprehd ///
    diuretic antiht ado insulin lipid epo pobinder kbinder vitdanalog bicarbonate

display _newline "  Variable                     N non-manquant   N manquant   % manquant"
display      "  ----------------------------------------------------------------------"
foreach v of local auditvars {
    capture confirm variable `v'
    if _rc {
        display "  " %-28s "`v'" "   [variable absente]"
    }
    else {
        quietly count if !missing(`v')
        local nok = r(N)
        local nmiss = _N - `nok'
        local pmiss = 100*`nmiss'/_N
        display "  " %-28s "`v'" %10.0f `nok' %13.0f `nmiss' %12.1f `pmiss'
    }
}

* Vue d'ensemble misstable (variables présentes uniquement)
display _newline "  --- misstable summarize (variables avec au moins 1 manquant) ---"
misstable summarize `auditvars'

* ###########################################################################
* SECTION B — SCATTERPLOTS biomarqueur vs KRU mesuré (R3 #5)
* ###########################################################################
display _newline(2) "########################################################"
display              "  B. Association biomarqueur ~ KRU mesuré (non-anuriques)"
display              "########################################################"

quietly count if kru_pos==1 & !missing(umod)
display "  N non-anuriques avec UMOD  : " r(N)
quietly count if kru_pos==1 & !missing(labb2mprehd)
display "  N non-anuriques avec β2M   : " r(N)

display _newline "  --- Spearman (non-anuriques) ---"
spearman kru_daugirdas_35 umod        if kru_pos==1, stats(rho p)
spearman kru_daugirdas_35 labb2mprehd if kru_pos==1, stats(rho p)

* --- Scatter UMOD vs KRU (annotation Spearman, coin haut-gauche) ---
quietly spearman kru_daugirdas_35 umod if kru_pos==1
local rho_u : display %4.2f r(rho)
quietly summarize kru_daugirdas_35 if kru_pos==1 & !missing(umod)
local xu = r(min) + 0.02*(r(max)-r(min))
quietly summarize umod if kru_pos==1
local yu = 0.97*r(max)
twoway (scatter umod kru_daugirdas_35 if kru_pos==1, mcolor(navy%60) msize(small)) ///
       (lowess umod kru_daugirdas_35 if kru_pos==1, lcolor(cranberry) lwidth(medthick)), ///
    xtitle("Measured KRU (mL/min/35 L)") ytitle("Serum uromodulin (ng/mL)") ///
    title("Serum uromodulin vs measured KRU") ///
    subtitle("Non-anuric patients") legend(off) scheme(s1mono) ///
    text(`yu' `xu' "Spearman {&rho} = `rho_u' (P < 0.001)", place(e) size(medsmall) color(black)) ///
    name(umod_kru, replace)
graph save   umod_kru "`gpath'\FigS_umod_vs_kru.gph", replace
graph export "`gpath'\FigS_umod_vs_kru.png", replace width(2000)

* --- Scatter β2M vs KRU (annotation Spearman, coin haut-droit) ---
quietly spearman kru_daugirdas_35 labb2mprehd if kru_pos==1
local rho_b : display %4.2f r(rho)
quietly summarize kru_daugirdas_35 if kru_pos==1 & !missing(labb2mprehd)
local xb = r(max) - 0.02*(r(max)-r(min))
quietly summarize labb2mprehd if kru_pos==1
local yb = 0.97*r(max)
twoway (scatter labb2mprehd kru_daugirdas_35 if kru_pos==1, mcolor(navy%60) msize(small)) ///
       (lowess labb2mprehd kru_daugirdas_35 if kru_pos==1, lcolor(cranberry) lwidth(medthick)), ///
    xtitle("Measured KRU (mL/min/35 L)") ytitle("Serum {&beta}2-microglobulin (mg/L)") ///
    title("Serum {&beta}2-microglobulin vs measured KRU") ///
    subtitle("Non-anuric patients") legend(off) scheme(s1mono) ///
    text(`yb' `xb' "Spearman {&rho} = `rho_b' (P < 0.001)", place(w) size(medsmall) color(black)) ///
    name(b2m_kru, replace)
graph save   b2m_kru "`gpath'\FigS_b2m_vs_kru.gph", replace
graph export "`gpath'\FigS_b2m_vs_kru.png", replace width(2000)

display _newline "  Graphes exportés : FigS_umod_vs_kru.png , FigS_b2m_vs_kru.png"

* ###########################################################################
* SECTION C — DISTRIBUTION DE L'UROMODULINE  (R3 #6)
* ###########################################################################
display _newline(2) "########################################################"
display              "  C. Distribution de l'uromoduline (N=151)"
display              "########################################################"

display _newline "  --- Résumé UMOD ---"
summarize umod, detail

quietly count if umod == 0
display _newline "  UMOD = 0 (undetectable)      : " r(N)
quietly count if umod <= 2.0 & !missing(umod)
display "  UMOD ≤ 2.0 ng/mL (≤ LOD)     : " r(N)
quietly count if umod < 2.0 & umod > 0
display "  UMOD entre 0 (excl) et 2.0   : " r(N)
quietly count if !missing(umod)
display "  UMOD non manquant (total)    : " r(N)

* Répartition ≤2.0 par statut anurique
display _newline "  --- UMOD ≤ 2.0 par statut anurique ---"
gen byte umod_le2 = (umod <= 2.0) if !missing(umod)
label define le2lbl 0 ">2.0" 1 "≤2.0 (≤LOD)", replace
label values umod_le2 le2lbl
tab umod_le2 kru_pos, col

* --- Histogramme avec ligne verticale à la limite de détection (2.0) ---
histogram umod, width(2) frequency ///
    fcolor(navy%50) lcolor(navy) ///
    xline(2.0, lcolor(cranberry) lpattern(dash) lwidth(medthick)) ///
    xtitle("Serum uromodulin (ng/mL)") ytitle("Number of patients") ///
    title("Distribution of serum uromodulin") ///
    subtitle("Dashed line = 2.0 ng/mL detection limit") ///
    scheme(s1mono) name(umod_hist, replace)
graph save   umod_hist "`gpath'\FigS_umod_distribution.gph", replace
graph export "`gpath'\FigS_umod_distribution.png", replace width(2000)

display _newline "  Graphe exporté : FigS_umod_distribution.png"

display _newline(2) "=== FIN do-file 16 ==="
