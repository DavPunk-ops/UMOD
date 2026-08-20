# UMOD Rebuttal — État de reprise (compact)

## Projet
Révision majeure CJASN : « Estimating Residual Kidney Function in HD Patients Using
Serum Uromodulin » (N=151, 2 centres). Réponse point-par-point à AE + 3 reviewers.
Auteur : David A. Jaques (néphrologue, HUG Genève).

## Consignes permanentes
- Branche git : `claude/check-umod-access-QEC1O`. Commit+push pré-autorisés (non-critiques).
- NE PAS modifier do-files originales 01–07 (commentaires d'en-tête OK si demandé).
- Nouvelles analyses = nouveaux do-files (08+).
- Unité KRU : **/35 L partout** (harmoniser le manuscrit ; PAS /1.73 m²).
- Réponses : concises, rigoureuses, chaque chiffre = sortie Stata réelle, honnêteté sur limites.
- Le TEXTE des réponses vit dans le .docx de l'utilisateur (PAS dans git).

## Do-files rebuttal (08–14)
- 08_anuria400_sensitivity — 400 mL (13 reclassés, tous KRU<2, AUC inchangée)
- 09_nonanuric_logit — Table S1 (non-anuriques : applied+refit)
- 10_subgroups_adpkd — Table S2 (sexe/âge/diabète/Charlson/ADPKD) + compare UMOD/β2M/KRU par ADPKD
- 11_clinical_thresholds — équation, calibration, Tables 4/5, S3 (+§5b breakdown seuils stricts)
- 12_nonanuric_robustness — §A refit+bootstrap, §B bandes/borderline, §C exclure UMOD≤2, §D <400mL, §E calib
- 13_b2m_dialysis — R3#7 (corrélations + modèle ajusté dialyse ; sessiontime+regimen)
- 14_kru_normalization — R3#5 (Table S5 : brut/35L/BSA)
- 15_table1_dialysis_frequency — R3 m3 (ligne fréquence Table 1 ; catégoriel + collapsé 2×/≥3×)
- 16_missing_scatter_distribution — R3 m4/m5/m6 (audit missing + scatterplots UMOD/β2M~KRU + histogramme UMOD)
- (originales : 01 merge, 02 KRU, 03 Table1, 04 univar, 05 multivar, 06 quantitatif continu, 07 temporel)

## Chiffres clés (vérifiés Stata)
- Cohorte : 151 analysés (89 non-anuriques dont 87 données β2M complètes ; 62 anuriques).
- AUC full cohort : UMOD 0.891, β2M 0.877, combiné **0.923**.
- Non-anuriques (Table S1) : UMOD 0.773 (0.668–0.879), β2M 0.802 (0.705–0.898),
  combiné **0.833** (0.746–0.920), optimism-corr **0.821**. Refit OR UMOD 1.10 (p=0.023), β2M 0.86 (p=0.003).
- Borderline KRU 1–3 (N=33, do-12§B) : combiné AUC 0.748 (0.57–0.93), UMOD 0.70, β2M 0.71.
  Spearman non-anur : UMOD ρ=0.45, β2M ρ=−0.57 (p<0.001). Continu = Figure 4 (déjà publiée).
- Équation : logit(P) = 2.453 + 0.157·UMOD − 0.175·β2M ; OR 1.17 / 0.84.
- Calibration : HL p=0.53 ; pente optimism-corr 0.96 ; Table S4 (déciles).
- Table 4 (full, N=148) : rule-out P̂<0.24 76(51.4%) NPV93.4% 5faux ; grey 22(14.9%) ;
  rule-in P̂≥0.55 50(33.8%) PPV82.0% 9faux ; 85.1% classés.
- Table 5 (non-anur, N=87) : rule-out 21(24.1%) NPV76.2%(54.9–89.4) 5faux ;
  grey 17(19.5%) ; rule-in 49(56.3%) PPV83.7%(71.0–91.5) 8faux ; 80.5% classés.
- Table S3 étendue (seuils conservateurs, non-anur, rule-out fixe P̂<0.24) :
  Sp≥90% P̂0.55 rule-in49(56.3) PPV83.7 8faux grey17(19.5) classé80.5% ;
  Sp≥95% P̂0.69 rule-in39(44.8) PPV92.3(79.7–97.3) 3faux grey27(31.0) classé69.0% ;
  Sp≥97.5% P̂0.82 rule-in33(37.9) PPV93.9(80.4–98.3) 2faux grey33(37.9) classé62.1%.
- R3#7 β2M×dialyse : corr vintage+0.48 UF+0.45 spKt/V+0.18 session+0.26 ; HDF>HD p=0.014 ;
  3×/sem β2M 29.9 vs 2×/sem 22.5 (KW p=0.0001) ; modèle pleinement ajusté OR β2M 0.79 (0.67–0.92) p=0.004 N=133.
- R3#5 normalisation (S5) : ρ 0.99–1.00 entre échelles ; ≤2/151 reclassés ; AUC combiné brut0.921/35L0.923/BSA0.920.
  Watson sur-estime volume urée (Kloppenburg KI 2001) → d'où sensibilité. Réfs : Casino&Basile NDT 2017 ; Casino NDT 2026 ; Watson 1980.
- ADPKD (n=16) : UMOD 11.5 vs 7.1 (p=0.51), β2M 24.2 vs 27.3 (p=0.56), KRU 0.83 vs 0.66 (p=0.92) → NS.
- Assay UMOD : dosé EN DUPLICAT (moyenne rendue ; β2M unique). LOD EUROIMMUN 2.0 ng/mL (blanc+3SD, verbatim).
  Calibrateurs 0/25/50/100/200/400 (min non-nul 25). CV fabricant intra≤3.2% inter≤7.8% à 30–228 ng/mL.
  Cohorte médiane UMOD 7.4 (IQR 2.6–16.3). CV>20% observé : 19/151 (12.6%), médiane 2.6 vs remaining 8.3.
  Below detection : 10 undetectable → codés 0, TOUS anuriques ; +17 <2.0 → 27 total <2.0.
  Croisement ≤2.0 (do-16) : anuriques 25, KRU<2 26, KRU≥2 1. Indétectables(=0) : 10, tous anuriques/KRU<2.
  ⚠ Patient discordant 128-T0 : UMOD 1.117 (FIABLE, pas de flag CV>20%), KRU/35L 2.82, urine 1800 mL → KRU≥2.
    = vraie discordance biologique (faux négatif genuine), PAS artéfact de seuil ni de mesure. Autre non-anur ≤2.0 = 041-T0 (KRU<2, concordant).
  CORRECTION R3#6 : le surclaim « below-LOD values correctly classified KRU<2 irrespective » était FAUX (26/27, pas 27/27).
    Phrasé corrigé (retenu) : « 10 zéros tous anuriques → assignation non biaisée ; exclure ≤2.0 → AUC 0.833→0.841 ». Sans revendiquer l'universalité.
    Le patient discordant est disclosé UNE fois (minor 6 : « only a single patient… KRU≥2 »), R3#6 ne le contredit plus.
  Robustesse exclure UMOD≤2 (do-12§C) : combiné 0.833→0.841 (exclut justement 128-T0 → AUC monte, cohérent).
  Processing : tube sérum (Hemogard jaune, clot-activator), coagulation 30 min RT, centrif 1300g×10min RT,
  aliquots Sarstedt, congélation −20°C(1–4h) puis −80°C, décongélation UNIQUE.
- Excel dosages_UMOD.xlsx ↔ .dta : correspondance VÉRIFIÉE ligne-à-ligne 151/151 (0 mismatch).
  2 exclus = 126-T0, 143-T0 (collecte incomplète, UMOD valides). Documenté dans en-tête do-01.
- Fréquence par KRU (do-13/regimen, N=148) : KRU<2 (n=89) 83 3×/sem 4 2×/sem 1 4× 1 6× ;
  KRU≥2 (n=59) 14 3×/sem 45 2×/sem.
- Flowchart (Fig S1) : 164 évalués → 11 exclus (4 refus, 3 quitté centre, 2 consent impossible,
  1 transplant fonctionnel, 1 macrohématurie) → 153 inclus → 2 exclus (collecte incomplète) → 151.
- Anurie : définie par auto-report patient (« ≥200 mL/j ? »), PAS mesure volumétrique. Collecte unique au 1er essai.
- R3 m3 fréquence (do-15, N=148, 3 manquants) : 2×/sem KRU<2 4(4.5%)/KRU≥2 45(76.3%) ; ≥3×/sem 85(95.5%)/14(23.7%) ; Fisher p<0.001.
- R3 m4 audit missing (do-16) : tout ≥94% complet ; 100% pour umod/kru/age/sex/race/charlson/dm/Vwatson ;
  β2M 148(98.0%) ; regimen 148 ; mode 147(97.4%) ; spktv 142(94.0%, le moins complet) ;
  urinevolume 89 = manquant STRUCTUREL (62 anuriques sans collecte). Aucune imputation, complete-case.
- R3 m5 scatter (do-16) : Spearman non-anur UMOD ρ=0.45 (n=89), β2M ρ=−0.57 (n=87), p<0.001. PNG FigS_umod/b2m_vs_kru.
- R3 m6 distribution UMOD (do-16) : médiane 7.4 (IQR 2.6–16.3), range 0–46 ; 27/151(17.9%) ≤2.0 ng/mL dont 10 =0 ;
  répartition ≤2.0 : 25 anuriques / 2 non-anuriques (≤LOD ↔ quasi excl. anurie). PNG FigS_umod_distribution.

## Numérotation figures supplémentaires (à harmoniser dans le .docx)
Figure S1 = flowchart (minor 2) ; Figure S2a/S2b = scatterplots UMOD/β2M~KRU (minor 5 + R3#4) ; Figure S3 = histogramme UMOD (minor 6).
Tables S numérotées séparément (S1 non-anur, S2 sous-groupes, S3 seuils, S4 calib, S5 normalisation).

## Statut des réponses
- ✅ FAITS : AE (tout) ; R1#1 (ADPKD+variabilité) ; R3 majors #1,#2,#3,#4,#5,#6,#7,#8 ; réfs Casino ×3.
- ✅ Minors R3 rédigés+data (do-15/16 poussés) : m1 (cross-sectional), m2 (screened/declined),
  m3 (fréquence/KRU + ligne Table 1), m4 (audit missing + Table SX), m5 (scatterplots), m6 (histogramme UMOD), m7 (PPV/prévalence).
- ✅ R3#4 (anurie cutoff 200 mL) FAIT : 4 bullets couvrant les 5 sous-questions —
  (a) déf+refs 1-3 (l'utilisateur a ses refs) ; (b) zéro/two-part (do-06 logit+OLS ; concession honnête de la séparation + réfutation empirique via non-anur Table S1) ;
  (c) collecte 1er essai unique (pas plus grand volume 1sem/1mois) ; (d) <200 mL non mesuré ; (e) sensibilité 400 mL (do-08 : 13 reclassés, 0 changement classif → AUC identique 0.891/0.877/0.923) + <100 infaisable.
  Two-part do-06 = part1 logit P(KRU>0) + part2 OLS E[KRU|KRU>0] sur non-anur (β2M en (β2M/10)^-2, SE robustes).
- ⏳ RESTE :
  - **R2 ×2** : (1) KRU pas seul critère HD incrémentale ; (2) faciliter adoption HD incrémentale. Pure Discussion.
  - **R1 #3-5** : footnote Table 1 (collecte chronométrée) ; légendes figures (grey zone/axe y/N) ; énoncé décisionnel pratique.
  - **Édits manuscrit** : harmoniser /35 L ; Methods anurie (patient-reported) ; demi-phrase CV bas UMOD ;
    Methods LOD garder « detection limit 2.0 » (verbatim EUROIMMUN, correct) ; Discussion (screening, longitudinal, PPV/prévalence) ; Abstract.

## Missing data (CONFIRMÉ par audit do-16 §A, N=151)
100% : age sex race V_watson charlson dm kru_daugirdas_35 umod.
99.3%(150) : bmi bsa kidneydisease posthdweight diuretic ado vitdanalog.
98.7%(149) : uf vintage sessiontime prehdsbp prehddbp antiht insulin lipid pobinder bicarbonate.
98.0%(148) : regimen labb2mprehd epo kbinder.  97.4%(147) : mode.  94.0%(142) : spktv (le moins complet).
urinevolume 89 (58.9%) = STRUCTUREL (62 anuriques sans collecte, KRU=0). Aucune imputation.
→ Supplementary Table SX (missing par variable) rédigée, à déposer dans le .docx supp.

## Git
Branche `claude/check-umod-access-QEC1O`. do-15 + do-16 (scatterplots annotés Spearman + .gph + croisements below-LOD) poussés. Arbre propre.
