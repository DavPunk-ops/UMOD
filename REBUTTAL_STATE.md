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
  Robustesse exclure UMOD≤2 (do-12§C) : combiné 0.833→0.841.
  Processing : tube sérum (Hemogard jaune, clot-activator), coagulation 30 min RT, centrif 1300g×10min RT,
  aliquots Sarstedt, congélation −20°C(1–4h) puis −80°C, décongélation UNIQUE.
- Excel dosages_UMOD.xlsx ↔ .dta : correspondance VÉRIFIÉE ligne-à-ligne 151/151 (0 mismatch).
  2 exclus = 126-T0, 143-T0 (collecte incomplète, UMOD valides). Documenté dans en-tête do-01.
- Fréquence par KRU (do-13/regimen, N=148) : KRU<2 (n=89) 83 3×/sem 4 2×/sem 1 4× 1 6× ;
  KRU≥2 (n=59) 14 3×/sem 45 2×/sem.
- Flowchart (Fig S1) : 164 évalués → 11 exclus (4 refus, 3 quitté centre, 2 consent impossible,
  1 transplant fonctionnel, 1 macrohématurie) → 153 inclus → 2 exclus (collecte incomplète) → 151.
- Anurie : définie par auto-report patient (« ≥200 mL/j ? »), PAS mesure volumétrique. Collecte unique au 1er essai.

## Statut des réponses
- ✅ FAITS : AE (tout) ; R1#1 (ADPKD+variabilité) ; R3 #1,#2,#3,#5,#6,#7,#8 ; réfs Casino ×3.
- ✅ Minors R3 rédigés (texte prêt, données OK) : m1 (cross-sectional), m2 (screened/declined),
  m3 (fréquence/KRU), m6 (distribution UMOD — texte), m7 (PPV/prévalence).
- ⏳ RESTE :
  - **do-15 à créer** : m5 scatterplots UMOD/β2M vs KRU (non-anur) + m6 histogramme UMOD + m4 audit missing.
  - **R3#4** : justif. cutoff 200 mL + RÉFÉRENCE (l'utilisateur cherche) + collecte unique 1er essai + <200 non mesuré (sensibilité 400 mL do-08 déjà faite).
  - **R2 ×2** : (1) KRU pas seul critère HD incrémentale ; (2) faciliter adoption HD incrémentale. Pure Discussion.
  - **R1 #3-5** : footnote Table 1 (collecte chronométrée) ; légendes figures (grey zone/axe y/N) ; énoncé décisionnel pratique.
  - **Édits manuscrit** : harmoniser /35 L ; Methods anurie (patient-reported) ; demi-phrase CV bas UMOD ;
    Methods LOD garder « detection limit 2.0 » (verbatim EUROIMMUN, correct) ; Discussion (screening, longitudinal, PPV/prévalence) ; Abstract.

## Missing data (partiel, à confirmer par audit do-15)
umod 151 · β2M 148 (3 manquants) · spktv 142 · mode 147 · uf 149 · vintage 149 · regimen 148 ·
posthdweight/bmi/bsa 150 · prehdsbp/dbp 149 · kidneydisease 150 · diuretic 150 · reste 151.

## Git
Dernier commit `5102f74` (do-11 §5b + Table S3 étendue). Arbre propre.
