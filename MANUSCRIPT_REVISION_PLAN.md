# Plan de révision du manuscrit — UMOD/CJASN

Guide des édits à porter dans le .docx du manuscrit. Principe : **ajouts ciblés,
groupés, renvois > prose**. Ne pas dénaturer le texte existant.

---

## STRUCTURE RESULTS (figée)
1. Description of the study cohort
2. Diagnostic performance of serum uromodulin for predicting KRU ≥2 mL/min/35L
3. Diagnostic performance of combined serum uromodulin and β2-microglobulin
4. **⭐ Discrimination restricted to non-anuric patients with measured KRU** (NOUVELLE, après la 3)
5. Quantitative relationship between serum uromodulin, β2-microglobulin, and KRU (ex-4)
6. **⭐ Sensitivity and robustness analyses** (NOUVELLE ; inclut §β2M-dialyse en paragraphe dédié)

---

## NUMÉROTATION SUPPLÉMENTAIRE (figée, ordre d'apparition)

### FIGURES
| N° | Contenu | Apparition | Source |
|----|---------|------------|--------|
| Fig S1 | Flowchart screening | Début Results | (existe) |
| Fig S2 | Distribution UMOD (ligne LOD 2.0) | R1 Cohort | do-16 FigS_umod_distribution |
| Fig S3a | Scatterplot UMOD vs KRU mesuré (non-anur, ρ=0.45) | R4 | do-16 FigS_umod_vs_kru |
| Fig S3b | Scatterplot β2M vs KRU mesuré (non-anur, ρ=−0.57) | R4 | do-16 FigS_b2m_vs_kru |

### TABLES
| N° | Contenu | Apparition | Source |
|----|---------|------------|--------|
| Table S1 | Missing data (completeness) | R1 Cohort | do-16 §A |
| Table S2 | Calibration (déciles obs/attendu) | R3 | do-11/05 |
| Table S3 | Non-anuric discrimination (AUC, OR) | R4 | do-09 |
| Table S4 | Conservative rule-in thresholds | R4 | do-11 §5 |
| Table S5 | Subgroups (sexe/âge/comorbidité/ADPKD) | R6 | do-10 |
| Table S6 | Normalization (brut/35L/BSA) | R6 | do-14 |

### Mapping ancien→nouveau (find-replace rebuttal + manuscrit)
FIGURES : flowchart S1→S1 · distribution S3→S2 · scatter UMOD S2a→S3a · scatter β2M S2b→S3b
TABLES : missing SX→S1 · calibration S4→S2 · non-anur S1→S3 · conservative S3→S4 · subgroups S2→S5 · normalization S5→S6
⚠ Faire le find-replace via étiquettes temporaires pour éviter les collisions en cascade.

---

## ÉDITS PAR SECTION

### INTRODUCTION
- Corriger typo `mL/min/1.73 m²` → **`mL/min/35 L`**.
- Re-sourcer « ~85 kDa » (réf biochimie/PM, PAS Zakrocka). Garder Zakrocka(1) pour « detectable UMOD ↔ preserved RKF » + « poorly cleared by conventional dialysis ».

### METHODS
- **Study population** : phrase « consecutively screened… » + renvoi Fig S1.
- **KRU/RKF** : anurie = auto-report patient (« ≥200 mL/j ? »), collecte chronométrée unique au 1er essai, KRU=0 sans collecte, une minorité de volumes collectés peut finir <200 mL. KRU /35 L (Watson) partout ; 1 phrase limites Watson.
- **Assays** : UMOD en duplicat (moyenne ; CV>20% flaggé) ; LOD 2.0 ng/mL (blanc+3SD) ; undetectable→0. Pré-analytique (tube clot-activator, coag 30 min, centrif 1300g×10min, congélation −20→−80°C, décongélation unique).
- **Statistics** : two-part (logit+OLS) ; bootstrap Harrell B=1000 ; calibration (HL+pente) ; two-cutoff (Se≥90%/Sp≥90% a priori + cibles Sp conservatrices) ; sensitivity analyses (nommer) ; missing data (complete-case, pas d'imputation).

### RESULTS
**1. Study population** (Table 1)
- + ligne fréquence de dialyse + footnote ; renvoi Fig S1 ; renvoi Table S1 (missing).
- 1 phrase distribution UMOD / below-detection (27 ≤2.0 dont 10 undetectable) + renvoi Fig S2.

**2. UMOD alone** (Fig 1, Fig 2) — inchangé sur le fond.

**3. Combined** (Fig 3, Table 4)
- + équation complète de prédiction : `logit(P) = 2.453 + 0.157·UMOD − 0.175·β2M`.
- + calibration : HL p=0.53, pente optimism-corrected 0.96, renvoi Table S2.
- + mal-classés Table 4 (révisée) + colonne bootstrap-corrected (NPV 92.3% [87.7–96.8], PPV 80.9% [75.0–90.0]).

**4. ⭐ Discrimination restricted to non-anuric patients with measured KRU** (NOUVELLE — le cœur)
- N=89 (87 biomarqueurs complets). Combiné AUC 0.833 (0.746–0.920 ; opt-corr 0.820) ; UMOD 0.773 ; β2M 0.802. OR refit : UMOD 1.10 (p=0.023), β2M 0.86 (p=0.003). Renvoi Table S3.
- Two-cutoff non-anur (Table 5) : rule-out 24.1% (NPV 76.2%), grey 19.5%, rule-in 56.3% (PPV 83.7%), 80.5% classés.
- Seuils conservateurs (Table S4).
- Bande borderline KRU 1–3 (N=33 ; AUC 0.748 [0.57–0.93] ; UMOD 0.70 ; β2M 0.71).
- Scatterplots UMOD/β2M vs KRU (Fig S3a/S3b ; ρ 0.45 / −0.57).

**5. Quantitative relationship (two-part)** (Fig 4) — inchangé.

**6. ⭐ Sensitivity and robustness analyses** (NOUVELLE)
- *Robustesse discrimination* : anurie 400 mL (13 reclassés, AUC inchangée 0.891/0.877/0.923) ; sous-groupes/ADPKD AUC 0.91–0.94 (Table S5) ; normalisation brut/35L/BSA AUC 0.921/0.923/0.920 (Table S6) ; exclusion UMOD≤2 (0.833→0.841).
- *Specificity of the β2-microglobulin signal* (paragraphe dédié) : corrélations positives β2M×dialyse (vintage+0.48, UF+0.45, spKt/V+0.18, session+0.26) ; plus haute en HDF (p=0.014) et 3×/sem (29.9 vs 22.5, p<0.001) → sens opposé à l'effet clairance ; OR ajusté 0.79 (0.67–0.92, p=0.004) survit au sur-ajustement (dont regimen, downstream) ; validation interne pente 0.96.

### DISCUSSION (ajouts groupés)
- Cadrage screening (estimer RKF, pas remplacer la collecte ; complément clinique).
- KRU = un seul composant de la décision de régime (R2#1).
- Usage two-threshold en routine (rule-out→éviter / grey→confirmer / rule-in→candidat).
- Cross-sectional vs longitudinal : monitoring sériel = perspective future (variabilité intra-personnelle + Δuromoduline vs ΔKRU).
- PPV dépend de la prévalence / case-mix (centres incrémentaux, patients incidents).
- Validation externe indispensable.
- Positionnement Zakrocka (précédent appuyant la biologie).
- Limitations (1 paragraphe) : mesure unique, précision d'assay à basse concentration, absence validation externe, N modeste, ADPKD peu nombreux.

### TABLES
- Table 1 : + ligne fréquence + footnote « KRU based on interdialytic timed collection ; anuric self-reported <200 mL assigned 0 » + note « frequency available for 148/151 ».
- Table 4 révisée (mal-classés + colonne bootstrap-corrected).
- Table 5 nouvelle (non-anur two-cutoff).
- Supplément : S1–S6 + Fig S1–S3 (voir numérotation).

### ABSTRACT
- /35 L (corriger si /1.73 présent).
- Cadrage screening dans la conclusion.
- (optionnel) robustesse non-anur AUC 0.833.

### CORRECTIONS PONCTUELLES
- optimism-corrected non-anur = **0.820** (pas 0.821) si présent dans le manuscrit.
- « 85 kDa » : re-sourcer (pas Zakrocka).

---

## CHIFFRES CLÉS VÉRIFIÉS (rappel — tous audités niveau 1)
Full cohort : UMOD 0.891 (0.838–0.944) · β2M 0.877 (0.821–0.933) · combiné 0.923 (0.883–0.964).
OR combiné full : UMOD 1.17 (1.09–1.26) · β2M 0.84 (0.77–0.92). Équation : 2.453+0.157·UMOD−0.175·β2M.
Non-anur (N=87) : UMOD 0.773 · β2M 0.802 · combiné 0.833 (0.746–0.920) opt-corr **0.820** ; OR refit 1.10 (p=0.023)/0.86 (p=0.003).
Calibration : HL p=0.53 · pente 0.96. Borderline KRU1–3 (N=33) : combiné 0.748.
Table 4 (N=148) : rule-out 76(51.4%) NPV93.4/corr92.3(87.7–96.8) 71/5 · grey 22(14.9%) · rule-in 50(33.8%) PPV82.0/corr80.9(75.0–90.0) 41/9 · 85.1% classés.
Table 5 (N=87) : rule-out 21(24.1%) NPV76.2(54.9–89.4) 16/5 · grey 17(19.5%) · rule-in 49(56.3%) PPV83.7(71.0–91.5) 41/8 · 80.5%.
Conservative (non-anur) : Sp90% 0.55/PPV83.7/8faux/grey19.5%/80.5% · Sp95% 0.69/92.3(79.7–97.3)/3/grey31.0%/69.0% · Sp97.5% 0.82/93.9(80.4–98.3)/2/grey37.9%/62.1%.
Subgroups AUC 0.906–0.944 ; ADPKD n=16, excl. 0.916 ; UMOD 11.5vs7.1 p=0.51, β2M 24.2vs27.3 p=0.56, KRU 0.83vs0.66 p=0.92.
Normalization : brut 0.921 / 35L 0.923 / BSA 0.920 ; ρ 0.99–1.00 ; ≤2/151 reclassés.
β2M×dialyse : corr +0.48/+0.45/+0.18/+0.26 ; HDF p=0.014 ; 3×29.9 vs 2×22.5 p<0.001 ; OR ajusté 0.79 (0.67–0.92) p=0.004.
Assay : LOD 2.0 ; 27 ≤2.0 (10 undetectable→0) ; CV>20% 19/151 (12.6%), médiane 2.6 vs 8.3 ; cohorte médiane UMOD 7.4 (IQR 2.6–16.3).
Frequency (N=148) : 2×/sem KRU<2 4(4.5%)/KRU≥2 45(76.3%) ; Fisher p<0.001.
Flowchart : 164→11(4/3/2/1/1)→153→2→151.
