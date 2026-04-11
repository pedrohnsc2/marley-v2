# Thermodynamic Landscape — Extended Analysis Report

**Generated:** 2026-04-10 23:51 UTC
**Runtime:** 0.05 seconds

---

## 1. Full Mutational Fitness Landscape

- **Matrix dimensions:** 25 positions x 4 bases = 100 variants
- **Wildtype dG:** -27.97 kcal/mol
- **Wildtype Tm:** 68.48 C
- **Global minimum dG:** -29.6200 kcal/mol (position 13, base G)
- **Wildtype is global minimum:** NO
- **Fitness gap (1st - 2nd):** 0.0500 kcal/mol

### Most Sensitive Positions (highest dG range)

| Position | WT Base | dG Range (kcal/mol) | dG Min | dG Max |
|----------|---------|---------------------|--------|--------|
|       23 |       C |              1.7900 | -27.97 | -26.18 |
|       13 |       A |              1.6500 | -29.62 | -27.97 |
|        5 |       A |              1.6400 | -29.33 | -27.69 |
|       21 |       A |              1.6000 | -29.57 | -27.97 |
|        6 |       A |              1.4500 | -28.86 | -27.41 |

## 2. Double Mutant Analysis & Epistasis

- **Total double mutants evaluated:** 2700
- **Position pairs tested:** 300
- **Better than wildtype:** 1769

### Epistasis Statistics

- **Mean epistasis:** 0.0052 kcal/mol
- **Std epistasis:** 0.1022 kcal/mol
- **Range:** [-0.7400, 0.6000] kcal/mol
- **Synergistic pairs:** 14 (0.5%)
- **Antagonistic pairs:** 21
- **Additive pairs:** 2665

**Best double mutant:** A13G+A21C (dG = -31.2200 kcal/mol, epistasis = 0.0000 kcal/mol)

## 3. Length Optimization (18-30 nt)

- **Lengths tested:** 18-30 nt
- **Current ASO length:** 25 nt
- **Optimal length:** 30 nt (dG = -33.3400 kcal/mol)
- **Current length is optimal:** NO

### Best dG by Length

| Length | Best dG (kcal/mol) | Tm (C) | GC | Off-target Risk | Current? |
|--------|-------------------|--------|------|-----------------|----------|
|     18 |          -19.6300 |  60.13 | 0.3333 | 3.49e-02 |       no |
|     19 |          -20.6400 |  61.30 | 0.3158 | 8.73e-03 |       no |
|     20 |          -21.6600 |  62.36 | 0.3000 | 2.18e-03 |       no |
|     21 |          -22.4600 |  62.73 | 0.2857 | 5.46e-04 |       no |
|     22 |          -23.7700 |  64.23 | 0.3182 | 1.36e-04 |       no |
|     23 |          -25.0700 |  65.60 | 0.3478 | 3.41e-05 |       no |
|     24 |          -26.5200 |  67.10 | 0.3333 | 8.53e-06 |       no |
|     25 |          -27.9700 |  68.48 | 0.3200 | 2.13e-06 |  YES <-- |
|     26 |          -28.9900 |  69.02 | 0.3077 | 5.30e-07 |       no |
|     27 |          -30.0200 |  69.59 | 0.3333 | 1.30e-07 |       no |
|     28 |          -31.3100 |  70.55 | 0.3214 | 3.00e-08 |       no |
|     29 |          -32.3200 |  70.96 | 0.3103 | 1.00e-08 |       no |
|     30 |          -33.3400 |  71.35 | 0.3000 | 0.00e+00 |       no |

## 4. Positional Sliding Window (25-mer)

- **SL RNA length:** 39 nt
- **Window length:** 25 nt
- **Total windows:** 15

- **Optimal window:** position 5-30 (dG = -27.9700 kcal/mol)
- **Current ASO covers optimal:** YES
- **dG deviation from optimal:** 0.0000 kcal/mol

### All Windows

| Start | End | dG (kcal/mol) | Tm (C) | GC | Current? |
|-------|-----|---------------|--------|------|----------|
|     0 |  25 |      -26.8100 |  66.49 | 0.2800 |       no |
|     1 |  26 |      -26.8100 |  66.49 | 0.2800 |       no |
|     2 |  27 |      -26.6700 |  66.25 | 0.3200 |       no |
|     3 |  28 |      -26.6700 |  66.25 | 0.2800 |       no |
|     4 |  29 |      -27.5300 |  67.72 | 0.3200 |       no |
|     5 |  30 |      -27.9700 |  68.48 | 0.3200 |  YES <-- |
|     6 |  31 |      -27.1100 |  67.02 | 0.3200 |       no |
|     7 |  32 |      -26.4000 |  66.00 | 0.3200 |       no |
|     8 |  33 |      -25.4600 |  64.48 | 0.2800 |       no |
|     9 |  34 |      -25.1800 |  63.93 | 0.2400 |       no |
|    10 |  35 |      -25.1800 |  63.93 | 0.2400 |       no |
|    11 |  36 |      -25.1800 |  63.93 | 0.2400 |       no |
|    12 |  37 |      -25.1800 |  63.93 | 0.2400 |       no |
|    13 |  38 |      -25.1800 |  63.93 | 0.2400 |       no |
|    14 |  39 |      -26.0500 |  65.41 | 0.2800 |       no |

## 5. Summary Statistics

- **Rank among single mutants:** 45/76
- **Binding energy percentile:** 41.33%
- **Variants within 2.0 kcal/mol of optimum:** 57
- **Robustness score (mean |ddG|):** 0.8524 +/- 0.4877 kcal/mol
- **Single mutants better than WT:** 44
- **Double mutants better than WT:** 1769

## Conclusion

MRL-ASO-001 is **NOT** the global optimum. Findings: 44 single mutant(s) with better dG; 1769 double mutant(s) with better dG; optimal length is 30 nt, not 25 nt.

---
*Analysis completed in 0.05 seconds.*
