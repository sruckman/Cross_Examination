# Cross Examination

An interactive R Shiny app for teaching genetic mapping. Students explore three approaches to detecting genetic variants: multi-parent populations (MPP/BSA), biparental F2 crosses, and GWAS.

---

## Run it in your browser

**https://sruckman.shinyapps.io/cross_examination/**

---

## Run it locally in R

```r
install.packages(c("shiny", "ggplot2", "gridExtra"))
shiny::runGitHub("Cross_Examination", "sruckman")
```

---

## The three tabs

### MPP / BSA
Simulates a multi-parent population (MAGIC or hub-and-spoke design). Students choose the number of founders, generations of recombination, and pool size. The app selects the top-phenotype individuals as cases and random individuals as controls, then runs a bulk segregant analysis (BSA) scan using a G-test (1 replicate) or Cochran-Mantel-Haenszel test (multiple replicates). Plots show founder haplotype mosaics, cases vs. controls, the LOD scan, and per-founder allele frequency differences between pools.

### Biparental QTL
Simulates a classic two-parent cross: P1 x P2 -> F1 -> F2. Each F2 chromosome is a mosaic of P1 (red) and P2 (blue) segments created by recombination. Students choose the number of F2 individuals and the QTL architecture, then scan using additive regression (F-test) at each of 202 positions across a 100 cM chromosome.

### GWAS
Simulates a case-control GWAS in an outbred population using a MAGIC-style mosaic model with 12 founder haplotypes. Students control sample size (up to 100k) and average LD block size. Larger blocks represent bottlenecked or domesticated species. Smaller blocks represent diverse outbred populations with long recombination history. The app runs a chi-square association scan at 1000 SNPs (10 SNPs/cM) and displays a haplotype mosaic of cases vs. controls alongside a Manhattan plot.

---

## Key features

- **Teaching mode**: produces identical results every run (fixed seed + fixed QTL positions), so instructors can use consistent scenarios for lectures. Turn it off for student exploration and homework.
- **Save plot**: exports the current haplotype and scan plots as a PNG with all parameter settings in the caption
- **LD block size control**: GWAS tab lets you vary block size to show how recombination history shapes the association signal
- **Black background**: optimized for projector display
- **Fixed Bonferroni thresholds**: LOD 3.6 for MPP and biparental (202 tests, α = 0.05), LOD 4.3 for GWAS (1000 tests, α = 0.05)
