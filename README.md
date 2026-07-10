# Cross Examination

An interactive R Shiny app for teaching genetic mapping. Students explore three approaches to detecting genetic variants: multi-parent populations (MPP), biparental F2 crosses, and GWAS.



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

## Teaching materials

**Using the Simulator** ([docx](Using%20the%20Simulator.docx) | [PDF](Using%20the%20Simulator.pdf)) — a visual guide to every tab and control, with annotated screenshots. Start here if you are new to the app.

**Student Worksheet** ([docx](Cross_Examination_Student_Worksheet.docx) | [PDF](Cross_Examination_Student_Worksheet.pdf)) — six scenario prompts with workspace for students to record their reasoning and simulator results. Each scenario targets a different mapping design and set of trade-offs.

**[Teaching Demo](teaching_demo.html)** ([docx](Cross_Examination_Teaching_Demo_Drosophila.docx) | [PDF](Cross_Examination_Teaching_Demo_Drosophila.pdf)) — a fully worked example using a *Drosophila* aggression scenario. Covers parameter choice, budget math, and interpretation of results under single-QTL and polygenic architectures. Includes instructor discussion notes.

---

## Using the simulator

Cross Examination has three tabs: **MPP**, **Biparental QTL**, and **GWAS**. Each tab simulates a different study design for mapping the genetic basis of a quantitative trait. The controls on the left set the experimental parameters; the plots on the right show the output after you click Simulate.

Turn **Teaching Mode** ON before running to fix the true QTL positions and random seed, so every student sees the same result.

### MPP tab controls

The MPP tab simulates a multiparent population mapped with bulk segregant analysis (BSA). Cases are the top N individuals by phenotype; controls are a random N drawn from the non-extreme individuals. The CMH test combines signal across replicates when more than one replicate is run. Key parameters: crossing design (fully intercrossed or hub-and-spoke), number of founders, generations of recombination, pool size, and number of replicates.

### MPP tab output

Three panels appear after clicking Simulate. The top panel shows the RIL haplotype mosaic for cases and controls. The middle panel shows the LOD scan with the Bonferroni threshold and true QTL position(s) marked. The bottom panel shows per-founder allele frequency differences between pools, which identifies which ancestral line carries the effect allele.

### Biparental QTL tab

Two inbred parents are crossed to make F1, then intercrossed to make F2 offspring. Each F2 genome is a mosaic of the two parental haplotypes shown in red and blue. A regression test at each position detects linkage between genotype and phenotype. Choose the number of F2 individuals and the QTL architecture, then click Simulate.

### GWAS tab

Simulates a genome-wide association study in an outbred population. Individual genotypes are tested for allele frequency differences between cases and controls at 1,000 SNPs. Control average LD block size to simulate populations with different recombination histories — larger blocks represent bottlenecked or domesticated species, smaller blocks represent diverse outbred populations.

See **[Using the Simulator](Using%20the%20Simulator.pdf)** for annotated screenshots of each tab.

---

## How to use the app

**Teaching Mode** is a blue checkbox at the top of each tab's sidebar. When on, every simulation run produces identical results (fixed random seed and fixed QTL positions). Use this for lectures and guided exercises. Turn it off for student exploration and homework, where results change each run.

**Running a simulation**: set parameters in the sidebar and click Simulate. GWAS with large sample sizes may take 10-30 seconds.

**Saving plots**: click Save Plot after running. A PNG downloads with a timestamp in the filename. The caption at the bottom of the image shows all parameter settings. Save multiple runs to compare across conditions.

**Significance thresholds**: a blue dashed line marks the Bonferroni threshold on every LOD plot. LOD 3.6 for MPP and Biparental QTL (202 tests, α = 0.05). LOD 4.3 for GWAS (1000 tests, α = 0.05). White vertical lines mark the true QTL position.

---

## The three tabs

### MPP
Simulates a multi-parent population (fully intercrossed or hub-and-spoke design). Students choose the number of founders, generations of recombination, and pool size. The app selects the top-phenotype individuals as cases and random individuals as controls, then runs a bulk segregant analysis scan using a G-test (1 replicate) or Cochran-Mantel-Haenszel test (multiple replicates). Plots show founder haplotype mosaics, cases vs. controls, the LOD scan, and per-founder allele frequency differences between pools.

### Biparental QTL
Simulates a classic two-parent cross: P1 x P2 → F1 → F2. Each F2 chromosome is a mosaic of P1 (red) and P2 (blue) segments created by recombination. Students choose the number of F2 individuals and the QTL architecture, then scan using additive regression (F-test) at each of 202 positions across a 100 cM chromosome.

### GWAS
Simulates a case-control GWAS in an outbred population using a MAGIC-style mosaic model with 12 founder haplotypes. Students control sample size (up to 100k) and average LD block size. Larger blocks represent bottlenecked or domesticated species. Smaller blocks represent diverse outbred populations with long recombination history. The app runs a chi-square association scan at 1000 SNPs (10 SNPs/cM) and displays a haplotype mosaic of cases vs. controls alongside a Manhattan plot.

---

## Key features

- **Teaching mode**: produces identical results every run (fixed seed + fixed QTL positions), so instructors can use consistent scenarios for lectures. Turn it off for student exploration and homework.
- **Save plot**: exports the current haplotype and scan plots as a PNG with all parameter settings in the caption
- **Light and dark mode**: dark background optimized for projector display; light mode available for print and screen use
- **LD block size control**: GWAS tab lets you vary block size to show how recombination history shapes the association signal
- **Fixed Bonferroni thresholds**: LOD 3.6 for MPP and Biparental QTL (202 tests, α = 0.05), LOD 4.3 for GWAS (1000 tests, α = 0.05)
