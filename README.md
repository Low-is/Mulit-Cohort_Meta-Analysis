# Pipeline Development - Multi Cohort Meta-analysis
## Study level effect size and pooled effect size estimations, inverse variance weighting, and Fisher p-value combinataion.


**Standard Deviation (SD)** - measures how spread out individual data points are around the mean.

**Standard Error** - measures how much the sample mean would vary if you repeated the experiment many times. Standard deviation of that distribution of sample means.

**Cohen's d** - used for comparing two groups, calculated as the differences between group means divided by the pooled standard deviaton.
  - Example:
      d = (M2 - M1)/ SDpooled; M1 and M2 are the mean of the first and second groups.
  - Hedges' g correction:
      cf = (1-3)/(4 * (n1+n2) - 9 )
      g = cf * (mean difference) / (pooled standard deviation)

**Pearson's r** - measures the strength of a linear relationship between two variables, with an effect size ranging from -1 to + 1.

**Pooled Standard Deviation** - a weighted average of standard deviation from two or more groups, used when assuming the groups have equal variances. It combines individual group varainces to provide a single, more precise estimate of the common variability across the populations.
  - Example:
      Pooled SD = sqrt((SD1^2 + SD2^2) / 2)

**Effect Size** - a statistical measure quantifying the magnitude of a relationship or difference between groups, helping to determine a finding's practical significance beyond statistical significance.
  - Example:
      Cohen's d, where a weight loss intervention group loses 10 kg, and the control group loses 5 kg, resulting in a Cohen's d approximately 0.67, indicating a medium effect size that is a noticeable difference in the weight loss between the groups.
    - Interpreting Cohen's d:
        - **0.2:** Small effect
        - **0.5:** Medium effect
        - **0.8:** Large effect

**Summary Effect Size** - A combined or aggregated effect size across multiple studies or measurements. Often used in meta-analysis to produce a single estimate that represents the overall effect. 

**Bioinformatics Meta-Analysis Toolkit for Gene Expression** 

This toolkit provides a robust framework for computing effect sizes, pooling meta-analytic estimates, and performing statistical analyses on gene expression datasets. It is designed to identify genes or gene sets with consistent differential expression across multiple studies, supporting biomarker discovery and predictive modeling.



## Key Features 🔑:
1. # Data Normalization and Batch Correction
  - Mitigates platform-specific variations using batch correction methods, median-centering, and rank-based gene transformations.
  - Ensures comparability of expression measures across studies and platforms.
2. # Gene-Level Effect Size Computation
  - Calculates study-level effects for individual genes using metrics such as Hedge's g and Cohen's d, quantifying the magnitude of expression differences between conditions.
  - Accounts for study-specfic variance using inverse-variance weigthing to manage heterogeneity across datasets. 
3. # Meta-Analysis Across Studies
  - Pools study-level effect sizes on a gene-wise basis.
  - Combines p-values across studies and adjusts them for multiple testing to control the false discovery rate (FDR).
  - Produces a list of "meta-genes" consistently associated with the condition or treatment of interest.
4. # Statistical Testing and Visualization
  - Performs t-tests while accommodating variances assumptions to determine statistical signficance. 
  - Generates comprehensive outputs, including:
      - Tables with effects sizes, summary effect sizes, p-values,and FDR-adjusted q-values.
      - Forest plots summarizing gene-level and pooled effect sizes.
5. # Predictive Modeling Integration
- Meta-genes can be used to build early diagnostic or prognostic models.
- Currently supports the ML algorithm, random forest, with customizable training control and optimized hyperparameter tuning grids.

## Applications 🖥️:
  - Identification of robust biomarkers across independent studies.
  - Construction of predictive models for early diagnosis or prognosis.
  - Insight into consistent molecular responses to conditions or treatments despite study-specific variability.
  - Combinging studies improves statistical power in selecting gene signatures.


---

## Below is an example of the output:

<!-- Forest plots side by side -->
<table>
  <tr>
    <td><img src="BOLA1.jpg" alt="BOLA1" width="250">></td>
    <td><img src="CYP4F3.jpg" alt="CYP4F3" width="250"></td>
    <td><img src="VEGFA.jpg" alt="VEGFA" width="250"></td>
  </tr>
</table>

---

<!-- Prediction models side by side below -->
<table>
  <tr>
    <td><img src="classification_models/DNA_internal.png" alt="Prediction Model 1" width="250"><br><center>Prediction Model 1</center></td>
    <td><img src="classification_models/DNA_internal_14.png" alt="Prediction Model 2" width="250"><br><center>Prediction Model 2</center></td>
    <td><img src="classification_models/internal_validation_gse69686.png" alt="Prediction Model 3" width="250"><br><center>Prediction Model 3</center></td>
  </tr>
</table>

---

## Future Directions 🚀:
  - Include mulitple feature selection methods, deafult setting is selecting shared genes across platforms before performing meta-analysis.
  - Additional machine learning algorithms optomized for diverse data types.
  - User-defined parameters will allow filtering by effect size thresholds and pooled effect size cutoffs.
  - Planned support for more data soures and processing methods, including ArrayExpress, GEO datasets with multiple platforms per study, and single-cell data.



# Citations: 
Zheng, H., Rao, A.M., Ganesan, A. et al. Multi-cohort analysis identifies a blood-based immune transcriptomic signature for early lung cancer detection. npj Precis. Onc. 9, 246 (2025).[https://doi.org/10.1038/s41698-025-01043-z](https://doi.org/10.1038/s41698-025-01043-z)

Walsh CJ, Batt J, Herridge MS, Mathur S, Bader GD, Hu P, Khatri P, Dos Santos CC. Comprehensive multi-cohort transcriptional meta-analysis of muscle diseases identifies a signature of disease severity. Sci Rep. 2022 Jul 4;12(1):11260.[https://www.nature.com/articles/s41598-022-15003-1](https://www.nature.com/articles/s41598-022-15003-1)

Rashid NU, Li Q, Yeh JJ, Ibrahim JG. Modeling Between-Study Heterogeneity for Improved Replicability in Gene Signature Selection and Clinical Prediction. J Am Stat Assoc. 2020;115(531):1125-1138.[https://pmc.ncbi.nlm.nih.gov/articles/PMC7528965/](https://pmc.ncbi.nlm.nih.gov/articles/PMC7528965/)

Pollard, Katherine S.; Dudoit, Sandrine; and van der Laan, Mark J., "Multiple Testing Procedures: R multtest Package and Applications to Genomics" (December 2004). U.C. Berkeley Division of Biostatistics Working Paper Series. Working Paper 164.[https://biostats.bepress.com/ucbbiostat/paper164](https://biostats.bepress.com/ucbbiostat/paper164)

Sweeney TE, Shidham A, Wong HR, Khatri P. A comprehensive time-course-based multicohort analysis of sepsis and sterile inflammation reveals a robust diagnostic gene set. Sci Transl Med. 2015 May 13;7(287):287ra71.[https://www.science.org/doi/10.1126/scitranslmed.aaa5993]([https://doi.org/10.1126/scitranslmed.aaa5993](https://www.science.org/doi/10.1126/scitranslmed.aaa5993))
