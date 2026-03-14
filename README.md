# A Haplotype-Resolved Transcriptome-Wide Association Framework

**HapTWAS** provides a unified haplotype-resolved transcriptome-wide association framework, employing haplotypes to capture the concerted *cis*-regulatory effects of local genetic architectures for genetically regulated expression modeling, effectively integrating the LD structure and phase information that are often fragmented in single-variant analyses. Specifically, it is engineered to be computationally universal, supporting both individual-level genotypes and summary-level statistics through a novel SNP-to-haplotype transformation strategy. It can also extend to other SNP-incorporated strategies for conditional analyses, such as integrating biological priors to divide specific functional regions.

Our results establish **HapTWAS** as a powerful and scalable strategy, not only enhancing the predictive accuracy of gene expression but also achieving a refined resolution for fine-mapping causal genes. The framework is publicly available via an interactive website (https://twas.farmgtex.org/), providing a statistically rigorous and versatile toolset to resolve an even greater fraction of the missing heritability and fully bridge the gap from genetic variation to biological function. This platform democratizes the use of our advanced models, empowering the research community without specialized computational expertise to conveniently apply HapTWAS on their own data, thereby accelerating discovery and interpretation across complex traits and species.

## Schematic overview

<p align="center">
  <img src="https://github.com/chengliu-zju/HapTWAS/blob/main/Schematic_overview.png" width="80%" height="80%" alt="Schematic_overview">
</p>

**HapTWAS** is a unified statistical pipeline designed to enhance the resolution of TWAS by incorporating local haplotype structures. The workflow is divided into three main stages.

1. Haplotype construction: combined SNPs into haplotypes by the methods of linkage disequilibrium or sliding window.
  
2. Prediction models learning: haplotype matrix `X_haplo` and tissue-specific gene expression `G` were regressed with elastic-net model to obtain tissue-specific haplotype weight `w` and reference *cis*-haplotype covariance `L`.

3. Gene-trait associations using prediction models and three types of input:

   a. individual-haplotype, directly predict gene expression `G^` and analyze correlation with phenotype `Y`;
   
   b. summary-haplotype: calculate TWAS results `Z_g` through `Z_haplo` from haplotype GWAS results;
   
   c. summary-SNP: transform the summary statistics of SNP level `Z_SNP` into haplotype level `Z^_haplo` with the dosage design matrix `H` of `k` SNPs to constitute `2^k` haplotypes firstly, and then calculate gene-level statistics.

The result will output the predicted gene expression (only individual input) and gene-level summary statistics. Together, these derivations provide a coherent methodological foundation that enables haplotype-level association discovery, bridging the gap between conventional SNP-based TWAS and higher-order genomic architecture.

