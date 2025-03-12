# MetaDICT-manuscript-code
This is the repository archiving code and data for Microbiome Data Integration via Shared Dictionary Learning.

# File Introduction
The "data" folder contains datasets we use in the manuscript. 

The "code" folder contains the RScript and Rmarkdown files for the algorithm and results in the manuscript. 
The file *preprocessing.Rmd* contains code for processing He’s data, which serves as the foundation for generating simulated data.
The file function.R includes the code for MetaDICT. 

We evaluated the performance of eight methods in the manuscript. DEBIAS-M and scANVI were implemented in Python, while the remaining methods were implemented in R. For the simulation experiments, first, generate the data and save the simulated count and metadata tables. Next, run DEBIAS-M and scANVI in Python (examples for executing these methods can be found in *DEBIAS-M_implementation.ipynb* and *scANVI_implementation.ipynb*). Finally, perform downstream analysis using the R Markdown file.

For experiments implemented in the R script, execute the script 500 times using different random seeds.

# Reference
Bo Yuan, Shulei Wang,
<b>Microbiome Data Integration via Shared Dictionary Learning</b>
(2024).
[<a href=https://www.biorxiv.org/content/10.1101/2024.10.04.616752v1>link</a>]

# Session info
R version 4.4.2 (2024-10-31)
Platform: aarch64-apple-darwin20
Running under: macOS Monterey 12.5

Matrix products: default
BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

time zone: America/Chicago
tzcode source: internal

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] TreeSummarizedExperiment_2.14.0 Biostrings_2.74.1               XVector_0.46.0                  SingleCellExperiment_1.28.1     SummarizedExperiment_1.36.0    
 [6] Biobase_2.66.0                  GenomicRanges_1.58.0            GenomeInfoDb_1.42.1             IRanges_2.40.1                  MatrixGenerics_1.18.0          
[11] matrixStats_1.5.0               S4Vectors_0.44.0                BiocGenerics_0.52.0             ANCOMBC_2.8.0                   Maaslin2_1.15.1                
[16] RDB_0.0.1                       mclust_6.1.1                    ggforce_0.4.2                   tidygraph_1.3.1                 igraph_2.1.2                   
[21] vegan_2.6-10                    permute_0.9-7                   ecodist_2.1.3                   PLSDAbatch_1.2.0                sva_3.54.0                     
[26] BiocParallel_1.40.0             genefilter_1.88.0               mgcv_1.9-1                      nlme_3.1-166                    MMUPHin_1.17.0                 
[31] doParallel_1.0.17               iterators_1.0.14                foreach_1.5.2                   ConQuR_2.0                      pROC_1.18.5                    
[36] caret_7.0-1                     lattice_0.22-6                  randomForest_4.7-1.2            reshape2_1.4.4                  MicrobiomeStat_1.2             
[41] MBESS_4.9.3                     effsize_0.8.1                   viridis_0.6.5                   viridisLite_0.4.2               scales_1.3.0                   
[46] stringr_1.5.1                   ggraph_2.2.1                    dplyr_1.1.4                     ggpubr_0.6.0                    ggplot2_3.5.1                  

loaded via a namespace (and not attached):
  [1] gld_2.6.6                rARPACK_0.11-0           nnet_7.3-19              TH.data_1.1-2            vctrs_0.6.5              energy_1.7-12           
  [7] proxy_0.4-27             digest_0.6.37            png_0.1-8                corpcor_1.6.10           shape_1.4.6.1            Exact_3.3               
 [13] pcaPP_2.0-5              ggrepel_0.9.6            deldir_2.0-4             parallelly_1.41.0        MASS_7.3-61              withr_3.0.2             
 [19] xfun_0.49                doRNG_1.8.6              survival_3.8-3           memoise_2.0.1            MatrixModels_0.5-3       gmp_0.7-5               
 [25] zoo_1.8-12               tidytree_0.4.6           gtools_3.9.5             DEoptimR_1.1-3-1         Formula_1.2-5            ellipse_0.5.0           
 [31] KEGGREST_1.46.0          httr_1.4.7               rstatix_0.7.2            globals_0.16.3           rhdf5filters_1.18.0      kmer_1.1.2              
 [37] rhdf5_2.50.1             rstudioapi_0.17.1        UCSC.utils_1.2.0         generics_0.1.3           base64enc_0.1-3          zlibbioc_1.52.0         
 [43] polyclip_1.10-7          statip_0.2.3             GenomeInfoDbData_1.2.13  SparseArray_1.6.0        xtable_1.8-4             ade4_1.7-22             
 [49] evaluate_1.0.1           S4Arrays_1.6.0           hms_1.1.3                glmnet_4.1-8             colorspace_2.1-1         getopt_1.20.4           
 [55] ROCR_1.0-11              readxl_1.4.3             magrittr_2.0.3           readr_2.1.5              future.apply_1.11.3      robustbase_0.99-4-1     
 [61] SparseM_1.84-2           DECIPHER_3.2.0           XML_3.99-0.17            Hmisc_5.2-1              class_7.3-22             pillar_1.10.1           
 [67] performance_0.12.4       pwalign_1.2.0            caTools_1.18.3           compiler_4.4.2           RSpectra_0.16-2          stringi_1.8.4           
 [73] biomformat_1.34.0        DescTools_0.99.58        gower_1.0.2              stabledist_0.7-2         minqa_1.2.8              lubridate_1.9.4         
 [79] GenomicAlignments_1.42.0 plyr_1.8.9               crayon_1.5.3             abind_1.4-8              timeSeries_4041.111      mixOmics_6.30.0         
 [85] cqrReg_1.2.1             haven_2.5.4              locfit_1.5-9.10          graphlayouts_1.2.1       bit_4.5.0.1              sandwich_3.1-1          
 [91] rootSolve_1.8.2.4        biglm_0.9-3              multcomp_1.4-26          codetools_0.2-20         recipes_1.1.0            e1071_1.7-16            
 [97] lmom_3.2                 phyloseq_1.50.0          multtest_2.62.0          splines_4.4.2            Rcpp_1.0.13-1            fastDummies_1.7.4       
[103] quantreg_5.99.1          cellranger_1.1.0         interp_1.1-6             ATE_0.4.0                knitr_1.49               blob_1.2.4              
[109] clue_0.3-66              lme4_1.1-35.5            fBasics_4041.97          fs_1.6.5                 checkmate_2.3.2          listenv_0.9.1           
[115] Rdpack_2.6.2             expm_1.0-0               gsl_2.1-8                ggsignif_0.6.4           tibble_3.2.1             Matrix_1.7-1            
[121] statmod_1.5.0            tzdb_0.4.0               tweenr_2.0.3             pkgconfig_2.0.3          pheatmap_1.0.12          tools_4.4.2             
[127] cachem_1.1.0             rbibutils_2.3            RSQLite_2.3.9            DBI_1.2.3                numDeriv_2016.8-1.1      phylogram_2.1.0         
[133] rmutil_1.1.10            fastmap_1.2.0            rmarkdown_2.29           grid_4.4.2               Rsamtools_2.22.0         broom_1.0.7             
[139] stable_1.1.6             BiocManager_1.30.25      insight_1.0.0            carData_3.0-5            rpart_4.1.23             farver_2.1.2            
[145] yaml_2.3.10              foreign_0.8-87           latticeExtra_0.6-30      spatial_7.3-17           bayesm_3.1-6             cli_3.6.3               
[151] purrr_1.0.2              lifecycle_1.0.4          mvtnorm_1.3-2            lava_1.8.0               backports_1.5.0          modeest_2.4.0           
[157] annotate_1.84.0          timechange_0.3.0         gtable_0.3.6             ape_5.8-1                limma_3.62.1             CVXR_1.0-15             
[163] jsonlite_1.8.9           edgeR_4.4.1              bitops_1.0-9             bit64_4.5.2              Rtsne_0.17               yulab.utils_0.1.8       
[169] RcppParallel_5.1.9       timeDate_4041.110        lazyeval_0.2.2           htmltools_0.5.8.1        glue_1.8.0               optparse_1.7.5          
[175] treeio_1.30.0            jpeg_0.1-10              gridExtra_2.3            boot_1.3-31              R6_2.5.1                 tidyr_1.3.1             
[181] gplots_3.2.0             forcats_1.0.0            Rmpfr_1.0-0              rngtools_1.5.2           cluster_2.1.8            Rhdf5lib_1.28.0         
[187] ipred_0.9-15             nloptr_2.1.1             compositions_2.0-8       DelayedArray_0.32.0      tidyselect_1.2.1         htmlTable_2.4.3         
[193] microbiome_1.28.0        tensorA_0.36.2.1         inline_0.3.20            car_3.1-3                AnnotationDbi_1.68.0     future_1.34.0           
[199] GUniFrac_1.8             ModelMetrics_1.2.2.2     munsell_0.5.1            KernSmooth_2.23-24       BiocStyle_2.34.0         data.table_1.16.4       
[205] htmlwidgets_1.6.4        RColorBrewer_1.1-3       hwriter_1.3.2.1          rlang_1.1.4              lmerTest_3.1-3           ShortRead_1.64.0        
[211] dada2_1.34.0             hardhat_1.4.0            prodlim_2024.06.25   
