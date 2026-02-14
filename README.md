# **topography_MDD**

This repository provides core code and relevant toolboxes for data analysis in the article entitled "Personalized Functional Topography-Based Brain Age Prediction Modeling Reveals Divergent Neurodevelopment Patterns in Major Depression ".

## Overview

Content includes standalone software, source code, and data. The project is structured into four parts corresponding to the major analyses in the article, including brain age prediction, brain, clinical and gene expression association analysis. Due to size limitation, the script relavent data can be found in [https://pan.bnu.edu.cn/l/A1NI3c](https://pan.bnu.edu.cn/l/91maR3).

## Toolboxes

SPM12, https://www.fil.ion.ucl.ac.uk/spm/software/spm12/

SeeCAT, a custom developed toolbox, available in /script_relavent_files/seecat.zip in https://pan.bnu.edu.cn/l/YFiIwE

ComBatHarmonization, ver. 20180620, https://github.com/Jfortin1/ComBatHarmonization

AHBAprocessing, ver. 20181025, https://github.com/BMHLab/AHBAprocessing

myPLS, https://github.com/MIPLabCH/myPLS

LIBSVM (3.25), https://www.csie.ntu.edu.tw/~cjlin/libsvmtools/

GAMLSS package v5.4-3, https://www.gamlss.com/

Scikit-learn, https://scikit-learn.org/stable/index.html

Scripts from  Cui et al.,  https://github.com/ZaixuCui/pncsinglefuncparcel_psychopathology

Scripts from Cui and Gong, https://github.com/ZaixuCui/Pattern_Regression_Clean

Brainnet Viewer, ver. 1.7, http://www.nitrc.org/projects/bnv/

Other software and web-based tools include:

BrainSMASH, https://brainsmash.readthedocs.io/en/latest/index.html

Metascape, http://metascape.org/gp/index.html/

We thank the authors and developers for providing these wonderful tools for data analysis.

## Predict_brain_age

This part was primarily carried out using our custom script for MATLAB with some functions fulfilled by  individualized topography generation (https://github.com/ZaixuCui/pncsinglefuncparcel_psychopathology) and ComBatHarmonization (https://github.com/Jfortin1/ComBatHarmonization). The analysis details include:

1. Generate individualized topography for Dataset 2 using NMF following Cui's code.
2. Correct the center effect of the topography by using ComBatHarmonization.
3. Calculate the whole brain age and system brain age based on SVM using topography as input feature.
4. Classify the patients into two groups based on the brain age gap (BAG).

## Analyse_subtypes_brain

1. Calculate the correlation between topography and BAG in two groups.
2. Characterize the BAG-related lifespan trajectories of participants.
3. Permut the trajectories across different groups.

## Analyse_subtypes_clinical

We evaluated the association between topography and clinical scores using PLS following cui and gong's codes (https://github.com/ZaixuCui/Pattern_Regression_Clean), and analyzed the contribution weight using myPLS toolbox (https://github.com/MIPLabCH/myPLS).

1. Evaluate the association between topography and HDRS-17 scores using 10-folds PLS in two groups. Significance of the PLS component was determined by permutation test.
2. Analyze the contribution weight of brain and clinical scores using PLS. The weights of PLS components were determined by the bootstrap method .
3. Characterized the change of topography and BAG after treatment.

## gene

In this analysis, we used the revised script from Xia et al. Molecular Psychitary 2022, ([mingruixia/MDD_ConnectomeGradient](https://github.com/mingruixia/MDD_ConnectomeGradient)). 

1. The Gene expression data from the Allen Institute for Brain Science was first preprocessed by using AHBAprocessing (https://github.com/BMHLab/AHBAprocessing), obtaining the gene expression profile for the Glasser-360 atlas.
2. Perform PLS analysis to examine the association between the BAG-related map and gene expression profiles.
3. Significance of the PLS component was determined by permutation test in which the spatial autocorrelations were corrected by generative modeling (gen_surrogate_map_for g1z.py)
4. The weights of PLS components were determined by the bootstrap method .
5. Both the descending order and ascending order of PLS weighted genes were submitted to  for enrichment analysis (http://metascape.org/gp/index.html/).
6. Calculation of transcription level differences for 4 neurodevelopmental gene sets between BAG-related and non-BAG-related  regions
