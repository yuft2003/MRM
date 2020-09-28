# MRM
A methylation imputation method using mixture of regression models (MRM) of radial basis functions (RBF)

## Instructions
DNA methylation is an important heritable epigenetic mark that plays a crucial role in transcriptional regulation and the pathogenesis of various human disorders. The commonly used DNA methylation measurement approaches, e.g., llumina Infinium HumanMethylat.ion-27 and -450 BeadChip arrays (27K and 450K arrays) and reduced representation bisulfite sequencing (RRBS), only cover a small proportion of the total CpG sites in the human genome, which considerably limited the scope of the genome-DNA methylation analysis.  We proposed a new computational strategy to impute the methylation value at the unmeasured CpG sites using the mixture of regression model (MRM) of radial basis functions (RBFs), integrating information of neighboring CpGs and the similarities in local methylation patterns across subjects and across multiple genomic regions. 

The in-house R script used for data analysis of our manuscript titled ‘A Novel Computational Strategy for DNA Methylation Imputation Using Mixture Regression Model (MRM)’ is shared in this repository.

## Example
Simulating data by MRM model

    N = 20           #----number of subject
    noise0.var = 0.2 #----variance of noise level
    pctMi = 0.2      #----missing rate
    I=50             #----number of CpGs per region
    M=100            #----number of regions
    K=4              # number of clusters of subjects
    N_mu = 10        #---- # of rfb centers 
    mu = seq(-1,1,length.out=N_mu) #---- rfb centers 
    gamma = -10
    noise0.var=0.2
    
    source(data_simulation.R)

Loading real WGBS data

    





Runing MRM method

    source ('YOUR_PATH/MRM_main.R')


Output
    
      > head(MRM_out)
    [[1]]
      POS         y          x subj regionID train     y_hat   y_hat_2   y_hat_3
    1:   1 0.6193507 -0.9822720    1        1     1 0.7165322 0.6930279 0.7056523
    2:   2 0.6266632 -0.9627019    1        1     0 0.7356034 0.6939539 0.7238293
    3:   3 0.8150562 -0.7539769    1        1     1 0.8896063 0.6747196 0.8696888
    4:   4 0.7947416 -0.7513913    1        1     1 0.8897192 0.6743983 0.8697851
    5:   5 1.0000000 -0.6866021    1        1     1 0.8692463 0.6675884 0.8500715

y_hat: imputed value by region model 
y_hat_2： imputed value by subject model
y_hat_3: imputed value by stack model

    









## Citation
TBD

