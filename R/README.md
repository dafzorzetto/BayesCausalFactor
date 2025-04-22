# Code files for BayesCausalFactor

### General functions (in **`src`** folder):
- **`BayesCausalFactor_Gibbs.R`**:

   Gibbs sampler for our proposed model: Causal Bayesian Regression Factor Model
- **`BART_BCF.R`**:
   
    code for competitor: Bart (Hill 2011) and BCF (Hahn at al 2020), with indipendence between components of the outcome
- **`StandardFactor_Gibbs.R`**:
   
    code for competitor: causal factor model with standard prior for factor scores
- **`code_parallelized.R`**:
   
    function to parallelize the code for replicated (useful for simulation study)
- **`simulation_functions.R`**:
   
    functions to simulate part of the generating process
- **`plots_functions.R`**:
   
    functions for some visualizations
    
### Code to reproduce the simulation study (in **`simulation_study`** folder):
 - **`simulation_study.Rmd`**:
   
    RMarkdown file to reproduce all the simulation study
