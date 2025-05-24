# Repeat-Aware_mutation_rate_estimator
This repository builds upon the ideas from [this work](https://github.com/medvedevgroup/Repeat-Aware_Substitution_Rate_Estimator/tree/main), with significant optimizations and extensions. It enables faster computation of mutation rates (genomic distances) between two whole genomes, leveraging KMC and other mature tools.



# Usage

## Requirement 

```
conda create -n <env_name> python=3.10 --file requirements.txt -c conda-forge -c bioconda -y
conda activate <env_name>
```

