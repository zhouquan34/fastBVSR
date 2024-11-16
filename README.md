"examples" gives some example commands for running fastBVSR in terminal. 

"test_data" contains example files for the input. The design matrix can be provided as a matrix file (see test.mat), a mean-genotype file (see test.mg) or a plink file (see test.plink). In the previous version, there was a bug that assumed all entries of the design matrix were integers for mean-genotype or plink inputs. 

If you have C++ GSL library installed, you can try compiling the code yourself. An example Makefile is provided. Some may find it convenient to run the program in R. An example is given in the folder R-mac. It is only for Mac OS. 

Thanks for your interest. My personal email: zhouquan.stat@gmail.com

**Reference:** Quan Zhou and Yongtao Guan. Fast model-fitting of Bayesian variable selection regression using the iterative complex factorization. Bayesian Analysis, 2018.  

**Argument List**

| Option             | Value Type          | Default Value    | Description                                |
|--------------------|---------------------|------------------|--------------------------------------------|
| `-g`     | `<filename>`        |  | design matrix in mean genotype format (same format as in piMASS)   |
| `-a`     | `<filename>`        |  | design matrix in PLINK 1.9 A-transpose format |
| `-m`     | `<filename>`        |  | design matrix in matrix format | 
| `-p`     | `<filename>`        |  | response vector (i.e., phenotype; same format as in piMASS) | 
| `-o`     | `<string>`          |  | prefix for output files | 
| `-w`     | `<integer>`         | 0    | number of burn-in iterations  | 
| `-s`     | `<integer>`         | 1000 | number of sampling iterations  | 
| `-R`     | `<integer>`         | 1000 | do Rao-Blackwellization once every R iterations | 
| `-r`     | `<integer>`         | 0    | random seed for GSL | 
| `--max-size` | `<integer>`  | 4000 | maximum number of predictors selected; an alias for `-smax`  |
| `--min-size` | `<integer>`  | 1    | minimum number of predictors selected (must be at least 1); an alias for `-smin` |
| `--pi-alpha` | `<double>`   | 0    | alpha parameter in the (truncated) beta prior of pi |
| `--pi-beta`  | `<double>`   | 1    | beta parameter in the (truncated) beta prior of pi |
| `--pi-max`   | `<double>`   | 1    | maximum value of pi |
| `--pi-min`   | `<double>`   | 0    | minimum value of pi |
| `--h2-max`   | `<double>`   | 0.999   | maximum value of heritability; an alias for `-hmax` |
| `--h2-min`   | `<double>`   | 1e-6    | minimum value of heritability; an alias for `-hin` |
| `--lunif-h2` |              | OFF  | use a uniform prior on log-heritability  | 
| `--start`    | `<integer>`  | 1    | initial model size (must be at least 1); an alias for `-nstart` | 
| `--add-unif` | `<double>`      | 0.3  | probability of using a uniform proposal when adding a predictor | 
| `--long` | `<double>`          | 0.1  | probability of a long-range proposal | 
| `--jump` | `<integer>`         | 5    | max size (i.e., number of single flips) of a long-range proposal | 
| `--exact`|                     | OFF  | exact calculation of marginal likelihood without using ICF and exchange algorithm | 
| `--prec` | `<integer>`         | 5000 | number of predictors used in the gram matrix pre-computation     |



