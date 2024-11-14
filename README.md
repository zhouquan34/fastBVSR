"examples" gives some example commands for running fastBVSR in terminal. 

"test_data" contains example files for the input. The design matrix can be provided as a matrix file (see test.mat), a mean-genotype file (see test.mg) or a plink file (see test.plink). In the previous version, there was a bug that assumed all entries of the design matrix were integers for mean-genotype or plink inputs. 

If you have C++ GSL library installed, you can try compiling the code yourself. An example Makefile is provided. 

Some may find it convenient to run the program in R. An example is given in the folder R-mac. It is only for Mac OS. 

Thanks for your interest. My personal email: zhouquan.stat@gmail.com

**Reference:** Quan Zhou and Yongtao Guan. Fast model-fitting of Bayesian variable selection regression using the iterative complex factorization. Bayesian Analysis, 2018.  

