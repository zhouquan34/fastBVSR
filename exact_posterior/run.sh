#! /bin/sh

Rscript exact_posterior.R >truth.txt
../fastBVSR-mac -g toy.txt -p toy.ph -w 10000 -s 200000 --h2-min 0.01 --h2-max 0.99 -o out

awk '{print $1"\t"$3}' out.beta.txt >mcmc.txt
echo "" >>mcmc.txt
awk 'NR > 1 {print $2}' out.path.txt >tmp_path
awk 'NR > 1 {print $1}' out.model.txt >tmp_model
paste tmp_path tmp_model >tmp
Rscript parse.R tmp  | sort -k1n >>mcmc.txt
rm tmp*
rm out*

