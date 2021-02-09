# PitfallsWPPA
PitfallsWPPA corrects the pitfall in GWAS when using WPPA to investigate the genetic association. All the code for the simulation study is shown here. The real dataset is available at http://www.ricediversity.org (Oryza sativa data).
## BayesC.jl
Details of conducting GWAS using WPPA are shown in BayesC.jl. All analyses were performed using the JWAS package, an open-source package for whole-genome analyses based on Bayesian regression methods. The modified GWAS function (using window specific threshold values) is in GWAS_fun.jl.
## Freq.R
A frequentist method of conduction GWAS based on genomic windows, EMMAX, was conducted in R. All the calculations were based on Chen et al. (2017. Genome-wide association analyses based on broadly different specifications for prior distributions, genomic windows, and estimation methods).
