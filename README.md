# Split-kl and PAC-Bayes-Split-kl Inequalities
This repository contains the code needed to replicate the experiments:
- Numerical studies of split-kl in Sec. 2.3 and Appendix D in the preprint [1].
- Empirical studies of PAC-Bayes-split-kl on linear classifiers in Sec. 4.1 in the preprint [1].

The details on numerical studies of split-kl are provided in `simulation/`.

## Experiment Environment
The implementation of PAC-Bayes-split-kl on linear classifiers is in R, and is tested in Windows.

## Basic usage
Before running the experiments, change the path in the script to the local repository.
Also, change the [data_option] to one of the options listed behind. The data sets are from [3]. Output files will be created in directory `out/`.

## Plots
After running the experiments, the plots of the paper can be generated by `cd`ing to the `plots` directory and running:`make plot`. Output figures will be generated in the subdirectories `figure`.

## Acknowledgements
Some of the implementations in `script`, `gen-data`, and `PBUB` are based on the implementation from [2].

## References
\[1\] Wu and Seldin: Split-kl and PAC-Bayes-split-kl Inequalities

\[2\] [Mhammedi, Grünwald and Guedj: PAC-Bayes Un-Expected Bernstein Inequality (NeurIPS 2019)](https://github.com/bguedj/PAC-Bayesian-Un-Expected-Bernstein-Inequality)

\[3\] [The UCI Repository](https://archive.ics.uci.edu/ml/index.php)
