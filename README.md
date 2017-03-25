# ELS-HMM — a method to detect ancient sweeps based on a signal of extended lineage sorting
ELS-HMM is a hidden Markov model that detects long regions in which an outgroup falls outside of the variation of a large number of individuals from another group. The model has been used to detect ancient selective sweeps on the human lineage after the split from archaics by using the Neandertal and Denisovan genome as outgroup genomes.

## Input Format
To configurate the emission probabilities that the outgroup individual share the derived allele at a site segregating in the other population, a tab-separated file is required specifying on each line the number of derived alleles in the population followed by the probability:

Number of derived alleles | Probability of the outgroup to be derived
------------------------- | -----------------------------------------
1 | 0.0168473
2 | 0.0188793
3 | 0.0225382
... | ...
n-1 | 0.825273

Probabilities of sharing at fixed sites are provided through the command line with the other parameters of the model.

The input file is also a tab-separated file with the following format:

Chromosome | Physical position | Information field | Total number of chromosomes in the population | Number of derived alleles | State of the outgroup individual: A for Ancestral and D for Derived | Genetic distance from the previous site
---------- | ----------------- | ----------------- | --------------------------------------------- | ------------------------- | ------------------------------------------------------------------- | ---------------------------------------
1 | 100000 | Neandertal | 370 | 100 | D | 0.0001
1 | 220000 | site2 | 370 | 20 | A | 0.0012

Note that this version allow you to run the model on only one chromosome at a time, do not combine data from different chromosomes into the same input file. 

## Examples of command line

For the 3-state model (with selection and parameter estimation):
```
hmm -e config-file -L 100000 -l 15000 -S 200000 -r 0.1 -F 0.99 -f 0.95 -N 0.0001 (-o log-file) input > output
```
The optional switch -o allows you to define the name of the log file which follows the steps of the hill-climbing algorithm, returning the log-likelihood and the associated values of the parameters. Otherwise, if you estimate parameters, this file will be named convergence.log.

For the 2-state model:
```
hmm -e config-file -L 100000 -l 15000 -F 0.99 -f 0.95 -B 0.0001 input > output
```
You can also run the HMM without any parameter estimation:
```
hmm -e config-file -L 100000 -l 15000 -S 200000 -r 0.1 -F 0.99 -f 0.95 input > output
```
```
hmm -e config-file -L 100000 -l 15000 -F 0.99 -f 0.95 input > output
```
