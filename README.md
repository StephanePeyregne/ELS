# ELS-hmm — a method to detect ancient sweeps based on a signal of extended lineage sorting
ELS-hmm is a hidden Markov model that detects long regions in which an outgroup falls outside of the variation of a large number of individuals from another group. The model has been used to detect ancient selective sweeps on the human lineage after the split from archaics by using the Neandertal and Denisovan genomes as outgroup genomes. ELS stands for Extended Lineage Sorting.

## Installation
To compile the software, we provide a Makefile that you can run from the command line by typing:
```make```
The software uses [nlopt](http://ab-initio.mit.edu/wiki/index.php/NLopt) for maximizing the likelihood and a development version needs to be installed to compile ELS-hmm. The software has been tested on Linux.

## Input Format
To configurate the probabilities (emissions) that the outgroup individual share the derived allele at a polymorphic site in the other population, a tab-separated file is required specifying on each line the number of derived alleles in the population followed by the probability:

Number of derived alleles | Probability of the outgroup to be derived
------------------------- | -----------------------------------------
1 | 0.0168473
2 | 0.0188793
3 | 0.0225382
... | ...
n-1 | 0.825273

Probabilities of sharing at fixed sites are provided through the command line with the other parameters of the model.

The input file is also a tab-separated file with the following format:

Chromosome | Physical position | Information field | Total number of chromosomes in the population | Number of derived alleles | State of the outgroup individual: A for Ancestral and D for Derived | distance from the previous site (genetic or physical distance)
---------- | ----------------- | ----------------- | --------------------------------------------- | ------------------------- | ------------------------------------------------------------------- | ---------------------------------------
1 | 10000 | Neandertal | 370 | 100 | D | 10000
1 | 10056 | site2 | 370 | 20 | A | 56
1 | 10100 | site3 | 370 | 55 | A | 44

Note that this version allows you to run the model on only one chromosome at a time, do not combine data from different chromosomes into the same input file.

## List of parameters

* Required parameters:
```

-e [path to the configuration file]

-L [average length of internal regions]

-l [average length of external regions]
```

Lengths should have the same unit as the distances provided in the last column of the input file. If you use the genetic distance and try to estimate parameters, we did not test the estimation if those parameters were small floating points: please use a unit of genetic distance that gives orders of magnitude similar to the physical distance (we used nM for humans).

```
-F [probability of the outgroup being Derived in an internal region at a fixed derived site]
```
In theory, for a region where the outgroup is internal to the population variation, this should be 1 as all chromosomes are derived at this site (if the outgroup was ancestral it would mean it is external to the population variation...) but, in practice, sequencing errors and back mutations could make the outgroup ancestral again so this probability is lower than 1.
```
-f [probability of the outgroup being derived in an external region at a fixed derived site]
```

* Optional parameters:
```
-S [average length of ELS regions]

-r [proportion of external regions that could represent events of ELS]
```
Those two parameters can be used if you want to model explicitely the second type of external regions that are very long and could represent events of selection. They have to be used together if you do so. ALso, both parameters influence each other and, when estimating those parameters, one should not interpret them too much. Under neutral simulations, it happens that the estimated r is not different from 0 if S is similar to the length of neutral external regions.

An advantage of this 3-state model is to be able to use a likelihood ratio test to compare to results obtained with a 2-state model. One can get this likelihood from the log file when estimating parameters (the log likelihood is in the last row, first column of this file).


* To estimate parameters:
```
-N followed by the log-likelihood maxima difference under which the estimation converged (this uses the COBYLA algorithm from the nlopt library). This parameter can, for now, only be used for the 3-state model.

-o [file name] : to provide a name to the log file allowing to follow the progression of the parameter estimation. Otherwise, this file is called convergence.log

-B : similar to -N but makes use of a very naive random search algorithm (can be used for both the 3-state and 2-state models)
```

## Examples of command line

For the 3-state model (with selection and parameter estimation):
```
hmm -e config-file -L 10000 -l 1500 -S 20000 -r 0.1 -F 0.99 -f 0.95 -N 0.0001 (-o log-file) input > output
```
The optional switch -o allows you to define the name of the log file which follows the steps of the hill-climbing algorithm, returning the log-likelihood and the associated values of the parameters. Otherwise, if you estimate parameters, this file will be named convergence.log.

For the 2-state model:
```
hmm -e config-file -L 10000 -l 1500 -F 0.99 -f 0.95 -B 0.0001 input > output
```
You can also run the HMM without any parameter estimation:
```
hmm -e config-file -L 10000 -l 1500 -S 20000 -r 0.1 -F 0.99 -f 0.95 input > output
```
```
hmm -e config-file -L 10000 -l 1500 -F 0.99 -f 0.95 input > output
```
### Output format

The software returns the input file with 2 or 3 additional columns (depending on the model used). The first additional column takes values of 0,1 and 2 if the probability of the external state or ELS state is above 0.8 at this site (represented as 1 and 2 respectively, 0 otherwise). The two following columns report the probability of the external state and ELS state at this site. The probability of the internal state can be calculated by substracting those probabilities to 1.
* Format of the log file:

Each line correspond to a new set of parameters. The first column corresponds to the log likelihood of the data given the model parameters. The following columns are: the proportion of ELS regions (if 3-state model), the average length of internal regions, the average length of external regions, the average length of ELS regions (if 3-state model), the probability of the outgroup being derived in an internal region at a fixed derived site (F parameter), the probability of the outgroup being derived in an external region at a fixed derived site (f parameter), the probability of the outgroup being derived in external regions at a polymorphic site (this should be 0 but it is hard coded as 0.01 to accomodate errors), followed by the probabilities of the outgroup being derived in internal regions at polymorphic sites depending on the number of derived sites in the population from 1 to n-1 (probabilities from the configuration file).

### Citation

The method is described in the following publication: Peyrégne, S., Boyle, M. J., Dannemann, M., & Prüfer, K. (2017). Detecting ancient positive selection in humans using extended lineage sorting. Genome research, 27(9), 1563–1572. https://doi.org/10.1101/gr.219493.116
