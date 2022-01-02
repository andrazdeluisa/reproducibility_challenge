# Reproducibility challenge 2021

This work aims to reproduce the experiments from the paper Thompson Sampling for Bandits with Clustered Arms (Carlsson et al. 2021, available on [link](https://www.ijcai.org/proceedings/2021/0305.pdf)), presented at the IJCAI 2021.

## Repository structure

Each *bandits_\*.R* script reproduces its part of the experiments, while *utils.R* contains functions that builds the required synthetic datasets and the implementations of all used algorithms.

In the folder *figures*, the outputs of the scripts (i.e. results of experiments in the form of visualisations) are saved.

## Dependencies

The code for the project is written in R (version 4.1.1.). To run all the scripts you need to have installed the following libraries (available on CRAN):

- ggplot2 (version 3.3.4),
- mvtnorm (version 1.1.3) (needed just for the contextual bandits).
