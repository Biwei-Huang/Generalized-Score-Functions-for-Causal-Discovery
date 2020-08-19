# Generalized-Score-Functions-for-Causal-Discovery

Copyright (c) 2017-2018 Biwei Huang

Causal structure learning by greedy equivalence search with generalized score functions (which is applicable to mixed continuous and discrete data, data with Gaussian or non-Gaussian distributions, linear or nonlinear causal mechanisms, and variables with multi-dimensionalities.)

### IMPORTANT FUNCTIONS:

function [Record] = GES(X,score_type,multi_sign,maxP,parameters)

* INPUT:
  * X: Data with T*D dimensions
  * score_type: the score function you want to use
     *     score_type = 1: cross-validated likelihood
     *     score_type = 2: marginal likelihood
  * multi_sign: 1: if are multi-dimensional variables; 0: otherwise
  * maxP: allowed maximum number of parents when searching the graph
  * parameters: 
     *   when using CV likelihood, 
     *      parameters.kfold: k-fold cross validation
     *      parameters.lambda: regularization parameter
     *   parameters.dlabel: for variables with multi-dimensions, indicate which dimensions belong to the i-th variable.


* OUTPUT:
  * Record.G: learned causal graph
      *      G(i,j) = 1: i->j (the edge is from i to j)
      *      G(i,j) = -1: i-j (the direction between i and j is not determined)
      *      G(i,j) = 0:  i j (no causal edge between i and j)
  * Record.update1: each update (Insert operator) in the forward step
  * Record.update2: each update (Delete operator) in the backward step
  * Record.G_step1: learned graph at each step in the forward step
  * Record.G_step2: learned graph at each step in the backward step
  * Record.score: the score of the learned graph


### EXAMPLE:
see example1.m


### CITATION:
	
Biwei Huang ,Kun Zhang , Yizhu Lin, Bernhard Scholkopf, Clark Glymour. Generalized Score Functions for Causal Discovery. KDD, 2018.
