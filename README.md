# Code for the paper: Good and Fast Row-Sparse ah-Symmetric Reflexive Generalized Inverses

@ Gabriel Ponte (gabponte@umich.edu), Marcia Fampa (fampa@cos.ufrj.br), Jon Lee (jonxlee@umich.edu), Luze Xu (xuluze@ust.hk).

This project includes the data and code for: "Good and Fast Row-Sparse ah-Symmetric Reflexive Generalized Inverses," Open Journal of Mathematical Optimization.

### Data File ###

All instances are stored in ``.csv`` format in the ``Instances`` directory. 


### Experimental Setting ###

Experiments are executed through ``runExperiments.jl``, where you can enable or disable each method via Boolean flags:
1. ``run_ls``: Local search method based on determinant
2. ``run_admm1n``: ADMM minimizing the 1-norm of H
3. ``run_admm21n``: ADMM minimizing the 2,1-norm of H
4. ``run_20n``: ADMM targeting 2,0-norm sparsity
5. ``run_2120n``: ADMM minimizing the 2,1-norm of H and targeting 2,0-norm sparsity
6. ``run_grb``: Gurobi solver for minimizing the 1-norm of H
7. ``run_msk``: Mosek solver for minimizing the 2,1-norm of H
8. ``fixed_tol``: If true, ADMM termination uses primal and dual residuals. If false, dynamic tolerances are used,  proposed by [Distributed Optimization and Statistical Learning via the Alternating Direction Method of Multipliers](https://stanford.edu/~boyd/papers/pdf/admm_distr_stats.pdf)

To reproduce the experiments involving ADMM with target 2,0-norm, set ``fixed_tol = true``.

### Results ###

Our results are stored as ``.csv`` in the ``ResultsCSV`` directory. 

#### Column Definitions ###

| Col | Description                            |
| --- | -------------------------------------- |
| A   | # rows of matrix **A**                 |
| B   | # columns of matrix **A**              |
| C   | rank of **A**                          |
| D   | objective value (depends on algorithm) |
| E   | 1-norm of **H**                        |
| F   | 0-norm of **H**                        |
| G   | 2,1-norm of **H**                      |
| H   | 2,0-norm of **H** (# nonzero rows)     |
| I   | Frobenius-norm of **H**                |
| J   | elapsed time (sec)                     |
| K   | # iterations                           |
