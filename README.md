# FER-KL-L1R
Sequential safe feature elimination rule for $L_{1}$-regularized regression with Kullback-Leibler divergence
=====================================================

Author: Hongmei Wang

This is a **Matlab** code corresponds to the following paper:

[1] C. F. Dantas, E. Soubies and C. Févotte  “Expanding Boundaries of GAP Safe Screening,” submitted to JMLR 2021.

[2] C. F. Dantas, E. Soubies and C. Févotte  “Safe screening for sparse regression with the Kullback-Leibler divergence,” submitted to ICASSP 2021.

[3] H. Wang, K. Jiang and Y. X "Sequential safe feature elimination rule for $L_{1}$-regularized regression with Kullback-Leibler divergence"  submitted to Neural Networks.

-----------
Disclaimers
-----------

- Datasets are not provided and need to be downloaded by the user and placed in a subfolder ./datasets (AllBooks, NIPSpaper, 20Newsgroup,  Encyclopedia, TasteProfile and MNIST). See file ./datasets/load_dataset.m and paper [3] for further instructions.
- Synthetic experiments can be performed directly.

-----------------
Files description
-----------------

- main.m: main script that launches all experiments. Simulation parameters can be set in this file. 
- path_KL_solvers.m: this script is called by the main.m script.
- Compute_final_result.m: this script is used to compute the final result used in our article.
- ./screening_rules/: contains the functions performing the proposed screening tests (as well as the required precalculations).
- ./solvers/: contains all tested solvers along with their corresponding version using screening:

Some plots are automatically generated, only if a single regularization value is chosen.
