# Extragradient Method: O(1/K) Last-Iterate Convergence for Monotone Variational Inequalities and Connections With Cocoercivity
Code for the paper "Extragradient Method: O(1/K) Last-Iterate Convergence for Monotone Variational Inequalities and Connections With Cocoercivity" by Eduard Gorbunov, Nicolas Loizou, Gauthier Gidel. The paper is accepted to AISTATS 2022. The code comes jointly with the paper.

ArXiv: https://arxiv.org/abs/2110.04261

## How to run the code?
**Packages.** To run our MATLAB code one should install Performance Estimation Toolbox (PESTO) https://github.com/AdrienTaylor/Performance-Estimation-Toolbox, SEDUMI https://yalmip.github.io/solver/sedumi/, YALMIP https://yalmip.github.io/ and add the to the path when executing our code. To tun the Python code one needs to install CVXPY https://www.cvxpy.org/.

**Folders.** Before running the code, one needs to create two folders in the directory with the code: "dump" and "plot".

## Files

- *EG_showing_non_cocoercivity.m*: the code for numerically solving problem (15) from the paper. This is needed to show non-cocoercivity of Extragradient operator.
- *EG_showing_non_cocoercivity.ipynb*: visualization of the results obtained via *EG_showing_non_cocoercivity.m* and constructing the example for Theorem 3.2 from the paper (this code generates Figures 1 and 2 from the paper).
- *EG_last_iter_squared_norm.m*: the code for computing the worst-case performance of Extragradient method in terms of reducing ||F(x^k)||^2.
- *EG_PEP_for_squared_norm_visualization.ipynb*: the code for visualizing the worst-case performance of Extragradient method in terms of reducing ||F(x^k)||^2 (this code generates Figure 3 from the paper). Before running this code, one needs to execute *EG_last_iter_squared_norm.m*.
- *EG_guessing_the_proof_last_iter_convergence.ipynb*: in this file, we solve numerically PEP for the worst-case ||F(x^1)||^2 - ||F(x^0)||^2, where x^1 is produced by Extragradient method from the point x^0. We also print out the solution of the dual problem for different values of the stepsize and conclude that they are always close to certain numerical values that we use in the regorous proof of Lemma 3.2 from our paper.
- *EG_norm_of_EG_operator_can_grow.m*: this code illustrates the remark about problem (21) from our paper (the norm of Extragradient operator can grow).
- *EG_different_stepsizes_break_the_proof.m*: this code illustrates the remark about problem (22) from our paper (\gamma_1 = \gamma_2 is crucial for having ||F(x^{k+1})||^2 <= ||F(x^k)||^2).

## Licence
If you want to use our code, please, cite our work:
> @article{gorbunov2021extragradient,
> title={Extragradient Method: $ O (1/K) $ Last-Iterate Convergence for Monotone Variational Inequalities and Connections With Cocoercivity},
> author={Gorbunov, Eduard and Loizou, Nicolas and Gidel, Gauthier},
> journal={arXiv preprint arXiv:2110.04261},
> year={2021}
> }
