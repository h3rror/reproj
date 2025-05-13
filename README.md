# Exact operator inference with minimal data

This is the Matlab implementation of the numerical experiments reported in: 

[arXiv citation]

Call the scripts <tt>go_chafee_infante.m</tt>, <tt>go_eight2.m</tt> and to <tt>go_burgers_periodic.m</tt> to perform the numerical experiments. Call <tt>plot_operator_errors.m</tt> to generate the operator error and condition number plots.



# Sampling low-dimensional Markovian dynamics for pre-asymptotically recovering reduced models from data with operator inference.

This is a Matlab implementation of data sampling with re-projection for recovering reduced models from data as described in:

[1] Peherstorfer, B. [Sampling low-dimensional Markovian dynamics for pre-asymptotically recovering reduced models from data with operator inference](https://arxiv.org/pdf/1908.11233).
arXiv:1908.11233, 2019. ([Download preprint](https://arxiv.org/pdf/1908.11233))<details><summary>BibTeX</summary><pre>
@article{P19ReProj,
title = {Sampling low-dimensional Markovian dynamics for pre-asymptotically recovering reduced models from data with operator inference},
author = {Peherstorfer, B.},
journal = {arXiv:1908.11233},
year = {2019},
}</pre></details>


Call the script <tt>reproj.m</tt> to generate the re-projected trajectory of the example discussed in Section 2.4 and Section 5.1 in [1]. The figure that is generated is Figure 2 in [1]. The matrix used to generate Figure 2 in [1] is stored in <tt>tcA.mat</tt> and <tt>tcA.txt</tt>.

See [https://cims.nyu.edu/~pehersto/](https://cims.nyu.edu/~pehersto/) for other publications on [nonintrusive model reduction and learning reduced models from data](https://cims.nyu.edu/~pehersto/).

---
Other references for nonintrusive model reduction with operator inference:

[2] Peherstorfer, B. and Willcox, K.
[Data-driven operator inference for non-intrusive projection-based model reduction.](https://www.sciencedirect.com/science/article/pii/S0045782516301104)
Computer Methods in Applied Mechanics and Engineering, 306:196-215, 2016.
([Download preprint](https://cims.nyu.edu/~pehersto/preprints/Non-intrusive-model-reduction-Peherstorfer-Willcox.pdf))

[3] Qian, E., Kramer, B., Marques, A., and Willcox, K.
[Transform & Learn: A data-driven approach to nonlinear model reduction](https://arc.aiaa.org/doi/10.2514/6.2019-3707).
In the AIAA Aviation 2019 Forum, June 17-21, Dallas, TX. ([Download preprint](https://www.dropbox.com/s/5znea6z1vntby3d/QKMW_aviation19.pdf?dl=0))

[4] Qian, E., Kramer, B., Peherstorfer, B. and Willcox, K. [Lift & Learn: Physics-informed machine learning for large-scale nonlinear dynamical systems](https://www.sciencedirect.com/science/article/abs/pii/S0167278919307651).
Physica D: Nonlinear Phenomena, 2020.([Download preprint](https://arxiv.org/pdf/1912.08177))

[5] Swischuk, R., Kramer, B., Huang, C. and Willcox, K. [Learning Physics-Based Reduced-Order Models for a Single-Injector Combustion Process](https://arc.aiaa.org/doi/abs/10.2514/1.J058943). AIAA Journal, 2020. ([Download preprint](https://arxiv.org/pdf/1908.03620))
