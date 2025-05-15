# Exact operator inference with minimal data

This is the Matlab implementation of the numerical experiments reported in: 

[arXiv citation]

[Bibtex]

Call the scripts <tt>go_chafee_infante.m</tt>, <tt>go_eight2.m</tt> and to <tt>go_burgers_periodic.m</tt> to perform the numerical experiments.

<tt>go_chafee_infante.m</tt> and <tt>go_burgers_periodic.m</tt> only take a few seconds.

<tt>go_eight2.m</tt> takes approx. 5 minutes if the boolean <tt>generatePODdata</tt> is set to  <tt>true</tt>, so the POD snapshot data is generated from scratch.
Otherwise, if the boolean <tt>generatePODdata</tt> is set to  <tt>false</tt>, the POD snapshot data is loaded and <tt>go_eight2.m</tt> takes approx. 3 minutes.

Call <tt>plot_operator_errors.m</tt> to generate the operator error and condition number plots.



