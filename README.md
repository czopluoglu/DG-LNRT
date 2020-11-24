# DG-LNRT

A repository for the code being used to run additional analysis and simulations 
related to DG-LNRT paper.

https://psyarxiv.com/bqa3t/

Below is a list of folders with explanations for each folder:

- **/dglnrt_v1**: This folder includes the code for fitting DG-LNRT to the Real
Dataset 1 by Toton & Maynes. It also has codes for running sensitivity analyses with
different priors.

- **/parameter_recovery_2**: This folder includes the code for a simulation study
to examine the quality of item parameter estimates relevant to Real Dataset 1 
when sample size is 93 and number of items is 25, and when the responses for 
60 examinees are missing for the same 12items. This mimicking the Real Dataset 1 
with no contamination of item preknowledge). The exact same set of $\beta$, 
$\alpha$, and $\tau$ parameters were being used and fixed across 100 replications.

- **/parameter_recovery**: This folder includes the code for a similar 
simulation study in **/parameter recovery**. It has the exact same setting. The 
only difference is that a different set of $\beta$, $\alpha$, and $\tau$ 
parameters were being generated for each replication. The distributions of 
$\beta$, $\alpha$, and $\tau$ were slightly different than the ones observed in
the Real Dataset 1.

- **/dglnrt_v1_simulation**: This folder includes the code for a simulation study where
there is item preknowledge. The exact same setting of Real Dataset 1 is replicated.
The $\beta$, $\alpha$, and $\tau_t$, and $\tau_c$ parameters were being generated 
using similar distributions obtained from the Real Dataset 1 by fitting a 
multigroup lognormal response time model with a gating mechanism.

- **/dglnrt_v1_null**: This folder includes the code for a simulaiton study where
there is no item preknowledge. The $\beta$, $\alpha$, and $\tau$ parameters
were being generated using similar distributions obtained from the Real Dataset 1
by fitting a multigroup lognormal response time model with a gating mechanism.

- **/dglnrt_v2**:


