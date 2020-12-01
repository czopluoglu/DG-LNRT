# DG-LNRT

A repository for the code being used to run additional analysis and simulations 
related to DG-LNRT paper.

https://psyarxiv.com/bqa3t/

Below is a list of folders with explanations for each folder:

- **/dglnrt_v1**: This folder includes the code for fitting DG-LNRT to the first real
dataset (Toton & Maynes, 2019). It also has codes for running sensitivity analyses with
different priors.

- **/dglnrt_v1_null**: This folder includes the code for a simulaiton study where
there is no item preknowledge. The $\beta$, $\alpha$, and $\tau$ parameters
were being generated using similar distributions obtained from the first real dataset
by fitting a multigroup lognormal response time model with a gating mechanism.

- **/dglnrt_v1_simulation**: This folder includes the code for a simulation study where
there is item preknowledge. The exact same setting of the first real dataset was replicated.
The $\beta$, $\alpha$, and $\tau_t$, and $\tau_c$ parameters were being generated 
using similar distributions obtained from the first real dataset by fitting a 
multigroup lognormal response time model with a gating mechanism.

- **/parameter_recovery_2**: This folder includes the code for a simulation study
to examine the quality of item parameter estimates relevant to first real dataset 
when sample size is 93 and number of items is 25, and when the responses for 
60 examinees are missing for the same 12 items. This mimics the first real dataset
with no contamination of item preknowledge. The exact same set of $\beta$, 
$\alpha$, and $\tau$ parameters were being used and fixed across 100 replications.

- **/parameter recovery**: This folder includes the code for a similar 
simulation study in **/parameter recovery**. It has the exact same setting. The 
only difference is that a different set of $\beta$, $\alpha$, and $\tau$ 
parameters were being generated for each replication. The distributions of 
$\beta$, $\alpha$, and $\tau$ were slightly different than the ones observed in
the first real dataset. This is on purpose to create a mismatch between true item
parameter distributions and prior distributions used in estimation.

- **/dglnrt_v2**: This folder includes the code for fitting DG-LNRT to the second real
dataset (Cizek & Wollack, 2016).

- **/dglnrt_v2_null**: This folder includes the code for a simulaiton study where
there is no item preknowledge. The $\beta$, $\alpha$, and $\tau$ parameters
were being generated using similar distributions obtained from the second real dataset
by fitting a multigroup lognormal response time model with a gating mechanism.

- **/dglnrt_v2_simulation**: This folder includes the code for a simulation study where
there is item preknowledge. The exact same setting of the second real dataset was replicated.
The $\beta$, $\alpha$, and $\tau_t$, and $\tau_c$ parameters were being generated 
using similar distributions obtained from the first real dataset by fitting a 
multigroup lognormal response time model with a gating mechanism.

