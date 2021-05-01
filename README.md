# DG-LNRT

A repository for the code being used to run additional analysis and simulations 
related to DG-LNRT paper.

Below is a list of folders with explanations for each folder:

- **/dglnrt_v1**: This folder includes the code for fitting DG-LNRT to the first real
dataset (Toton & Maynes, 2019). 

- **/dglnrt_v1_null**: This folder includes the code for a simulation study where
there is no item preknowledge. The $\beta$, $\alpha$, and $\tau$ parameters
were being generated using similar distributions obtained from the first real dataset
by fitting a multigroup lognormal response time model with a gating mechanism.

- **/dglnrt_v1_simulation**: This folder includes the code for a simulation study where
there is item preknowledge. The exact same setting of the first real dataset was replicated.
The $\beta$, $\alpha$, and $\tau_t$, and $\tau_c$ parameters were being generated 
using similar distributions obtained from the first real dataset by fitting a 
multigroup lognormal response time model with a gating mechanism.

- **/dglnrt_v2**: This folder includes the code for fitting DG-LNRT to the second real
dataset (Cizek & Wollack, 2016).

- **/dglnrt_v2_null**: This folder includes the code for a simulation study where
there is no item preknowledge. The $\beta$, $\alpha$, and $\tau$ parameters
were being generated using similar distributions obtained from the second real dataset
by fitting a multigroup lognormal response time model with a gating mechanism.

- **/dglnrt_v2_simulation**: This folder includes the code for a simulation study where
there is item preknowledge. The exact same setting of the second real dataset was replicated.
The $\beta$, $\alpha$, and $\tau_t$, and $\tau_c$ parameters were being generated 
using similar distributions obtained from the first real dataset by fitting a 
multigroup lognormal response time model with a gating mechanism.

- **/modelfit**: This folder includes the code for examining the standardized residuals for both datasets after fitting the DG-LNRT model.

- **/sinharay2020**: This folder includes the code to compute the frequentis probabilities
based on the Z-statistic provided by Sinaharay (2020) for both datasets and follow-up simulations

