# DG-LNRT

A repository for the code being used to run additional analysis and simulations 
related to DG-LNRT paper.

We would like to thank the Research Advanced Computing Services that maintain Talapas, The University of Oregon’s performance computing cluster, for providing resources to run extensive simulations that contributed to the research results reported within this paper. URL: https://racs.uoregon.edu/talapas

Below is a list of folders with explanations:

- **dglnrt_v1**: This folder includes the code for fitting DG-LNRT to the first real
dataset (Toton & Maynes, 2019). 

- **dglnrt_v1_null**: This folder includes the code for a simulation study where
there is no item preknowledge. The $\beta$, $\alpha$, and $\tau$ parameters
were being generated using similar distributions obtained from the first real dataset
by fitting a multigroup lognormal response time model with a gating mechanism [(Zopluoglu et al., 2001).](https://onlinelibrary.wiley.com/doi/full/10.1111/emip.12428)

- **dglnrt_v1_simulation**: This folder includes the code for a simulation study where
there is item preknowledge. The exact same setting of the first real dataset was replicated.
The $\beta$, $\alpha$, and $\tau_t$, and $\tau_c$ parameters were being generated 
using similar distributions obtained from the first real dataset by fitting a 
multigroup lognormal response time model with a gating mechanism [(Zopluoglu et al., 2001).](https://onlinelibrary.wiley.com/doi/full/10.1111/emip.12428). It is assumed that there was a perfect recovery of compromised items, and the item compromise status was correctly specified for all 12 compromised items during the model fitting process.

- **dglnrt_v1_simulation_partially_identified**: This folder includes the code for a simulation study where there is item preknowledge with partially identified compromised items. The setting is exactly identical to the **/dglnrt_v1_simulation** with one important difference. It is assumed that there was an imperfect recovery of compromised items, and the item compromise status was correctly specified for only six out of 12 compromised items during the model fitting process. The remaining six compromised items were treated as not compromised during the model fitting process.

- **dglnrt_v1_simulation_misidentified**: This folder includes the code for a simulation study where
there is item preknowledge with partially identified compromised items and misidentified uncompromised items. The setting is exactly identical to the **/dglnrt_v1_simulation** with one important difference. It is assumed that the item compromise status was correctly specified for only six out of 12 compromised items during the model fitting process. The remaining six compromised items were treated as not compromised during the model fitting process. In addition, six out of 13 uncompromised items were assumed to be incorrectly identified as compromised, and they were treated as compromised items during the model fitting process.

- **dglnrt_v2**: This folder includes the code for fitting DG-LNRT to the second real
dataset (FormA, Cizek & Wollack, 2016).

- **dglnrt_v2_null**: This folder includes the code for a simulation study where
there is no item preknowledge. The $\beta$, $\alpha$, and $\tau$ parameters
were being generated using similar distributions obtained from the second real dataset
by fitting a multigroup lognormal response time model with a gating mechanism [(Zopluoglu et al., 2001).](https://onlinelibrary.wiley.com/doi/full/10.1111/emip.12428).

- **dglnrt_v2_simulation**: This folder includes the code for a simulation study where
there is item preknowledge. The exact same setting of the second real dataset was replicated.
The $\beta$, $\alpha$, and $\tau_t$, and $\tau_c$ parameters were being generated 
using similar distributions obtained from the first real dataset by fitting a 
multigroup lognormal response time model with a gating mechanism [(Zopluoglu et al., 2001).](https://onlinelibrary.wiley.com/doi/full/10.1111/emip.12428). It is assumed that there was a perfect recovery of compromised items, and the item compromise status was correctly specified for all 64 compromised items during the model fitting process.

- **dglnrt_v2_simulation_partially_identified**: This folder includes the code for a simulation study where there is item preknowledge with partially identified compromised items. The setting is exactly identical to the **/dglnrt_v2_simulation** with one important difference. It is assumed that there was an imperfect recovery of compromised items, and the item compromise status was correctly specified for 32 out of 64 compromised items during the model fitting process. The remaining 32 compromised items were treated as not compromised during the model fitting process.

- **dglnrt_v2_simulation_misidentified**: This folder includes the code for a simulation study where
there is item preknowledge with partially identified compromised items and misidentified uncompromised items. The setting is exactly identical to the **/dglnrt_v2_simulation** with one important difference. It is assumed that the item compromise status was correctly specified for only 32 out of 64 compromised items during the model fitting process. The remaining 32 compromised items were treated as not compromised during the model fitting process. In addition, 32 out of 106 uncompromised items were assumed to be incorrectly identified as compromised, and they were treated as compromised items during the model fitting process.

- **modelfit**: This folder includes the code for examining the standardized residuals for both datasets after fitting the DG-LNRT model.

- **sinharay2020**: This folder includes the code to compute the frequentis probabilities
based on the Z-statistic provided by Sinaharay (2020) for both datasets and follow-up simulations

- **dgirt_v1**: This folder inclused the code to fit the Deterministic Gated IRT Model using raw item responses for the first dataset by Toton & Maynes.

- **dgirt_v2**: This folder inclused the code to fit the Deterministic Gated IRT Model using raw item responses for the second dataset by Cizek and Wollack (FormA).
