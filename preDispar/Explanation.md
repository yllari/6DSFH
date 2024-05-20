Here, all data preprocessing for the disPar/dirSFH pipeline is explained.
An extinction, absolute magnitude and color is assigned with AG_generation. 
After that, two datasets are generated with QSHAG_selection

## Full
A subset of data feeded to disPar with:
- A filter on absolute magnitude MG < 5.5.
- Distance cut for quality purposes '1/new_parallax' < 1.2 (kpc)

## QSHAG
A subset of data feeded to dirSFH with the following cuts:
- Same brigthness cut as Full 'MG' <5.5
- Distance cut for quality purposes '1/new_parallax' < 1.2 (kpc)
  
Additionally

- Error in velocity cut 'radial_velocity_error' < 20 km/s (to be selected as criteria. Must be consistent with disPar)
- Error in parallax 'parallax_over_error' > 5
- Cut in extinction 'AG' < 0.5
- Cut in excess color for quality purposes
  '(0.001+0.039*(bp_rp) < log10(phot_bp_rp_excess_factor)) & (log10(phot_bp_rp_excess_factor) < 0.12 + 0.039*(bp_rp)'

First subdataset will be given to disPar to generate the dispersed completeness simulation of the mother CMD.

Second dataset will be the target CMD for dirSFH.
