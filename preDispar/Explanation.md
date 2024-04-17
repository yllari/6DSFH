Here is explained all data preprocessing for the disPar/dirSFH pipeline.
## Full
A subset of data feeded to disPar with JUST a filter on absolute magnitude MG < 5.5.
Distance cut for quality purposes '1/new_parallax' < 1.2 (kpc)

## QSHAG
A subset of data feeded to dirSFH with the following cuts:
- Same brigthness cut as Full 'MG' <5.5
- Distance cut for quality purposes '1/new_parallax' < 1.2 (kpc)
  
Additionally

- Error in velocity cut 'radial_velocity_error' < 20 km/s
- Error in parallax 'parallax_over_error' > 5
- Cut in extinction 'AG' < 0.5
- Cut in excess color for quality purposes
  '(0.001+0.039*(bp_rp) < log10(phot_bp_rp_excess_factor)) & (log10(phot_bp_rp_excess_factor) < 0.12 + 0.039*(bp_rp)'
Nothing is being done to RUWE, disPar cannot simulate RUWE so it is not included

First subdataset will be given to disPar to generate the dispersed, completeness simulation of the mother CMD
Second dataset will be the target CMD for dirSFH.
