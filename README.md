# Prefrontal_sensory_GC

Current theories of consciousness can be categorized to some extent by their predictions about the role of the prefrontal cortex (PFC) in conscious perception. One family of theories propose that the PFC is necessary for conscious perception. The other family postulates that the PFC is not necessary and that other areas (e.g., posterior cortical areas) are more important for conscious perception. 

No-report paradigms could potentially arbitrate the debate as they have been proposed to distinguish the role of the PFC in task reporting from conscious perception. While previous no-report paradigms tend to point to a reduction of PFC activity, they have not examined the critical role of the PFC in “monitoring” the patterns of activity in the sensory cortex to generate conscious perception. 

To address this, we reanalysed EEG data from a no-report inattentional blindness paradigm (Shafto & Pitts, 2015) to examine the role of feedforward input patterns to the PFC from sensory cortices using nonparametric spectral Granger causality and quantified the amount of information that reflects the contents of consciousness using multivariate classifiers.

##########################################

Below we include the scripts used to preprocess the raw EEG files to obtain to GC stimates from sensory to PFC ROIs before applying classification analysis to determine whether these input patterns can differentiate between FACE and RANDOM stimuli.

List of scripts:
s01a_EEG_preProcessing.m
s01b_rmLineNoise.m
s01c_final_EEG_preProcessing.m
s02_sourceInv.m
s03_extractSW.m
s04_topCoord.m
s05_applyGC_setupDecode.m
s06_decode_GC_F.m


# Purpose of each of the scripts
*s01a_EEG_preProcessing.m

Here we take the raw EEG files per participant and assign the channels and trial type triggers.
We apply basic preprocessing steps of artefact detecton usin eyeblinks and noisy channels (>100 uv)
Next, we split the files into either FACE or RANDOM and epoch based on the trigger times

Inputs:
Outputs:


********************************
*s01b_rmLineNoise.m

Due to the presence of 60 Hz line noise in some participant's data (which was interferring with the decoding analysis), we removed the line noise per participant and FACE or RANDOM condition before implementing the final pre-processing script

Inputs:
Outputs:


********************************
*s01c_final_EEG_preProcessing.m

The final preprcoessing step of baseline correct is applied to the datefiles with line noise removed (rmline)

Inputs:
Outputs:


********************************
*s02_sourceInv.m

We can now take the single-trial ERPs and apply source conversion to determine cortical activity across the entire 3D brain during each trial (total of 8,196 voxels). For this, we rely on the SPM in-built source reconstruction function that utilises the Mutliple Sparse Priors greedy search method.

Intputs:
Outputs:


********************************
*s03_extractSW.m

Determine which voxel coordinates across the brain belong to each of the 8 defined ROIs: left and right Occipital cortex, Fusiform Region, dorsal-rostral prefrontal cortex and ventro-orbital prefrontal cortex.

Using this information, extract each of the single-trial source reconstructed estimates from EACH of the ROIs and save them in a file that will be used for the Granger Causality estimation in the next script.


********************************
*s04_topCoord.m

Using the extracted source waveforms from each voxel, we now want to determine (using the data from Phase 3), which (single) voxel within each ROI has the 'greatest difference between the FACE and RANDOM trails' (as evidenced by the highest t-statistic at any point within the 0 to 500 ms time window).
