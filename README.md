# Prefrontal_sensory_GC
Supporting code for Rowe, Garrido & Tsuchiya (2023)

Current theories of consciousness can be categorized to some extent by their predictions about the role of the prefrontal cortex (PFC) in conscious perception. One family of theories propose that the PFC is necessary for conscious perception. The other family postulates that the PFC is not necessary and that other areas (e.g., posterior cortical areas) are more important for conscious perception. 

No-report paradigms could potentially arbitrate the debate as they have been proposed to distinguish the role of the PFC in task reporting from conscious perception. While previous no-report paradigms tend to point to a reduction of PFC activity, they have not examined the critical role of the PFC in “monitoring” the patterns of activity in the sensory cortex to generate conscious perception. 

To address this, we reanalysed EEG data from a no-report inattentional blindness paradigm (Shafto & Pitts, 2015) to examine the role of feedforward input patterns to the PFC from sensory cortices using nonparametric spectral Granger causality and quantified the amount of information that reflects the contents of consciousness using multivariate classifiers.

##########################################

Below we include the scripts used to preprocess the raw EEG files to obtain to GC stimates from sensory to PFC ROIs before applying classification analysis to determine whether these input patterns can differentiate between FACE and RANDOM stimuli. <br />

List of scripts:  <br />
s01a_EEG_preProcessing.m  <br />
s01b_rmLineNoise.m  <br />
s01c_final_EEG_preProcessing.m  <br />
s02_sourceInv.m  <br />
s03_extractSW.m  <br />
s04_topCoord.m  <br />
s05_applyGC_setupDecode.m  <br />
s06_decode_GC_F.m  <br />
s07_decode_SWaveGrouped.m  <br />


# Description of each script
s01a_EEG_preProcessing.m <br />

Here we take the raw EEG files per participant and assign the channels and trial type triggers. <br />
We apply basic preprocessing steps of artefact detecton usin eyeblinks and noisy channels (>100 uv)  <br />
Next, we split the files into either FACE or RANDOM and epoch based on the trigger times  <br />

Inputs: Raw EEG files in EDF format: for example, 'P0002_Move_Markers.edf'  <br />
        Fudicials used for source-allocaton 'chanFudicials_ShaftoFINAL.mat'  <br />
 <br />
Outputs: Partially pre-processed EEG data split by FACE or RANDOM trials: 'e_random_aEBn_noHP_spmeeg_P0002_Move_Markers.mat'  <br />


****************************************************************
s01b_rmLineNoise.m  <br />

Due to the presence of 60 Hz line noise in some participant's data (which was interferring with the decoding analysis), we removed the line noise per participant and FACE or RANDOM condition before implementing the final pre-processing script  <br />

Inputs: ae_random_aEBn_noHP_spmeeg_P0002_Move_Markers.mat <br />
        ae_faces_aEBn_noHP_spmeeg_P0002_Move_Markers.mat  <br />
Outputs: ae_rm_final_random_aEBn_noHP_spmeeg_P0002_Move_Markers.mat  <br />
         ae_rm_final_faces_aEBn_noHP_spmeeg_P0002_Move_Markers  <br />


****************************************************************
s01c_final_EEG_preProcessing.m  <br />

The final preprcoessing step of baseline correct is applied to the datefiles with line noise removed (rmline) after we split the participant FACE and RANDOM trial files into INDIVIDUAL files per trial (we also remove any noisy trials or those with a button response).  <br />

Inputs: ae_rm_final_random_aEBn_noHP_spmeeg_P0002_Move_Markers.mat <br />
        ae_rm_final_faces_aEBn_noHP_spmeeg_P0002_Move_Markers  <br />
        
Outputs:  <br />
bm_Sae_rm_final_faces_aEBn_noHP_P0002_Phase1_-100to500ms_trial_X.mat (where X will equal from 1 to the number of trials)  <br />
bm_Sae_rm_final_random_aEBn_noHP_P0002_Phase1_-100to500ms_trial_X.mat <br />

****************************************************************
s02_sourceInv.m  <br />

We can now take the single-trial ERPs and apply source conversion to determine cortical activity across the entire 3D brain during each trial (total of 8,196 voxels). For this, we rely on the SPM in-built source reconstruction function that utilises the Mutliple Sparse Priors greedy search method.  <br />

Intputs: bm_Sae_rm_final_faces_aEBn_noHP_P0002_Phase1_-100to500ms_trial_X.mat  <br />
         bm_Sae_rm_final_random_aEBn_noHP_P0002_Phase1_-100to500ms_trial_X.mat  <br />
Outputs: As above (but now we have a source-level matrix assigned to these files and sourcewaves images)  <br />


****************************************************************
s03_extractSW.m  <br />

Determine which voxel coordinates across the brain belong to each of the 8 defined ROIs: left and right Occipital cortex, Fusiform Region, dorsal-rostral prefrontal cortex (dr-PFC) and ventro-lateral prefrontal cortex (vo-PFC).  <br />

Using this information, extract each of the single-trial source reconstructed estimates from EACH of the ROIs and save them in a file that will be used for the Granger Causality estimation in the next script.  <br />

Inputs:  bm_Sae_rm_final_faces_aEBn_noHP_P0002_Phase1_-100to500ms_trial_X.mat  <br />
         bm_Sae_rm_final_random_aEBn_noHP_P0002_Phase1_-100to500ms_trial_X.mat  <br />
          <br />
Outputs: Final_P0002_PhaseX_rmline_0to500ms_faces_Source_Waveform_at_for_ALL_trials_Coords_from_1_to_8196.mat (where X can equal 1 to 3)  <br />
Final_P0002_PhaseX_rmline_0to500ms_random_Source_Waveform_at_for_ALL_trials_Coords_from_1_to_8196.mat  <br />
CoordIdxs_for_P0002_all8196_Phase_X.mat  <br />

****************************************************************
s04_topCoord.m  <br />

Using the extracted source waveforms from each voxel, we now want to determine (using the data from Phase 3), which (single) voxel within each ROI has the 'greatest difference between the FACE and RANDOM trails' (as evidenced by the highest t-statistic at any point within the 0 to 500 ms time window).  <br />

Once we have determined this 'TOP VOXEL" within each of the 8 ROIs, we use this voxel coordinate and extract the data from Phase 1 and Phase 2 for use in our GC analysis.  <br />

Inputs: Coordinates_8196_voxels_byROI.mat (this determines which voxel coordinates below to what ROI)  <br />
        Final_P0002_PhaseX_rmline_0to500ms_faces_Source_Waveform_at_for_ALL_trials_Coords_from_1_to_8196.mat (where X = 1 to 3)  <br />
        Final_P0002_PhaseX_rmline_0to500ms_random_Source_Waveform_at_for_ALL_trials_Coords_from_1_to_8196.mat (where X = 1 to 3)  <br />
 <br />
Outputs: CoordSorted_rmline_P0002_PhaseX_EBRem_XXX_by_TstatOverTime_0_to_500ms.mat (where X = 1 to 3 and XXX = one of the 8 ROIs) <br />


****************************************************************
s05_applyGC_setupDecode.m <br />

For one participant, one phase and one conditions (i.e., face or random) at a time, extract the source waveform (over time) from pairs of the ROIs (i.e., 4 PFC and 4 sensory ROIs = 16 ROI combinations) and estimate the Granger Causality in the FORWARDS direction from the sensory node to the PFC node. Save each of the GC results by COMBINING all GC estimates from the 16 pairs and in a format ready for classification analysis in the next script. <br />

Inputs: for sensory node = CoordSorted_rmline_P0002_PhaseX_EBRem_XXX_by_TstatOverTime_0_to_500ms.mat (XXX = one of 4 sensory ROIs)  <br />
for PFC node = CoordSorted_rmline_P0002_PhaseX_EBRem_XXX_by_TstatOverTime_0_to_500ms.mat (XXX = one of the 4 PFC ROIs)  <br />
 <br />
Outputs: GC_P0002_PhaseX_CVmethod_thisRep3_zscore_1_bw_XX_and_XXX_by_TstatOverTime_usingTrainProp_0.7_and_0_to_500ms.m (X = 1 to 3, XX = one of the 4 sensory ROIs and XXX = one of the PFC ROIs) <br />


****************************************************************
s06_decode_GC_F.m  <br />

Using the GC estimates from the previous step (i.e., from the 16 combinatons of ROI pairs), run SVM classification (70/30 train/test split) using 10 random splits of the data (cross-validation) to determine if we can classify between a FACE or RANDOM trial using the GC estimates. The output results reflect the classifcation accuracy (10 repetitions) for this participant and phase.  <br />

Inputs: GC_P0002_PhaseX_CVmethod_thisRep3_zscore_1_bw_XX_and_XXX_by_TstatOverTime_usingTrainProp_0.7_and_0_to_500ms.m  <br />
Outputs: results_GC_toPFC_ForwardsFrom_Sensory_rmline_P0002_PhaseX_55Hz_MultiROI_ALL_PFC_ALL_FFA_ALL_OCC_0to500ms_10Reps.m  <br />


****************************************************************
s07_decode_SWaveGrouped.m  <br />

For one participant and one phase, use the source waveform (over time) from either sensory or PFC ROIs (i.e., 4 PFC or 4 sensory ROIs).  Run SVM classification (70/30 train/test split) using 10 random splits of the data (cross-validation) to determine if we can classify between a FACE or RANDOM trial using the GC estimates. The output results reflect the classifcation accuracy (10 repetitions) for this participant and phase.  <br />

Inputs: for sensory node = CoordSorted_rmline_P0002_PhaseX_EBRem_XXX_by_TstatOverTime_0_to_500ms.mat (XXX = one of 4 sensory ROIs)  <br />
for PFC node = CoordSorted_rmline_P0002_PhaseX_EBRem_XXX_by_TstatOverTime_0_to_500ms.mat (XXX = one of the 4 PFC ROIs)  <br />
 <br />
Outputs = results_SourceWaveform_allPFC_rmline_P0002_PhaseX_55Hz_MultiROI_0to500ms_10Reps.mat (X = Phase 1 - 3)  <br />
          results_SourceWaveform_allSensory_rmline_P0002_PhaseX_55Hz_MultiROI_0to500ms_10Reps.mat  <br />

****************************************************************
****************************************************************
# Supplementary code
The following scripts enable the analyses discussed in the Supplementary Material <br />

(a) supp_s08_decode_APz.m  <br />
(b) supp_s08_decode_GC_B.m  <br />
(c) supp_s08_GC_sensory_B.m  <br />
(d) supp_s08_GC_sensory_F.m  <br />
(e) supp_s08_decode_Oz.m  <br />
(f) supp_s08_decode_PWR.m  <br />
(g) supp_s08_decode_SWave.m  <br />

****************************************************************
****************************************************************

supp_s08_decode_APz.m = Classification analysis using scalp-EEG from elecltrode APz  <br />

Input -> NEW_ALL_APZ_XXXXX_TrialData_P0002_All_Phases.mat (where XXXX = faces or random) <br />
Output -> results_APZ_decode_P0002_PhaseX_CVmethod_10REPS__AT_APZ__0to500ms_thisRepXX.mat (where X = Phase 1 to 3, XX = 1 to 10 repetitions) <br />

****************************************************************

supp_s08_decode_GC_B.m = Classification analysis using GC estimates BACKWARDS from PFC to sensory cortices  <br />

Input ->
Output ->

****************************************************************

supp_s08_decode_sensory_B.m = Classification analysis using GC estimates BACKWARDS from Fusiform Region to Occipital (sensory only)  <br />

Input ->
Output ->

****************************************************************

supp_s08_decode_sensory_F.m = Classification analysis using GC estimates FORWARDS from Occipital to Fusiform Region (sensory only)  <br />

Input ->
Output ->

****************************************************************

supp_s08_decode_Oz.m = Classification analysis using scalp-EEG from elecltrode Oz  <br />

Input -> NEW_ALL_OZ_XXXXX_TrialData_P0002_All_Phases.mat (where XXXX = faces or random) <br />
Output -> results_OZ_decode_P0002_PhaseX_CVmethod_10REPS__AT_APZ__0to500ms_thisRepXX.mat (where X = Phase 1 to 3, XX = 1 to 10 repetitions) <br />


****************************************************************

supp_s08_decode_PWR.m = Classificaton analysis using the PWR at each ROI (or grouped by PFC and sensory ROIs)  <br />

Input ->
Output ->

****************************************************************

supp_s08_decode_SWave.m = Classificaton analysis using the source-leve waveforms at each ROI (individually)  <br />

Input ->
Output ->

****************************************************************
