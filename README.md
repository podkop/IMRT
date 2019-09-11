# IMRT
Python programs for research in Intensity Modulated Radiotherapy Treatment planning

File "dose_vol_hist.py" contains functions for producing dose-volume histograms (DVHs) as Excel files
The usage scenario: given several alternative IMRT plans for a patient, we want to produce DVHs
  a) showing dose distributions in all regions of interests (ROIs) for each plan,
  b) comparing dose distributions in specific ROIs for different plans

The main function **dvh2excel** takes as input data of dose distributions in all ROIs for all treatment plans (simply an array \[plans x rois x doses_in_voxels\]) and creates an *.xlsx* file with 
  * data for DVHs (each treatment plan -> a separate worksheet),
  * DVHs grouped by treatment plans (for each treatment plan, DVHs of all ROIs in a separate chart in the corresponding worksheet),
  * in the case of more than one plan, DVHs grouped by ROIs in the "Comparison" worksheet (for each ROI, DVHs resulted from different plans in a separate chart)

*I hope the comments in the code are explanatory*
