# dLGN-IN-imaging

Source code to analyze calcium imaging data 
as described in Muellner & Roska 2023, https://doi.org/10.1101/2023.03.22.533751

Please cite our manuscript if you use this code for your own publication.

The scripts can be run in the following order on your raw calcium imaging data:

```
zigzagcorr.m

% __________________________________
% Inputs:
%
% image: raw data
% zdim: dimension containing zigzag artifacts
% rdim: channel used for correction, if image has more than one channel
% __________________________________
% Outputs:
%
% bdshift: average shift per column (to double-check that average across columns is adequate)
% corrx1: shift in forward lines
% corrx2: shift in backward lines
% __________________________________

```
Correct zigzag artifacts from scanning.

```
rigidmotioncorr.m

%%
% NOTE: this function can be considerably sped up by restricting the
% calculation of the normalized cross correlation to the area of interest.
% To achieve this, e.g. add the following two lines to the Matlab function
% normxcorr2.m at line #84, following the definition of variable 'mn':
%     [ma, na] = size(A);
%     xcorr_TA = xcorr_TA(m:ma,n:na);
% If modified, uncomment the lines following normxcorr2 below.
%%
% __________________________________
% Inputs:
% 
% image: raw data; 1st/2nd dimension space, 3rd dimension time
% options.downs: factor for downsampling in time
% options.ups: factor for upsampling in space
% options.maxoff: maximally allowed offset in x or y, default 30
% __________________________________
% Outputs:
% 
% corrx: shift in 1st dimension
% corry: shift in 2nd dimension
% corrv: correlation coefficient with reference image

```
Performs rigid motion correction of the time series.

```
tembackg.m
```
Determine darkest pixels within a given radius of each pixel, which will be used for temporal background estimation;
apply for example to 95th percentile projection of your data: myprctile2D(data,95,dim).

```
ROIdetect.m

% Inputs:
%
% Cn: local correlation matrix to seed ROIs (see localcorr.m); optional 3rd dimension provides corr per time intervals (default: per 90 seconds)
% maxdata: projected data to add seeds for optimized detection of sparsely responding cells (default: 95th percentile projection, see myprctile2D.m)
% data: motion corrected data; dimension 1: time, 2/3: space
% chunks: data intervals to calculate correlations (default: per 90 seconds)
% rem: binary vector indicating datapoints to treat as NaN
% bgt: matrix for temporal background calculation (see tempbackg.m)
% maxradius: maximum ROI radius (default: 15 um)
% maxnroi: maximum number of ROIs detected (default: 10000)
% filters: filter length for running average (default: 0.5 seconds)
% file_name: filename to save code
% root: folder to save code 
% __________________________________
% Outputs:
%
% roi: spatial ROIs 
% sur: ROI surrounds
% orig: ROI origins
% C: temporal footprints
% Cfilt: filtered temporal footprints
% Cbg: temporal background
% M: origin coordinates
% nroi: number of detected ROIs
% copyInd: priority matrix for initialization of ROIs
% Ind: remaining priority matrix (with intialized ROIs removed)
% mCn: matrix of local correlations (projected across time intervals)
% corrthr: adaptive correlation threshold per ROI

```
Detects regions of interest (ROIs) based on correlation-based seeding, adaptive thresholding and expansion.

```
ResponseMovGrat.m

% Inputs:
%
% data: temporal footprint of the ROI
% proi: binary mask of the ROI
% plotdata: e.g. maximum projection of the full dataset
% screenrot: screen rotation with respect to horizontal axis
% time_first: onset of stimuli
% time_last: offset of stimuli
% orientation_val: stimuli orientation in degree 
% speed_val: stimuli speed 
% data_frame_start_ms: start of each imaging frame in ms
% save_dir: directory to save Results
% file_name: filename to save Results
% filters: filter length
% plotsw: 1 = plot results, 0 = do not plot results 
% roin: ROI number to save Results
% savesw: 1 = save results, 0 = do not save results 
% __________________________________
% Outputs:
%
% saved as 'Results' structure

```
Analyzes responses to moving grating stimuli.

```
ExtractReceptiveFields.m

% Inputs:
%
% timec: temporal footprint, 1st dim: time, 2nd dim: trials
% angle: visual angle, 1st dim: trials, 2nd dim: space axes
% reftime: start of each imaging frame in ms
% time_first: onset of stimuli
% time_last: offset of stimuli
% pre: interval considered baseline in ms
% ontime: flash length in ms
% offtime: time between flashes in ms
% filters: filter length
% nrep: shuffling repetitions
% __________________________________
% Outputs:
%
% onmat: matrix of ON response
% offmat: matrix of OFF responses
% onmat_shuf: shuffling control matrices for ON responses
% offmat_shuf: shuffling control matrices for OFF responses 
% on_est: 3x3 average filtered ON responses
% off_est: 3x3 average filtered OFF responses 
% onoff_est: sumed on_est+off_est
% x_on,y_on: ON recpective field center
% x_off,y_off: OFF receptive field center
% x_onoff,y_onoff: ON-OFF receptive field center

```
Analyzes calcium responses to receptive field sparse noise stimuli (white flashes on black background).

```
utils
```
Contains several functions used in the scripts above.



