# dLGN-IN-imaging

Source code to analyze calcium imaging data 
as described in Muellner & Roska 2023, https://doi.org/10.1101/2023.03.22.533751

ROI detection is optimized for irregular and dense structures and therefore works well for both dendrites and somata. 

Please cite our manuscript if you use this code for your own publication.

The functions can be run in the following order on your raw calcium imaging data. 
Inputs and outputs are described in the file headers.
The folder 'utils' contains several sub-functions used in the main functions below. 



```
zigzagcorr.m
```
Correct zigzag artifacts from scanning.

```
rigidmotioncorr.m
```
Performs rigid motion correction of the time series. Please note the comment on computation speed in the file header.

```
tembackg.m
```
Determine darkest pixels within a given radius of each pixel, which will be used for temporal background estimation;
apply for example to 95th percentile projection of your data: myprctile2D(data,95,dim).

```
ROIdetect.m
```
Detects regions of interest (ROIs) based on correlation-based seeding, adaptive thresholding and expansion.
No assumptions are made regarding the temporal or spatial footprints of the ROIs, but ROIs of size <= 1 pixel are disregarded (comment the last paragraph to switch this feature off).

```
ResponseMovGrat.m
```
Analyzes responses to moving grating stimuli.

```
ExtractReceptiveFields.m
```
Analyzes calcium responses to receptive field sparse noise stimuli (white flashes on black background).


