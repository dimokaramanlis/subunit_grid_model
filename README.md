# subunit_grid_model
Functions to fit subunit grid models to spiking data from retinal ganglion cells. 

This repository is a supplement to the manuscript ["Nonlinear receptive fields evoke redundant retinal coding of natural scenes"](https://www.biorxiv.org/content/10.1101/2023.01.10.523412v1?rss=1).

## Requirements and installation

- MATLAB >= R2016b
- MATLAB Toolboxes:
    - Parallel Computing Toolbox
    - Signal Processing Toolbox
    - Statistics and Machine Learning Toolbox
- An NVIDIA GPU and CUDA capabilities

Model optimization happens on the GPU for performance reasons, but a GPU is not strictly required. Replacing `gpuArray` arrays with standard `single` or `double` arrays will do the job.

To "install" the repository, just clone it or download and unpack a ZIP version of it (<1 min).

## Model fitting and predictions using data from the paper

Data can be downloaded from our [G-node repository](https://gin.g-node.org/gollischlab/Karamanlis_Gollisch_2023_RGC_spiketrains_natural_movies_and_subunit_models). First you need to download and unpack the datasets in the same path as the repository. Demos are expected to produce informative figures at their end and are expected to run for up to a few minutes for subunit model fitting (on a normal desktop computer with a GPU). Fitting is slower for flickering gratings.

### Flashed gratings
The script `demo_flashed_gratings.m` will guide you through data loading and model fitting. This script also uses the fitted model to generate predictions for natural images. 

You will need datasets for which `gratingflashes_data.m` is available. You can change the path of the dataset at the beginning of the script. Model components can be visualized during the fitting process by setting the parameter `options.showfig` to `True`.

### Flickering gratings
Similar to the flashed gratings, the script `demo_flickering_gratings.m` will guide you through data loading, model fitting, and generating model predictions for natural movies. 

You will need datasets for which `gratingflicker_data.m` is available. You can change the path of the dataset at the beginning of the script.

## Simulated data
Under construction...

## Using your own data
You can easily reuse the code for your own recordings, as long as the data are formatted in a similar manner. 
```
gfdata = load('gratingflashes_data')
icell  = 1 
xx     = stiminfo(presentOrder,:); %stimulus values
yy     = trialCounts(icell,:)';

```

