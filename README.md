# subunit_grid_model
Functions to fit subunit grid models to spiking data from retinal ganglion cells.

## Requirements

## Data for fitting

### Data from the paper

These can be downloaded from the repository https://gin.g-node.org/gollischlab/Karamanlis_Gollisch_2023_RGC_spiketrains_natural_movies_and_subunit_models 


#### Flashed gratings
```

gfdata = load(datapath)
icell  = 1 
xx     = stiminfo(presentOrder,:); %stimulus values
yy     = trialCounts(icell,:)';

```
#### Flickering gratings


### Your own data

### Run simulation

## Model fitting

Specify fitting parameters and run the ADAM-based optimizers. Model components can be visualized during the fitting process by setting the parameter `options.showfig` to `True`.

### Flashed gratings


```
let message = 'Hello world';
alert(message);
```


### Flickering gratings


```
let message = 'Hello world';
alert(message);
```


## Model predictions

### Natural images

### Natural movies

