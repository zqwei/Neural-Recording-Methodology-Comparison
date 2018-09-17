# Neural Recording Methodology Comparison

## Dataset list

### Simultaneous ephys and imaging
* GCaMP6s/f - AAV
* GCaMP6s/f - Transgenic

### Short delay
* Spiking neurons (99 sessions; from Li et al. Nature, 2015)
* GCaMP6s Ca++ imaging __AAV__
* GCaMP6s Ca++ imaging __Transgenic GP4.3__
* GCaMP6f Ca++ imaging (1 session) __Transgenic GP5.17__
* GCaMP6f Ca++ imaging (22 sessions) __Transgenic GP5.17__
* Spiking neurons
* Whole-cell intracellular recording

### Long delay
* Spiking neurons
* GCaMP6f Ca++ imaging (3 sessions) __Transgenic GP5.17__
* GCaMP6s Ca++ imaging (3 sessions) __Transgenic GP4.3__

### Preprocessing
* removing the units with low number of trials (*< 20*) in either trial type condition
* spiking dataset -- removing low firing units; removing low ROC units
* spiking dataset -- maximum intersection after SD prepocessing
* Ca++ imaging dataset -- removing low ROC units


## Models
### Spike to image models
* Linear model
* Hill function model
* Sigmoid function model

### Image to spike models
* C2S Nonnegative Wienner Filter
* C2S FOOPSI
* C2S Finite Rate of Innovation
* C2S Constrained OOPSI AR(n)
* C2S Constrained OOPSI MCMC
* C2S Peel Linear

## Analyses
### Single neuron analyses
* dynamical selectivity
    * non-selective
    * monophaisc selective
    * multiphasic selective
* Selectivity conservation
    * ramp-up dynamics
    * ramp-down dynamics
* Peakiness
    * peak activity

### Population analyses
* PCA contents
    * trial type
    * time
    * other
* Coding directions
    * LDA -- trial type
    * LDA -- behavioral epochs

## Future analyses
### Single neuron analyses
* Activity distribution across trial periods (pre-sample, sample, delay, response)
* Dynamics of neuronal modulation
* k-means
* Selectivity over time (two sample z-score no sorting; sorting)
* Single neuron choice probability index distribution (ROC), across trial periods
* Mixed selectivity (ANOVA)
* Kurtosis distribution analysis
* Fano factors
* Peakiness
    * center of mass
    * saturation point

### Population analyses
* Noise correlations
* Decoding which trial period one is in, how fast can you tell the data that the trial period switched
    * three LDA decoders for adjacent epochs
    * single LDA decoder with time taper
* Variance-Information directions
    * LDA-LDA
    * PCA-PCA
    * LDA-PCA
* dPCAa -- rotation of the PCA components to capture the variances exclusively
for temporary dynamics or choice selection (code from *C.K. Machens*)
* GPFA -- Time-series analysis (GPFA code from *B. Yu*)
* Tensor decomposition

## Citations
### Papers

### Website
http://im-phys.org/

## Questions, comments, issues
