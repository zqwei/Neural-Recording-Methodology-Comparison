# Neural Recording Methodology Comparison
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3979786.svg)](https://doi.org/10.5281/zenodo.3979786)


## Dataset list
Beside the existing datasets from Svoboda Lab (see description for details), the 6f-TG imaging data were recorded by Kayvon Daie and their rois were mannually selected by Ziqiang Wei. 

### Dataset description and resources
A description of the data can be found at http://im-phys.org/data and documented breifly in file https://ndownloader.figshare.com/files/24190289 

### Preprocessing
* removing the units with low number of trials (*< 20*) in either trial type condition
* spiking dataset -- removing low firing units; removing low ROC units
* spiking dataset -- maximum intersection after SD prepocessing
* Ca++ imaging dataset -- removing low ROC units

### Precompiled data
A precompiled version of full datasets can be downloaded from https://ndownloader.figshare.com/articles/12786296/versions/1 , which can be used directly in the codes. In order to make the code functionally, please copy the data to folder `TempDat` in the root folder of the project.

Citation: Wei, Ziqiang; Lin, Bei-Jung; Chen, Tsai-Wen; Daie, Kayvon; Svoboda, Karel; Druckmann, Shaul (2020): Precompiled data for "A comparison of neuronal population dynamics measured with calcium imaging and electrophysiology". figshare. Dataset. https://doi.org/10.6084/m9.figshare.12786296.v1

### Raw data
The existing data can be found from CRCNS.org, while we release the new 6f-TG raw data at https://doi.org/10.6084/m9.figshare.12792587

Citation: Wei, Ziqiang; Lin, Bei-Jung; Chen, Tsai-Wen; Daie, Kayvon; Svoboda, Karel; Druckmann, Shaul (2020): Raw data for "A comparison of neuronal population dynamics measured with calcium imaging and electrophysiology". figshare. Dataset. https://doi.org/10.6084/m9.figshare.12792587

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

## Analyses in Wei et al., 2020
### Single neuron analyses
* dynamical selectivity (Figure 2A-I)
    * non-selective
    * monophaisc selective
    * multiphasic selective
    * code is named as `switchSelectivityDistByTime.m`
* Selectivity conservation (Figure 2J, 4C-E)
    * ramp-up dynamics
    * ramp-down dynamics
    * code could be sent upon request (raw functions for the plots can be found in `Func`, but not compiled in api)
* Peakiness (Figure 7)
    * peak activity
    * code is named as `rawActivityPositivePeak.m` (Figure 7 raster plots)
    * code is named as `peakLocationDistrPosNeuron.m` (Figure 7 statistical plots; histogram)

### Population analyses
* PCA contents (Figure 5)
    * trial type
    * time
    * other
    * code is named as `meanPCA.m`
* Coding directions (Figure 6)
    * LDA -- trial type
    * LDA -- behavioral epochs
    * code for trial type is named as `trialTypeLDANumExampleTrial.m`
    * code for behavioral epochs could be sent upon request (raw functions for the plots can be found in `Func`, but not compiled in api)

## Analyses beyond Wei et al., 2020
All analyses listed here can be sent out upon request.

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

### Contact
weiz (at) janelia (dot) hhmi (dot) org
