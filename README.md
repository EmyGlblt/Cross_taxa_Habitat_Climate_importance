# Cross_taxa_Habitat_Climate_importance

Code to implement the Cross taxa analyses of terrestrial animal communities changes in Finland (Hierarchical Modelling of Species Communities (HMSC) and adapted variance partition summaries). The code is linked to a paper in review.

We used Bayesian joint species distribution models and variance partition to investigate the occurrence and abundance patterns of 503 animal species, including birds, large mammals, small mammals (rodents and shrews), butterflies, and moths across a ~1,100 km latitudinal gradient in Finland. For each species group, spatial context (boreal region, within and between region), and functional group (species grouped by life history traits) we evaluate the relative and shared importance of 3 main environmental drivers: climate, habitat composition, landscape configuration.

_Disclaimer_: The code in this repository represents one version of the code developed for the project and may yet undergo changes and revisions.
_Authors_: The Hmsc code was developed through collaboration between Emy Guilbault and Mirkka Jones. The variance partitioning code was developed through collaboration between Emy Guilbault and Jarno Vanhatalo.
<br>

## Data files

The code and data used in the current analyses are available in GitHub [code available here] and Zenodo [link to be provided upon acceptance]. Only processed data for the analyses (without coordinates) are available. 
Data have been collected in systematic and organized monitoring programs lead by the Natural Resources Institute Finland and Finnish Environment Institute. Species specific monitoring programs have followed their own sampling strategy, which are described in the Materials and Methods section. Given the substantial number of sampled sites and to have manageable computational running times, we used stratified sampling to subsample representative sites in time and space for birds and mammals. Finally, we retained species with at least 40 observations in this period in our models; thus, rare species were omitted to prevent model overfitting.
<br>


### Species data

_Bird_ monitoring data is curated by the Finnish Museum of Natural History, and available through the Finnish Biodiversity Information Facility, FinBIF (https://laji.fi).

_Butterfly_ monitoring data is curated by the Finnish Environment Institute (SYKE), data can be accessed through the European butterfly transect dataset through eBMS after formally completing a license form (https://butterfly-monitoring.net/ebms-dataaccess).

_Moth_ monitoring data is curated by the Finnish Environment Institute (SYKE), and available through the Finnish Biodiversity Information Facility, FinBIF (https://laji.fi/en/observation/list?sourceId=KE.1501).

_Large mammals_ and _small mammals_ datasets belong to the Natural Resources Institute Finland (Luke), any data requests should be sent to kirjaamo@luke.fi.
<br>

### Environmental data

The *climatic* variables are available from the Finnish Meteorological Institute (https://etsin.fairdata.fi/datasets/fmi?keys=Finnish%20Meteorological%20Insitute&terms=organization_name_en.keyword&p=1&sort=best). 
Habitat and landscape variables were derived from Corine Land Cover (CLC) database (https://land.copernicus.eu/en/products/corine-land-cover). 
The method section described the preparation of environmental covariates used in the models.


## R analyses
<br>
The scripts provided in this repository cover the following:
_________________________________________________________________________
<br>

1. Model
   - Hmsc_models
     - occ
     - ab
   - Performances             
     - occ
     - ab
   - Sensitivity               

2. Postprocessing     
   - Script_VP_extract
   - VP_joint_interp
   - VP_plotting

X. Extra functions   

<br>
_Description_

<br>

1. Model ---------------------------------------------------------------------------------------
   - Hmsc_models
> [!NOTE]
> _Hmsc models run for both occurrence and abundance data for all taxa in a separate R file_
     - occ
     - ab
   - Performances              
> [!NOTE]
> _Hmsc models run for both occurrence and abundance data for all taxa in a separate R file_
     - occ
     - ab
   - Sensitivity               
> [!NOTE]
> _Contains one script related to buffer size sensitivity results investigation (Check_sensitivityResch4.R) and 2 scripts looking at the variance partitioning (VP) grouping for occurrence and abundance data._

<br>

2. Postprocessing ---------------------------------------------------------------------------------------
> [!NOTE]
> _Present the new variance partition summaries code and how to extract them from Hmsc models, plot results and perform joint variance partition interpretation._
   
   - Script_VP_extract
> [!IMPORTANT]
> _*New* variance partition summaries code (computeVarPartSummaries_111125.R), linear term correlation function (Corr_function.R), script to extract variance partition summaries adapted from Schulz et al. (2025) (Extract_NMCDP.R), conditional variance partition (Extract_VCPcond.R) and for functional traits (Extract_VP_traits.R)._
   
   - VP_joint_interp
> [!IMPORTANT]
> _*New* code to interpret the joint variation between environmental drivers provided for both occurrence and abundance data on separate scripts._
   
   - VP_plotting
> [!IMPORTANT]
> _ Provide all script to recreate the figures in the main and supplementary part of the manuscript Euler plots (figure 2), termplots (figure 2, 3 and 4), joint effects of environemntal frivers (figure 5) and variance partition summaries from the supplementary: density plots (VP_val_density.R) and driver dominance (Driver_dominance.R and Driver_dominance_cond.R)

<br>

X. Extra functions ---------------------------------------------------------------------------------------
> [!NOTE]
> _Contains scripts for data summarizing (TableCLCdata.R, TableCLCprop.R), data processing (cHmsc.R from the HMSC package), data analysing (Corr_function.R) and data plotting (Finbr_niceplot.R and function_ggplot_strip_color.R from https://stackoverflow.com/questions/53455092/r-ggplot2-change-colour-of-font-and-background-in-facet-strip)_

<br>
Analyses were performed in R software (>4.4). Specifically, we used Hmsc software (3.3).

<br>
<br>

_References_

Ovaskainen, O., & Abrego, N. (2020). Joint species distribution modelling: With applications in R. Cambridge University Press.
Schulz, T., Saastamoinen, M. and Vanhatalo, J. (2025). “ Model-Based Variance Partitioning for Statistical Ecology.” Ecological Monographs 95(1): e1646. https://doi.org/10.1002/ecm.1646

