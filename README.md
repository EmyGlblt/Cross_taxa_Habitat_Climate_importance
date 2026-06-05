# Cross_taxa_Habitat_Climate_importance

Code to implement the cross taxa analyses of terrestrial animal communities changes in Finland (Hierarchical Modelling of Species Communities (HMSC) and adapted variance partition summaries). The code is linked to a paper in review.

We used Bayesian joint species distribution models and variance partition to investigate the occurrence and abundance patterns of 503 animal species, including birds, large mammals, small mammals (rodents and shrews), butterflies, and moths across a ~1,100 km latitudinal gradient in Finland. For each species group, spatial context (boreal region, within and between region), and functional group (species grouped by life history traits) we evaluate the relative and shared importance of 3 main environmental drivers: climate, habitat composition, landscape configuration.

_Disclaimer_: The code in this repository represents one version of the code developed for the project and may yet undergo changes and revisions.
_Authors_: The Hmsc code was developed through collaboration between Emy Guilbault and Mirkka Jones. The variance partitioning code was developed through collaboration between Emy Guilbault and Jarno Vanhatalo.
<br>

## Data files

The code and data used in the current analyses are available in GitHub [code available here] and Zenodo [https://doi.org/10.5281/zenodo.17952060]. Only processed data for the analyses (without coordinates) are available. 
Data have been collected in systematic and organized monitoring programs lead by the Natural Resources Institute Finland and Finnish Environment Institute. Species specific monitoring programs have followed their own sampling strategy, which are described in the Materials and Methods section. Given the substantial number of sampled sites and to have manageable computational running times, we used stratified sampling to subsample representative sites in time and space for birds and mammals. Finally, we retained species with at least 40 observations in this period in our models; thus, rare species were omitted to prevent model overfitting.
<br>
Location:  Finland, Northern Europe.

Time period: 1999-2019


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

### Description

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
> _*New* code to interpret the joint variation between environmental drivers provided for both occurrence and abundance data._
   
   - VP_plotting
> [!IMPORTANT]
> _ Provide all script to recreate the figures in the main and supplementary part of the manuscript Euler plots (figure 2), termplots (figure 2, 3 and 4), joint effects of environemntal frivers (figure 5) and variance partition summaries from the supplementary: density plots (VP_val_density.R) and driver dominance (Driver_dominance.R and Driver_dominance_cond.R)._

<br>

X. Extra functions ---------------------------------------------------------------------------------------
> [!NOTE]
> _Contains scripts for data summarizing (TableCLCdata.R, TableCLCprop.R), data processing (cHmsc.R from the HMSC package), data analysing (Corr_function.R) and data plotting (Finbr_niceplot.R and function_ggplot_strip_color.R from https://stackoverflow.com/questions/53455092/r-ggplot2-change-colour-of-font-and-background-in-facet-strip)._

<br>
<br>

_References_

Ovaskainen, O., & Abrego, N. (2020). Joint species distribution modelling: With applications in R. Cambridge University Press.

Schulz, T., Saastamoinen, M. and Vanhatalo, J. (2025). “ Model-Based Variance Partitioning for Statistical Ecology.” Ecological Monographs 95(1): e1646. https://doi.org/10.1002/ecm.1646


## R analyses

<br>
Hmsc models were run on HPC (csc):
<br>

| Module name (R version) | CRAN package dating | Bioconductor version | RStudio Server version | oneMKL version | Cmdstan version |
|:-----------------------:|---------------------|:--------------------:|:----------------------:|:--------------:|:---------------:|
|        r-env/440        |     May 15 2024     |         3.19         |      2024.04.0-735     |    2024.1.0    |      2.35.0     |
|        r-env/432        |     Jan 15 2024     |         3.18         |      2023.12.0-369     |    2024.0.0    |      2.34.1     |

<br>

Post analyses were performed in R software (>4.4). 


In addition, we used the following R packages: abind v. 1.4.8 (Plate and Heiberger 2024), class v. 7.3.23 (Venables and Ripley 2002), colorblindr v. 0.1.0 (McWhite and Wilke 2026), corrplot v. 0.95 (Wei and Simko 2024), cowplot v. 1.2.0 (Wilke 2025a), eulerr v. 7.0.4 (Larsson and Gustafsson 2018), formattable v. 0.2.1 (Ren and Russell 2021), geomtextpath v. 0.2.0 (Cameron and van den Brand 2025), ggeasy v. 0.1.6 (Carroll, Schep, and Sidi 2025), ggExtra v. 0.11.0 (Attali and Baker 2025), ggpubr v. 0.6.3 (Kassambara 2026), ggridges v. 0.5.7 (Wilke 2025b), ggtern v. 4.0.0 (Hamilton and Ferry 2018), grDevices v. 4.5.3 (R Core Team 2026a), grid v. 4.5.3 (R Core Team 2026b), gridExtra v. 2.3 (Auguie 2017), gtable v. 0.3.6 (Wickham and Pedersen 2024), Hmsc v. 3.3.7 (Tikhonov et al. 2025), kableExtra v. 1.4.0 (Zhu 2024), knitr v. 1.51 (Xie 2014, 2015, 2025), mapdata v. 2.3.1 (Richard A. Becker and Ray Brownrigg. 2022), maps v. 3.4.3 (Becker et al. 2025), mapsFinland v. 0.2.2 (Haukka, 2020), maptools v. 1.1-8 (Bivand, 2003) - before archived, patchwork v. 1.3.2 (Pedersen 2025), pdftools v. 3.9.0 (Ooms 2026), png v. 0.1.9 (Urbanek 2026), quantmod v. 0.4.28 (Ryan and Ulrich 2025), RColorBrewer v. 1.1.3 (Neuwirth 2022), remotes v. 2.5.0 (Csárdi et al. 2024), rempsyc v. 0.2.0 (Thériault 2023), reshape2 v. 1.4.5 (Wickham 2007), scales v. 1.4.0 (Wickham, Pedersen, and Seidel 2025), sf v. 1.1.0 (E. Pebesma 2018; E. Pebesma and Bivand 2023), sjPlot v. 2.9.0 (Lüdecke 2025), sp v. 2.2.1 (E. J. Pebesma and Bivand 2005; Bivand, Pebesma, and Gomez-Rubio 2013), tidyverse v. 2.0.0 (Wickham et al. 2019), wesanderson v. 0.3.7 (Ram and Wickham 2023).

<br>

_Package citations_


> Attali, Dean, and Christopher Baker. 2025. ggExtra: Add Marginal Histograms to “ggplot2,” and More “ggplot2” Enhancements. https://doi.org/10.32614/CRAN.package.ggExtra.
>
> 
> Auguie, Baptiste. 2017. gridExtra: Miscellaneous Functions for “Grid” Graphics. https://doi.org/10.32614/CRAN.package.gridExtra.
> 
> Becker, Richard A., Allan R. Wilks, Ray Brownrigg, Thomas P. Minka, and Alex Deckmyn. 2025. maps: Draw Geographical Maps. https://doi.org/10.32614/CRAN.package.maps.
>
> 
> Bivand, Roger S., Edzer Pebesma, and Virgilio Gomez-Rubio. 2013. Applied Spatial Data Analysis with R, Second Edition. Springer, NY. https://asdar-book.org/.
> Bivand R. (2003). maptools: Tools for Handling Spatial Objects. R package version 1.1-8, https://cran.r-project.org/web/packages/maptools.
> Cameron, Allan, and Teun van den Brand. 2025. geomtextpath: Curved Text in “ggplot2”. https://doi.org/10.32614/CRAN.package.geomtextpath.
> Carroll, Jonathan, Alicia Schep, and Jonathan Sidi. 2025. ggeasy: Easy Access to “ggplot2” Commands. https://doi.org/10.32614/CRAN.package.ggeasy.
> Csárdi, Gábor, Jim Hester, Hadley Wickham, Winston Chang, Martin Morgan, and Dan Tenenbaum. 2024. remotes: R Package Installation from Remote Repositories, Including “GitHub”. https://doi.org/10.32614/CRAN.package.remotes.
> Hamilton, Nicholas E., and Michael Ferry. 2018. “ggtern: Ternary Diagrams Using ggplot2.” Journal of Statistical Software, Code Snippets 87 (3): 1–17. https://doi.org/10.18637/jss.v087.c03.
> Haukka, J. (2020). mapsFinland: Maps of Finland. R package version 0.2.2, https://cran.r-project.org/web/packages/mapsFinland.
> Kassambara, Alboukadel. 2026. ggpubr: “ggplot2” Based Publication Ready Plots. https://doi.org/10.32614/CRAN.package.ggpubr.
> Larsson, Johan, and Peter Gustafsson. 2018. “A Case Study in Fitting Area-Proportional Euler Diagrams with Ellipses Using Eulerr.” In Proceedings of International Workshop on Set Visualization and Reasoning, 2116:84–91. Edinburgh, United > Kingdom: CEUR Workshop Proceedings. https://ceur-ws.org/Vol-2116/paper7.pdf.
> Lüdecke, Daniel. 2025. sjPlot: Data Visualization for Statistics in Social Science. https://CRAN.R-project.org/package=sjPlot.
> McWhite, Claire D., and Claus O. Wilke. 2026. colorblindr: Simulate Colorblindness in r Figures. https://github.com/clauswilke/colorblindr.
> Neuwirth, Erich. 2022. RColorBrewer: ColorBrewer Palettes. https://doi.org/10.32614/CRAN.package.RColorBrewer.
> Ooms, Jeroen. 2016. pdftools: Text Extraction, Rendering and Converting of PDF Documents. https://doi.org/10.32614/CRAN.package.pdftools.
> Pebesma, Edzer. 2018. “Simple Features for R: Standardized Support for Spatial Vector Data.” The R Journal 10 (1): 439–46. https://doi.org/10.32614/RJ-2018-009.
> Pebesma, Edzer J., and Roger Bivand. 2005. “Classes and Methods for Spatial Data in R.” R News 5 (2): 9–13. https://CRAN.R-project.org/doc/Rnews/.
> Pebesma, Edzer, and Roger Bivand. 2023. Spatial Data Science: With applications in R. Chapman and Hall/CRC. https://doi.org/10.1201/9780429459016.
> Pedersen, Thomas Lin. 2025. patchwork: The Composer of Plots. https://doi.org/10.32614/CRAN.package.patchwork.
> Plate, Tony, and Richard Heiberger. 2024. abind: Combine Multidimensional Arrays. https://doi.org/10.32614/CRAN.package.abind.
> R Core Team. 2026. R: A Language and Environment for Statistical Computing. Vienna, Austria: R Foundation for Statistical Computing. https://www.R-project.org/.
> Ram, Karthik, and Hadley Wickham. 2023. wesanderson: A Wes Anderson Palette Generator. https://doi.org/10.32614/CRAN.package.wesanderson.
> Ren, Kun, and Kenton Russell. 2021. formattable: Create “Formattable” Data Structures. https://doi.org/10.32614/CRAN.package.formattable.
> Richard A. Becker, Original S code by, and Allan R. Wilks. R version by Ray Brownrigg. 2022. mapdata: Extra Map Databases. https://doi.org/10.32614/CRAN.package.mapdata.
> Ryan, Jeffrey A., and Joshua M. Ulrich. 2007. quantmod: Quantitative Financial Modelling Framework. https://doi.org/10.32614/CRAN.package.quantmod.
> Thériault, Rémi. 2023. “rempsyc: Convenience Functions for Psychology.” Journal of Open Source Software 8 (87): 5466. https://doi.org/10.21105/joss.05466.
> Tikhonov, Gleb, Otso Ovaskainen, Jari Oksanen, Melinda de Jonge, Oystein Opedal, and Tad Dallas. 2025. Hmsc: Hierarchical Model of Species Communities. https://doi.org/10.32614/CRAN.package.Hmsc.
> Urbanek, Simon. 2010. png: Read and Write PNG Images. https://doi.org/10.32614/CRAN.package.png.
> Venables, W. N., and B. D. Ripley. 2002. Modern Applied Statistics with s. Fourth. New York: Springer. https://www.stats.ox.ac.uk/pub/MASS4/.
> Wei, Taiyun, and Viliam Simko. 2024. R Package “corrplot”: Visualization of a Correlation Matrix. https://github.com/taiyun/corrplot.
> Wickham, Hadley. 2007. “Reshaping Data with the reshape Package.” Journal of Statistical Software 21 (12): 1–20. https://www.jstatsoft.org/v21/i12/.
> Wickham, Hadley, Mara Averick, Jennifer Bryan, Winston Chang, Lucy D’Agostino McGowan, Romain François, Garrett Grolemund, et al. 2019. “Welcome to the tidyverse.” Journal of Open Source Software 4 (43): 1686. https://doi.org/10.21105/joss.01686.
> Wickham, Hadley, and Thomas Lin Pedersen. 2024. gtable: Arrange “Grobs” in Tables. https://doi.org/10.32614/CRAN.package.gtable.
> Wickham, Hadley, Thomas Lin Pedersen, and Dana Seidel. 2025. scales: Scale Functions for Visualization. https://doi.org/10.32614/CRAN.package.scales.
> Wilke, Claus O. 2015. cowplot: Streamlined Plot Theme and Plot Annotations for “ggplot2”. https://doi.org/10.32614/CRAN.package.cowplot.
> Wilke, Claus O. 2017. ggridges: Ridgeline Plots in “ggplot2”. https://doi.org/10.32614/CRAN.package.ggridges.
> Xie, Yihui. 2014. “knitr: A Comprehensive Tool for Reproducible Research in R.” In Implementing Reproducible Computational Research, edited by Victoria Stodden, Friedrich Leisch, and Roger D. Peng. Chapman; Hall/CRC.
> Xie, Yihui. 2015. Dynamic Documents with R and Knitr. 2nd ed. Boca Raton, Florida: Chapman; Hall/CRC. https://yihui.org/knitr/.
> Xie, Yihui. 2025. knitr: A General-Purpose Package for Dynamic Report Generation in R. https://yihui.org/knitr/.
> Zhu, Hao. 2024. kableExtra: Construct Complex Table with “kable” and Pipe Syntax. https://doi.org/10.32614/CRAN.package.kableExtra.
