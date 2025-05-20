# Prev_surface

This repository presents a methodological framework for conducting Bayesian geostatistical analysis of disease prevalence using simulated data. It outlines the key analytical steps required to integrate spatial covariates, explore their relationships with a health outcome, and fit a geostatistical model that accounts for spatial dependence. Although the data used here are synthetic and simplified, the structure is designed to be broadly applicable to epidemiological research questions where spatial heterogeneity in disease risk is hypothesized.

**What the project does**

Illustrates how to derive spatial covariates—such as distance to land use features—and incorporate them into a geostatistical modeling framework.

Applies exploratory analysis, empirical semivariogram modeling, and Bayesian inference to estimate spatial structure and generate predictions across a continuous surface.
Includes full spatial prediction with visualization tools to support interpretation of spatial patterns in prevalence.

**Why the project is useful**

Serves as a flexible template for researchers seeking to apply geostatistical methods to their own prevalence data and spatial covariates.

Enables investigation of how spatially explicit environmental factors correlate with disease risk while accounting for spatial autocorrelation.

Promotes reproducibility through transparent and adaptable code using standard tools in the R statistical environment.

**Important note on simulated data**

The analysis is based on simplified simulated data to clarify methodological steps. However, this synthetic structure lacks the complexity typically found in real-world epidemiological datasets.

Users applying these methods to empirical data—particularly datasets with meaningful spatial autocorrelation—can expect more robust and informative model outputs. The model’s capacity to capture spatial structure is constrained by the simplicity of the dummy data, and real data will likely yield more realistic and insightful spatial predictions.

**NOTE ON DATA**

I have created simulated data for this work through. ideally for the final linear.model.MLE() model (line 439) the input variables (partial sill, range) are taken from the variogram completed upon the residuals of the glm() on line 391-401 which. i beleive this error is due to the simplicity of my simulated spatial data and not something one is likely to encounter in real world datasets. However i am still trying to sort it and create a simulated dataset that will allow this analysis. I insteasd use a partial sill and range that allowed the final model to work and used to produce a surface later. one will need to amend these accordingly to work in other analysis. 
