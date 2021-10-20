# Unmeasured Confounding in High-Dimensional Data

Unmeasured confounding is a strong limitation to drawing causal inferences and, in general, a key threat to observational studies.

This repository contains an [essay](https://github.com/bglbrt/UCHDD/blob/main/essay.pdf), code and some example data which explore multiple novel approaches to overcome this problem.

## Data

All data was either simulated or sourced on the Urban Statistics section of the [European Data Portal](https://ec.europa.eu/eurostat/). For the most part, we only use data on German cities and urban areas that were collected in 2011 and which record many characteristics ranging over multiple themes covering demographics, economic activity, environment, urban planning, etc.

## Required libraries

The algorithms and examples detailed in the essay require the following libraries:
 - in Python
	 - [NumPy](https://numpy.org)
	 - [Pandas](https://pandas.pydata.org)
 - in R
	 - [MASS](https://cran.r-project.org/web/packages/MASS/index.html)
	 - [glmnet](https://cran.r-project.org/web/packages/glmnet/index.html)
	 - [corpcor](https://cran.r-project.org/web/packages/corpcor/index.html)
	 - [cate](https://cran.r-project.org/web/packages/cate/index.html)
