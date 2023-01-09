# Growth modeling applied to several fruit species

Daniel Jacob (INRAE BFP 2022)

<br/>

**Contents**

* **R scripts** for growth modeling on 9 fruit species as part of the ANR FRIMOUSS project (see below). 2 models based on sigmoids (single and double) were chosen. Optimization of model parameters uses the [R package minpack.lm](https://cran.r-project.org/web/packages/minpack.lm/index.html) which implements the Levenberg-Marquart nonlinear least-squares algorithm.
   * odam.R : retrieves data directly in the FRIMOUSS data collection from an ODAM server using the API (see [R package Rodam](https://github.com/inrae/Rodam/))
   * fitmodels.R : general routines for model fitting
   * growth.R : routines to interface growth modeling with FRIMOUSS collection data

* 2 **Jupyter notebooks** implementing growth modeling
   * Growth_model_nb1.ipynb : comparison of the growth modeling based on both models for one species
   * Growth_model_nb2.ipynb : comparison of the growth modeling based on the second model for all species 

---

**Frimouss project**: FRuit Integrative MOdelling for a Unified Selection System

  * **ANR Project ID**: [ANR-15-CE20-0009](http://www.agence-nationale-recherche.fr/Project-ANR-15-CE20-0009)

  * **Publication** : Léa Roch, Sylvain Prigent, Holger Klose, Coffi-Belmys Cakpo, Bertrand Beauvoit, et al.. Biomass composition explains fruit relative growth rate and discriminates climacteric from non-climacteric species. Journal of Experimental Botany, Oxford University Press (OUP), 2020, 71 (19), pp.5823-5836. [doi:10.1093/jxb/eraa302](https://academic.oup.com/jxb/article/71/19/5823/5864020)

**Frimouss dataset interfaced by ODAM**

  * All the data concerning each fruit species have been organized and managed in the same way using the ODAM approach and the associated tools.

  * ODAM Dataexplorer : https://pmb-bordeaux.fr/dataexplorer/?dc=Frimouss

  * ODAM online documentation: https://inrae.github.io/ODAM/


