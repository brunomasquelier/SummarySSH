# A new method for estimating recent adult mortality from summary sibling histories
In low- and middle-income countries with limited death registration statistics, adult mortality rates are commonly estimated through sibling survival histories (SSH). In full SSH, respondents are asked about either the age, or the age and time of death, of each of their siblings in turn. Full SSH allow direct mortality estimation but can be time-consuming to collect. We provide here some R code used to develop a new indirect estimation method using summary SSH, requiring only a limited set of questions to produce recent mortality estimates. The code requires a local installation of SOCSIM, available here: https://lab.demog.berkeley.edu/socsim/. It can also be revised to use the rsocsim package in R, available here: https://github.com/MPIDR/rsocsim.

The method is detailed in 
> Bruno Masquelier, Ashira Menashe-Oren, Georges Reniers et al. A new method for estimating recent adult mortality from summary sibling histories, 07 June 2024, PREPRINT (Version 1) available at Research Square [https://doi.org/10.21203/rs.3.rs-4460057/v1]

This repository contains the following folders:
- Code: R Code to re-create the set of 192 microsimulations and to reconstruct summary SSH within these simulations as if they had been collected from adults aged 15-49. For each age group of respondents, the script allows re-calculating the coefficients that convert the proportion of adult siblings who died in the previous 5 years into age-specific mortality rates. Another file is used to evaluate the performance of this new method with real data, using 154 Demographic and Health Surveys. Users interested in applying the method as detailed in the paper without re-creating the simulations should use this second file.
- Data: All the necessary input data to create the simulations. DHS data used to evaluate the method need to be downloaded by the users, using the R package rdhs.  
- Results: CSV data files with age-specific mortality rates and estimates of 35q15 from the direct and indirect estimation.
 
