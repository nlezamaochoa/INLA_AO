# INLA_AO

This project describes the distribution of Mobula mobular bycatch species from the European observer database (AZTI, IEO, IRD) in the tropical tuna purse-seine fishery in the eastern Atlantic Ocean using a Bayesian INLA-SPDE approach during the years 2005-2015.

A Generalized Additive Model was used as the algorith to run the Bayesians approach, using a binomial family distribution and a logit function. The presence of M. mobular (473) and absences (116928) were included in the model as presence-absences.

Non-relationships were considered in the model, using a cc for Month.

The environmental variables considered in the final model were obtained from the EU Copernicus database:

Type of set, month, Chl, Ni, O2, SSH

