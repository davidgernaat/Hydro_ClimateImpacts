# Hydro_ClimateImpacts
Scripts used to generate results as published in the journal Nature Climate Change titled 'Climate change impacts on renewable energy supply'.



ISIMIP calculation steps with scripts

1.	Downloaded data from ISIMIP server

---- data_prep/ISIMIP

2.	Convert_RCPx_mms_to_m3s(_2050)(_hist): converts it to m3/s averages for 30yrs
3.	downscaleR_monthly_ISIMIP(_2050).m: downscales global 0.5x0.5deg runoff data to 15s continental data per month
4.	DownscaleQ_ISIMIP: routes the runoff to discharge maps

---- Analysis/ISIMIP

5.	ISIMIP_GetNewQs(_Exis).m: Saves the new discharge values based on the lat/lon locations from a Hydrus run and the downscaled discharge maps
6.	ISIMIP_recalc(_Exis): Takes the new discharge values and recalculates the energy potential
