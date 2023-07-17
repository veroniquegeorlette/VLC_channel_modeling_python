# VLC channel simulator in Python
[![DOI](https://zenodo.org/badge/667326582.svg)](https://zenodo.org/badge/latestdoi/667326582)

This folder contains a large part of the channel simulator for a VLC (Visible Light Communication) channel simulator.
The various source codes linked to the case studies:
- Bus shelter
- Street lamp
- Industrial lighting
- Impulse response 
- 4-transmitter system
The "BUS_STATION_simul_python" folder contains 4 parts of code 
- Main_Bus_Smoke" contains the main code defining the various scenario parameters.
- fAttenuationsBus_Smoke" are the functions used to calculate the presence of smoke or fog and return the parameters to "Main".
- Plot_levelcurve_Bus_Smoke" contains all the functions for displaying the results. There are also save and display functions. 
The "ON_HOUSE_simul_python" folder contains 5 pieces of code 
- On_house_pattern" displays the IES file, which has been processed beforehand and is found in the form of the Z vector in the code.
- Main_House_Smoke" contains the main code defining the various scenario parameters.
- fAttenuationsHouse_Smoke" are the functions used to calculate the presence of smoke or fog and return the parameters to "Main".
- Plot_levelcurve_House_Smoke" contains all the functions for displaying the results. There are also save and display functions. 
The "Pattern_industry" folder contains 2 files.
- Industry_pattern" processes the test_indus.txt file to display it and calculate the optical power distribution in a plane. I recommend playing the code by section.
- test_indus.txt uses the coefficients present in the original IES file.
The "Rep_imp_bus" file contains 2 sections of code 
- Rep_imp_bus_station" introduces the scenario's main parameters and displays the curves and impulse response of the system under study. It separates and joins the contributions from smoke and fog.
- fAttenuations_Gl" is the set of functions used to calculate the attenuation for fog and smoke. 
The "4_emitters_single_receiver_LOS" code displays the optical power distribution when 4 transmitters are present in the room.