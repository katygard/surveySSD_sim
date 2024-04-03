# surveySSD_sim
Public repository of the code and saved data outputs for Gardner and Hawkins 2024.

simulation_SSD.R contains all the code to generate the simulated macroinvertebrate community and run the SSDs for part 1 of manuscript. 
NRSA_analyses_SSD.R contains downloading the NRSA data files from the USEPA's website, as well as subsetting it and running the analyses found in the manuscript.
rarify.R can be sourced from John Van Sickle, USEPA, but is saved here for ease of access.
post_bootstrap.RData is a saved data file after line 373 in simulation_SSD.R. This can be used to save time, as well as replicate Figure 3, which needs the exact same species names to be present to work.
