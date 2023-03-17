Harmonized data and R code for "Plankton response to global warming is characterized by non-uniform shifts in assemblage
 composition since the last ice age" by Anne Strack, Lukas Jonkers, Marina C. Rillo, Helmut Hillebrand and Michal Kucera.

Analyse planktonic foraminifera species assemblages from the North Atlantic Ocean over the past 24,000 years.

Scripts written by Anne Strack

DATA SOURCES
* WOA18: Locarnini, R. A. et al. World Ocean Atlas 2018, Volume 1: Temperature. A. Mishonov, Technical Editor. NOAA Atlas NESDIS 81, 52 (2019).
* LGMR: Osman, M. B. et al. Globally resolved surface temperatures since the Last Glacial Maximum. Nature 599, 239-244, doi:10.1038/s41586-021-03984-4 (2021).
* MARGO: Kucera, M., Rosell-Melé, A., Schneider, R., Waelbroeck, C. & Weinelt, M. Multiproxy approach for the reconstruction of the glacial ocean surface (MARGO). Quat. Sci. Rev. 24, 813-819, doi:10.1016/j.quascirev.2004.07.017 (2005). Kucera, M. et al. Reconstruction of sea-surface temperatures from assemblages of planktonic foraminifera: multi-technique approach based on geographically constrained calibration data sets and its application to glacial Atlantic and Pacific Oceans. Quat. Sci. Rev. 24, 951-998, doi:10.1016/j.quascirev.2004.07.014 (2005).
* planktonic foraminifera assemblage data: individual citations provided in CoreList_PlanktonicForaminifera.csv

DATA
1. Harmonized assemblage data*: FullDataTable_PF_harmonized.txt
2. Core list with additional information to time series: CoreList_PlanktonicForaminifera.csv
3. Reference list for PF names: ReferenceList_PlanktonicForaminifera.csv

CODE
1. 01_DataAnalysis_PCA.R: principal component analysis on assemblage data of individual time series as well as on whole dissimilarity matrix (results shown in Fig. 1 and 2)
2. 02_DataAnalysis_LocalBiodiversityChange.R: local biodiversity change analysis of individual time series (results shown in Fig. 3 and Extended Data Fig. 1); also recalculates resolution of time-series
3. 03_DataAnalysis_NoAnalogueAssemblages.R: calculates compositional dissimilarity to the nearest LGM sample to analyse existence of no-analogues (results shown in Fig. 4, as well as Extended Data Fig. 3 and 4)
4. 04_DataAnalysis_LDG_LGMresiduals.R: visualises latitudinal diversity gradient through time and the difference between richness and Shannon diversity to their respective LGM mean values (results shown in Fig. 5)

*Assemblage data of individual time series were manually downloaded, checked and harmonized following the taxonomy of
Siccha and Kucera (2017) and combined into one data file. Species not reported in the time series data were assumed to
be absent (i.e., zero abundance). We merged Globigerinoides ruber ruber and Globigerinoides ruber albus, because some 
studies only reported them together as Globigerinoides ruber. Also, P/D intergrades (an informal category of morphological
intermediates between Neogloboquadrina incompta and Neogloboquadrina dutertrei) were merged with Neogloboquadrina incompta.
In total, 41 species of planktonic foraminifera were included in this study.

Siccha, M. & Kucera, M. ForCenS, a curated database of planktonic foraminifera census counts in marine surface sediment samples. Sci. Data 4, 170109, doi:10.1038/sdata.2017.109 (2017).

