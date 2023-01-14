# grouping EPILOG data
library(data.table)
library(tidyverse)

## read symbiosis table
source("code/data_clean/read_symbiosis_table.R")

north_atlantic <- read_csv("~/Desktop/margo_modern_NorthAtlantic.csv") %>% filter(`%(1) or Raw(2)` == 1)

north_atlantic <- north_atlantic %>% pivot_longer(cols=c(`Orbulina universa`:`Globigerinita uvula`), names_to = "Species", values_to="Relative_Abundance") %>%
  select(!c(`Coring device`, `Other ID`, `Other UnID`, `%(1) or Raw(2)`))


north_atlantic <- north_atlantic %>% mutate(Species=recode(Species, 
                                           "Globigerinoides ruber (pink)" =  "Globigerinoides ruber",
                                           "Globigerinoides ruber (white)" =  "Globigerinoides ruber",
                                           "Globigerinoides ruber total" = "Globigerinoides ruber",
                                           "Globigerinoides sacculifer w/o sac" = "Globigerinoides sacculifer",
                                           "Globigerinoides sacculifer with sac" = "Globigerinoides sacculifer",
                                           "Globigerinoides sacc total"  = "Globigerinoides sacculifer",       
                                           "Neogloboquadrina pachyderma L" = "Neogloboquadrina pachyderma",
                                           "Neogloboquadrina pachyderma R" = "Neogloboquadrina incompta",
                                           "Globorotalia truncatulinoides L" = "Globorotalia truncatulinoides",
                                           "Globorotalia truncatulinoides R" = "Globorotalia truncatulinoides",
                                           "Globorotalia truncatulinoides total" = "Globorotalia truncatulinoides",
                                          "P/D integrade + N. pachyderma R" = "N. pachyderma",
                                          "Globorotalia truncatulinoides L" = "Globorotalia truncatulinoides",
                                          "Globorotalia truncatulinoides R" = "Globorotalia truncatulinoides",
                                          "Globorotalia truncatulinoides total" = "Globorotalia truncatulinoides",
                                          "Globorotalia menardii + tumida"= "Globorotalia menardii & Globorotalia tumida",
                                          "Globorotalia menardii flexuosa"  ="Globorotalia flexuosa",    
                                          "Dentagloborotalia anfracta" = "Dentigloborotalia anfracta",
                                          "Berggrenia pumilio + T. humilis"="Turborotalita humilis & Berggrenia pumilio",
                                          "Globigerinella siphonifera (=aequilateralis)" = "Globigerinella siphonifera",
                                          "Globorotalia crassula" = "Others"
                                          )) %>%
  mutate(Species=recode(Species,   
                        "Beela digitala" = "Beella digitata",
                        "N. pachyderma" = "Neogloboquadrina pachyderma",
#                        "N. pachyderma" = "Neogloboquadrina pachyderma"
))


find_sp(north_atlantic, "Species")


north_atlantic <- merge(north_atlantic, symbiosis_tbl, by.x="Species", by.y = "Species") %>%
  select(!c("Sample depth - lower (m)","Sample depth - upper (m)", "sedimentation rate (cm/ky)", "Publication", "Total Planktics", 
            "Chronostratigraphic quality (chronozone level)", "Ocean", "Species"))

north_atlantic <- north_atlantic %>% group_by(Core, Latitude, Longitude, `Water depth (m)`, Symbiosis, Spinose) %>% 
  summarise_all(.funs = sum, na.rm=T) %>% ungroup() %>% filter(Symbiosis != "Undetermined" & Spinose != 'Undetermined')

north_atlantic <- north_atlantic %>% mutate(Relative_Abundance = Relative_Abundance/100)

write_csv(north_atlantic, "data/modern_data/margo_modern_groupped(manual).csv")

