# MultAnalysis.R
### load packages ####
ld_pkgs <- c("tidyverse","mvabund")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE)
rm(ld_pkgs)

source("R/datfol.R")
ggplot2::theme_set(ggthemes::theme_few())
ppi <- 300
cbPalette <- c("#999999", #153/153/153
               "#E69F00",#230/159/000
               "#56B4E9",#086/180/233
               "#CC79A7", #204/121/167
               "#009E73",#000/158/115
               "#F0E442",#240/228/066
               "#0072B2",#000/114/178
               "#D55E00",#213/094/000
               
               "#444444", 
               "#C34D55",
               "#33A2C4",
               "#554C31",
               "#C5C221",
               "#5531A1",
               "#B32C55",
               "#BB3593" 
)

dir(datfol)

# load taxon data ####
df_tx_l <- readxl::read_xlsx(paste0(datfol,"inf_ts_longRAW_USE.xlsx"),
                             sheet="dat_all")

## widen and throw away silly variables.  Calculate mean across replicates
df_tx_l %>% 
  select(.,year:abundance, taxonUSE) %>%
  select(., -rep, -zone2.1, -zone2.2, -yr.trn.sh.meth.rep,-taxonReported,
         -yr.trn,-yr.trn.sh,-yr.trn.sh.meth) %>% 
  group_by(year, transect, shore, zone1, mesh,taxonUSE) %>% 
  # group_by(year, transect, shore, zone1,mesh,yr.trn,
  #          yr.trn.sh,yr.trn.sh.mesh, taxonUSE) %>% 
  summarise(abundance = mean(abundance, na.omit=TRUE),
            .groups = "drop") %>% 
  pivot_wider(names_from = taxonUSE, values_from = abundance,
              values_fill = 0) %>% 
  ungroup() -> df_tx_w
View(df_tx_w)

# load sediment data ####
df_sed_l <- readxl::read_xlsx(paste0(datfol,"sed.psa.bulkWIP_use.xlsx"),
                             sheet=1)

df_sed_l %>% 
  select(.,c(year:zone1,"SAMPLE TYPE":CLAY_perc))
         