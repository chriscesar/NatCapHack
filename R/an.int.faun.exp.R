# an.int.faun.exp.R ####
# initial exploration of intertidal infauna data

# load packages ####
ld_pkgs <- c("tidyverse", "tictoc", "gllvm", "Hmsc","vegan")
vapply(ld_pkgs, library, logical(1L),
       character.only = TRUE, logical.return = TRUE);rm(ld_pkgs)

# set universals ####
tictoc::tic.clearlog();tic("Set universals");print("set universals")
perms <- 9999 ### number of permutations to run for multivariate analyses
theme_set(ggthemes::theme_few())###set theme for all ggplot objects
nit <- 200 #number of iterations
ppi <- 300 #image resolution
# colourblind friendly colour palette (RGB values also commented)
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

cbPalette2 <- c("#646464", #100/100/100
                "#B46D00",#180/109/0
                "#2482BA",#036/130/186
                "#006C41",#000/108/065
                "#BEB210",#190/178/016
                "#004080",#000/064/128
                "#A32C00",#163/044/000
                "#9A4775"#154/071/117
)

toc(log=TRUE)

# load data ####
tic("load data");print("load data")
df_abd_0 <- as_tibble(read.csv("data/abund_all.csv", header = TRUE)) # faunal abundance data
df_traits_0 <- as_tibble(read.csv("data/species_all_traits.csv", header=TRUE)) # sp. traits
df_sites_0 <- as_tibble(read.csv("data/sites_all_folkupdated.csv", header = TRUE)) # sites data
df_env_0 <- as_tibble(readxl::read_xlsx("data/env_psa_all_NEW.xlsx",
                                        sheet = "env_psa_all",guess_max = 10000
                                        )) # environmental data
toc(log=TRUE)

# generate taxon list
# df_abd_0 %>% 
#   dplyr::select(.,TaxonName) %>% distinct(.) -> taxa
# write.csv(taxa, file = "data/taxa_raw.csv",row.names = FALSE)

# Format data ####
# Append abundance & trait data 
# df_abd_0 %>% 
#   ## remove eggs and non-counted 'bits' of taxa
#   dplyr::filter(., !Observations %in% c(
#     "Eggs", "egg", "eggs", "Epitoke", "epitoke", "zoea", "Parts", "Megalopa"
#     )) %>% 
#   left_join(x=.,y=df_traits_0, by = "TaxonName") -> dfl_abdTrt
# 
# write.csv(dfl_abdTrt, file = "data/abd_Traits.csv",row.names = FALSE)

#combine raw & 'tweaked' taxon data

taxdat <- readxl::read_xlsx("data/taxa_matched.xlsx", sheet=1)

df_abd_0 %>%
  left_join(., taxdat, by = "TaxonName") %>% 
  filter(., is.na(NOTES)) %>% 
  ## Convert to presnce/absence
  mutate(.,AbundUse = 1) %>% 
  dplyr::select(., c(Site:Mesh_Size_mm,ScientificName_accepted:Species,
                     TaxonNameUSE,AbundUse, tmp)) %>% 
  dplyr::select(.,-ScientificName_accepted) %>% 
  group_by(across(!AbundUse)) %>% 
  summarise(., AbundUse = sum(AbundUse),.groups = "drop") %>% 
  ungroup(.) -> df_abd_1

## quick ordination

# df_abd_1 %>% #names()
#   dplyr::select(., Site, Station, Date, Mesh_Size_mm,
#                 TaxonNameUSE, AbundUse) %>% 
#   group_by(across(!AbundUse)) %>% 
#   summarise(., AbundUse = sum(AbundUse),.groups = "drop") %>% 
#   pivot_wider(.,
#               names_from = TaxonNameUSE,
#               values_from = AbundUse,
#               values_fill = list(AbundUse = 0)) %>%
#   dplyr::select(.,!c(Site:Mesh_Size_mm)) -> df_tmp_w
# 
# ord1 <- metaMDS(df_tmp_w)
########################
# Look at Plymouth data only ####
## 01 create temporary variables
# df_abd_1$tmp <- paste0(df_abd_1$Site,df_abd_1$Station,df_abd_1$Date)
# df_env_0$tmp <- paste0(df_env_0$Site,df_env_0$Station,df_env_0$Date)
# df_sites_0$tmp <- paste0(df_sites_0$Site,df_sites_0$Station,df_sites_0$Date)

## 02 Pivot Wider
df_abd_1 %>%
  dplyr::select(., -c(Kingdom:Species)) %>% 
  group_by(across(!AbundUse)) %>% 
  summarise(AbundUse = sum(AbundUse),.groups = "drop") %>%
  ungroup() %>% 
  pivot_wider(., names_from = TaxonNameUSE, values_from = AbundUse,
              values_fill = 0) -> df_abd_1_w
  
df_env_0 %>% 
  dplyr::select(., tmp, PSA_Variable, Value) %>% 
group_by(across(!Value)) %>% 
  pivot_wider(.,values_from = PSA_Variable, values_from = Value)

## 03 Append BSH codes to abundance data ####
df_abd_1 %>% 
  left_join(.,df_sites_0, by = "tmp") -> dfall
  # left_join(.,df_env_0, by = "tmp") %>% 
  # dplyr::select(., -c(Site.y,Station.y,Date.y,Site,Station,
  #                     Easting.y,Northing.y, Water_Depth_m.y, Date)) %>% #names(.)
  # rename(.,Site = Site.x, Station = Station.x, Easting = Easting.x,
  #        Northing = Northing.x,Water_Depth_m = Water_Depth_m.x,
  #        Date = Date.x)

write.csv(dfall, file="data/dfall.csv", row.names = FALSE)

## 02 Extract only Plymouth data ####
dfall %>% #names(.)
  dplyr::filter(., Site == "Plymouth") %>% 
  dplyr::filter(.,Mud_Sand_Gravel)
