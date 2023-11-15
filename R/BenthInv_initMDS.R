### traits explore ###
#load packages
require(tidyverse)
require(vegan)

ggplot2::theme_set(ggthemes::theme_few())

# load data

df0 <- readxl::read_xlsx("data/BENTH_OPEN_DATA_TAXA_2023-11-15.xlsx",sheet = "out")
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

### remove trait data and widen
wbs <- c("SOUTHAMPTON WATER", "TAMAR",
         "CARRICK ROADS","POOLE HARBOUR")

dfw <- df0 %>% 
  filter(!str_detect(TAXON_GROUP_NAME, "^insect")) %>% #remove insects
  dplyr::select(c(AGENCY_AREA:SIEVE_SIZE,PREFERRED_TAXON_NAME_Broad,
                  NUMBER_FOUND)) %>% 
  filter(!is.na(NUMBER_FOUND)) %>%
  filter(., WATER_BODY %in% wbs) %>%
  group_by(AGENCY_AREA,REPORTING_AREA,
          SEA_AREA, WATERBODY_TYPE_DESCRIPTION, 
          WATER_BODY, SITE_ID,                    
          SITE_VERSION, SITE_NGR_PREFIX,            
          SITE_EASTING, SITE_NORTHING,              
          SITE_NGR_10_FIG, SITE_FULL_EASTING,          
          SITE_FULL_NORTHING, WIMS_SITE_ID,               
          WFD_WATERBODY_ID, SAMPLE_ID,                  
          SAMPLE_VERSION, REPLICATE_CODE,             
          SURVEY_CODE, SAMPLE_NGR_PREFIX,          
          SAMPLE_EASTING, SAMPLE_NORTHING,            
          SAMPLE_NGR_10_FIG, SAMPLE_FULL_EASTING,        
          SAMPLE_FULL_NORTHING, SAMPLE_DATE,                
          SAMPLE_TYPE_DESCRIPTION, SAMPLE_METHOD_DESCRIPTION,  
          SAMPLE_REASON, BENT_GRAB_DEPTH,            
          BENT_PSA, BENT_DEPTH_RPD_LAYER,       
          ANALYSIS_ID, DATE_OF_ANALYSIS,           
          ANALYSIS_TYPE_DESCRIPTION, ANALYSIS_METHOD_DESCRIPTION,
          SIEVE_SIZE, PREFERRED_TAXON_NAME_Broad) %>% 
  summarise(NUMBER_FOUND=sum(NUMBER_FOUND), .groups = "drop") %>% 
  pivot_wider(names_from = PREFERRED_TAXON_NAME_Broad,
              values_from = NUMBER_FOUND,
              values_fill = 0)

kp <- vegan::specnumber(dfw[,-c(1:37)])>5 #remove samples with <5 taxa
dfw <- dfw[kp,]

## initial MDS
# chop off metadata
dfw_ord <- dfw[,-c(1:37)]

ptm <- Sys.time()###
set.seed(pi);ord <- vegan::metaMDS(dfw_ord,
                                   autotransform = TRUE,
                                   trymax = 100)
saveRDS(ord, file = "data/out/ord3d.Rdata")
Sys.time() - ptm;rm(ptm)
plot(ord)

#### extract ordination axes ####
scores_site <- dfw %>% 
  dplyr::select(c(1:37))
tmp_sites <- as_tibble(as.data.frame(scores(ord,display = "site")))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
scores_site$NMDS1 <- tmp_sites$NMDS1 #location of individual samples in NMDS space
scores_site$NMDS2 <- tmp_sites$NMDS2 #location of individual samples in NMDS space

# Using the scores function from vegan to extract the species scores and
# convert to a data.frame
scores_species <- as.data.frame(scores(ord,display = "species"))
scores_species$lbfull <-  row.names(scores_species)
scores_species$lb <-  make.cepnames(row.names(scores_species))#shorten names

#### generate mean centroids by WB ####
scores_site %>% 
  group_by(WATER_BODY) %>%
  summarise(mn_ax1_WB=mean(NMDS1),mn_ax2_WB=mean(NMDS2)) %>%
  ungroup() -> centr

scores_site <- left_join(scores_site,centr,by="WATER_BODY");rm(centr)

ggplot()+
  geom_hline(colour="grey",yintercept = 0, lty=2)+
  geom_vline(colour="grey",xintercept = 0, lty=2)+
  geom_text(data=scores_species, aes(x = NMDS1, y=NMDS2, label=lb),
            size=3,
            alpha=0.2)+
geom_segment(data=scores_site,aes(x=NMDS1,y=NMDS2,
                                  colour=WATER_BODY,
                                  xend=mn_ax1_WB,yend=mn_ax2_WB),
             show.legend = FALSE)+
  geom_point(data=scores_site, show.legend=TRUE,
             aes(x=NMDS1, y=NMDS2,
                 fill = WATER_BODY,
                 shape = WATER_BODY),
             size=3)+
  scale_fill_manual(values=c(cbPalette))+
  scale_colour_manual(values=c(cbPalette))+
  scale_shape_manual(values = rep(c(24,25,23),each=2))+
  coord_equal()+
  theme(legend.title = element_blank(),
        axis.title = element_text(face="bold"))
