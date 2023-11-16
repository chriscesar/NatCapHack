### traits explore ###
#load packages
require(tidyverse)
require(vegan)
require(ggtext)
require(patchwork)

ggplot2::theme_set(ggthemes::theme_few())

# load data
ppi <- 300
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

# ptm <- Sys.time()###
# set.seed(pi);ordAbund <- vegan::metaMDS(dfw_ord,
#                                    autotransform = TRUE,
#                                    trymax = 20)
# saveRDS(ordAbund, file = "data/out/ordAbund.Rdata")
# Sys.time() - ptm;rm(ptm)
ordAbund <- readRDS("data/out/ordAbund.Rdata")
plot(ordAbund)

#### extract ordination axes ####
scores_siteAbund <- dfw %>% 
  dplyr::select(c(1:37))
tmp_sites <- as_tibble(as.data.frame(scores(ordAbund,display = "site")))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
scores_siteAbund$NMDS1 <- tmp_sites$NMDS1 #location of individual samples in NMDS space
scores_siteAbund$NMDS2 <- tmp_sites$NMDS2 #location of individual samples in NMDS space

# Using the scores function from vegan to extract the species scores and
# convert to a data.frame
scores_speciesAbund <- as.data.frame(scores(ordAbund,display = "species"))
scores_speciesAbund$lbfull <-  row.names(scores_speciesAbund)
scores_speciesAbund$lb <-  make.cepnames(row.names(scores_speciesAbund))#shorten names

#### generate mean centroids by WB ####
scores_siteAbund %>% 
  group_by(WATER_BODY) %>%
  summarise(mn_ax1_WB=mean(NMDS1),mn_ax2_WB=mean(NMDS2)) %>%
  ungroup() -> centrAbund

scores_siteAbund <- left_join(scores_siteAbund,centrAbund,
                              by="WATER_BODY");rm(centrAbund)

scores_siteAbund$WBlbl <- ifelse(scores_siteAbund$WATER_BODY=="CARRICK ROADS",
                            "Carrick Rd",
                            ifelse(scores_siteAbund$WATER_BODY=="TAMAR",
                                   "Tamar",
                                   ifelse(scores_siteAbund$WATER_BODY=="SOUTHAMPTON WATER",
                                          "Soton W",
                                          ifelse(scores_siteAbund$WATER_BODY=="POOLE HARBOUR",
                                                 "Poole",NA))))

ggplot()+
  geom_hline(colour="grey",yintercept = 0, lty=2)+
  geom_vline(colour="grey",xintercept = 0, lty=2)+
  geom_text(data=scores_speciesAbund,
            aes(x = NMDS1, y=NMDS2, label=lb),
            size=3,
            alpha=0.2)+
geom_segment(data=scores_siteAbund,aes(x=NMDS1,y=NMDS2,
                                  colour=WATER_BODY,
                                  xend=mn_ax1_WB,yend=mn_ax2_WB),
             show.legend = FALSE)+
  geom_point(data=scores_siteAbund, show.legend=TRUE,
             aes(x=NMDS1, y=NMDS2,
                 fill = WATER_BODY,
                 shape = WATER_BODY),
             size=3)+
  scale_fill_manual(values=c(cbPalette))+
  scale_colour_manual(values=c(cbPalette))+
  scale_shape_manual(values = rep(c(24,25,23),each=2))+
  geom_textbox(size=3,data=scores_siteAbund,aes(x=mn_ax1_WB,
                                           y=mn_ax2_WB,
                                           label=WBlbl,
                                           fill=WATER_BODY),
               width = unit(0.15, "npc"),
               inherit.aes = FALSE,show.legend = FALSE,hjust=0.5)+
  coord_equal()+
  labs(title = "NMDS: Taxon data")+
  theme(legend.title = element_blank(),
        legend.position="none",
        axis.title = element_text(face="bold")) -> plabund
ggsave(filename = "figs/MDS_by_tax_abund.pdf",width = 12,height = 12,
       units = "in",
       plot=plabund)

# weight trait data by abundance ####

#steps#
# remove taxon names
# convert to long
# multiply abundance by trait values
# group by samples and sum traits

dfTrt <- df0 %>% 
  filter(!str_detect(TAXON_GROUP_NAME, "^insect")) %>% #remove insects
  dplyr::select(c(AGENCY_AREA:SIEVE_SIZE,
                  NUMBER_FOUND, sr_Less_than_10:b_None)) %>% 
  filter(!is.na(NUMBER_FOUND)) %>%
  filter(., WATER_BODY %in% wbs) %>%
  pivot_longer(cols = sr_Less_than_10:b_None,
                names_to = "trait", values_to = "score") %>% 
  filter(!is.na(score)) %>% 
  mutate(wt_trt = NUMBER_FOUND*score) %>% 
  select(-NUMBER_FOUND, -score) %>% 
  group_by(
   AGENCY_AREA, REPORTING_AREA,             
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
   SIEVE_SIZE, trait
   ) %>% 
  summarise(wt_trt=sum(wt_trt), .groups = "drop") %>% 
  pivot_wider(names_from = trait,
              values_from = wt_trt,
              values_fill = 0)

## remove traitless rows
sms <- rowSums(dfTrt[,-c(1:37)])!=0
dfTrt <- dfTrt[sms,]

## initial MDS
# chop off metadata
dfTrt_ord <- dfTrt[,-c(1:37)]

# ptm <- Sys.time()###
# set.seed(pi);ordTrt <- vegan::metaMDS(dfTrt_ord,
#                                    autotransform = TRUE,
#                                    trymax = 20)
# saveRDS(ordTrt, file = "data/out/ordTrt.Rdata")
# Sys.time() - ptm;rm(ptm)
ordTrt <- readRDS("data/out/ordTrt.Rdata")
plot(ordTrt)

#### extract ordination axes ####
scores_site <- dfTrt %>% 
  dplyr::select(c(1:37))
tmp_sites <- as_tibble(as.data.frame(scores(ordTrt,display = "site")))  #Using the scores function from vegan to extract the site scores and convert to a data.frame
scores_site$NMDS1 <- tmp_sites$NMDS1 #location of individual samples in NMDS space
scores_site$NMDS2 <- tmp_sites$NMDS2 #location of individual samples in NMDS space

# Using the scores function from vegan to extract the species scores and
# convert to a data.frame
scores_species <- as.data.frame(scores(ordTrt,display = "species"))
scores_species$lbfull <-  row.names(scores_species)
scores_species$lb <-  make.cepnames(row.names(scores_species))#shorten names

#### generate mean centroids by WB ####
scores_site %>% 
  group_by(WATER_BODY) %>%
  summarise(mn_ax1_WB=mean(NMDS1),mn_ax2_WB=mean(NMDS2)) %>%
  ungroup() -> centr

scores_site <- left_join(scores_site,centr,by="WATER_BODY");rm(centr)

scores_site$WBlbl <- ifelse(scores_site$WATER_BODY=="CARRICK ROADS",
                            "Carrick Rd",
                            ifelse(scores_site$WATER_BODY=="TAMAR",
                                   "Tamar",
                                   ifelse(scores_site$WATER_BODY=="SOUTHAMPTON WATER",
                                          "Soton W",
                                          ifelse(scores_site$WATER_BODY=="POOLE HARBOUR",
                                                 "Poole",NA))))

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
  geom_textbox(size=3,data=scores_site,aes(x=mn_ax1_WB,
                                           y=mn_ax2_WB,
                                           label=WBlbl,
                                           fill=WATER_BODY),
               width = unit(0.15, "npc"),
               inherit.aes = FALSE,show.legend = FALSE,hjust=0.5)+
  coord_equal()+
  labs(title = "NMDS: Trait data")+
  theme(legend.title = element_blank(),
        legend.position="none",
        axis.title = element_text(face="bold")) -> pltrt
ggsave(filename = "figs/MDS_by_traits.pdf",width = 12,height = 12,
       units = "in",
       plot=pltrt)

png(file = "figs/MDS.png",
    width=12*ppi, height=6*ppi, res=ppi)
plabund+pltrt
dev.off()
