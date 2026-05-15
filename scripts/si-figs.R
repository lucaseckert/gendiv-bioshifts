#### SI FIGURES ####

## This script begins with a copy of 0_ModelFit.R from Oliveira et al.
## To reproduce data filtering and weights
## The additions to the code are clearly marked

gc();rm(list=ls())

################################################################################
#required packages
list.of.packages <- c(
    "doParallel", "parallel","foreach",
    "pbapply","dplyr", "tidyr", "data.table",
    "scales","effects","psych", 
    "glmmTMB", "lme4", "lmerTest","here","rlist","partR2",
    "ggtext","gridExtra","grid","lattice","viridis","performance","MuMIn","glmm.hp") 


new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages)) install.packages(new.packages)

sapply(list.of.packages, require, character.only = TRUE)

################################################################################
#define the data repository
# if(!dir.exists(here("Output/full_model"))){
#     dir.create(here("Output/full_model"),recursive = TRUE)
# }
# 
# dir.in=here("Data") #to change accordingly to the location of the data
# dir.out=here("Output/full_model") #to change. It's the repository where the results are saved
# 
# # Load data
# mydataset <- read.csv2(here(dir.in,"gen_data_final_fonseca.csv"),
#                        sep=",",dec=".",h=T) 


mydataset<-read.csv("data/gen_data_final_fonseca.csv")

#Data selection
## Latitude data
mydatatogo <- mydataset  %>%
    dplyr::filter(Type == "LAT",
                  shift_vel_sign == "pospos" | shift_vel_sign == "negneg", # Select only shifts in the same direction of velocity
                  SHIFT != 0, # remove non-significant shifts
                  # Nucleotide_diversity > 0 # select only GD values > 0
    ) %>% 
    mutate(
        # Climate velocity
        vel = as.numeric(velocity),
        vel_abs = abs(vel),
        vel_abs_log = log(vel_abs),
        vel_abs_log1p = log1p(vel_abs),
        
        # Shift
        SHIFT = SHIFT, # to deal with zero shift
        SHIFT_abs = abs(SHIFT),
        SHIFT_abs_log = log(SHIFT_abs),
        SHIFT_abs_log1p = log1p(SHIFT_abs),
        
        # Genetic diversity
        GD = Nucleotide_diversity, 
        GD_sqrt = sqrt(GD),
        GD_log = log(GD),
        GD_log1p = log1p(GD),
        
        # Methods
        Lat = abs(Latitude),
        Lat_band = round(Lat,0),
        ID.area = scale(ID.area),
        DUR = scale(DUR),
        LogExtent = log(Extent),
        START = scale(START),
        Param = factor(Param),
        Group = factor(Group),
        spp = factor(spp), 
        ExtentF = cut(Extent,
                      ordered_result = TRUE,
                      breaks=seq(min(Extent), max(Extent), length.out=10),
                      include.lowest=TRUE),
        NtempUnitsF = cut(NtempUnits,
                          ordered_result = TRUE,
                          breaks=seq(min(NtempUnits), max(NtempUnits), length.out=10),
                          include.lowest=TRUE)) %>%
    
    dplyr::select(
        # Genetic diversity
        GD, GD_log, GD_log1p, GD_sqrt, TajimasD,
        # Shift
        SHIFT, SHIFT_abs, SHIFT_abs_log, SHIFT_abs_log1p, 
        # SHIFT_cor, SHIFT_cor_abs, SHIFT_cor_abs_log, SHIFT_cor_raw, SHIFT_abs_log_scale,
        # Velocity
        vel, vel_abs, vel_abs_log, vel_abs_log1p, 
        trend.mean,
        # Methods + Taxonomy
        Article_ID, 
        Hemisphere,
        shift_vel_sign,
        Lat, # latitudinal position where GD was collected
        Lat_band,
        DUR, Nperiodes, LogNtempUnits, NtempUnits, Extent, LogExtent, ContinuousGrain, Quality, PrAb, ExtentF, NtempUnitsF,
        Param, Group, spp, Class, Order, Family, Genus, 
        ECO, Uncertainty_Parameter, Uncertainty_Distribution, Grain_size, Data, Article_ID
    ) 

# transform continuous variables
cont_vars <- c(1:14,18, 20:26)

mydatatogo[,cont_vars] <- lapply(mydatatogo[,cont_vars], as.numeric)
mydatatogo[,-cont_vars] <- lapply(mydatatogo[,-cont_vars], function(x) factor(x, levels = unique(x)))

##eliminating species with uncertain obs
mydatatogo=subset(mydatatogo,spp!="Agrilus planipennis")
mydatatogo=subset(mydatatogo,spp!="Chrysodeixis eriosoma")
mydatatogo=droplevels(mydatatogo)


#### ADDITION: TAJIMAS D HIST ####

tajima_df<-mydatatogo %>% select(spp, TajimasD) %>% distinct()

d_fig<-ggplot(tajima_df, aes(x=TajimasD))+
  #geom_vline(xintercept = 0, linetype="dashed", color = "grey", linewidth=0.75)+
  geom_density(fill="orange", alpha=0.5)+
  geom_vline(xintercept = -2, linetype="dashed", color = "grey40", linewidth=0.75)+
  geom_vline(xintercept = 2, linetype="dashed", color = "grey40", linewidth=0.75)+
  labs(x="Tajimas D", y="Density")+
  theme_bw()+
  theme(panel.grid = element_blank())

## proportion of species with significant Tajimas D
nrow(filter(tajima_df, TajimasD<=-2))/nrow(tajima_df)

## saving
ggsave(d_fig, 
       filename = "figures/tajimas-d.tiff",
       device = "tiff",
       width = 3, height = 2,
       units = "in", dpi = 300,
       scale = 2)

## Filter Classes with at least 5 species per Param
n_sps = 10

test <- mydatatogo %>%
    group_by(Class,Param) %>%
    summarise(N_spp = length(unique(spp))) %>% # how many species per parameter?
    dplyr::filter(N_spp >= n_sps) # select classes with > n_sps per param

test <- mutate(test, Class_Param = paste(Class, Param))

mydatatogo <- mydatatogo %>%
    mutate(Class_Param = paste(Class, Param)) %>%
    dplyr::filter(Class_Param %in% test$Class_Param) %>%
    dplyr::select(-Class_Param)

mydatatogo[,-cont_vars] <- lapply(mydatatogo[,-cont_vars], function(x) factor(x, levels = unique(x)))

## Extra fixes
# Set the reference param level to the centroid of species obs
mydatatogo$Param <- relevel(mydatatogo$Param, ref = "O") 

# Taxonomy data
taxatable <- mydatatogo %>%
    group_by(Class, Order, Family) %>%
    summarise(Species = length(unique(spp)),
              Shift = n())

# write.csv(taxatable,
#           here("Output/taxa_table.csv"),
#           row.names = FALSE)

################################################################################
##############################Variance partitioning ############################
################################################################################


mydatatogo2=mydatatogo
mydatatogoCE=subset(mydatatogo2,Param=="O")
x1=data.frame(table(mydatatogoCE$spp))
x1=subset(x1,Freq>0)
x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
mydatatogoCE=merge(mydatatogoCE,x1[,c(1,3)],by.x="spp",by.y="Var1")

mydatatogoLE=subset(mydatatogo2,Param=="LE")
x1=data.frame(table(mydatatogoLE$spp))
x1=subset(x1,Freq>0)
x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
mydatatogoLE=merge(mydatatogoLE,x1[,c(1,3)],by.x="spp",by.y="Var1")

mydatatogoTE=subset(mydatatogo2,Param=="TE")
x1=data.frame(table(mydatatogoTE$spp))
x1=subset(x1,Freq>0)
x1$weight_obs_spp=(1/x1$Freq)*(1/nrow(x1))
mydatatogoTE=merge(mydatatogoTE,x1[,c(1,3)],by.x="spp",by.y="Var1")

mydatatogoTE$weight_obs_spp=mydatatogoTE$weight_obs_spp*(1/3)
mydatatogoCE$weight_obs_spp=mydatatogoCE$weight_obs_spp*(1/3)
mydatatogoLE$weight_obs_spp=mydatatogoLE$weight_obs_spp*(1/3)

mydatatogo2=rbind(mydatatogoTE,mydatatogoCE,mydatatogoLE)

#### ADDITION: CHECK WEGHTS ####

mydatatogo2 %>% arrange(weight_obs_spp) %>% head()
mydatatogo2 %>% arrange(-weight_obs_spp) %>% head()

param_weights<-mydatatogo2 %>% group_by(Param) %>% 
  summarize(mean_weight=mean(weight_obs_spp),
            median_weight=median(weight_obs_spp),
            min_weight=min(weight_obs_spp),
            max_weight=max(weight_obs_spp),
            n=n())

range_weights<-ggplot(mydatatogo2, aes(x=Param, y=weight_obs_spp))+
  geom_boxplot(fill = "orange", alpha = 0.5) +
  geom_text(data=param_weights,
            aes(label = paste0("n=",n), x = Param, y = max_weight+0.00005), 
            color="grey35",
            hjust = 0.5,      
            vjust = 0.5,          
            size = 3.5, 
            check_overlap = TRUE) +
  scale_x_discrete(labels = c("O" = "CE", "LE" = "LE", "TE" = "TE"))+
  labs(y="Weight", x="Range Position")+
  theme_bw()+
  theme(panel.grid = element_blank())

weight_data<-mydatatogo2 %>% group_by(Class, Param) %>% 
  summarize(mean_weight=mean(weight_obs_spp),
            median_weight=median(weight_obs_spp),
            min_weight=min(weight_obs_spp),
            max_weight=max(weight_obs_spp),
            n=n()) %>% 
  mutate(col = if_else(mean_weight > 0.00075, "white", "black"))

tiles<-ggplot(weight_data, aes(x=Param, y=Class, fill=mean_weight))+
  geom_tile()+
  geom_text(data=weight_data,
            aes(label = paste0("n=",n), x=Param, y=Class, colour = col), 
            hjust = 0.5,      
            vjust = 0.5,          
            size = 3.5, 
            check_overlap = TRUE) +
  scale_color_manual(values=c("white","black"))+
  scale_fill_viridis_c(option="inferno")+
  scale_x_discrete(labels = c("O" = "CE", "LE" = "LE", "TE" = "TE"))+
  labs(y="Class", x="Range Position", fill="Mean weight")+
  guides(colour = "none")+ 
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(angle=45, hjust=1),
        legend.position="right",
        legend.title = element_text(vjust=0.75))

weight_plot<-cowplot::plot_grid(range_weights, tiles,
                     nrow = 1,
                     rel_widths = c(0.4, 0.6),
                     labels = c('A', 'B'),
                     label_size = 14,
                     #align = 'vh',
                     axis = 'tb')

ggsave(weight_plot, 
       filename = "figures/weights.tiff",
       device = "tiff",
       width = 4, height = 2,
       units = "in", dpi = 300,
       scale = 2)

############################################################


