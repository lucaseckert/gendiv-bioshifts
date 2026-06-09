#### LRT BOOTSTRAP ####

## packages
pkgs<-c("tidyverse", "DHARMa", "here", "glmmTMB")
invisible(lapply(pkgs, library, character.only = TRUE))

#### DATA PREP ####
data<-read.csv(here("data", "gen_data_final_fonseca.csv"))

mydatatogo <- data  %>%
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

mydatatogo2<- mydatatogo
mydatatogo2$vel_abs_s <- scale(mydatatogo2$vel_abs)
mydatatogo2$GD_s <- scale(mydatatogo2$GD)

#### LRT BOOTSTRAP ####

## full model, with species random effect
full_model<-glmmTMB(SHIFT_abs ~ vel_abs_s * GD_s * Param + 
                  LogNtempUnits + LogExtent + ContinuousGrain + PrAb + Quality +
                  (1|Class/spp),
                family = Gamma(link = "log"), REML = FALSE, data = mydatatogo2)

## null model without any GD terms
null_gd<-glmmTMB(SHIFT_abs ~ vel_abs_s * Param +
                    LogNtempUnits + LogExtent + ContinuousGrain + PrAb + Quality +
                    (1|Class/spp),
                  family = Gamma(link = "log"), REML = FALSE, data = mydatatogo2)

## observed LR
obsLR <- as.numeric(2 * (logLik(full_model) - logLik(null_gd)))

## lrt bootstrap
lrt_res<-simulateLRT(null_gd, full_model, n = 1000)
saveRDS(lrt_res, here("data", "lrt_bootstrap.rds"))

## results
lrt_res


