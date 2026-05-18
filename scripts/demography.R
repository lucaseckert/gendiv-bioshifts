#### DEMOGRAPHY ####

## This script explored demographic correlates of genetic diversity
## uses pi data from Fonsecca and demographic data from Avonet

## packages
pkgs<-c("tidyverse", "lme4", "lmerTest", "ggeffects", "cowplot", "parameters")
invisible(lapply(pkgs, library, character.only = TRUE))

#### DATA PREP ####

## bioshifts data
bioshifts<-read.csv("data/gen_data_final_fonseca.csv") %>% 
  filter(Class=="Aves",
         shift_vel_sign %in% c("pospos","negneg"),
         Type=="LAT") %>% 
  select(sp=sp_name_std_v1, pi=Nucleotide_diversity) %>% 
  distinct()

## avonet data
avonet<-read.csv("data/avonet-data.csv") %>% 
  select(sp=Species1, order=Order1, Range.Size, Centroid.Latitude, Migration, Hand.Wing.Index, Mass, Wing.Length)

## mmismatched species - just fixing manually
bioshifts %>% filter(!sp %in% avonet$sp) %>% pull(sp)
bioshifts[bioshifts$sp=="Coloeus monedula",]$sp<-"Corvus monedula"
bioshifts[bioshifts$sp=="Phylloscopus sibillatrix",]$sp<-"Phylloscopus sibilatrix"
bioshifts[bioshifts$sp=="Charadrius morinellus",]$sp<-"Eudromias morinellus"
bioshifts[bioshifts$sp=="Chroicocephalus ridibundus",]$sp<-"Larus ridibundus"
bioshifts[bioshifts$sp=="Cuculus optatus",]$sp<-"Cuculus horsfieldi" # no data in avonet
bioshifts[bioshifts$sp=="Luscinia svecica",]$sp<-"Cyanecula svecica"
bioshifts[bioshifts$sp=="Spizelloides arborea",]$sp<-"Passerella arborea"
bioshifts[bioshifts$sp=="Tetrastes bonasia",]$sp<-"Bonasa bonasia"

## combine
data<-bioshifts %>%
  filter(sp %in% avonet$sp) %>% 
  # losing that one cuckoo with no range data
  left_join(avonet, by=c("sp")) %>%
  rename(mass=Mass, hand_wing=Hand.Wing.Index, wing_length=Wing.Length) %>%
  mutate(range_size=as.numeric(Range.Size),
         lat=abs(as.numeric(Centroid.Latitude)),
         order=factor(order),
         pi_sqrt=sqrt(pi),
         pi_log=log(pi+1e-4),
         migration=factor(case_when(
           Migration == "1" ~ 1,
           Migration %in% c("2", "3") ~ 2)))

## summary - n=341 species across 18 orders
summary(data)

#### QUICK PLOTS ####

## pi ~ range size
ggplot(data, aes(x=log(range_size), y=pi_log))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()

## pi ~ latitude
ggplot(data, aes(x=lat, y=pi_log))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()

## pi ~ order
ggplot(data, aes(x=order, y=pi_log))+
  geom_boxplot()+
  theme_bw()

#### LINEAR MIXED MODELS ####

## latitude lmm
lat_model<-lmer(pi_log ~ lat + (1 | order), data=data)
summary(lat_model)
model_parameters(lat_model, standardize = "refit")

## range size lmm with order as random effect
model<-lmer(pi_log~ log(range_size) + (1 | order), 
            data=data)
summary(model)
model_parameters(model, standardize = "refit")

## with latitude
model_lat<-lmer(pi_log ~ log(range_size) + lat + (1 | order), 
                      data=data)
summary(model_lat)
model_parameters(model_lat, standardize = "refit")

#### EFFECT OF RANGE SIZE ####

## get marginal effects
pred_data<-ggpredict(model, terms = "range_size [n=20]", back_transform = FALSE)

## plotting
model_label<-"log(pi) %~% log('range size') + (1 * '|' * order)"
model_plot<-ggplot(filter(data, range_size > exp(12)), aes(x = log(range_size), y = pi_log)) +
  geom_point(alpha = 0.75, color="grey") +
  geom_line(data = filter(pred_data, x>exp(12)), aes(x = log(x), y = predicted), color = "orange", linewidth = 1) +
  geom_ribbon(data = filter(pred_data, x>exp(12)), aes(x = log(x), ymin = conf.low, ymax = conf.high), 
              alpha = 0.2, fill = "orange", inherit.aes = F) +
  annotate("text", x = 15, y = -8.8, label = model_label, size = 3.5, parse=T, color="grey35") +
  labs(
    x = expression(log(range~size)),
    y = expression(log(pi)))+
  theme_bw() +
  theme(panel.grid = element_blank())

#### EFFECT OF LATITUDE ####

## predict latitude effect
lat_pred<-ggpredict(lat_model, terms = "lat [n=20]", back_transform = FALSE)

## plotting
lat_plot<-ggplot(data, aes(x=lat, y=pi_log))+
  geom_point(alpha = 0.75, color="grey") +
  geom_line(data = lat_pred, aes(x = x, y = predicted), color = "orange", linewidth = 1) +
  geom_ribbon(data = lat_pred, aes(x = x, ymin = conf.low, ymax = conf.high), 
              alpha = 0.2, fill = "orange", inherit.aes = F) +
  annotate("text",
           x = 37.5, y = -8.8, 
           label = "log(pi) %~% 'abs(latitude)' + (1 * '|' * order)", 
           size = 3.5, parse=T, color="grey35") +
  labs(x = "abs(latitude)",
       y = expression(log(pi)))+
  theme_bw() +
  theme(panel.grid = element_blank())

## saving
ggsave(lat_plot, 
       filename = "figures/latitude.tiff",
       device = "tiff",
       width = 3, height = 2,
       units = "in", dpi = 300,
       scale = 2)

#### EFFECT OF PHYLOGENY ####

## median pi by order
order_means<-data %>% 
  group_by(order) %>% 
  summarise(med_pi=median(pi),
            n=n()) %>% 
  filter(n>2) %>% 
  arrange(med_pi)

## data for figure
order_data<-filter(data, order %in% count(data,order)[count(data,order)$n>2,]$order) %>% 
  mutate(order=factor(order, levels = order_means$order))

y_top<-max(order_data$pi_log)
y_bottom<-min(order_data$pi_log)
order_means<-order_means %>% mutate(pos=c(rep(y_top, 6), rep(y_bottom, 6)))

## plotting
order_plot<-ggplot(order_data, aes(x = order, y = pi_log)) +
  geom_boxplot(fill = "orange", alpha = 0.5) +
  geom_text(data=order_means,
            aes(label = order, x = order, y = pos), 
            color="grey35",
            angle = 90,           
            hjust = c(rep(0.9, 6), rep(0.1, 6)),         
            vjust = 1.5,          
            size = 3.5, 
            check_overlap = TRUE) +
  labs(y = expression(log(pi)),
       x = "Order") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color="white"),
        axis.ticks.x = element_blank())

#### SAVING FIGURE ####

## combine
main_plot<-plot_grid(order_plot, model_plot,
          nrow = 1,
          labels = c('A', 'B'),
          label_size = 14,
          align = 'vh',
          axis = 'tblr')

## saving
ggsave(main_plot, 
       filename = "figures/demography.tiff",
       device = "tiff",
       width = 4, height = 2,
       units = "in", dpi = 300,
       scale = 2)
