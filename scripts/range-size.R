#### RANGE SIZE BIRDS ####

## packages
library(tidyverse)
library(terra)
library(lme4)
library(lmerTest)
library(ggeffects)
library(cowplot)

## range shift data
data<-read.csv("data/gen_data_final_fonseca.csv") %>% 
  filter(Class=="Aves",
         shift_vel_sign %in% c("pospos","negneg")) %>% 
  select(sp=sp_name_std_v1, pi=Nucleotide_diversity) %>% 
  distinct()

## avonet data
avonet<-readxl::read_excel("data/avonet-data.xlsx", sheet = 2) %>% 
  right_join(data, by=c("Species1"="sp")) %>% 
  mutate(range_size=as.numeric(Range.Size),
         lat=abs(as.numeric(Centroid.Latitude))) %>% 
  filter(pi>0, !is.na(Order1)) %>% 
  mutate(order=factor(Order1),
         migration=factor(case_when(
           Migration == "1" ~ 1,
           Migration %in% c("2", "3") ~ 2,
           TRUE ~ NA_real_)))

## pi ~ range size
ggplot(avonet, aes(x=log(range_size), y=log(pi)))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()

## pi ~ latitude
ggplot(avonet, aes(x=lat, y=log(pi)))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_bw()

## latitude lm
lm(log(pi) ~ lat, data=avonet) %>% summary()
lat_model<-lmer(log(pi) ~ lat + (1 | order), data=avonet)

## LMM with order as random effect
model<-lmer(log(pi) ~ log(range_size) + lat + (1 | order), 
            data=avonet)
summary(model)

## migration should maybe be in there
model_migration<-lmer(log(pi) ~ log(range_size) + lat + migration + (1 | order), 
                      data=avonet)
summary(model_migration)

## without latitude
model_nolat<-lmer(log(pi) ~ log(range_size) + (1 | migration) + (1 | order), 
                      data=avonet)
summary(model_nolat)

## get marginal effects
pred_data<-ggpredict(model, terms = "range_size [n=20]", back_transform = FALSE)

## plotting
model_label<-"log(pi) %~% log('range size') + 'abs(latitude)' + (1 * '|' * order)"
model_plot<-ggplot(filter(avonet, range_size > exp(10)), aes(x = log(range_size), y = log(pi))) +
  geom_point(alpha = 0.75, color="grey") +
  geom_line(data = filter(pred_data, x>exp(10)), aes(x = log(x), y = predicted), color = "orange", linewidth = 1) +
  geom_ribbon(data = filter(pred_data, x>exp(10)), aes(x = log(x), ymin = conf.low, ymax = conf.high), 
              alpha = 0.2, fill = "orange") +
  annotate("text", x = 14, y = -8.85, label = model_label, size = 3.5, parse=T, color="grey35") +
  labs(
    x = expression(log(range~size)),
    y = expression(log(pi)))+
  theme_bw() +
  theme(panel.grid = element_blank())

## pi ~ latitude nice
lat_pred<-ggpredict(lat_model, terms = "lat [n=20]", back_transform = FALSE)
lat_label<-"log(pi) %~% 'abs(latitude)' + (1 * '|' * order)"
lat_plot<-ggplot(avonet, aes(x=lat, y=log(pi)))+
  geom_point(alpha = 0.75, color="grey") +
  geom_line(data = lat_pred, aes(x = x, y = predicted), color = "orange", linewidth = 1) +
  geom_ribbon(data = lat_pred, aes(x = x, ymin = conf.low, ymax = conf.high), 
              alpha = 0.2, fill = "orange") +
  annotate("text", x = 40, y = -8.85, label = lat_label, size = 3.5, parse=T, color="grey35") +
  labs(
    x = "abs(latitude)",
    y = expression(log(pi)))+
  theme_bw() +
  theme(panel.grid = element_blank())

## pi ~ order
y_top<-log(max(order_data$pi))
y_bottom<-log(min(order_data$pi))

order_means<-avonet %>% 
  group_by(order) %>% 
  summarise(med_pi=median(pi),
            n=n()) %>% 
  filter(n>2) %>% 
  arrange(med_pi) %>% 
  mutate(pos=c(rep(y_top, 6), rep(y_bottom, 6)))

h<-c(rep(0.9, 6), rep(0.1, 6))

order_data<-filter(avonet, order %in% count(avonet,order)[count(avonet,order)$n>2,]$order) %>% 
  mutate(order=factor(order, levels = order_means$order))

order_plot<-ggplot(order_data, aes(x = order, y = log(pi))) +
  geom_boxplot(fill = "orange", alpha = 0.5) +
  geom_text(data=order_means,
            aes(label = order, x = order, y = pos), 
            angle = 90,           
            hjust = h,         
            vjust = 1.5,          
            size = 3.5, 
            check_overlap = TRUE) +
  labs(y = expression(log(pi)),
       x = "Order") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color="white"),
        axis.ticks.x = element_blank())

## combine
main_plot<-plot_grid(order_plot, lat_plot, model_plot, 
          nrow = 1,
          labels = c('A', 'B', 'C'), 
          label_size = 14,
          align = 'vh',       
          axis = 'tblr')

ggsave(main_plot, 
       filename = "figures/main-plot.tiff",
       device = "tiff",
       width = 6, height = 2,
       units = "in", dpi = 300,
       scale = 2)

## migration
ggplot(avonet, aes(x=migration, y=log(pi)))+
  geom_boxplot(fill = "orange", alpha = 0.5) +
  labs(y = expression(log(pi)),
       x = "Migration") +
  theme_bw() +
  theme(panel.grid = element_blank())


