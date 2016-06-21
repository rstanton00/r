#update packages as necessary
require(utils)
update.packages(ask=FALSE)

#verify that necessary packages are installed
list.of.packages <- c("ggplot2", "Rcpp", "Rmisc", "car", "asbio", "lsmeans", "dplyr", "gridExtra", "cowplot")
list.of.packages <- c("ggplot2", "Rcpp", "Rmisc", "car", "asbio", "lsmeans", "dplyr", "gridExtra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(Rmisc)
library(dplyr)
library(car)
library(lsmeans)
library(asbio)
library(ggplot2)
library(gridExtra)
#library(cowplot)
require(stats)
require(graphics)

rm(list=ls())

#linux
setwd('/home/rstanton/Documents/biology/writings/thesis/thesis_stats_and_comments/csv_data')

#mac
#setwd('/Users/rstanton/Documents/pers/thesis/N orbicollis data/Data Analysis/csv')

diet_data <- read.csv('no_diet_lab_and_field.csv', head=TRUE, stringsAsFactors=FALSE)

#tell R that certain data are categorical
diet_data$Sex <- factor(diet_data$Sex, levels = c("M", "F"), labels = c("Male", "Female"))
diet_data$Treatment <- factor(diet_data$Treatment, levels = c("a", "b", "c", "TL_Davis_2010", "Washingon_Co_2011", "TL_Davis_2011"))
diet_data$Measured <- factor(diet_data$Measured, levels = c("PO_P", "Melaniz"))

#remove known outliers
#assign rows to data, using dataset data, where column BeetleId in data != "delete"
diet_data <- diet_data[diet_data$BeetleID != 'delete',]
diet_data <- diet_data[diet_data$BeetleID != 'bF004',]
#diet_data <- diet_data[diet_data$Population != 'TL_Davis_2010',]
# does this work? diet_data <- filter(diet_data, !grepl('delete', BeetleID))

#add columns of transformed data
diet_data["ln_premass_mg"] <- c(log2(diet_data$PreMass_mg))
diet_data["ln_postmass_mg"] <- c(log2(diet_data$PostMass_mg))
diet_data["ln_elytra"] <- c(log2(diet_data$ElytraLength_mm))
diet_data["ln_melanization"] <- c(log2(diet_data$z_AGV))
diet_data["ln_po"] <- c(log2(diet_data$zPO_avg_abs_min))
diet_data["sqrt_p"] <- c(log2(diet_data$z_Protein_avg_adjusted_mg_ml))

#create subsets of data for when we do analysis only on lab or field data
lab_data <- subset(diet_data, diet_data$Treatment != 'TL_Davis_2010' & 
                     diet_data$Treatment != 'TL_Davis_2011' &
                     diet_data$Treatment != 'Washingon_Co_2011')
#re-declare treatment levels for lab data as we only had a,b,c treatments in lab
lab_data$Treatment <- factor(lab_data$Treatment, levels = c("a", "b", "c"))

field_data <- subset(diet_data, diet_data$Treatment != 'a' & diet_data$Treatment != 'b' & 
                       diet_data$Treatment != 'c')
field_data$Treatment <- factor(field_data$Treatment, levels = c("TL_Davis_2010", "TL_Davis_2011", "Washingon_Co_2011"))

#set options for all linear models/anovas
options(contrasts=c("contr.sum","contr.poly"))

######################################
# lab analysis section
######################################

#test for equal start condions for beetle age, body mass, elytra length, and body condition */
#All data together measured stands for "PO and Protein measurements" and "Filaments= Melanization" respectively */
#age
model <- lm(BeetleAge_Start_days ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data=lab_data, na.action=na.omit)
age.anova <- Anova(model, type=c(3))
age.rg <- ref.grid(model)
lsmeans(age.rg, "Measured")

#pre body mass
model <- lm(ln_premass_mg ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data=lab_data, na.action=na.omit)
prebodymass.anova <- Anova(model, type=c(3))
prebodymass.rg <- ref.grid(model)
lsmeans(prebodymass.rg, "Measured")

#elytra length
model <- lm(ln_elytra ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data=lab_data, na.action=na.omit)
elytra.anova <- Anova(model, type=c(3))
elytra.rg <- ref.grid(model)
lsmeans(elytra.rg, "Measured")

#create a labeller function for renaming the treatment facets for facet wrapping
to_string_trt <- as_labeller(c('PO_P'="Phenoloxidase", 'Melaniz'="Melanization"))

#plot elytra length point plots with SE bars
#create SE measurements
dataSum <- summarySE(lab_data, measurevar="ElytraLength_mm", groupvars=c("Measured", "Treatment", "Sex"))

ggplot(dataSum, aes(x=factor(Treatment), y=ElytraLength_mm, pch=Sex,
                  ymax=ElytraLength_mm + dataSum$se,
                  ymin=ElytraLength_mm - dataSum$se)) +
  geom_point(position=position_dodge(width=0.25), size=2.5) +
  geom_errorbar(position=position_dodge(width=0.25), width=0.25) +
  scale_x_discrete(breaks = c("a", "b", "c"), labels=c("Ad lib.", "3 days", "5 days")) +
  facet_wrap(~ Measured, labeller = to_string_trt) +
  ylab("Elytra Length (mm)") + xlab("Treatment") +
  theme_bw() +
  theme(text = element_text(size=11),
        strip.text.x=element_text(size=11),
        strip.text.y=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank())

#alternately, six different plots, though "limits" doesn't work as of 2016/06/15
#ggplot(elSum, aes(x=factor(Measured), y=ElytraLength_mm), group=Sex) +
#geom_point() + geom_errorbar(limits, width=0.1) + facet_wrap(~ Treatment + Sex)

#pre-condition calculated
model <- lm(resid_premass_calculated ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = lab_data, na.action=na.omit)
residpremass.anova <- Anova(model, type=c(3))
residpremass.rg <- ref.grid(model)
lsmeans(residpremass.rg, "Sex")
lsmeans(residpremass.rg, "Measured")
lsmeans(residpremass.rg, "Treatment")
#produces similar results, slight rounding differences between SAS and R output
#conduct pairwise comparisons between treatments using the scheffe test
pairw.anova(y=lab_data$resid_premass_calculated, x=lab_data$Treatment, method="scheffe")
#plot(pairw.anova(y=lab_data$resid_premass_calculated, x=lab_data$Treatment, method="scheffe"))

#plot pre condition point plots with SE bars
#create SE measurements
dataSum <- summarySE(lab_data, measurevar="resid_premass_calculated", groupvars=c("Measured", "Treatment", "Sex"))

ggplot(dataSum, aes(x=factor(Treatment), y=resid_premass_calculated, pch=Sex,
                    ymax=resid_premass_calculated + dataSum$se,
                    ymin=resid_premass_calculated - dataSum$se)) +
  geom_point(position=position_dodge(width=0.25), size=2.5) +
  geom_errorbar(position=position_dodge(width=0.25), width=0.25) +
  scale_x_discrete(breaks = c("a", "b", "c"), labels=c("Ad lib.", "3 days", "5 days")) +
  facet_wrap(~ Measured, labeller = to_string_trt) +
  ylab("Pre-Experiment Body Condition (mg)") + xlab("Treatment") +
  theme_bw() +
  theme(text = element_text(size=11),
        strip.text.x=element_text(size=11),
        strip.text.y=element_text(size=11),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank())