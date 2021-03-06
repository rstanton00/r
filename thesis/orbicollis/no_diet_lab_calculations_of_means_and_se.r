#update packages as necessary
require(utils)
update.packages(ask=FALSE)

#verify that necessary packages are installed
#list.of.packages <- c("ggplot2", "Rcpp", "Rmisc", "car", "asbio", "lsmeans", "dplyr", "gridExtra", "cowplot")
list.of.packages <- c("agricolae", "ggplot2", "Rcpp", "Rmisc", "car", "asbio", "lsmeans", "dplyr", "gridExtra")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(Rmisc)
library(dplyr)
library(car)
library(lsmeans)
library(asbio)
library(ggplot2)
library(gridExtra)
library(agricolae)
#library(cowplot)
require(stats)
require(graphics)

rm(list=ls())

#linux
#setwd('/home/rls/Documents/biology/writings/thesis/thesis_stats_and_comments/csv_data')

#mac
setwd('/Users/rstanton/Documents/pers/thesis/N orbicollis data/Data Analysis/csv')

diet_data <- read.csv('no_diet_lab_and_field.csv', head=TRUE, stringsAsFactors=FALSE)

#remove known outliers
#assign rows to data, using dataset data, where column BeetleId in data != "delete"
diet_data <- diet_data[diet_data$BeetleID != 'delete',]
diet_data <- diet_data[diet_data$BeetleID != 'bF004',]
#remove TL Davis population from 2010, as elytra were not measured
diet_data <- diet_data[diet_data$Treatment != 'TL_Davis_2010',]
# does this work? diet_data <- filter(diet_data, !grepl('delete', BeetleID))

#tell R that certain data are categorical
diet_data$Sex <- factor(diet_data$Sex, levels = c("M", "F"), labels = c("Male", "Female"))
diet_data$Treatment <- factor(diet_data$Treatment, levels = c("a", "b", "c", "Washington_Co_2011", "TL_Davis_2011"))
#diet_data$Treatment <- factor(diet_data$Treatment, levels = c("a", "b", "c", "TL_Davis_2010", "Washington_Co_2011", "TL_Davis_2011"))
diet_data$Measured <- factor(diet_data$Measured, levels = c("PO_P", "Melaniz"))

#add columns of transformed data
diet_data["ln_premass_mg"] <- c(log2(diet_data$PreMass_mg))
diet_data["ln_postmass_mg"] <- c(log2(diet_data$PostMass_mg))
diet_data["ln_elytra"] <- c(log2(diet_data$ElytraLength_mm))
diet_data["ln_melanization"] <- c(log2(diet_data$z_AGV))
diet_data["ln_po"] <- c(log2(diet_data$zPO_avg_abs_min))
diet_data["sqrt_p"] <- c(log2(diet_data$z_Protein_avg_adjusted_mg_ml))

#create subsets of data for when we do analysis only on lab or field data
lab_data <- subset(diet_data, diet_data$Treatment != 'TL_Davis_2011' & diet_data$Treatment != 'Washington_Co_2011')
#re-declare treatment levels for lab data as we only had a,b,c treatments in lab
lab_data$Treatment <- factor(lab_data$Treatment, levels = c("a", "b", "c"))

field_data <- subset(diet_data, diet_data$Treatment != 'a' & diet_data$Treatment != 'b' & 
                       diet_data$Treatment != 'c')
field_data$Treatment <- factor(field_data$Treatment, levels = c("TL_Davis_2010", "TL_Davis_2011", "Washington_Co_2011"))

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
to_string_trt <- as_labeller(c('PO_P'="A", 'Melaniz'="B"))

#get summary of means and std. errors for analysis table
dataSum <- summarySE(lab_data, measurevar="ElytraLength_mm", groupvars=c("Treatment"))

#plot elytra length point plots with SE bars
#create SE measurements
dataSum <- summarySE(lab_data, measurevar="ElytraLength_mm", groupvars=c("Measured", "Treatment", "Sex"))

aovElytra <- aov(ln_elytra ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data=lab_data, na.action=na.omit)
summary(aovElytra)
HSD.test(aovElytra, "ln_elytra")
#TukeyHSD(aovElytra, "Treatment")

#field site 1 = WaCo 2011
#field site 2 = TLDavis 2011
#field site 3 = TLDavis 2010
ggplot(dataSum, aes(x=factor(Treatment), y=ElytraLength_mm, pch=Sex,
                  ymax=ElytraLength_mm + dataSum$se,
                  ymin=ElytraLength_mm - dataSum$se)) +
  xlab("Lab Starvation Treatment") +
  ylab("Elytra Length (mm)") +
  geom_point(position=position_dodge(width=0.3), size=2.1) +
  geom_errorbar(position=position_dodge(width=0.3), width=0.21) +
  scale_x_discrete(breaks = c("a", "b", "c"), labels=c("Fed AL", "Fed 2x", "Fed 0x")) +
  facet_wrap(~ Measured, labeller = to_string_trt) +
  theme_bw() +
  theme(text = element_text(size=10),
        #strip.text.x = element_text(size=10),
        strip.text.x = element_blank(),
        strip.text.y = element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        #legend.key = element_blank(),
        legend.position = "bottom",
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "black") )

#effken uno thesis:
#width = 4.23"
#height = 3.44"

#BOXPLOT
ggplot(lab_data, aes(x=factor(Treatment), y=ElytraLength_mm, pch=Sex, fill=Sex)) +
  geom_boxplot() +
  scale_x_discrete(breaks = c("a", "b", "c"), labels=c("Fed AdLib.", "Fed 2x", "Fed 0x")) +
  facet_wrap(~ Measured, labeller = to_string_trt) +
  ylab("Elytra Length (mm)") + xlab("Lab Starvation Treatment") +
  theme_bw() +
  theme(text = element_text(size=10),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
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

#field site 1 = WaCo 2011
#field site 2 = TLDavis 2011
#field site 3 = TLDavis 2010
ggplot(dataSum, aes(x=factor(Treatment), y=resid_premass_calculated, pch=Sex,
                    ymax=resid_premass_calculated + dataSum$se,
                    ymin=resid_premass_calculated - dataSum$se)) +
  xlab("Lab Starvation Treatment") +
  ylab("Pre-Experiment Body Condition (residuals)") + 
  geom_point(position=position_dodge(width=0.3), size=2.1) +
  geom_errorbar(position=position_dodge(width=0.3), width=0.21) +
  scale_x_discrete(breaks = c("a", "b", "c"), labels=c("Fed AL", "Fed 2x", "Fed 0x")) +
  facet_wrap(~ Measured, labeller = to_string_trt) +
  theme_bw() +
  theme(text = element_text(size=10),
        #strip.text.x=element_text(size=10),
        strip.text.x = element_blank(),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        #legend.key = element_blank())
        legend.position = "bottom",
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "black") )
        
#BOXPLOT
ggplot(lab_data, aes(x=factor(Treatment), y=diet_data$resid_premass_calculated, pch=Sex, fill=Sex)) +
  geom_boxplot() +
  scale_x_discrete(breaks = c("a", "b", "c"), labels=c("Fed AdLib.", "Fed 2x", "Fed 0x")) +
  facet_wrap(~ Measured, labeller = to_string_trt) +
  ylab("Pre-Experiment Body Condition (residuals)") + xlab("Lab Starvation Treatment") +
  theme_bw() +
  theme(text = element_text(size=10),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank())

#post mass
model <- lm(ln_postmass_mg ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = lab_data, na.action=na.omit)
postmass.anova <- Anova(model, type=c(3))
postmass.rg <- ref.grid(model)
lsmeans(postmass.rg, "Sex")
lsmeans(postmass.rg, "Measured")
lsmeans(postmass.rg, "Treatment")
#conduct pairwise comparisons between treatments using the scheffe test
pairw.anova(y=lab_data$ln_postmass, x=lab_data$Treatment, method="scheffe")
#plot(pairw.anova(y=lab_data$ln_postmass, x=lab_data$Treatment, method="scheffe"))

#post condition calculated
model <- lm(resid_postmass_calculated ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = lab_data, na.action=na.omit)
residpostmass.anova <- Anova(model, type=c(3))
residpostmass.rg <- ref.grid(model)
lsmeans(residpostmass.rg, "Sex")
lsmeans(residpostmass.rg, "Measured")
lsmeans(residpostmass.rg, "Treatment")
#conduct pairwise comparisons between treatments using the scheffe test
pairw.anova(y=lab_data$resid_postmass_calculated, x=lab_data$Treatment, method="scheffe")
#plot(pairw.anova(y=lab_data$resid_postmass_calculated, x=lab_data$Treatment, method="scheffe"))

#show analysis for results table
dataSum <- summarySE(lab_data, measurevar="resid_postmass_calculated", groupvars=c("Treatment"))

#plot LAB ONLY post condition point plots with SE bars
#create SE measurements
dataSum <- summarySE(lab_data, measurevar="resid_postmass_calculated", groupvars=c("Measured", "Treatment", "Sex"))

#field site 1 = WaCo 2011
#field site 2 = TLDavis 2011
#field site 3 = TLDavis 2010
ggplot(dataSum, aes(x=factor(Treatment), y=resid_postmass_calculated, pch=Sex,
                    ymax=resid_postmass_calculated + dataSum$se,
                    ymin=resid_postmass_calculated - dataSum$se)) +
  xlab("Lab Starvation Treatment") +
  ylab("Post-Experiment Body Condition (residuals)") +
  geom_point(position=position_dodge(width=0.3), size=2.1) +
  geom_errorbar(position=position_dodge(width=0.3), width=0.21) +
  scale_x_discrete(breaks = c("a", "b", "c"), labels=c("Fed AL", "Fed 2x", "Fed 0x")) +
  facet_wrap(~ Measured, labeller = to_string_trt) +
  theme_bw() +
  theme(text = element_text(size=10),
        #strip.text.x=element_text(size=10),
        strip.text.x = element_blank(),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        #legend.key = element_blank())
        legend.position = "bottom",
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "black") )

#BOXPLOT
ggplot(lab_data, aes(x=factor(Treatment), y=resid_postmass_calculated, pch=Sex, fill=Sex)) +
  geom_boxplot() +
  scale_x_discrete(breaks = c("a", "b", "c"), labels=c("Fed AdLib.", "Fed 2x", "Fed 0x")) +
  facet_wrap(~ Measured, labeller = to_string_trt) +
  ylab("Post-Experiment Body Condition (residuals)") + xlab("Lab Starvation Treatment") +
  theme_bw() +
  theme(text = element_text(size=10),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank())

#plot ALL DATA post condition point plots with SE bars
#create SE measurements
dataSum <- summarySE(diet_data, measurevar="resid_postmass_calculated", groupvars=c("Treatment", "Sex"))

#field site 1 = WaCo 2011
#field site 2 = TLDavis 2011
#field site 3 = TLDavis 2010
ggplot(dataSum, aes(x=factor(Treatment), y=resid_postmass_calculated, pch=Sex,
                    ymax=resid_postmass_calculated + dataSum$se,
                    ymin=resid_postmass_calculated - dataSum$se)) +
  geom_point(position=position_dodge(width=0.3), size=2.1) +
  geom_errorbar(position=position_dodge(width=0.3), width=0.21) +
  scale_x_discrete(breaks = c("a", "b", "c", "Washington_Co_2011", "TL_Davis_2011", "TL_Davis_2010"), labels=c("Fed AL", "Fed 2x", "Fed 0x", "FS1", "FS2", "FS3")) +
  xlab("Lab Starvation Treatment or Field Site") +
  ylab("Post-Experiment Body Condition (residuals)") + 
  theme_bw() +
  theme(text = element_text(size=10),
        #strip.text.x=element_text(size=10),
        strip.text.x = element_blank(),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        #legend.key = element_blank())
        legend.position = "bottom",
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "black") )

#BOXPLOT
ggplot(diet_data, aes(x=factor(Treatment), y=resid_postmass_calculated, pch=Sex, fill=Sex)) +
  geom_boxplot() +
  scale_x_discrete(breaks = c("a", "b", "c", "Washington_Co_2011", "TL_Davis_2011", "TL_Davis_2010"), labels=c("Fed AdLib.", "Fed 2x", "Fed 0x", "FS1", "FS2", "FS3")) +
  facet_wrap(~ Measured, labeller = to_string_trt) +
  ylab("Post-Experiment Body Condition (residuals)") + xlab("Lab Starvation Treatment") +
  theme_bw() +
  theme(text = element_text(size=10),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank())

#test for effect of starvation treatment
df <- subset(lab_data, select = c("BeetleID", "Sex", "Treatment", "Measured", "ln_premass_mg", "ln_postmass_mg"))
dfL <- reshape(df, varying = list(5:6), idvar = "BeetleID", direction = "long", v.names = "ln_mass_mg")
dfL$time <- factor(dfL$time, levels = 1:2, labels = c("pre", "post"))
head(dfL); tail(dfL)
aovMod <- aov(ln_mass_mg ~ Sex * Treatment * Measured * time + Error(BeetleID/time), data = dfL)

#repeated measures ANOVA for condition calculated
df <- subset(lab_data, select = c("BeetleID", "Sex", "Treatment", "Measured", "resid_premass_calculated", "resid_postmass_calculated"))
dfL <- reshape(df, varying = list(5:6), idvar = "BeetleID", direction = "long", v.names = "resid_mass_calc")
dfL$time <- factor(dfL$time, levels = 1:2, labels = c("pre", "post"))
head(dfL); tail(dfL)
aovMod <- aov(resid_mass_calc ~ Sex * Treatment * Measured * time + Error(BeetleID/time), data = dfL)

#plot protein avg
#create SE measurements
dataSum <- summarySE(diet_data, measurevar="z_Protein_avg_adjusted_mg_ml", groupvars=c("Measured", "Treatment", "Sex"))
levels(dataSum$Measured)[levels(dataSum$Measured)=="PO_P"] <- "Phenoloxidase"
levels(dataSum$Measured)[levels(dataSum$Measured)=="Melaniz"] <- "Melanization"

#create protein plot
#field site 1 = WaCo 2011
#field site 2 = TLDavis 2011
#field site 3 = TLDavis 2010
ggplot(dataSum, aes(x=factor(Treatment), y=z_Protein_avg_adjusted_mg_ml, pch=Sex,
                    ymax=z_Protein_avg_adjusted_mg_ml + dataSum$se,
                    ymin=z_Protein_avg_adjusted_mg_ml - dataSum$se)) +
  xlab("Lab Starvation Treatment or Field Site") +
  ylab("Protein Concentration in Hemolymph (mg/ml)") + 
  geom_point(position=position_dodge(width=0.3), size=2.1) +
  geom_errorbar(position=position_dodge(width=0.3), width=0.21) +
  scale_x_discrete(breaks = c("a", "b", "c", "Washington_Co_2011", "TL_Davis_2011", "TL_Davis_2010"), labels=c("Fed AL", "Fed 2x", "Fed 0x", "FS1", "FS2", "FS3")) +
  theme_bw() +
  theme(text = element_text(size=10),
        #strip.text.x=element_text(size=10),
        strip.text.x = element_blank(),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        #legend.key = element_blank())
        legend.position = "bottom",
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "black") )

#BOXPLOT
ggplot(diet_data, aes(x=factor(Treatment), y=z_Protein_avg_adjusted_mg_ml, pch=Sex, fill=Sex)) +
  geom_boxplot() +
  scale_x_discrete(breaks = c("a", "b", "c", "Washington_Co_2011", "TL_Davis_2011", "TL_Davis_2010"), labels=c("Fed AdLib.", "Fed 2x", "Fed 0x", "FS1", "FS2", "FS3")) +
  ylab("Protein Concentration in Hemolymph (mg/ml)") + xlab("Lab Starvation Treatment or Field Site") +
  theme_bw() +
  theme(text = element_text(size=10),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank())

#set up new dataset that only examines phenoloxidase/protein data for analysis
newData <- lab_data[ which(lab_data$Measured=='PO_P'), ]

#protein
model <- lm(sqrt_p ~ Sex + Treatment + Sex:Treatment, data = newData, na.action=na.omit)
labProteinAnova <- Anova(model, type=c(3))
protein.rg <- ref.grid(model)
lsmeans(protein.rg, "Sex")
lsmeans(protein.rg, "Treatment")
#produces similar results, slight rounding differences between SAS and R output
#conduct pairwise comparisons between treatments using the scheffe test
pairw.anova(y=newData$sqrt_p, x=newData$Treatment, method="scheffe")
#plot(pairw.anova(y=newData$sqrt_p, x=newData$Treatment, method="scheffe"), main="Treatment Effect on Protein",
#     xlab="Treatment", ylab="Protein mg/ml")

#plot phenoloxidase avg data
#create SE measurements
dataSum <- summarySE(diet_data, measurevar="zPO_avg_abs_min", groupvars=c("Measured", "Treatment", "Sex"))
levels(dataSum$Measured)[levels(dataSum$Measured)=="PO_P"] <- "Phenoloxidase"
levels(dataSum$Measured)[levels(dataSum$Measured)=="Melaniz"] <- "Melanization"

#create phenoloxidase plot
#field site 1 = WaCo 2011
#field site 2 = TLDavis 2011
#field site 3 = TLDavis 2010
ggplot(dataSum, aes(x=factor(Treatment), y=zPO_avg_abs_min, pch=Sex,
                    ymax=zPO_avg_abs_min + dataSum$se,
                    ymin=zPO_avg_abs_min - dataSum$se)) +
  xlab("Lab Starvation Treatment or Field Site") +
  ylab("Phenoloxidase (abs/min)") + 
  geom_point(position=position_dodge(width=0.3), size=2.1, na.rm = TRUE) +
  geom_errorbar(position=position_dodge(width=0.3), width=0.21) +
  scale_x_discrete(breaks = c("a", "b", "c", "Washington_Co_2011", "TL_Davis_2011", "TL_Davis_2010"), labels=c("Fed AL", "Fed 2x", "Fed 0x", "FS1", "FS2", "FS3")) +
  theme_bw() +
  theme(text = element_text(size=10),
        #strip.text.x=element_text(size=10),
        strip.text.x=element_blank(),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        #legend.key = element_blank())
        legend.position = "bottom",
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "black") )

#BOXPLOT
ggplot(diet_data, aes(x=factor(Treatment), y=zPO_avg_abs_min, pch=Sex, fill=Sex)) +
  geom_boxplot() +
  scale_x_discrete(breaks = c("a", "b", "c", "Washington_Co_2011", "TL_Davis_2011", "TL_Davis_2010"), labels=c("Fed AdLib.", "Fed 2x", "Fed 0x", "FS1", "FS2", "FS3")) +
  ylab("Phenoloxidase (abs/min)") + xlab("Lab Starvation Treatment or Field Site") +
  theme_bw() +
  theme(text = element_text(size=10),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank())

#phenoloxidase analysis
model <- lm(ln_po ~ Sex + Treatment + Sex:Treatment, data = newData, na.action=na.omit)
protein <- Anova(model, type=c(3))
po.rg <- ref.grid(model)
lsmeans(po.rg, "Treatment")
#produces similar results, slight rounding differences between SAS and R output
#conduct pairwise comparisons between treatments using the scheffe test
pairw.anova(y=newData$ln_po, x=newData$Treatment, method="scheffe")
#plot(pairw.anova(y=newData$ln_po, x=newData$Treatment, method="scheffe"), main="Treatment Effect on Phenoloxidase",
#     xlab="Treatment", ylab="Phenoloxidase")

##############
# Melanization
##############

#set up dataset that uses only melanization data
newData <- lab_data[ which(lab_data$Measured=='Melaniz'), ]

model <- lm(z_AGV ~ Sex + Treatment + Sex:Treatment, data = newData, na.action=na.omit)
labMelanAnova <- Anova(model, type=c(3))
labMelan.rg <- ref.grid(model)
lsmeans(labMelan.rg, "Sex")
lsmeans(labMelan.rg, "Treatment")
#produces similar results, slight rounding differences between SAS and R output
#conduct pairwise comparisons between treatments using the scheffe test
pairw.anova(y=newData$z_AGV, x=newData$Treatment, method="scheffe")

#plot melanization avg data
#create SE measurements
dataSum <- summarySE(diet_data, measurevar="z_AGV", groupvars=c("Measured", "Treatment", "Sex"))


levels(dataSum$Measured)[levels(dataSum$Measured)=="PO_P"] <- "Phenoloxidase"
levels(dataSum$Measured)[levels(dataSum$Measured)=="Melaniz"] <- "Melanization"
#create plot
#field site 1 = WaCo 2011
#field site 2 = TLDavis 2011
#field site 3 = TLDavis 2010
ggplot(dataSum, aes(x=factor(Treatment), y=z_AGV, pch=Sex,
                  ymax=z_AGV + dataSum$se,
                  ymin=z_AGV - dataSum$se)) +
  xlab("Lab Starvation Treatment or Field Site") +
  ylab("Average Grey Value") + 
  geom_point(position=position_dodge(width=0.3), size = 2.1, na.rm=TRUE) +
  geom_errorbar(position=position_dodge(width=0.3), width=0.21) +
  scale_x_discrete(breaks = c("a", "b", "c", "Washington_Co_2011", "TL_Davis_2011", "TL_Davis_2010"), labels=c("Fed AL", "Fed 2x", "Fed 0x", "FS1", "FS2", "FS3")) +
  scale_y_continuous(trans = "reverse") +
  theme_bw() +
  theme(text = element_text(size=10),
        #strip.text.x=element_text(size=10),
        strip.text.x=element_blank(),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        #legend.key = element_blank())
        legend.position = "bottom",
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "black") )

ggplot(dataSum, aes(x=factor(Treatment), y=z_Protein_avg_adjusted_mg_ml, pch=Sex,
                    ymax=z_Protein_avg_adjusted_mg_ml + dataSum$se,
                    ymin=z_Protein_avg_adjusted_mg_ml - dataSum$se)) +
  xlab("Lab Starvation Treatment or Field Site") +
  ylab("Protein Concentration in Hemolymph (mg/ml)") + 
  geom_point(position=position_dodge(width=0.3), size=2.1) +
  geom_errorbar(position=position_dodge(width=0.3), width=0.21) +
  scale_x_discrete(breaks = c("a", "b", "c", "Washington_Co_2011", "TL_Davis_2011", "TL_Davis_2010"), labels=c("Fed AL", "Fed 2x", "Fed 0x", "FS1", "FS2", "FS3")) +
  theme_bw() +
  theme(text = element_text(size=10),
        #strip.text.x=element_text(size=10),
        strip.text.x = element_blank(),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        #legend.key = element_blank())
        legend.position = "bottom",
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "black") )


#BOXPLOT
ggplot(diet_data, aes(x=factor(Treatment), y=z_AGV, pch=Sex, fill=Sex)) +
  geom_boxplot() +
  scale_x_discrete(breaks = c("a", "b", "c", "Washington_Co_2011", "TL_Davis_2011", "TL_Davis_2010"), labels=c("Fed AdLib.", "Fed 2x", "Fed 0x", "FS1", "FS2", "FS3")) +
  ylab("Average Grey Value") + xlab("Lab Starvation Treatment or Field Site") +
  theme_bw() +
  theme(text = element_text(size=10),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank())

#model <- lm(cbind(ln_premass_mg, ln_postmass_mg) ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = diet_lab_data, na.action=na.omit)
#effect_of_treatment <- Anova(model, type=c(3))

#df <- subset(diet_lab_data, select = c("BeetleID", "Sex", "Treatment", "Measured", "ln_premass_mg", "ln_postmass_mg"))
#TRY TO TELL R THAT PRE & POST MASS ARE CORRELATED VIA LME?  idea from learning r for sas/spss section

#t.test(diet_lab_data$ln_postmass_mg, diet_lab_data$ln_premass_mg, paired=TRUE, data=diet_lab_data, na.action=na.omit)
#summary(lme(ln_postmass_mg ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = diet_lab_data, na.action=na.omit))


######################################
# field analysis section
######################################
#need to create new dataset not including tl davis 2010 in order to get some models working
field_data_no_tl2010 <- field_data[field_data$Treatment != 'TL_Davis_2010',]
field_data_no_tl2010$Treatment <- factor(field_data_no_tl2010$Treatment, levels = c("TL_Davis_2011", "Washington_Co_2011"))

#body mass
model <- lm(ln_postmass_mg ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = field_data_no_tl2010, na.action=na.omit)
fieldbodymass.anova <- Anova(model, type=c(3))
fieldbodymass.rg <- ref.grid(model)
lsmeans(fieldbodymass.rg, "Measured")

#elytra length
model <- lm(ln_elytra ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = field_data_no_tl2010, na.action=na.omit)
fieldelytra.anova <- Anova(model, type=c(3))
fieldelytra.rg <- ref.grid(model)
lsmeans(fieldelytra.rg, "Measured")

#create a labeller function for renaming the treatment facets for facet wrapping
#to_string_trt <- as_labeller(c('PO_P'="Phenoloxidase", 'Melaniz'="Melanization"))

#get summary of means and std. errors for analysis table
dataSum <- summarySE(field_data_no_tl2010, measurevar="ElytraLength_mm", groupvars=c("Treatment"))

#plot elytra length point plots with SE bars
#create SE measurements
dataSum <- summarySE(field_data_no_tl2010, measurevar="ElytraLength_mm", groupvars=c("Measured", "Treatment", "Sex"))

aovFieldElytra <- aov(ln_elytra ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data=field_data_no_tl2010, na.action=na.omit)
summary(aovFieldElytra)
HSD.test(aovFieldElytra, "ln_elytra")
#TukeyHSD(aovElytra, "Treatment")

#field site 1 = WaCo 2011
#field site 2 = TLDavis 2011
#field site 3 = TLDavis 2010
ggplot(dataSum, aes(x=factor(Treatment), y=ElytraLength_mm, pch=Sex,
                    ymax=ElytraLength_mm + dataSum$se,
                    ymin=ElytraLength_mm - dataSum$se)) +
  geom_point(position=position_dodge(width=0.25), size=2.5) +
  geom_errorbar(position=position_dodge(width=0.25), width=0.25) +
  scale_x_discrete(breaks = c("TL_Davis_2011", "Washington_Co_2011"), labels=c("FS1", "FS2")) +
  facet_wrap(~ Measured, labeller = to_string_trt) +
  ylab("Elytra Length (mm)") + xlab("Field Site") +
  theme_bw() +
  theme(text = element_text(size=10),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank())

#BOXPLOT
ggplot(field_data_no_tl2010, aes(x=factor(Treatment), y=ElytraLength_mm, pch=Sex, fill=Sex)) +
  geom_boxplot() +
  scale_x_discrete(breaks = c("TL_Davis_2011", "Washington_Co_2011"), labels=c("FS1", "FS2")) +
  facet_wrap(~ Measured, labeller = to_string_trt) +
  ylab("Elytra Length (mm)") + xlab("Field Site") +
  theme_bw() +
  theme(text = element_text(size=10),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank())


#condition
model <- lm(resid_postmass_calculated ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = field_data_no_tl2010, na.action=na.omit)
fieldpostmass.anova <- Anova(model, type=c(3))
fieldpostmass.rg <- ref.grid(model)
lsmeans(fieldpostmass.rg, "Measured")

#get summary of means and std. errors for analysis table
dataSum <- summarySE(field_data_no_tl2010, measurevar="resid_postmass_calculated", groupvars=c("Measured", "Treatment", "Sex"))

aovFieldCondition <- aov(resid_postmass_calculated ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data=field_data_no_tl2010, na.action=na.omit)
summary(aovFieldCondition)
HSD.test(aovFieldCondition, "resid_postmass_calculated")
pairw.anova(y=field_data_no_tl2010$resid_postmass_calculated, x=field_data_no_tl2010$Treatment, method="scheffe")

#field site 1 = WaCo 2011
#field site 2 = TLDavis 2011
#field site 3 = TLDavis 2010
ggplot(dataSum, aes(x=factor(Treatment), y=resid_postmass_calculated, pch=Sex,
                    ymax=resid_postmass_calculated + dataSum$se,
                    ymin=resid_postmass_calculated - dataSum$se)) +
  xlab("Field Site") +
  ylab("Field Beetle Body Condition (residuals)") + 
  geom_point(position=position_dodge(width=0.3), size=2.1) +
  geom_errorbar(position=position_dodge(width=0.3), width=0.21) +
  scale_x_discrete(breaks = c("TL_Davis_2011", "Washington_Co_2011"), labels=c("FS1", "FS2")) +
  facet_wrap(~ Measured, labeller = to_string_trt) +
  theme_bw() +
  theme(text = element_text(size=10),
        #strip.text.x=element_text(size=10),
        strip.text.x=element_blank(),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        #legend.key = element_blank())
        legend.position = "bottom",
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "black") )

#BOXPLOT
ggplot(field_data_no_tl2010, aes(x=factor(Treatment), y=resid_postmass_calculated, pch=Sex, fill=Sex)) +
  geom_boxplot() +
  scale_x_discrete(breaks = c("TL_Davis_2011", "Washington_Co_2011"), labels=c("FS1", "FS2")) +
  facet_wrap(~ Measured, labeller = to_string_trt) +
  ylab("Field Beetle Body Condition (residuals)") + xlab("Field Site") +
  theme_bw() +
  theme(text = element_text(size=10),
        strip.text.x=element_text(size=10),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank())

###########
#Combined field and lab data for post-condition
###########
model <- lm(resid_postmass_calculated ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = diet_data, na.action=na.omit)
combinedpostmass.anova <- Anova(model, type=c(3))
combinedpostmass.rg <- ref.grid(model)
lsmeans(combinedpostmass.rg, "Measured")

#get summary of means and std. errors for analysis table
dataSum <- summarySE(diet_data, measurevar="resid_postmass_calculated", groupvars=c("Measured", "Treatment", "Sex"))

aovCombinedCondition <- aov(resid_postmass_calculated ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data=diet_data, na.action=na.omit)
summary(aovCombinedCondition)
HSD.test(aovCombinedCondition, "resid_postmass_calculated")
pairw.anova(y=diet_data$resid_postmass_calculated, x=diet_data$Treatment, method="scheffe")

#field site 1 = WaCo 2011
#field site 2 = TLDavis 2011
#field site 3 = TLDavis 2010
ggplot(dataSum, aes(x=factor(Treatment), y=resid_postmass_calculated, pch=Sex,
                    ymax=resid_postmass_calculated + dataSum$se,
                    ymin=resid_postmass_calculated - dataSum$se)) +
  xlab("Lab Starvation Treatment or Field Site") +
  ylab("Beetle Body Condition (residuals)") + 
  geom_point(position=position_dodge(width=0.3), size=2.1) +
  geom_errorbar(position=position_dodge(width=0.3), width=0.21) +
  scale_x_discrete(breaks = c("a", "b", "c", "TL_Davis_2011", "Washington_Co_2011"), labels=c("Fed AL", "Fed 2x", "Fed 0x", "FS1", "FS2")) +
  facet_wrap(~ Measured, labeller = to_string_trt) +
  theme_bw() +
  theme(text = element_text(size=10),
        #strip.text.x=element_text(size=10),
        strip.text.x=element_blank(),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        #legend.key = element_blank())
        legend.position = "bottom",
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "black") )


############
#Elytra
############
dataSum <- summarySE(diet_data, measurevar="ElytraLength_mm", groupvars=c("Measured", "Treatment", "Sex"))

#elytra length plot for both laboratory and field-caught beetles
ggplot(dataSum, aes(x=factor(Treatment), y=ElytraLength_mm, pch=Sex,
                    ymax=ElytraLength_mm + dataSum$se,
                    ymin=ElytraLength_mm - dataSum$se)) +
  xlab("Lab Starvation Treatment or Field Site") +
  ylab("Elytra length (mm)") + 
  geom_point(position=position_dodge(width=0.3), size=2.1) +
  geom_errorbar(position=position_dodge(width=0.3), width=0.21) +
  scale_x_discrete(breaks = c("a", "b", "c", "Washington_Co_2011", "TL_Davis_2011", "TL_Davis_2010"), labels=c("Fed AL", "Fed 2x", "Fed 0x", "FS1", "FS2", "FS3")) +
  facet_wrap(~ Measured, labeller = to_string_trt) +
  theme_bw() +
  theme(text = element_text(size=10),
        #strip.text.x=element_text(size=10),
        strip.text.x = element_blank(),
        strip.text.y=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        #legend.key = element_blank())
        legend.position = "bottom",
        legend.background = element_rect(size = 0.1, linetype = "solid", colour = "black") )

############
#Protein (no significant factors, so no plots)
############
newData <- field_data_no_tl2010[ which(field_data_no_tl2010$Measured=='PO_P'), ]
model <- lm(sqrt_p ~ Sex + Treatment + Sex:Treatment, data = newData, na.action=na.omit)
protein <- Anova(model, type=c(3))
protein.rg <- ref.grid(model)
lsmeans(protein.rg, "Sex")
lsmeans(protein.rg, "Treatment")
#produces similar results, slight rounding differences between SAS and R output
#conduct pairwise comparisons between treatments using the scheffe test
pairw.anova(y=newData$sqrt_p, x=newData$Treatment, method="scheffe")
pairw.anova(y=newData$sqrt_p, x=newData$Sex, method="scheffe")


############
#PO (no significant factors, so no plots)
############
model <- lm(ln_po ~ Sex + Treatment + Sex:Treatment, data = newData, na.action=na.omit)
po <- Anova(model, type=c(3))
po.rg <- ref.grid(model)
lsmeans(po.rg, "Sex")
lsmeans(po.rg, "Treatment")
#produces similar results, slight rounding differences between SAS and R output
#conduct pairwise comparisons between treatments using the scheffe test
pairw.anova(y=newData$ln_po, x=newData$Treatment, method="scheffe")
pairw.anova(y=newData$ln_po, x=newData$Sex, method="scheffe")


############
#melanizaiton (no significant factors, so no plots)
############
newData <- field_data_no_tl2010[ which(field_data_no_tl2010$Measured=='Melaniz'), ]
model <- lm(z_AGV ~ Sex + Treatment + Sex:Treatment, data = newData, na.action=na.omit)
melan <- Anova(model, type=c(3))
melan.rg <- ref.grid(model)
lsmeans(melan.rg, "Sex")
lsmeans(melan.rg, "Treatment")
#produces similar results, slight rounding differences between SAS and R output
#conduct pairwise comparisons between treatments using the scheffe test
pairw.anova(y=newData$z_AGV, x=newData$Treatment, method="scheffe")
pairw.anova(y=newData$z_AGV, x=newData$Sex, method="scheffe")



####
# Pairwise analysis between field & lab data for PROTEIN
####
newData <- diet_data[ which(diet_data$Measured=='PO_P'), ]
pairw.anova(y=newData$sqrt_p, x=newData$Treatment, method="scheffe")

#In R, how do I merge data from multiple groups into a new super group?
#I want all field data to be in a field_data group and all lab data and then compare the two against each other

####
# Pairwise analysis between field & lab data for PO
####
pairw.anova(y=newData$ln_po, x=newData$Treatment, method="scheffe")

#
# Pairwise for melanization
#
newData <- diet_data[ which(diet_data$Measured=='Melaniz'), ]
pairw.anova(y=newData$ln_melanization, x=newData$Treatment, method="scheffe")

newData <- newData[ which(newData$Sex=='Female'), ]
pairw.anova(y=newData$ln_melanization, x=newData$Treatment, method="scheffe")