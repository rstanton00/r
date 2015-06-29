#verify that necessary packages are installed
list.of.packages <- c("ggplot2", "Rmisc", "Rcpp", "car", "asbio", "lsmeans", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(dplyr)
library(car)
library(lsmeans)
library(ggplot2)
library(Rmisc)
require(graphics)
require(utils)
require(stats)
#library(asbio)

rm(list=ls())

#linux
#setwd('/home/rstanton/Documents/biology/writings/thesis/thesis_stats_and_comments/csv_data')

#mac
setwd('/Users/rstanton/Documents/pers/thesis/N orbicollis data/Data Analysis/csv')

diet_lab_data <- read.csv('no_diet_lab.csv', head=TRUE, stringsAsFactors=FALSE)

#tell R that certain data are categorical
diet_lab_data$Sex <- factor(diet_lab_data$Sex, levels = c("M", "F"), labels = c("Male", "Female"))
diet_lab_data$Treatment <- factor(diet_lab_data$Treatment, levels = c("a", "b", "c"))
diet_lab_data$Measured <- factor(diet_lab_data$Measured, levels = c("PO_P", "Melaniz"))

#remove known outliers
#assign rows to data, using dataset data, where column BeetleId in data != "delete"
diet_lab_data <- diet_lab_data[diet_lab_data$BeetleID != 'delete',]
diet_lab_data <- diet_lab_data[diet_lab_data$BeetleID != 'bF004',]
#diet_lab_data <- filter(diet_lab_data, !grepl('delete', BeetleID))

#add columns of transformed data
diet_lab_data["ln_premass_mg"] <- c(log2(diet_lab_data$PreMass_mg))
diet_lab_data["ln_postmass_mg"] <- c(log2(diet_lab_data$PostMass_mg))
diet_lab_data["ln_elytra"] <- c(log2(diet_lab_data$ElytraLength_mm))
diet_lab_data["ln_melanization"] <- c(log2(diet_lab_data$z_AGV))
diet_lab_data["ln_po"] <- c(log2(diet_lab_data$zPO_avg_abs_min))
diet_lab_data["sqrt_p"] <- c(log2(diet_lab_data$z_Protein_avg_adjusted_mg_ml))

#set options for all linear models/anovas
options(contrasts=c("contr.sum","contr.poly"))

#test for equal start condions for beetle age, body mass, elytra length, and body condition */
#All data together measured stands for "PO and Protein measurements" and "Filaments= Melanization" respectively */
#age
model <- lm(BeetleAge_Start_days ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = diet_lab_data, na.action=na.omit)
age.anova <- Anova(model, type=c(3))
age.rg <- ref.grid(model)
lsmeans(age.rg, "Measured")

#pre body mass
model <- lm(ln_premass_mg ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = diet_lab_data, na.action=na.omit)
prebodymass.anova <- Anova(model, type=c(3))
prebodymass.rg <- ref.grid(model)
lsmeans(prebodymass.rg, "Measured")

#elytra length
model <- lm(ln_elytra ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = diet_lab_data, na.action=na.omit)
elytra.anova <- Anova(model, type=c(3))
elytra.rg <- ref.grid(model)
lsmeans(elytra.rg, "Measured")

#plot elytra length box plots with SE bars
elSum <- summarySE(diet_lab_data, measurevar="ln_elytra", groupvars=c("Measured", "Treatment", "Sex"))


#ggplot(elytraSummary, aes(x=Measured, y=ln_elytra, fill=Sex)) +
#  geom_bar(position=position_dodge(), stat="identity") +
#  geom_errorbar(aes(ymin=ln_elytra-se, ymax=ln_elytra+se), width = .2, position=position_dodge(.9))

#elytraSummary <- summarySE(diet_lab_data, measurevar="ElytraLength_mm", groupvars=c("Sex", "Treatment", "Measured"))
#ggplot(elytraSummary, aes(x=c(Measured, Treatment), y=ElytraLength_mm, fill=Sex)) +
#  geom_bar(position=position_dodge(), stat="identity") +
#  geom_errorbar(aes(ymin=ElytraLength_mm-se, ymax=ElytraLength_mm+se), width = .2, position=position_dodge(.9))

#BEGIN EXAM HERE

#pre-condition calculated
model <- lm(resid_premass_calculated ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = diet_lab_data, na.action=na.omit)
residpremass.anova <- Anova(model, type=c(3))
residpremass.rg <- ref.grid(model)
lsmeans(residpremass.rg, "Sex")
lsmeans(residpremass.rg, "Measured")
lsmeans(residpremass.rg, "Treatment")
#produces similar results, slight rounding differences between SAS and R output
#conduct pairwise comparisons between treatments using the scheffe test
pairw.anova(y=diet_lab_data$resid_premass_calculated, x=diet_lab_data$Treatment, method="scheffe")
plot(pairw.anova(y=diet_lab_data$resid_premass_calculated, x=diet_lab_data$Treatment, method="scheffe"))

#post mass
model <- lm(ln_postmass_mg ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = diet_lab_data, na.action=na.omit)
postmass.anova <- Anova(model, type=c(3))
postmass.rg <- ref.grid(model)
lsmeans(postmass.rg, "Sex")
lsmeans(postmass.rg, "Measured")
lsmeans(postmass.rg, "Treatment")
#conduct pairwise comparisons between treatments using the scheffe test
pairw.anova(y=diet_lab_data$ln_postmass, x=diet_lab_data$Treatment, method="scheffe")
plot(pairw.anova(y=diet_lab_data$ln_postmass, x=diet_lab_data$Treatment, method="scheffe"))

#post condition calculated
model <- lm(resid_postmass_calculated ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = diet_lab_data, na.action=na.omit)
residpostmass.anova <- Anova(model, type=c(3))
residpostmass.rg <- ref.grid(model)
lsmeans(residpostmass.rg, "Sex")
lsmeans(residpostmass.rg, "Measured")
lsmeans(residpostmass.rg, "Treatment")
#conduct pairwise comparisons between treatments using the scheffe test
pairw.anova(y=diet_lab_data$resid_postmass_calculated, x=diet_lab_data$Treatment, method="scheffe")
plot(pairw.anova(y=diet_lab_data$resid_postmass_calculated, x=diet_lab_data$Treatment, method="scheffe"))

#test for effect of starvation treatment
df <- subset(diet_lab_data, select = c("BeetleID", "Sex", "Treatment", "Measured", "ln_premass_mg", "ln_postmass_mg"))
dfL <- reshape(df, varying = list(5:6), idvar = "BeetleID", direction = "long", v.names = "ln_mass_mg")
dfL$time <- factor(dfL$time, levels = 1:2, labels = c("pre", "post"))
head(dfL); tail(dfL)
aovMod <- aov(ln_mass_mg ~ Sex * Treatment * Measured * time + Error(BeetleID/time), data = dfL)

#repeated measures ANOVA for condition calculated
df <- subset(diet_lab_data, select = c("BeetleID", "Sex", "Treatment", "Measured", "resid_premass_calculated", "resid_postmass_calculated"))
dfL <- reshape(df, varying = list(5:6), idvar = "BeetleID", direction = "long", v.names = "resid_mass_calc")
dfL$time <- factor(dfL$time, levels = 1:2, labels = c("pre", "post"))
head(dfL); tail(dfL)
aovMod <- aov(resid_mass_calc ~ Sex * Treatment * Measured * time + Error(BeetleID/time), data = dfL)

#set up new dataset that only examines phenoloxidase/protein data for analysis
newData <- diet_lab_data[ which(diet_lab_data$Measured=='PO_P'), ]

#protein
model <- lm(sqrt_p ~ Sex + Treatment + Sex:Treatment, data = newData, na.action=na.omit)
protein <- Anova(model, type=c(3))
protein.rg <- ref.grid(model)
lsmeans(protein.rg, "Sex")
lsmeans(protein.rg, "Treatment")
#produces similar results, slight rounding differences between SAS and R output
#conduct pairwise comparisons between treatments using the scheffe test
pairw.anova(y=newData$sqrt_p, x=newData$Treatment, method="scheffe")
plot(pairw.anova(y=newData$sqrt_p, x=newData$Treatment, method="scheffe"))

#phenoloxidase
model <- lm(ln_po ~ Sex + Treatment + Sex:Treatment, data = newData, na.action=na.omit)
protein <- Anova(model, type=c(3))
po.rg <- ref.grid(model)
lsmeans(po.rg, "Treatment")
#produces similar results, slight rounding differences between SAS and R output
#conduct pairwise comparisons between treatments using the scheffe test
pairw.anova(y=newData$ln_po, x=newData$Treatment, method="scheffe")
plot(pairw.anova(y=newData$ln_po, x=newData$Treatment, method="scheffe"))

#model <- lm(cbind(ln_premass_mg, ln_postmass_mg) ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = diet_lab_data, na.action=na.omit)
#effect_of_treatment <- Anova(model, type=c(3))

#df <- subset(diet_lab_data, select = c("BeetleID", "Sex", "Treatment", "Measured", "ln_premass_mg", "ln_postmass_mg"))
#TRY TO TELL R THAT PRE & POST MASS ARE CORRELATED VIA LME?  idea from learning r for sas/spss section

#t.test(diet_lab_data$ln_postmass_mg, diet_lab_data$ln_premass_mg, paired=TRUE, data=diet_lab_data, na.action=na.omit)
#summary(lme(ln_postmass_mg ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = diet_lab_data, na.action=na.omit))