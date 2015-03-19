#if you haven't already
#install.packages("dplyr")
#library(dplyr)
#install.packages("car")
#install.packages("asbio")
library(car)
library(lsmeans)
library(asbio)
require(graphics)
require(utils)
require(stats)

rm(list=ls())

#linux
#setwd('/home/rstanton/Documents/biology/writings/thesis/thesis_stats_and_comments/csv_data')

#mac
setwd('/Users/rstanton/Documents/thesis/N orbicollis data/Data Analysis/csv')

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
