#update packages as necessary
require(utils)
update.packages(ask=FALSE)

#verify that necessary packages are installed
list.of.packages <- c("ggplot2", "Rcpp", "Rmisc", "car", "asbio", "lsmeans", "dplyr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(Rmisc)
library(dplyr)
library(car)
library(lsmeans)
library(asbio)
library(ggplot2)
require(stats)
require(graphics)

rm(list=ls())

#linux
setwd('/home/rstanton/Documents/biology/writings/thesis/thesis_stats_and_comments/csv_data')

#mac
#setwd('/Users/rstanton/Documents/pers/thesis/N orbicollis data/Data Analysis/csv')

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

#plot elytra length point plots with SE bars
#create SE measurements
elSum <- summarySE(diet_lab_data, measurevar="ElytraLength_mm", groupvars=c("Measured", "Treatment", "Sex"))

ggplot(elSum, aes(x=factor(Treatment), y=ElytraLength_mm, pch=Sex,
                  ymax=ElytraLength_mm + elSum$se,
                  ymin=ElytraLength_mm - elSum$se)) +
  geom_point(position=position_dodge(width=0.5), size = 5) +
  geom_errorbar(position=position_dodge(width=0.5)) +
  scale_x_discrete(breaks = c("a", "b", "c"), labels=c("Ad libitum", "3 days", "5 days")) +
  facet_wrap(~ Measured) +
  ylab("Elytra Length (mm)") + xlab("Treatment")

#ggplot(elSum, aes(x=factor(Measured), y=ElytraLength_mm, pch=Sex)) +
#  geom_point(position=position_dodge(width=0.5), size = 5) +
#  geom_errorbar(limits, position=position_dodge(width=0.5)) +
#  scale_x_discrete(breaks = c("PO_P", "Melaniz"), labels=c("PO", "Melanization")) +
#  facet_wrap(~ Treatment) +
#  ylab("Elytra Length (mm)") + xlab("")

#alternately, six different plots
#ggplot(elSum, aes(x=factor(Measured), y=ElytraLength_mm), group=Sex) +
#  geom_point() + geom_errorbar(limits, width=0.1) + facet_wrap(~ Treatment + Sex)

#plot elytra length box plots with SE bars
#elSum <- summarySE(diet_lab_data, measurevar="ln_elytra", groupvars=c("Measured", "Treatment", "Sex"))
#ggplot(diet_lab_data, aes(x=Measured, y=ln_elytra, col=Sex)) +
#  geom_boxplot() +
#  facet_wrap(~ Treatment)

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

#plot pre condition
#create SE measurements
elSum <- summarySE(diet_lab_data, measurevar="resid_premass_calculated", groupvars=c("Measured", "Treatment", "Sex"))
#create plot
ggplot(elSum, aes(x=factor(Treatment), y=resid_premass_calculated, pch=Sex,
                  ymax=resid_premass_calculated + elSum$se,
                  ymin=resid_premass_calculated - elSum$se)) +
  geom_point(position=position_dodge(width=0.5), size = 5) +
  geom_errorbar(position=position_dodge(width=0.5)) +
  scale_x_discrete(breaks = c("a", "b", "c"), labels=c("Ad libitum", "3 days", "5 days")) +
  facet_wrap(~ Measured) +
  ylab("Pre Body Mass (mg)") + xlab("Treatment")


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

#plot post condition
#create SE measurements
elSum <- summarySE(diet_lab_data, measurevar="resid_postmass_calculated", groupvars=c("Measured", "Treatment", "Sex"))
#create plot
ggplot(elSum, aes(x=factor(Treatment), y=resid_postmass_calculated, pch=Sex,
                  ymax=resid_postmass_calculated + elSum$se,
                  ymin=resid_postmass_calculated - elSum$se)) +
  geom_point(position=position_dodge(width=0.5), size = 5) +
  geom_errorbar(position=position_dodge(width=0.5)) +
  scale_x_discrete(breaks = c("a", "b", "c"), labels=c("Ad libitum", "3 days", "5 days")) +
  facet_wrap(~ Measured) +
  ylab("Post Body Mass (mg)") + xlab("Treatment")

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

#plot protein avg
#create SE measurements
elSum <- summarySE(diet_lab_data, measurevar="z_Protein_avg_adjusted_mg_ml", groupvars=c("Measured", "Treatment", "Sex"))
#create plot
ggplot(elSum, aes(x=factor(Treatment), y=z_Protein_avg_adjusted_mg_ml, pch=Sex,
                  ymax=z_Protein_avg_adjusted_mg_ml + elSum$se,
                  ymin=z_Protein_avg_adjusted_mg_ml - elSum$se)) +
  geom_point(position=position_dodge(width=0.5), size = 5, na.rm=TRUE) +
  geom_errorbar(position=position_dodge(width=0.5)) +
  scale_x_discrete(breaks = c("a", "b", "c"), labels=c("Ad libitum", "3 days", "5 days")) +
  ylab("Protein (mg/ml)") + xlab("Treatment")

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
plot(pairw.anova(y=newData$sqrt_p, x=newData$Treatment, method="scheffe"), main="Treatment Effect on Protein",
     xlab="Treatment", ylab="Protein mg/ml")

#plot phenoloxidase avg data
#create SE measurements
elSum <- summarySE(diet_lab_data, measurevar="zPO_avg_abs_min", groupvars=c("Measured", "Treatment", "Sex"))
#create plot
ggplot(elSum, aes(x=factor(Treatment), y=zPO_avg_abs_min, pch=Sex,
                  ymax=zPO_avg_abs_min + elSum$se,
                  ymin=zPO_avg_abs_min - elSum$se)) +
  geom_point(position=position_dodge(width=0.5), size = 5, na.rm=TRUE) +
  geom_errorbar(position=position_dodge(width=0.5)) +
  scale_x_discrete(breaks = c("a", "b", "c"), labels=c("Ad libitum", "3 days", "5 days")) +
  ylab("Phenoloxidase (abs/min)") + xlab("Treatment")

#phenoloxidase
model <- lm(ln_po ~ Sex + Treatment + Sex:Treatment, data = newData, na.action=na.omit)
protein <- Anova(model, type=c(3))
po.rg <- ref.grid(model)
lsmeans(po.rg, "Treatment")
#produces similar results, slight rounding differences between SAS and R output
#conduct pairwise comparisons between treatments using the scheffe test
pairw.anova(y=newData$ln_po, x=newData$Treatment, method="scheffe")
plot(pairw.anova(y=newData$ln_po, x=newData$Treatment, method="scheffe"), main="Treatment Effect on Phenoloxidase",
     xlab="Treatment", ylab="Phenoloxidase")

#plot melanization avg data
#create SE measurements
elSum <- summarySE(diet_lab_data, measurevar="z_AGV", groupvars=c("Measured", "Treatment", "Sex"))
#create plot
ggplot(elSum, aes(x=factor(Treatment), y=z_AGV, pch=Sex,
                  ymax=z_AGV + elSum$se,
                  ymin=z_AGV - elSum$se)) +
  geom_point(position=position_dodge(width=0.5), size = 5, na.rm=TRUE) +
  geom_errorbar(position=position_dodge(width=0.5)) +
  scale_x_discrete(breaks = c("a", "b", "c"), labels=c("Ad libitum", "3 days", "5 days")) +
  scale_y_continuous(trans = "reverse") +
  ylab("Average Grey Value") + xlab("Treatment")

#model <- lm(cbind(ln_premass_mg, ln_postmass_mg) ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = diet_lab_data, na.action=na.omit)
#effect_of_treatment <- Anova(model, type=c(3))

#df <- subset(diet_lab_data, select = c("BeetleID", "Sex", "Treatment", "Measured", "ln_premass_mg", "ln_postmass_mg"))
#TRY TO TELL R THAT PRE & POST MASS ARE CORRELATED VIA LME?  idea from learning r for sas/spss section

#t.test(diet_lab_data$ln_postmass_mg, diet_lab_data$ln_premass_mg, paired=TRUE, data=diet_lab_data, na.action=na.omit)
#summary(lme(ln_postmass_mg ~ Sex + Treatment + Measured + Sex:Treatment + Sex:Measured + Measured:Treatment + Sex:Treatment:Measured, data = diet_lab_data, na.action=na.omit))