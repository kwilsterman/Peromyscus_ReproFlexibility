
library(ggplot2)
library(ggbeeswarm)
library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)

b <-read.csv(file="LitterSize_fulldataset.csv",header=TRUE,sep=",")
b$Strain <- as.factor(b$Strain)
b$Species <- as.factor(b$Species)

# SUPPL TABLE 2 -----------------------------------------------------------
b_ls <- b[which(b$Parity <4),]

## parity impacts litter size within Peromyscus ##
full_model <- lmer(Npup ~ Parity + Species + (1|F_ID),
              data = b_ls)
qqnorm(resid(full_model))
qqline(resid(full_model))
anova(full_model)
summary(full_model)

## strain-specific parity effects on LS ##
b_ls$Parity <- as.factor(b_ls$Parity)

strain <- b_ls[which(b_ls$Strain=="GOS"),]
strain_specific_model <- lmer(Npup ~ Parity + (1|F_ID),
              data = strain)
qqnorm(resid(strain_specific_model))
qqline(resid(strain_specific_model))
anova(strain_specific_model)
pairs(emmeans(strain_specific_model, ~Parity))

## strains differ in mean littersize (maniculatus only) ##
pman <- b_ls[which(b_ls$Species=="PEMA"),]
pman_model <- lmer(Npup ~ Parity + Strain + (1|F_ID),
                              data = pman)
qqnorm(resid(pman_model))
qqline(resid(pman_model))
anova(pman_model)
pairs(emmeans(pman_model, ~Strain))


# SUPPL TABLE 3 -----------------------------------------------------------
env <- read.csv(file="Env_LitterSize_LabOnly.csv",header=TRUE,sep=",")

cor.test(env$LitterSize.Mean, env$TEMP_DegreeOfSeasonality)

# SUPPL TABLE 5 -----------------------------------------------------------

bIBI <- b[which(b$Parity <10),]
bIBI <- bIBI[which(bIBI$IBI < 45),]

bIBI_pema <- bIBI[which(bIBI$Species == "PEMA"),]
bIBI_pema <- bIBI_pema[-which(bIBI_pema$Strain == "AUS"),]
bIBI_pema <- bIBI_pema[-which(bIBI_pema$Strain == "SAT"),]
m3 <- lmer(IBI ~ Parity + Npup + PrevLitterSize + Strain + (1|F_ID),
           data = bIBI_pema)
summary(m3)
anova(m3)
pairs(emmeans(m3, ~Strain))

bIBI_pele <- bIBI[which(bIBI$Species == "PELE"),]
m4 <- lmer(IBI ~ Parity + Npup + PrevLitterSize + Strain + (1|F_ID),
           data = bIBI_pele)
summary(m4)
anova(m4)

bIBI_peca <- bIBI[which(bIBI$Species == "PECA"),]
m5 <- lmer(IBI ~ Parity + Npup + PrevLitterSize + (1|F_ID),
           data = bIBI_peca)
summary(m5)
anova(m5)

bIBI_pepo <- bIBI[which(bIBI$Species == "PEPO"),]
m6 <- lmer(IBI ~ Parity + Npup + PrevLitterSize  + Strain + (1|F_ID),
           data = bIBI_pepo)
summary(m6)
anova(m6)


# SUPPL TABLE 6 -----------------------------------------------------------
##stock center strains have << increase in littersize across parities compared with WD strains
b_ls_origin <- b_ls[-which(b_ls$Species == "PECA"),]
b_ls_origin <- b_ls_origin[-which(b_ls_origin$Species == "PMEL"),]
b_ls_origin <- b_ls_origin[-which(b_ls_origin$Species == "PEGO"),]

strainorigin_model <- lmer(Npup ~ StrainType + Parity + Species + (1|F_ID),
              data = b_ls_origin)
qqnorm(resid(strainorigin_model))
qqline(resid(strainorigin_model))
anova(strainorigin_model)
pairs(emmeans(strainorigin_model, ~ StrainType))


b_ibi_origin <- bIBI[-which(bIBI$Species == "PECA"),]
b_ibi_origin <- b_ibi_origin[-which(b_ibi_origin$Species == "PMEL"),]
b_ibi_origin <- b_ibi_origin[-which(b_ibi_origin$Species == "PEGO"),]

strainorigin_model <- lmer(IBI ~ StrainType + Species + Parity + Npup + PrevLitterSize  + (1|F_ID),
                           data = b_ibi_origin)
qqnorm(resid(strainorigin_model))
qqline(resid(strainorigin_model))
anova(strainorigin_model)
pairs(emmeans(strainorigin_model, ~ StrainType))


# SUPPL TABLE 7 -----------------------------------------------------------
##BW
strain <- b_ls[which(b_ls$Strain == "BW"),]
strain$LOC <- as.factor(strain$LOC)

model <- lmer(Npup ~ LOC + Parity + (1|F_ID), data = strain)
anova(model)

##LL
strain <- b[which(b$Strain == "LL"),]
strain$LOC <- as.factor(strain$LOC)

model <- lm(Npup ~ LOC, data = strain[which(strain$Parity == 1),])
anova(model)

