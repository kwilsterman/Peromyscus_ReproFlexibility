setwd("~/Documents/Research/Permyscus_ReproFlexibility/TO SUBMIT")
library(ggplot2)
library(ggbeeswarm)
library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)

b <-read.csv(file="LitterSize_fulldataset.csv",header=TRUE,sep=",")
b$Strain <- as.factor(b$Strain)
b$Species <- as.factor(b$Species)

b_ls <- b[which(b$Parity <4),]

##change parity for summary table
b_lw$Parity <- as.factor(b_ls$Parity)
LSSummaryTable <- b_ls %>% 
  group_by(Species, Strain, Parity) %>%
  summarize(NLitterSize = n(), 
            NDams = n_distinct(F_ID, na.rm=TRUE),
            medianLitterSize = median(Npup, na.rm = TRUE),
            averageLitterSize = mean(Npup, na.rm = TRUE),
            maxLitterSize = max(Npup, na.rm = TRUE),
            minLitterSize = min(Npup, na.rm = TRUE)
  )

ggplot() +
  geom_violin(data = b2wd, aes(x=as.factor(Parity), y = Npup, color = Strain), adjust = 2, scale = "count") +
  geom_point(data = LSSummaryTable, aes(x=as.factor(Parity), y = averageLitterSize, color = Strain)) + 
  facet_wrap(~Strain) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


model <- lmer(Npup ~ StrainType*Parity + (1|Strain),
              data = b2wd)
qqnorm(resid(model))
qqline(resid(model))
anova(model)
pairs(emmeans(model, ~ StrainType))

LL <- b2wd[which(b2wd$Strain == "LL"),]

##stock center strains have << inccrase in littersize across parities compared with WD strains

model <- lmer(Npup ~ Strain + Parity + (1|F_ID),
              data = b)
qqnorm(resid(model))
qqline(resid(model))
anova(model)
pairs(emmeans(model, ~Strain))

##var. amoing strains in litter sizes

strain <- b2wd[which(b2wd$Strain=="GOS"),]
model <- lmer(Npup ~ Parity + (1|F_ID),
              data = strain)
qqnorm(resid(model))
qqline(resid(model))
anova(model)
pairs(emmeans(model, ~Parity))

#strains w/ sig var in litter size among parities



b <-read.csv(file="LitterSize_sh_072221.csv",header=TRUE,sep=",")
bIBI <- b[which(b$Parity <10),]
bIBI <- bIBI[-which(bIBI$Strain == "HL"),]
bIBI$Strain <- as.factor(bIBI$Strain)
bIBI$Species <- as.factor(bIBI$Species)

bIBI <- bIBI[which(bIBI$IBI < 45),]

IBISummaryTable <- bIBI %>% 
  group_by(Species, Strain) %>%
  summarize(NIBI = n(),
            NDam = n_distinct(F_ID, na.rm=TRUE),
            medianIBI = median(IBI, na.rm = TRUE),
            averageIBI = mean(IBI, na.rm = TRUE),
            maxIBI = max(IBI, na.rm = TRUE),
            minIBI = min(IBI, na.rm = TRUE)
  )
IBISummaryTable_sub <- IBISummaryTable[-which(IBISummaryTable$NIBI < 10),]

ggplot() +
  geom_violin(data = bIBI, aes(x=Strain, y = IBI, color = Strain), adjust = 2) +
  geom_point(data = IBISummaryTable_sub, aes(x=Strain, y = medianIBI, color = Strain)) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


m2 <- lmer(IBI ~ Parity + Npup + PrevLitterSize + Species + (1|Strain/F_ID),
           data = bIBI)
summary(m2)
anova(m2)

bIBI_pema <- bIBI[which(bIBI$Species == "PEMA"),]
bIBI_pema <- bIBI_pema[-which(bIBI_pema$Strain == "SW"),]
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



pairs(emmeans(m2, ~Species))
plot(emmeans(m2, ~Parity))
plot(emmeans(m2, ~Npup))
plot(emmeans(m2, ~PrevLitterSize))



qqnorm(resid(m2))
qqline(resid(m2))


strain <- b2wd[which(b2wd$Strain == "LL"),]
strain$LOC <- as.factor(strain$LOC)

summary(strain$LOC[which(strain$Parity == 1)])

model <- lm(Npup ~ LOC, data = strain[which(strain$Parity == 1),])
anova(model)
plot(emmeans(model, ~LOC))

##we find no effect of location, where breeding conditions vary among labs, on litter size in LL, BW, ME or LN lines
#however, sample size is small for these analyses and restricted to first births, where we have the most
#comparable data aong labs. Better analyses would also look at time to breed after pairing, litter IBI, 
#and litter size across parities

me <- b2wd[which(b2wd$Strain == "ME"),]
me$LOC <- as.factor(me$LOC)
me <- me[-which(me$LOC == "IL"),]

plot(me$IBI ~ me$LOC)
t <- lmer(Npup ~ LOC + Parity + (1|F_ID), data = me)
anova(t)
pairs(emmeans(t, ~LOC))
