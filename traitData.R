setwd("~/Documents/PADI/traitdata/")
library(reshape2)
library(tidyverse)
library(MCMCglmm)
library(readxl)

# read in data -----
lipid0 <- read_excel("PADI_data_July2017.xlsx", sheet = 1, na = "NA")
carb0 <- read_excel("PADI_data_July2017.xlsx", sheet = 2, na = "NA")
prot0 <- read_excel("PADI_data_July2017.xlsx", sheet = 3, na = "NA")
color0 <- read_excel("PADI_data_July2017.xlsx", sheet = 5, na = "NA")
sa0 <- read_excel("PADI_data_July2017.xlsx", sheet = 6, na = "NA")
weight0 <- read_excel("PADI_data_July2017.xlsx", sheet = 7, na = "NA")
dryweight0 <- read_excel("PADI_data_July2017.xlsx", sheet = 8, na = "NA")
antibac0 <- read_excel("PADI_data_July2017.xlsx", sheet = 4, na = "NA")

#
antibac <- antibac0 %>% mutate(geno = factor(geno))
summary(antibac)

lipid <- lipid0 %>% select(c(geno,treat,lipidConc_mg_mL)) %>% mutate(geno = factor(geno))
summary(lipid)

carb <- carb0 %>% select(c(geno,treat,carbConc_mg_mL)) %>% mutate(geno = factor(geno))
summary(carb)

prot <- prot0 %>% filter(protExtraction=="1" & protZoox=="n") %>% select(c(geno,treat,protConc_mg_mL)) %>% mutate(geno = factor(geno))
summary(prot)

# for color, average the two sides
color <- color0 %>% mutate(geno = factor(geno)) %>% group_by(geno,treat) %>% summarize(meanColorChange = mean(colorChange))
summary(color)

sa <- sa0 %>% select(c(geno,treat,saCm2)) %>% mutate(geno = factor(geno))
summary(sa)

weight_3oct <- weight0 %>% filter(wetweight_volDate=="3oct") %>% select(c(geno,treat,wetweight_weightMucus_g, wetweight_volMucus_mL)) %>% rename(wetweightg_3oct = wetweight_weightMucus_g, wetvolumemL_3oct = wetweight_volMucus_mL) %>% mutate(geno = factor(geno))
summary(weight_3oct)

weight_9oct <- weight0 %>% filter(wetweight_volDate=="9oct") %>% select(c(geno,treat,wetweight_weightMucus_g, wetweight_volMucus_mL)) %>% rename(wetweightg_9oct = wetweight_weightMucus_g, wetvolumemL_9oct = wetweight_volMucus_mL) %>% mutate(geno = factor(geno))
summary(weight_9oct)

dryweight <- dryweight0 %>% mutate(dryweight = drymass_wetTubeWeight-drymass_afterSpeedvacWeight) %>% select(c(geno,treat,dryweight))
summary(dryweight)


# make big data table ------
dats <- merge(lipid,carb,by=c("geno","treat"), all = T)
dats <- merge(dats,prot,by=c("geno","treat"), all = T)
dats <- merge(dats,color,by=c("geno","treat"), all = T)
dats <- merge(dats,sa,by=c("geno","treat"), all = T)
dats <- merge(dats,weight_3oct,by=c("geno","treat"), all = T)
dats <- merge(dats,weight_9oct,by=c("geno","treat"), all = T)
dats <- merge(dats,dryweight,by=c("geno","treat"), all = T)
summary(dats)
head(dats)

# normalize lipid, carb, prot to weight and surface area
alldata <- dats %>% 
  mutate(lipid_ug_cm2 = (1000*lipidConc_mg_mL*wetvolumemL_3oct)/saCm2, 
         carb_ug_cm2 = (1000*carbConc_mg_mL*wetvolumemL_3oct)/saCm2, 
         prot_mg_cm2 = (protConc_mg_mL*wetvolumemL_3oct)/saCm2, 
         weight_g_cm2_3oct = wetweightg_3oct/saCm2, 
         weight_g_cm2_9oct = wetweightg_9oct/saCm2, 
         wetvolumemL_cm2_3oct = wetvolumemL_3oct/saCm2,
         wetvolumemL_cm2_9oct = wetvolumemL_9oct/saCm2) 

head(alldata)

# replace negative values with NA
summary(alldata)
alldata$lipid_ug_cm2[alldata$lipid_ug_cm2<0] <- NA
alldata$carb_ug_cm2[alldata$carb_ug_cm2<0] <- NA

# save data ------
# save(alldata, file="alldata.Rdata")
load("alldata.Rdata")
head(alldata)

# just get what I want ....
# alldata <- alldata %>% select(geno, treat, meanColorChange, lipid_ug_cm2, carb_ug_cm2, prot_mg_cm2, weight_g_cm2_3oct, weight_g_cm2_9oct)
# head(alldata)
# summary(alldata)

# Compare mucus production across days -----
head(alldata)

production_g <- alldata %>% select(geno,treat,weight_g_cm2_3oct, weight_g_cm2_9oct) %>% rename(initial=weight_g_cm2_3oct, recovery=weight_g_cm2_9oct)
head(production_g)

pg <- melt(production_g, id=c("geno","treat"))
head(pg)

production_mL <- alldata %>% select(geno,treat,wetvolumemL_cm2_3oct, wetvolumemL_cm2_9oct) %>% rename(initial=wetvolumemL_cm2_3oct, recovery=wetvolumemL_cm2_9oct)
head(production_mL)

perc_production_change <- production_mL %>% mutate(perc_production_change = (recovery/initial)*100)
mean(perc_production_change$perc_production_change, na.rm=T)
sd(perc_production_change$perc_production_change, na.rm=T)

pmL <- melt(production_mL, id=c("geno","treat"))
head(pmL)

pmL$time.treat <- interaction(pmL$treat, pmL$variable)
pg$time.treat <- interaction(pg$treat, pg$variable)

ggVolume = ggplot(pmL, aes(x=time.treat, y=value, fill=treat))+
  scale_fill_manual(values=c("brown","tan"))+
  geom_boxplot()+
  geom_jitter(width = 0.1)+
  ylab("Volume Mucus Produced (mL/ m2)")+
  theme_bw()

ggWeight = ggplot(pg, aes(x=time.treat, y=value, fill=treat))+
  scale_fill_manual(values=c("brown","tan"))+
  geom_boxplot()+
  geom_jitter(width = 0.1)+
  ylab("Weight Mucus Produced (g/ m2)")+
  theme_bw()

gridExtra::grid.arrange(ggVolume, ggWeight, ncol=1)

set.seed(1)
mcmcVolume <- MCMCglmm(value~treat*variable, random=~geno, data = pg)
summary(mcmcVolume)
#                           post.mean  l-95% CI  u-95% CI eff.samp  pMCMC    
# (Intercept)              0.188286  0.172894  0.202811   1000.0 <0.001 ***
# treath                   0.012247 -0.010498  0.031847    920.2  0.248    
# variablerecovery        -0.104986 -0.124818 -0.083486   1122.8 <0.001 ***
# treath:variablerecovery -0.002847 -0.032096  0.029113   1000.0  0.836 
set.seed(1)
mcmcWeight <- MCMCglmm(value~treat*variable, random=~geno, data = pmL)
summary(mcmcWeight)
#                         post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)               2.07945  1.85564  2.30824     1000 <0.001 ***
# treath                    0.13867 -0.10765  0.39218     1000  0.274    
# variablerecovery         -1.09658 -1.33777 -0.85875     1123 <0.001 ***
# treath:variablerecovery  -0.04156 -0.40605  0.31075     1000  0.790   

# Stats ------

hist(alldata$meanColorChange)
hist(alldata$lipid_ug_cm2)
hist(alldata$prot_mg_cm2)
hist(alldata$carb_ug_cm2)

# remove strangly large carb value
alldata$carb_ug_cm2[alldata$carb_ug_cm2>70] <- NA

# color change
mcmcColor <- MCMCglmm(meanColorChange~treat, random=~geno, data = alldata)
summary(mcmcColor)
#                 post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)    0.8356   0.6341   1.0466   1238.9 <0.001 ***
# treath        -1.1618  -1.4448  -0.8794    910.8 <0.001 ***

# protein
mcmcProtein = MCMCglmm(prot_mg_cm2~treat, random = ~geno, data = alldata)
summary(mcmcProtein)
#                 post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)     1.400    1.012    1.876     1000 <0.001 ***
# treath          2.059    1.490    2.665     1000 <0.001 ***

# carb
mcmcCarb = MCMCglmm(carb_ug_cm2~treat, random = ~geno, data = alldata)
summary(mcmcCarb)
#             post.mean l-95% CI u-95% CI eff.samp pMCMC
# (Intercept)    1.3386   0.6152   1.9909    756.9 <0.001 ***
# treath         0.6412  -0.1054   1.3495    367.7  0.102 

# lipid
mcmcLipid = MCMCglmm(lipid_ug_cm2~treat, random=~geno, data=alldata)
summary(mcmcLipid)
#               post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)    20.677   11.756   28.331     1000 <0.001 ***
# treath         15.715    3.028   27.992     1128  0.016 * 

# plot --------------

ggColor = ggplot(alldata, aes(x=treat, y=meanColorChange, fill=treat, group=treat))+
  scale_fill_manual(values=c("brown","tan"))+
  geom_boxplot()+
  geom_jitter(width = 0.1)+
  annotate("text", x = 1, y = -1.4, label = "treat p < 0.001", size = 3)+
  ylab("Color Change (Final - Initial)")+
  theme_bw()

ggProtein = ggplot(alldata, aes(x=treat, y=prot_mg_cm2, fill=treat, group=treat))+
  scale_fill_manual(values=c("brown","tan"))+
  geom_boxplot()+
  geom_jitter(width = 0.1)+
  annotate("text", x = 2, y = 0.5, label = "treat p < 0.001", size = 3)+
  ylab("Protein (mg/cm2)")+
  theme_bw()

ggCarb = ggplot(alldata, aes(x=treat, y=carb_ug_cm2, fill=treat, group=treat))+
  geom_boxplot()+
  scale_fill_manual(values=c("brown","tan"))+
  geom_jitter(width = 0.1)+
  annotate("text", x = 2, y = 0, label = "treat p = 0.100", size = 3)+
  ylab("Carb (ug/cm2)")+
  theme_bw()

ggLipid = ggplot(alldata, aes(x=treat, y=lipid_ug_cm2, fill=treat, group=treat))+
  geom_boxplot()+
  scale_fill_manual(values=c("brown","tan"))+
  geom_jitter(width = 0.1)+
  annotate("text", x = 2, y = 0, label = "treat p < 0.001", size = 3)+
  ylab("Lipid (ug/cm2)")+
  theme_bw()

quartz()
gridExtra::grid.arrange(ggColor, ggProtein, ggLipid, ggCarb, ncol=2)

# ANTIBAC ------
d0 <- melt(antibac, id=c("well", "geno", "treat", "rep"))
head(d0)
summary(d0)

# calculate coefficient of variance (sd/mean*100) ------
sd = dcast(d0, geno+treat ~ variable, sd)
means = dcast(d0, geno+treat ~ variable, mean)
names(sd)
cv = sd[c(3:27)]/means[c(3:27)]*100
head(cv)
table(cv>20) # how many samples have CV > 20%?
# FALSE 800

# remove the negative control (ASW- "0A M) and positive (LB+ "0L P") controls ----
head(d0)
d <- d0 %>% filter(treat==c("C", "H")) %>% filter(!geno=="0L")
summary(d)

# average replicates and plot --------
head(d)

df <- d %>% group_by(geno,treat,variable) %>% summarise(mean.value = mean(value), sd = sd(value))
head(df)

df.treat <- d %>% group_by(treat,variable) %>% summarise(mean.value = mean(value), sd = sd(value))
head(df.treat)
limits <- aes(ymax = mean.value + sd, ymin=mean.value - sd)

gg1 <- ggplot(df, aes(x=variable, y=mean.value, color=treat)) +
  geom_point() +
  geom_line(aes(group=treat))+
  geom_text(aes(label=geno), nudge_x = 0.0001)+
  theme_bw()
gg1

gg2 <- ggplot(df.treat, aes(x=variable, y=mean.value, color=treat)) +
  geom_point() +
  geom_line(aes(group=treat))+
  geom_errorbar(limits, width=0.2)+
  ylab("OD600")+
  xlab("time (hr)")+
  scale_color_manual(values=c("brown","tan", "black"))+
  theme_bw()
quartz()
gg2

# it looks like the wells were dry by hour 5. Remove all timepoints before then
good_times <- c("0", "0.5", "1", "1.5", "2", "2.5", "3", "3.5", "4", "4.5","5")
dat <- d %>% filter(variable %in% good_times )
summary(dat)

# stats for time series data ---------
summary(dat)
library(nlme)

model <- lme(value ~ treat*variable, 
             random= ~ 1|well, 
             method = "ML",
             data = dat)
summary(model)
                    # Value   Std.Error  DF   t-value p-value
# treatH:variable0.5 -0.01624609 0.002588309 410  -6.27672  0.0000
# treatH:variable1   -0.01709152 0.002588309 410  -6.60335  0.0000
# treatH:variable1.5 -0.01928565 0.002588309 410  -7.45106  0.0000
# treatH:variable2   -0.01676174 0.002588309 410  -6.47594  0.0000
# treatH:variable2.5 -0.01249326 0.002588309 410  -4.82680  0.0000
# treatH:variable3   -0.00968826 0.002588309 410  -3.74308  0.0002
# treatH:variable3.5 -0.00818870 0.002588309 410  -3.16372  0.0017
# treatH:variable4   -0.00711109 0.002588309 410  -2.74739  0.0063
# treatH:variable4.5 -0.00603935 0.002588309 410  -2.33332  0.0201
# treatH:variable5   -0.00619543 0.002588309 410  -2.39362  0.0171

anova(model)
#                 numDF denDF   F-value p-value
# (Intercept)        1   410 22106.302  <.0001
# treat              1    41     8.148  0.0067
# variable          10   410    67.439  <.0001
# treat:variable    10   410    10.802  <.0001

# qpcr --------
library(MCMC.qpcr)

qpcr0 <- read_excel("PADI_data_July2017.xlsx", sheet = 9, na = "NA")
head(qpcr0)

eff3 <- read.table("PrimerEfficienciesITS2AcerActin.txt",header=T,sep="\t") 
eff3
eff <- PrimEff(eff3)

mL_mucus <- 14/1000

d <- qpcr0
head(d)
co <- cq2counts(data=d, genecols=c(3:4), condcols=c(1:2), effic=eff)
head(co)
co$sample <- row.names(co)

# Summarize counts per mL mucus
co %>% ungroup() %>% group_by(gene,treat) %>% summarize(log(mean(count))/mL_mucus)

countsActin <- co %>% filter(gene=="actin")
head(countsActin)

countsITS2 <- co %>% filter(gene=="its2")
head(countsITS2)

# stats for qpcr-----
set.seed(1)
mActin <- MCMCglmm(count~treat,
                   random=~geno,
                   family="poisson",
                   nitt=20000,
                   data=countsActin)
summary(mActin)
#             post.mean l-95% CI u-95% CI eff.samp   pMCMC   
#(Intercept)   -0.3427  -1.5151   0.7446     1256 0.56353   
#treath         2.9740   1.6044   4.5663     1489 0.00118 **

set.seed(1)
mITS2 <- MCMCglmm(count~treat,
                   random=~geno,
                   family="poisson",
                   nitt=20000,
                   data=countsITS2)
summary(mITS2)
#                 post.mean l-95% CI u-95% CI eff.samp  pMCMC    
# (Intercept)   -15.411  -26.268   -4.881    1.471 <6e-04 ***
# treath         14.537    4.265   25.742    2.183 <6e-04 ***

# plot qpcr -----
ggActin = ggplot(countsActin, aes(x=treat, y=log(count,2)/mL_mucus, fill=treat, group=treat))+
  geom_boxplot()+
  scale_y_continuous(trans='log2')+
  scale_fill_manual(values=c("brown","tan"))+
  geom_jitter(width = 0.1)+
  annotate("text", x = 2, y = 8, label = "treat p = 0.001", size = 3)+
  ylab("Log2 Actin Counts")+
  ggtitle("Actin Counts per mL Mucus")+
  theme_bw()
ggActin

ggITS2 = ggplot(countsITS2, aes(x=treat, y=log(count,2)/mL_mucus, fill=treat, group=treat))+
  geom_boxplot()+
  scale_y_continuous(trans='log2')+
  scale_fill_manual(values=c("brown","tan"))+
  geom_jitter(width = 0.1)+
  annotate("text", x = 1, y = 4, label = "treat p < 0.001", size = 3)+
  ylab("Log2 ITS2 Counts")+
  ggtitle("ITS2 Counts per mL Mucus")+
  theme_bw()
ggITS2

gridExtra::grid.arrange(ggActin,ggITS2,ncol=2)