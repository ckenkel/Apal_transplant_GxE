#Morphological Plasticity Analysis
# by Holland Elder adapted for Acropora palmata from Million 2021

#install.packages("gdata")
#utils::install.packages("gap")
#utils::install.packages("lme4")#solving Matrix-glmer incompatibility for R v 4.3.2; say no to installing from source
library(tidyr)
library(dplyr)
options(dplyr.summarise.inform=F)
library(ggplot2)
library(lme4)
library(lmerTest)
library(multcomp)
library(pbkrtest)
library(ggfortify)
#library(gridExtra)
library(ggpubr)
library(lsmeans)
#library(tidyverse)
library(reshape2)
library(naniar)
library(piecewiseSEM)
library(foreach)
library(car)
library(agricolae)
library(DHARMa) 

#Start here with some checks
grow<-read.csv("~/Dropbox/CarlsLab/ResearchProjects/NOAA_CRCP/PALMATA/ApalmMorphologicalPlasticity_rready.csv", header = T,stringsAsFactors=TRUE)
str(grow)
summary(grow)

# Convert numbers into factors for these columns, and generate faactor for nested arrays within sites  
grow$Tag<-as.factor(grow$Tag)
grow$Array<-as.factor(grow$Array)
grow$FragID<-as.factor(grow$FragID)
grow$ArrayInSite<-with(grow,factor(Site:Array))

#grow$TissueSampledT0<-as.factor(grow$TissueSampledT0)
#grow$TissueSampledT24<-as.factor(grow$TissueSampledT24)


#relevel sites from East to West (can choose to order in anyway, eventually they will be ordered by survival)
grow= grow %>% 
  mutate(Site = factor(Site, levels=c("Eastern Sambo","Marker 32","Western Sambo","Big Pine","Dave's Ledge","Looe Key","Eastern Dry Rocks","Maryland Shoals","Bahia Honda")))
grow= grow %>%
  mutate(Genotype=factor(Genotype,levels=c("AP13-X5","AP13-X7","AP13-X9","AP13-XK","AP14-5")))


#how many measurements do we have for each time point
sum(!is.na(grow$T0_SA)) #135
sum(!is.na(grow$T12_SA)) #35
sum(!is.na(grow$T15_SA)) #107
sum(!is.na(grow$T18_SA)) #104
sum(!is.na(grow$T24_SA)) #92

#A common mistake when phenotyping is recording V as the SA of a mesh with closed holes. 
#If this is the case, V will be greater than SA but this shouldn't happen just based on these coral sizes.
#this should be 0, just as a check
sum(grow$T24_SA<grow$T24_V,na.rm = T) #looks good for all time points

# Imputation would not be valid in an Acropora palmata data set. Can test some correlations perhaps at a later date but not now.
# Breakage is also not relevant for Acropora palmata.
# So calculations for these steps have been removed.


###################################### Survival Analysis #############################
#using coxPH ratios
library(survival) # for survival analysis 
library(vegan) # for plotting ellipses in principal coordinates analysis plots
#library(corr) # for correlations, not available for R v 3.6.3
#library(MCMCglmm) # for mcmcglmm stats
library(nlme) #for lme
library(MASS) # for stepAIC
library(coxme) # for mixed effects cox model
library(ggridges) # for ridge plots
library(corrplot) # for correlations
#library(summarytools) # for dfSummary
library(cowplot)
library(RColorBrewer)
library(survminer)
library(psych)
library(gdata)


surv<-grow
str(surv)


# Classify this column as numeric 0=alive and 1=dead
surv$T24_TOD<-as.numeric(surv$T24_TOD_mod)

####Making a survival object where TOD is time of death event, and Dead is whether or not a fragment was dead
score<-surv[complete.cases(surv$T24_TOD),] #removing ramets Missing from outplant, note there are none

score$Site<-relevel(score$Site,ref="Eastern Sambo")## picked site with best survival as reference (want to know how much worse others are)
score$Genotype<-relevel(score$Genotype,ref="AP14-5")## picked genet with best survival as reference

died=Surv(score$T24_TOD, score$T24_Dead)

####stats


score$Site<-factor(score$Site,levels=c("Eastern Sambo","Western Sambo","Dave's Ledge","Looe Key","Maryland Shoals","Eastern Dry Rocks","Marker 32","Big Pine","Bahia Honda"))## picked site with best survival as reference (want to know how much worse others are)
scurveS=survfit(died~Site, data=score)
scurveS
#                        n events median 0.95LCL 0.95UCL
#Site=Eastern Sambo     15      0     NA      NA      NA
#Site=Western Sambo     15      1     NA      NA      NA
#Site=Dave's Ledge      15      1     NA      NA      NA
#Site=Looe Key          15      2     NA      NA      NA
#Site=Maryland Shoals   15      3     NA      NA      NA
#Site=Eastern Dry Rocks 15      5     NA       6      NA
#Site=Marker 32         15      6     NA      20      NA
#Site=Big Pine          15      9     20      20      NA
#Site=Bahia Honda       15     15     12      12      NA

score$Genotype<-factor(score$Genotype,levels=c("AP14-5","AP13-X5","AP13-X7","AP13-XK","AP13-X9"))## picked genet with best survival as reference
scurveG<- survfit(died~Genotype,data=score)
scurveG

#                  n events median 0.95LCL 0.95UCL
#Genotype=AP14-5  27      7     NA      NA      NA
#Genotype=AP13-X5 27      8     NA      NA      NA
#Genotype=AP13-X7 27      8     NA      NA      NA
#Genotype=AP13-XK 27      9     NA      NA      NA
#Genotype=AP13-X9 27     10     NA      20      NA

#to see which genotypes survive best at each site
scurveGandE<- survfit(died~Site+Genotype,data=score)
scurveGandE

##- identify model that accounts for full effects and tests for potential interaction
control<-coxme.control(iter.max=30)
m1<-coxme(died~Site*Genotype+(1|ArrayInSite)+(1|UniqueID), data=score,control=control) #site by genet effects controlling for array and repeated measures

#getting warning that coefficients don't converge - do have some very large coefficients
#and looks like risk estimate is not jiving with counts above - e.g.
#very high risk for Looe Key, on par with Bahia Honda
#general effects okay, but will have to use releveling and subsetting to pull out proper pairwise coefficients

summary(m1)
m1$loglik #Integrated
-2*-150.4977   #300.9954

m2<-coxme(died~Site+Genotype+(1|ArrayInSite)+(1|UniqueID), data=score,control=control) #site and genet effects controlling for array and repeated measures

summary(m2)
m2$loglik
-2*-166.6199 #333.2398 #marginal change in loglik

m3<-coxme(died~Site+(1|ArrayInSite)+(1|UniqueID), data=score) #still getting warning about infinite coefficient - will have to make WS reference for site coefficients

summary(m3)
m3$loglik
-2*-166.6556 #333.3112 genotype doesn't change log likelihood at all

m4<-coxme(died~Genotype+(1|ArrayInSite)+(1|UniqueID), data=score) #still getting warning about infinite coefficient - will have to make WS reference for site coefficients

summary(m4)
m4$loglik
-2*-179.6147 #359.2294 effect of site

m0<-coxme(died~1+(1|ArrayInSite)+(1|UniqueID), data=score) #null model
summary(m0)

anova(m0,m3) #big improvement in log lik with site
anova(m0,m4) #genet no effect
anova(m3,m2) #meh, adding genet doesn't help
anova(m2,m1) #nor does addition of interaction term

m1a<-coxme(died~Site+Site:Genotype+(1|ArrayInSite)+(1|UniqueID), data=score,control=control) #excluding fixed genet effect
summary(m1a)

anova(m3,m1a) #again no improvement with addition of interaction term

#lets try to get some resolution on site coefficients - need to relevel- put ES at end

score$Site<-factor(score$Site,levels=c("Western Sambo","Dave's Ledge","Looe Key","Maryland Shoals","Eastern Dry Rocks","Marker 32","Big Pine","Bahia Honda","Eastern Sambo"))## picked site with best survival as reference (want to know how much worse others are)

m3a<-coxme(died~Site+(1|ArrayInSite)+(1|UniqueID), data=score) #still getting warning about infinite coefficient - will have to make WS reference

summary(m3a)

#now try with Big Pine as reference

score$Site<-factor(score$Site,levels=c("Big Pine","Western Sambo","Dave's Ledge","Looe Key","Maryland Shoals","Eastern Dry Rocks","Marker 32","Bahia Honda","Eastern Sambo"))## picked site with best survival as reference (want to know how much worse others are)

m3b<-coxme(died~Site+(1|ArrayInSite)+(1|UniqueID), data=score) #still getting warning about infinite coefficient - will have to make WS reference

summary(m3b)
exp(m3b$coefficients)

###issue seems to be that effects are time-dependent - 
#high mortality at BP only happens in final time step...
#maybe try time-dependent estimate of site hazard?
head(surv)
surv$T24_TOD_mod_cat<-as.factor(surv$T24_TOD_mod) #create dummy categorical for time

surv$SiteTime<-with(surv,factor(Site:T24_TOD_mod_cat))
head(surv)

control<-coxph.control(iter.max=500)
m3t<-coxph(Surv(T24_TOD_mod,T24_Dead) ~ Site + SiteTime,data=surv)


# Visualize
#++++++++++++++++++++++++++++++++++++

#Using custom color palette
colorS <- c("plum4", "lightcoral", "lightsalmon", "#0099c6", "darkseagreen1", "rosybrown1", "mediumseagreen", "lightgoldenrod", "paleturquoise1")

colorG <- c("salmon","rosybrown", "paleturquoise4", "goldenrod","seagreen")

score$Site<-factor(score$Site,levels=c("Bahia Honda","Big Pine","Dave's Ledge","Eastern Sambo","Eastern Dry Rocks", "Looe Key","Marker 32","Maryland Shoals","Western Sambo"))## picked site with best survival as reference (want to know how much worse others are)
score$Genotype<-factor(score$Genotype,levels=c("AP14-5","AP13-X5","AP13-X7","AP13-XK","AP13-X9"))## picked genet with best survival as reference

scurveS=survfit(died~Site, data=score) #just for purposes of plotting
scurveS

ggsurvplot(scurveG, data = score,
           palette=colorG,fun = "pct", conf.int = FALSE,legend=c(.136,.361),legend.title="Genotype",
           legend.labs=c("AP14-5","AP13-X5","AP13-X7","AP13-XK","AP13-X9"),
           ggtheme = theme_bw(),ylim=c(0,100))

ggsurvplot(scurveS, data = score,
           palette=colorS,fun = "pct", conf.int = FALSE,legend=c(.2,.361),legend.title="Site",
           legend.labs=c("Bahia Honda","Big Pine","Dave's Ledge","Eastern Sambo","Eastern Dry Rocks", "Looe Key","Marker 32","Maryland Shoals","Western Sambo"),
           ggtheme = theme_bw(),ylim=c(0,100))

### Now, use correlogram to show correlations in genet performance across sites


scurveSG <- survfit(died~Site+Genotype, data=score)

res<-data.frame(summary(scurveSG)$table)

site<-str_split(rownames(res),",",simplify=TRUE)[,1]
site<-gsub('Site=','',site)
genet<-str_split(rownames(res),",",simplify=TRUE)[,2]
genet<-gsub('Genotype=','',genet)
df<-data.frame(res,site,genet,stringsAsFactors=TRUE)

df$n.surv<-df$n.start-df$events
df$PercSurv<-df$n.surv/df$n.start
rownames(df)<-seq_len(nrow(df))
df.sort<-df[order(df$genet),]
df.short<-df.sort[,c(10,11,13)]
df.short$site = factor(df.short$site, levels=c("Eastern Sambo","Western Sambo","Dave's Ledge","Looe Key","Maryland Shoals","Eastern Dry Rocks","Marker 32","Big Pine","Bahia Honda"))
# Geographic order "Bahia Honda","Big Pine","Looe Key","Dave's Ledge","Maryland Shoals","Eastern Sambo","Western Sambo","Marker 32", "Eastern Dry Rocks"
# Mean Survival order "Eastern Sambo","Western Sambo","Dave's Ledge","Looe Key","Maryland Shoals","Eastern Dry Rocks","Marker 32","Big Pine","Bahia Honda"
GxSMatrix<-spread(df.short,site,PercSurv)
GcorrBySite<-cor(GxSMatrix[,c(2:10)]) 

#note - no variance in either ES or Bahia Honda - therefore must restrict to sites 
#for which some variance in survival was observed

GcorrBySite2<-cor(GxSMatrix[,c(3:9)]) 

corrplot(GcorrBySite2,type="upper",tl.col="black",tl.srt=45,method='ellipse')

corrplot.mixed(GcorrBySite2,lower = "number",upper="ellipse",tl.col="black")


GcorrBySiteSp<-corr.test(GxSMatrix[,c(3:9)],method = "spearman")
GcorrBySiteSp[["p"]]

upperTriangle(GcorrBySiteSp[["p"]],diag=F,byrow = T)<-lowerTriangle(GcorrBySiteSp[["p"]],diag=F,byrow = F)

quartz()
corrplot(GcorrBySiteSp[["r"]],type="upper",tl.col="black",tl.srt=45,method='ellipse',p.mat = GcorrBySiteSp[["p"]],sig.level = c(0.1), insig = "label_sig",pch.col = 'grey',diag = T)

##### is there any correlation among site survival ranks in Wyatt's and Hollands datasets?

#ES, WS, DL, LK, MS, EDR, M32, BP, BH
holland=c(1,0.9333333,0.9333333,0.8666667,0.8,0.6666667,0.6,0.4,0) #mean percent survival by site
wyatt=c(0.9642857,0.8571429,0.8461538,0.76,0.7,0.76,0.9259259,0.8518519,0.6) #mean percent survival by site

plot(holland~wyatt)

hist(holland) #not normal, need nonparametric
hist(wyatt)

cor.test(holland, wyatt, method="kendall") #reasonable, but nonsig

### From Wyatt's R script
#                     n     events
#Site=E. Sambo        28      1     
#Site=Marker 32       27      2     
#Site=W. Sambo        28      4     
#Site=Big Pine        27      4    
#Site=Dave's Ledge    26      4     
#Site=Looe Key        25      6    
#Site=EDR             25      6    
#Site=Maryland Shoals 30      9     
#Site=Bahia Honda     30     12   

########End of survival analysis#####################################################

######## Moving on to Size & Growth ####################################################

summary(grow)
head(grow)

#First, calculate monthly growth rates over time-windows

grow$T12_SAgr<-(grow$T12_SA-grow$T0_SA)/12
grow$T12_Vgr<-(grow$T12_V-grow$T0_V)/12
grow$T12_ISgr<-(grow$T12_IS-grow$T0_IS)/12

grow$T15_SAgr<-(grow$T15_SA-grow$T12_SA)/3
grow$T15_Vgr<-(grow$T15_V-grow$T12_V)/3
grow$T15_ISgr<-(grow$T15_IS-grow$T12_IS)/3

grow$T18_SAgr<-(grow$T18_SA-grow$T15_SA)/3
grow$T18_Vgr<-(grow$T18_V-grow$T15_V)/3
grow$T18_ISgr<-(grow$T18_IS-grow$T15_IS)/3

grow$T24_SAgr<-(grow$T24_SA-grow$T18_SA)/6
grow$T24_Vgr<-(grow$T24_V-grow$T18_V)/6
grow$T24_ISgr<-(grow$T24_IS-grow$T18_IS)/6

#Second, generate shape-based traits

#calculating the shape data
grow<-grow %>% mutate(T12_SAtoV=T12_SA/T12_V,
                      T12_Sphericity=(4*pi*(((3*T12_V)/(4*pi))^(1/3))^2)/T12_SA,
                      T15_SAtoV=T15_SA/T15_V,
                      T15_Sphericity=(4*pi*(((3*T15_V)/(4*pi))^(1/3))^2)/T15_SA,
                      T18_SAtoV=T18_SA/T18_V,
                      T18_Sphericity=(4*pi*(((3*T18_V)/(4*pi))^(1/3))^2)/T18_SA,
                      T24_SAtoV=T24_SA/T24_V,
                      T24_Sphericity=(4*pi*(((3*T24_V)/(4*pi))^(1/3))^2)/T24_SA)

### Now some general data exploration

hist(grow$T0_SA) #very normal

hist(grow$T12_SA) #a little skewed, but not bad
hist(grow$T12_V) #more skewed
hist(grow$T12_IS) #even more skewed

hist(grow$T15_SA) #a little skewed, but not bad
hist(grow$T15_V) #more skewed
hist(grow$T15_IS) #even more skewed

hist(grow$T18_SA) #more skewed
hist(grow$T18_V) #very skewed
hist(grow$T18_IS) #very skewed

hist(grow$T24_SA) #not bad
hist(grow$T24_V) #very skewed
hist(grow$T24_IS) #very very skewed

#Quick look at overall averages for sites and genotype
grow %>% group_by(Genotype) %>% #seems like some diffs among genets in growth - X9 is outpacing
  summarise(mean(T24_SA,na.rm=T),
            mean(T24_V,na.rm=T),
            mean(T24_IS,na.rm=T))

grow %>% group_by(Site) %>%. #EDR is not great for survival, but growth is good
  summarise(mean(T24_SA,na.rm=T),
            mean(T24_V,na.rm=T),
            mean(T24_IS,na.rm=T))
  
#Does surface area vary among genets originally?
  
anova(lm(T0_SA~Genotype,data=grow)) #yes

plot(T0_SA~Genotype,data=grow) 
#yes, some genets larger initially (but note X9 is middle), 
#must therefore explore including as covariate
#no diffs in initial V or IS bc there is no 3D growth

####### Mixed models for absolute size 
####################

#1. Effects on absolute size in order to have largest data set
#   - basically considering all of things that happen to a coral, what impacts final size

#pivoting data set to have a column for trait value and time and trait, omitting Bahia Honda bc insufficient survival
grow$T0_SAdup<-grow$T0_SA #need to also keep T0_SA as covariate to test for now

size<-grow %>% 
  filter(Site!="Bahia Honda") %>%
  pivot_longer(cols=c(10,12,13,15,17,18,20,22,23,25,27,28,30,34,35,50:57),names_to = c('time','trait'),names_sep='_',values_to="size",values_drop_na=T)
size$time<-factor(size$time,ordered = T,levels=c("T0","T12","T15","T18","T24")) #making time ordinal

# This version does not have T12 because this time point is missing data 
# inferrence doesn't change however, so using full dataset version above
#size2<-grow %>% 
#  filter(Site!="Bahia Honda") %>%
#  pivot_longer(cols=c(10,34,35,17,18,20,22,23,25,27,28,30),names_to = c('time','trait'),names_sep='_',values_to="size",values_drop_na=T)
#size2$time<-factor(size2$time,ordered = T,levels=c("T0","T15","T18","T24")) #making time ordinal

### testing out transforms ###
dat=size[which(size$trait=="SA"),]

#dat=size2[which(size$trait=="SA"),]

hist(dat$size)
hist(log(dat$size)) #not quite far enough
hist(sqrt(dat$size)) #much improved
summary(dat$size)
head(dat$time)
#####

#Modeling SA

SA<-lmer(size~T0_SAdup+Genotype*Site*time+(1|Site:Array)+(1|UniqueID),data=size[which(size$trait=="SA"),],REML = T)

resSA <- simulateResiduals(fittedModel = SA) #check fits with DHARMa
plot(resSA) #not great
testDispersion(resSA) #significantly underdispersed 

#try transform
SAmod<-lmer((sqrt(size))~T0_SAdup+Genotype*Site*time+(1|Site:Array)+(1|UniqueID),data=size[which(size$trait=="SA"),],REML = T)

resSAmod <- simulateResiduals(fittedModel = SAmod) #check fits with DHARMa
plot(resSAmod) #still overfit
testDispersion(resSAmod) #still underdispersed 

#try simplification
SAmod2<-lmer(sqrt(size)~T0_SAdup+Genotype*Site+time+(1|Site:Array)+(1|UniqueID),data=size[which(size$trait=="SA"),],REML = T)

resSAmod2 <- simulateResiduals(fittedModel = SAmod2) #check fits with DHARMa
plot(resSAmod2) #still overfit, but better
testDispersion(resSAmod2) #no longer underdispersed 

#how much of an effect does initial SA have?
dat=size[which(size$trait=="SA"),]
par(mfrow=c(2,2))
sub=subset(dat,time=="T12")
plot(sqrt(size)~T0_SAdup,sub,xlab="T0 Surface Area",ylab="sqrt(Surface Area at T12)",main="T12")
abline(lm(sqrt(size)~T0_SAdup,sub)) 
summary(lm(sqrt(size)~T0_SAdup,sub)) #NS
sub=subset(dat,time=="T15")
plot(sqrt(size)~T0_SAdup,sub,xlab="T0 Surface Area",ylab="sqrt(Surface Area at T15)",main="T15")
abline(lm(sqrt(size)~T0_SAdup,sub)) 
summary(lm(sqrt(size)~T0_SAdup,sub)) #NS
sub=subset(dat,time=="T18")
plot(sqrt(size)~T0_SAdup,sub,xlab="T0 Surface Area",ylab="sqrt(Surface Area at T18)",main="T18")
abline(lm(sqrt(size)~T0_SAdup,sub)) 
summary(lm(sqrt(size)~T0_SAdup,sub)) #NS
sub=subset(dat,time=="T24")
plot(sqrt(size)~T0_SAdup,sub,xlab="T0 Surface Area",ylab="sqrt(Surface Area at T24)",main="T24")
abline(lm(sqrt(size)~T0_SAdup,sub)) 
summary(lm(sqrt(size)~T0_SAdup,sub)) #NS

#effect minimal, now without initial size term
SAmod3<-lmer((sqrt(size))~Genotype*Site+time+(1|Site:Array)+(1|UniqueID),data=size[which(size$trait=="SA"),],REML = T)

resSAmod3 <- simulateResiduals(fittedModel = SAmod3) #check fits with DHARMa
plot(resSAmod3) #still overfit, but no warnings, QQ is better
testDispersion(resSAmod3) #no longer underdispersed 

anova(SAmod3,ddf="Kenward-Roger") # additive effect of site, def time effect, no interaction
#Type III Analysis of Variance Table with Kenward-Roger's method
#              Sum Sq Mean Sq NumDF  DenDF  F value  Pr(>F)    
#Genotype        38.9    9.71     4  61.13   2.2529 0.07366 .  
#Site            91.4   13.06     7  15.94   3.0287 0.03161 *  
#time          4743.3 1185.82     4 343.51 274.9943 < 2e-16 ***
#Genotype:Site   74.1    2.65    28  59.96   0.6137 0.92090    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

summary(SAmod3) 
rand(SAmod3) # both significant
#ANOVA-like table for random-effects: Single term deletions
#Model:
#  (sqrt(size)) ~ Genotype + Site + time + (1 | Site:Array) + (1 | UniqueID) + Genotype:Site
#                 npar  logLik    AIC    LRT Df Pr(>Chisq)    
#<none>             47 -970.95 2035.9                         
#(1 | Site:Array)   46 -977.71 2047.4 13.523  1  0.0002356 ***
#  (1 | UniqueID)     46 -977.39 2046.8 12.888  1  0.0003307 ***
##  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

#need to do Tukey's post-hoc for site and genotype

SAs <- glht(SAmod3, linfct = mcp(Site = "Tukey"))
summary(SAs) #no sig pairwise difs post MTC

SAt <- glht(SAmod3, linfct = mcp(time = "Tukey"))
summary(SAt) #strong differences over time


hist(dat$T0_SAdup)

#plot of average genotype size over time at a site; need to re-incorporate T0 SA, col 10

colorG <- c("rosybrown", "paleturquoise4","seagreen", "goldenrod","salmon")
colorS <- c("#0099c6","mediumseagreen", "paleturquoise1", "lightcoral", "lightsalmon", "rosybrown1","darkseagreen1","lightgoldenrod")


#quartz()
size[which(size$trait=="SA"),] %>% #change trait as needed
  group_by(Genotype,Site,time) %>%
  summarise(size=mean(size)) %>%
  mutate(time=factor(time,levels = c("T0","T12","T15","T18","T24"))) %>%
  ggplot(aes(x=time,y=size,group=Genotype))+
  geom_point(aes(color=Genotype))+
  geom_line(aes(color=Genotype,group=Genotype))+
  scale_color_manual(values=c("rosybrown", "paleturquoise4","seagreen", "goldenrod","salmon"))+
  ylab(expression(paste("Surface Area (",cm^{2},")")))+
  xlab(label = NULL)+theme_bw()+
  facet_wrap(~Site,ncol=4,scales = "fixed") #free

#now average for sites with error bars
pd <- position_dodge(0.3)
size[which(size$trait=="SA"),] %>% #change trait as needed
  group_by(Site,time) %>%
  summarise(avg=mean(size),sd=sd(size),n=n(),se=sd/sqrt(n)) %>%
  mutate(time=factor(time,levels = c("T0","T12","T15","T18","T24"))) %>%
  ggplot(aes(x=time,y=avg,group=Site))+
  geom_point(aes(color=Site),position=pd)+
  geom_line(aes(color=Site,group=Site),position=pd)+
  geom_errorbar(aes(ymin=avg-se,ymax=avg+se,color=Site),lwd=0.4,width=0.3,position=pd)+
  scale_color_manual(values=c("#0099c6","mediumseagreen", "paleturquoise1", "lightcoral", "lightsalmon", "rosybrown1","darkseagreen1","lightgoldenrod"))+
  ylab(expression(paste("Surface Area (",cm^{2},")")))+
  xlab(label = NULL)+theme_bw()


# V 
### testing out transforms ###
dat=size[which(size$trait=="V"),]
summary(dat)

hist(dat$size) #very skewed
hist(log(dat$size+1)) #better
hist(sqrt(dat$size+1)) #still skew
hist((dat$size)^1/3)
summary(dat$size)
head(dat$time)
#####

#Modeling V

# T0_V is zero as the frags do not have volume cause they are flat
Vmod<-lmer((log(size+1))~Genotype*Site+time+(1|Site:Array)+(1|UniqueID),data=size[which(size$trait=="V"),],REML = T)

resVmod<- simulateResiduals(fittedModel = Vmod) #check fits with DHARMa
plot(resVmod) #overfit, but residuals vs predicted better than SA
testDispersion(resVmod) #barely not underdispersed 

qqPlot(residuals(Vmod),xlab="Theoretical Quantiles",ylab="Observed Quantiles",main = "Volume") #log(x+1) transform makes the data more normal - just lower tail off
plot(Vmod,which=1,main="Vmod") #want random scatter; no apparent trendline, looks better after log(x+1)


anova(Vmod,ddf="Kenward-Roger") 
#Type III Analysis of Variance Table with Kenward-Roger's method
#              Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
#Genotype       14.12   3.531     4  61.13  11.0558 8.231e-07 ***
#Site            7.91   1.130     7  15.92   3.5383   0.01735 *  
#time          578.84 144.711     4 343.73 453.1163 < 2.2e-16 ***
#Genotype:Site   5.01   0.179    28  59.97   0.5599   0.95275    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# significant factors: genotype, site, time, genotype by site, and site by time
summary(Vmod)
rand(Vmod) # both random effects significant

#pairwise tests
#need to do Tukey's post-hoc for site and genotype

Vg <- glht(Vmod, linfct = mcp(Genotype = "Tukey"))
summary(Vg) #no sig pairwise difs post MTC, closest is XK vs X9 and 14-5

Vs <- glht(Vmod, linfct = mcp(Site = "Tukey"))
summary(Vs) #no sig pairwise difs post MTC

Vt <- glht(Vmod, linfct = mcp(time = "Tukey"))
summary(Vt) #strong differences over time

## Plot
size[which(size$trait=="V"),] %>% 
  group_by(Genotype,Site,time) %>%
  summarise(size=mean(size)) %>%
  mutate(time=factor(time,levels = c("T0","T12","T15","T18","T24"))) %>%
  ggplot(aes(x=time,y=size,group=Genotype))+
  geom_point(aes(color=Genotype))+
  geom_line(aes(color=Genotype,group=Genotype))+
  scale_color_manual(values=c("rosybrown", "paleturquoise4","seagreen", "goldenrod","salmon"))+
  ylab(expression(paste("Volume (",cm^{3},")")))+
  xlab(label = NULL)+theme_bw()+
  facet_wrap(~Site,ncol=4,scales = "fixed") #free


#now average for genets with error bars
pd <- position_dodge(0.3)
size[which(size$trait=="V"),] %>% #change trait as needed
  group_by(Genotype,time) %>%
  summarise(avg=mean(size),sd=sd(size),n=n(),se=sd/sqrt(n)) %>%
  mutate(time=factor(time,levels = c("T0","T12","T15","T18","T24"))) %>%
  ggplot(aes(x=time,y=avg,group=Genotype))+
  geom_point(aes(color=Genotype),position=pd)+
  geom_line(aes(color=Genotype,group=Genotype),position=pd)+
  geom_errorbar(aes(ymin=avg-se,ymax=avg+se,color=Genotype),lwd=0.4,width=0.3,position=pd)+
  scale_color_manual(values=c("rosybrown", "paleturquoise4","seagreen", "goldenrod","salmon"))+
  ylab(expression(paste("Volume (",cm^{3},")")))+
  xlab(label = NULL)+theme_bw()


#now average for sites with error bars
pd <- position_dodge(0.3)
size[which(size$trait=="V"),] %>% #change trait as needed
  group_by(Site,time) %>%
  summarise(avg=mean(size),sd=sd(size),n=n(),se=sd/sqrt(n)) %>%
  mutate(time=factor(time,levels = c("T0","T12","T15","T18","T24"))) %>%
  ggplot(aes(x=time,y=avg,group=Site))+
  geom_point(aes(color=Site),position=pd)+
  geom_line(aes(color=Site,group=Site),position=pd)+
  geom_errorbar(aes(ymin=avg-se,ymax=avg+se,color=Site),lwd=0.4,width=0.3,position=pd)+
  scale_color_manual(values=c("#0099c6","mediumseagreen", "paleturquoise1", "lightcoral", "lightsalmon", "rosybrown1","darkseagreen1","lightgoldenrod"))+
  ylab(expression(paste("Volume (",cm^{3},")")))+
  xlab(label = NULL)+theme_bw()

### testing out transforms ###
dat=size[which(size$trait=="IS"),]
summary(dat)

hist(dat$size) #very skewed
hist(log(dat$size+1)) #better
hist(sqrt(dat$size+1)) #still skew
hist((dat$size)^1/3) #terrible
summary(dat$size)
head(dat$time)
#####

#Modeling IS
ISmod<-lmer(log(size+1)~Genotype*Site+time+(1|Site:Array)+(1|UniqueID),data=size[which(size$trait=="IS"),],REML = T)

resISmod<- simulateResiduals(fittedModel = ISmod) #check fits with DHARMa
plot(resISmod) #QQ plot looks good; residuals kind of wild
testDispersion(resISmod) #technically not underdispersed

# other viz
qqPlot(residuals(ISmod),xlab="Theoretical Quantiles",ylab="Observed Quantiles",main = "ISmod") #log(x+1) transform makes the data more normal

plot(ISmod,which=1,main="ISmod") #want random scatter; zeros skew at T0

anova(ISmod,ddf="Kenward-Roger")
#Type III Analysis of Variance Table with Kenward-Roger's method
#              Sum Sq Mean Sq NumDF  DenDF  F value    Pr(>F)    
#Genotype        9.88   2.469     4  62.06   6.8347 0.0001264 ***
#Site            6.62   0.946     7  15.92   2.6185 0.0527623 .  
#time          566.42 141.604     4 340.67 391.9976 < 2.2e-16 ***
#Genotype:Site   4.82   0.172    28  61.20   0.4765 0.9831983    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(ISmod) 
rand(ISmod) # Random effects significant

#pairwise tests
#need to do Tukey's post-hoc for site and genotype

ISg <- glht(ISmod, linfct = mcp(Genotype = "Tukey"))
summary(ISg) #no sig pairwise difs post MTC, closest is XK vs X9 and 14-5

ISt <- glht(ISmod, linfct = mcp(time = "Tukey"))
summary(ISt) #strong differences over time


## Plot
size[which(size$trait=="IS"),] %>% 
  group_by(Genotype,Site,time) %>%
  summarise(size=mean(size)) %>%
  mutate(time=factor(time,levels = c("T0","T12","T15","T18","T24"))) %>%
  ggplot(aes(x=time,y=size,group=Genotype))+
  geom_point(aes(color=Genotype))+
  geom_line(aes(color=Genotype,group=Genotype))+
  scale_color_manual(values=c("rosybrown", "paleturquoise4","seagreen", "goldenrod","salmon"))+
  ylab(expression(paste("Interstitial Space (",cm^{3},")")))+
  xlab(label = NULL)+theme_bw()+
  facet_wrap(~Site,ncol=4,scales = "fixed") #free


#now average for genets with error bars
pd <- position_dodge(0.3)
size[which(size$trait=="IS"),] %>% #change trait as needed
  group_by(Genotype,time) %>%
  summarise(avg=mean(size),sd=sd(size),n=n(),se=sd/sqrt(n)) %>%
  mutate(time=factor(time,levels = c("T0","T12","T15","T18","T24"))) %>%
  ggplot(aes(x=time,y=avg,group=Genotype))+
  geom_point(aes(color=Genotype),position=pd)+
  geom_line(aes(color=Genotype,group=Genotype),position=pd)+
  geom_errorbar(aes(ymin=avg-se,ymax=avg+se,color=Genotype),lwd=0.4,width=0.3,position=pd)+
  scale_color_manual(values=c("rosybrown", "paleturquoise4","seagreen", "goldenrod","salmon"))+
  ylab(expression(paste("Interstitial Space (",cm^{3},")")))+
  xlab(label = NULL)+theme_bw()

### testing out transforms ###
dat=size[which(size$trait=="SAtoV"),]
summary(dat)

hist(dat$size) #very skewed
hist(log(dat$size+1)) #better
hist(sqrt(dat$size+1)) #still skew
hist((dat$size)^1/3) #terrible
summary(dat$size)
head(dat$time)
#####

#Modeling SAtoV
SAtoVmod<-lmer(log(size+1)~Genotype*Site+time+(1|Site:Array)+(1|UniqueID),data=size[which(size$trait=="SAtoV"),],REML = T)

resSAtoVmod<- simulateResiduals(fittedModel = SAtoVmod) #check fits with DHARMa
plot(resSAtoVmod) #QQ plot looks better
testDispersion(resSAtoVmod) #underdispersed

# other viz
qqPlot(residuals(SAtoVmod),xlab="Theoretical Quantiles",ylab="Observed Quantiles",main = "SAtoVmod") #log(x+1) transform makes the data more normal; a few outliers

plot(SAtoVmod,which=1,main="SAtoVmod") #want random scatter; not too bad

anova(SAtoVmod,ddf="Kenward-Roger")
#Type III Analysis of Variance Table with Kenward-Roger's method
#              Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
#Genotype      9.3719 2.34298     4  54.253 35.9237 1.168e-14 ***
#Site          2.3671 0.33815     7  15.226  5.1830  0.003527 ** 
#time          0.0831 0.04157     2 192.807  0.6373  0.529817    
#Genotype:Site 2.3205 0.08288    28  52.627  1.2708  0.223145    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(SAtoVmod) 
rand(SAtoVmod) # only UniqueID significant

#pairwise tests
#need to do Tukey's post-hoc for site and genotype

SAtoVg <- glht(SAtoVmod, linfct = mcp(Genotype = "Tukey"))
summary(SAtoVg) #big diffs between genets - XK is outlier

SAtoVSs <- glht(SAtoVmod, linfct = mcp(Site = "Tukey"))
summary(SAtoVSs) #no sig pairwise difs post MTC


SAtoVt <- glht(SAtoVmod, linfct = mcp(time = "Tukey"))
summary(SAtoVt) #strong differences over time


## Plot
size[which(size$trait=="SAtoV"),] %>% 
  group_by(Genotype,Site,time) %>%
  summarise(size=mean(size)) %>%
  mutate(time=factor(time,levels = c("T0","T12","T15","T18","T24"))) %>%
  ggplot(aes(x=time,y=size,group=Genotype))+
  geom_point(aes(color=Genotype))+
  geom_line(aes(color=Genotype,group=Genotype))+
  scale_color_manual(values=c("rosybrown", "paleturquoise4","seagreen", "goldenrod","salmon"))+
  ylab("Surface Area to Volume Ratio")+
  xlab(label = NULL)+theme_bw()+
  facet_wrap(~Site,ncol=4,scales = "fixed") #free


#now average for genets with error bars
pd <- position_dodge(0.3)
size[which(size$trait=="SAtoV"),] %>% #change trait as needed
  group_by(Genotype,time) %>%
  summarise(avg=mean(size),sd=sd(size),n=n(),se=sd/sqrt(n)) %>%
  mutate(time=factor(time,levels = c("T0","T12","T15","T18","T24"))) %>%
  ggplot(aes(x=time,y=avg,group=Genotype))+
  geom_point(aes(color=Genotype),position=pd)+
  geom_line(aes(color=Genotype,group=Genotype),position=pd)+
  geom_errorbar(aes(ymin=avg-se,ymax=avg+se,color=Genotype),lwd=0.4,width=0.3,position=pd)+
  scale_color_manual(values=c("rosybrown", "paleturquoise4","seagreen", "goldenrod","salmon"))+
  ylab("Surface Area to Volume Ratio")+
  xlab(label = NULL)+theme_bw()

#now average for sites with error bars
pd <- position_dodge(0.3)
size[which(size$trait=="SAtoV"),] %>% #change trait as needed
  group_by(Site,time) %>%
  summarise(avg=mean(size),sd=sd(size),n=n(),se=sd/sqrt(n)) %>%
  mutate(time=factor(time,levels = c("T0","T12","T15","T18","T24"))) %>%
  ggplot(aes(x=time,y=avg,group=Site))+
  geom_point(aes(color=Site),position=pd)+
  geom_line(aes(color=Site,group=Site),position=pd)+
  geom_errorbar(aes(ymin=avg-se,ymax=avg+se,color=Site),lwd=0.4,width=0.3,position=pd)+
  scale_color_manual(values=c("#0099c6","mediumseagreen", "paleturquoise1", "lightcoral", "lightsalmon", "rosybrown1","darkseagreen1","lightgoldenrod"))+
  ylab("Surface Area to Volume Ratio")+
  xlab(label = NULL)+theme_bw()


### testing out transforms ###
dat=size[which(size$trait=="Sphericity"),]
summary(dat)

hist(dat$size) #pretty normal

head(dat$time)
#####

#Modeling Sphericity
SPmod<-lmer(size~Genotype*Site+time+(1|Site:Array)+(1|UniqueID),data=size[which(size$trait=="Sphericity"),],REML = T)

resSPmod<- simulateResiduals(fittedModel = SPmod) #check fits with DHARMa
plot(resSPmod) #not terrible
testDispersion(resSPmod) #just underdispersed

# other viz
qqPlot(residuals(SPmod),xlab="Theoretical Quantiles",ylab="Observed Quantiles",main = "SAtoVmod") #log(x+1) transform makes the data more normal; a few outliers

plot(SPmod,which=1,main="SAtoVmod") #get bigger over time...

anova(SPmod,ddf="Kenward-Roger")
#Type III Analysis of Variance Table with Kenward-Roger's method
#               Sum Sq  Mean Sq NumDF   DenDF F value    Pr(>F)    
#Genotype      0.69094 0.172735     4  53.079 11.0596 1.346e-06 ***
#Site          0.17970 0.025671     7  15.749  1.6436    0.1949    
#time          0.08314 0.027714     3 232.744  1.7744    0.1528    
#Genotype:Site 0.51060 0.018236    28  51.742  1.1676    0.3084    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(SPmod) 
rand(SPmod) # both significant

#pairwise tests
#need to do Tukey's post-hoc for site and genotype

SPg <- glht(SPmod, linfct = mcp(Genotype = "Tukey"))
summary(SPg) #big diffs between genets - XK is outlier


## Plot
size[which(size$trait=="Sphericity"),] %>% 
  group_by(Genotype,Site,time) %>%
  summarise(size=mean(size)) %>%
  mutate(time=factor(time,levels = c("T0","T12","T15","T18","T24"))) %>%
  ggplot(aes(x=time,y=size,group=Genotype))+
  geom_point(aes(color=Genotype))+
  geom_line(aes(color=Genotype,group=Genotype))+
  scale_color_manual(values=c("rosybrown", "paleturquoise4","seagreen", "goldenrod","salmon"))+
  ylab("Sphericity")+
  xlab(label = NULL)+theme_bw()+
  facet_wrap(~Site,ncol=4,scales = "fixed") #free




############## PCA of shape related traits by time, genet and site

#Wanted to do a PCA for morphology but using absolute size isn't meaningful because the dominant signal is size change over time so instead looking at shape in traits not dependent on size
#Evaluating invariant shape of colonies at end of the experiment at each site
#conclusion is: not very informative...no real separation by site or genet

shape2<-size %>% filter(time=="T24") 
shape2<-size[,c(1:6,36,37)]
shape2<-pivot_wider(shape2,names_from = trait,values_from = size)

#Morphology principal components 
pcaMor<-prcomp(shape2[,10:12],scale=T)
summary(pcaMor)

#if you want to color points by categories of good vs bad sites
shape2$Quality<-ifelse(shape2$Site %in% c("E. Sambo","Marker 32","W. Sambo"),"Top 3",ifelse(shape2$Site %in% c("EDR","Maryland Shoals","Bahia Honda"),"Bottom 3","Mid 3"))
shape2$Quality<-as.factor(shape2$Quality)

#change components of autoplot to color by site, quality, genet, etc
quartz()
autoplot(pcaMor,data=shape2,loadings=T,loadings.label=T,colour='Genotype',size=2,frame=T,frame.type=,
         loadings.label.size=3,loadings.label.repel=T, loading.label.color="darkgrey", loadings.colour='lightgrey',
         label.size=3)+
  #geom_text(vjust=1.5,label=envC[-6,1],size=3)+
  #geom_text_repel(label=envC[-6,1])+
  theme_classic()+
  theme(legend.position = "right")+
  scale_color_manual(values = colorG)+scale_fill_manual(values=colorG)

autoplot(pcaMor,data=shape2,loadings=T,loadings.label=T,colour='Site',size=2,frame=T,frame.type=,
         loadings.label.size=3,loadings.label.repel=T, loading.label.color="darkgrey", loadings.colour='lightgrey',
         label.size=3)+
  #geom_text(vjust=1.5,label=envC[-6,1],size=3)+
  #geom_text_repel(label=envC[-6,1])+
  theme_classic()+
  theme(legend.position = "right")+
  scale_color_manual(values = colorS)+scale_fill_manual(values=colorS)


#####
#End of mixed modeling for absolute size & shape

####Mixed models for growth rate

names(grow)

growthR<-grow %>% 
  filter(Site!="Bahia Honda") %>%
  pivot_longer(cols=c(37:48),names_to = c('time','trait'),names_sep='_',values_to="grow",values_drop_na=T)
growthR$time<-factor(growthR$time,ordered = T,levels=c("T12","T15","T18","T24")) #making time ordinal


#####
summary(growthR)
names(growthR)
head(growthR)




######Running the mixed effects models and apply transformations if needed

#SA
###looking at distribution
dat=growthR[which(growthR$trait=="SAgr"),]

hist(dat$grow) #pretty normal

hist(log(dat$grow+16)) #bad

hist(sqrt(dat$grow+16))

hist((dat$grow+16)^1/3)

#Modeling SA growth rate
SAGmod<-lmer(((grow))~Genotype*Site+time+(1|Site:Array)+(1|UniqueID),data=growthR[which(growthR$trait=="SAgr"),],REML = T)

resSAGmod<- simulateResiduals(fittedModel = SAGmod) #check fits with DHARMa
plot(resSAGmod) #issues
testDispersion(resSAGmod) #underdispersed

# other viz
qqPlot(residuals(SAGmod),xlab="Theoretical Quantiles",ylab="Observed Quantiles",main = "SAgr") 

plot(SAGmod,which=1,main="SAgr") #get bigger at tails

anova(SAGmod,ddf="Kenward-Roger")
#Type III Analysis of Variance Table with Kenward-Roger's method
#              Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)    
#Genotype       297.9   74.47     4  52.916  1.9333 0.118440    
#Site          1181.4  168.78     7  15.563  4.3815 0.007269 ** 
#time          1369.9  456.64     3 161.340 11.8538 4.67e-07 ***
#Genotype:Site 1006.1   35.93    28  49.589  0.9331 0.569074    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(SAGmod) 
rand(SAGmod) # both significant

#pairwise tests
#need to do Tukey's post-hoc for site and genotype

SAGs <- glht(SAGmod, linfct = mcp(Site = "Tukey"))
summary(SAGs) #no sig diff for site

SAGt <- glht(SAGmod, linfct = mcp(time = "Tukey"))
summary(SAGt) #gr increases over time


## Plot

pd <- position_dodge(0.3)
growthR[which(growthR$trait=="SAgr"),] %>% #change trait as needed
  group_by(Site,time) %>%
  summarise(avg=mean(grow),sd=sd(grow),n=n(),se=sd/sqrt(n)) %>%
  mutate(time=factor(time,levels = c("T0","T12","T15","T18","T24"))) %>%
  ggplot(aes(x=time,y=avg,group=Site))+
  geom_point(aes(color=Site),position=pd)+
  geom_line(aes(color=Site,group=Site),position=pd)+
  geom_errorbar(aes(ymin=avg-se,ymax=avg+se,color=Site),lwd=0.4,width=0.3,position=pd)+
  scale_color_manual(values=c("#0099c6","mediumseagreen", "paleturquoise1", "lightcoral", "lightsalmon", "rosybrown1","darkseagreen1","lightgoldenrod"))+
  ylab("Monthly Growth Rate of Surface Area")+
  xlab(label = NULL)+theme_bw()


#V
###looking at distribution
dat=growthR[which(growthR$trait=="Vgr"),]
summary(dat)

hist(dat$grow) #long tail

hist(log(dat$grow+5)) #bad

hist(sqrt(dat$grow+5)) #better

hist((dat$grow+5)^1/3) #not great

#Modeling SA growth rate
VGmod<-lmer((sqrt(grow+5))~Genotype*Site+time+(1|Site:Array)+(1|UniqueID),data=growthR[which(growthR$trait=="Vgr"),],REML = T)

resVGmod<- simulateResiduals(fittedModel = VGmod) #check fits with DHARMa
plot(resVGmod) #issues
testDispersion(resVGmod) #underdispersed

# other viz
qqPlot(residuals(VGmod),xlab="Theoretical Quantiles",ylab="Observed Quantiles",main = "Vgr") 

plot(VGmod,which=1,main="Vgr") #get bigger at tails

anova(VGmod,ddf="Kenward-Roger")
#Type III Analysis of Variance Table with Kenward-Roger's method
#              Sum Sq Mean Sq NumDF   DenDF F value   Pr(>F)    
#Genotype      2.1315 0.53287     4  53.391  5.3792 0.001028 ** 
#Site          3.4814 0.49734     7  15.416  5.0204 0.003974 ** 
#time          3.5786 1.19285     3 159.961 12.0410 3.78e-07 ***
#Genotype:Site 2.5891 0.09247    28  50.671  0.9336 0.568292    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(VGmod) 
rand(VGmod) # both significant

#pairwise tests
#need to do Tukey's post-hoc for site and genotype

VGs <- glht(VGmod, linfct = mcp(Site = "Tukey"))
summary(VGs) #no sig diff for site

VGg <- glht(VGmod, linfct = mcp(Genotype = "Tukey"))
summary(VGg) #XK is lower than 14-5

VGt <- glht(VGmod, linfct = mcp(time = "Tukey"))
summary(VGt) #gr increases over time


## Plot

pd <- position_dodge(0.3)
growthR[which(growthR$trait=="Vgr"),] %>% #change trait as needed
  group_by(Site,time) %>%
  summarise(avg=mean(grow),sd=sd(grow),n=n(),se=sd/sqrt(n)) %>%
  mutate(time=factor(time,levels = c("T0","T12","T15","T18","T24"))) %>%
  ggplot(aes(x=time,y=avg,group=Site))+
  geom_point(aes(color=Site),position=pd)+
  geom_line(aes(color=Site,group=Site),position=pd)+
  geom_errorbar(aes(ymin=avg-se,ymax=avg+se,color=Site),lwd=0.4,width=0.3,position=pd)+
  scale_color_manual(values=c("#0099c6","mediumseagreen", "paleturquoise1", "lightcoral", "lightsalmon", "rosybrown1","darkseagreen1","lightgoldenrod"))+
  ylab("Monthly Growth Rate of Volume")+
  xlab(label = NULL)+theme_bw()

pd <- position_dodge(0.3)
growthR[which(growthR$trait=="Vgr"),] %>% #change trait as needed
  group_by(Genotype,time) %>%
  summarise(avg=mean(grow),sd=sd(grow),n=n(),se=sd/sqrt(n)) %>%
  mutate(time=factor(time,levels = c("T0","T12","T15","T18","T24"))) %>%
  ggplot(aes(x=time,y=avg,group=Genotype))+
  geom_point(aes(color=Genotype),position=pd)+
  geom_line(aes(color=Genotype,group=Genotype),position=pd)+
  geom_errorbar(aes(ymin=avg-se,ymax=avg+se,color=Genotype),lwd=0.4,width=0.3,position=pd)+
  scale_color_manual(values=c("rosybrown", "paleturquoise4","seagreen", "goldenrod","salmon"))+
  ylab("Monthly Growth Rate of Volume")+
  xlab(label = NULL)+theme_bw()


#IS
###looking at distribution
dat=growthR[which(growthR$trait=="ISgr"),]
summary(dat)

hist(dat$grow) #long tail

hist(log(dat$grow+3)) #better

hist(sqrt(dat$grow+3)) #not great

hist((dat$grow+3)^1/3) #bad

#Modeling SA growth rate
ISGmod<-lmer((log(grow+3))~Genotype*Site+time+(1|Site:Array)+(1|UniqueID),data=growthR[which(growthR$trait=="ISgr"),],REML = T)

resISGmod<- simulateResiduals(fittedModel = ISGmod) #check fits with DHARMa
plot(resISGmod) #issues
testDispersion(resISGmod) #VERY underdispersed

# other viz
qqPlot(residuals(ISGmod),xlab="Theoretical Quantiles",ylab="Observed Quantiles",main = "ISgr") 

plot(ISGmod,which=1,main="ISgr") #get bigger at tails

anova(ISGmod,ddf="Kenward-Roger")
#Type III Analysis of Variance Table with Kenward-Roger's method
#              Sum Sq Mean Sq NumDF   DenDF F value    Pr(>F)    
#Genotype      2.5906 0.64764     4  53.831  4.3982  0.003776 ** 
#Site          4.6185 0.65979     7  15.283  4.4803  0.006838 ** 
#time          6.5312 2.17708     3 158.339 14.7842 1.557e-08 ***
#Genotype:Site 3.6752 0.13126    28  51.584  0.8915  0.621200    
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
summary(ISGmod) 
rand(ISGmod) # only UniqueID sig

#pairwise tests
#need to do Tukey's post-hoc for site and genotype

ISGs <- glht(ISGmod, linfct = mcp(Site = "Tukey"))
summary(ISGs) #no sig diff for site

ISGg <- glht(ISGmod, linfct = mcp(Genotype = "Tukey"))
summary(ISGg) #XK is almost outlier again

ISGt <- glht(ISGmod, linfct = mcp(time = "Tukey"))
summary(ISGt) #gr increases over time


## Plot

pd <- position_dodge(0.3)
growthR[which(growthR$trait=="ISgr"),] %>% #change trait as needed
  group_by(Site,time) %>%
  summarise(avg=mean(grow),sd=sd(grow),n=n(),se=sd/sqrt(n)) %>%
  mutate(time=factor(time,levels = c("T0","T12","T15","T18","T24"))) %>%
  ggplot(aes(x=time,y=avg,group=Site))+
  geom_point(aes(color=Site),position=pd)+
  geom_line(aes(color=Site,group=Site),position=pd)+
  geom_errorbar(aes(ymin=avg-se,ymax=avg+se,color=Site),lwd=0.4,width=0.3,position=pd)+
  scale_color_manual(values=c("#0099c6","mediumseagreen", "paleturquoise1", "lightcoral", "lightsalmon", "rosybrown1","darkseagreen1","lightgoldenrod"))+
  ylab("Monthly Growth Rate of Interstitial Space")+
  xlab(label = NULL)+theme_bw()

pd <- position_dodge(0.3)
growthR[which(growthR$trait=="ISgr"),] %>% #change trait as needed
  group_by(Genotype,time) %>%
  summarise(avg=mean(grow),sd=sd(grow),n=n(),se=sd/sqrt(n)) %>%
  mutate(time=factor(time,levels = c("T0","T12","T15","T18","T24"))) %>%
  ggplot(aes(x=time,y=avg,group=Genotype))+
  geom_point(aes(color=Genotype),position=pd)+
  geom_line(aes(color=Genotype,group=Genotype),position=pd)+
  geom_errorbar(aes(ymin=avg-se,ymax=avg+se,color=Genotype),lwd=0.4,width=0.3,position=pd)+
  scale_color_manual(values=c("rosybrown", "paleturquoise4","seagreen", "goldenrod","salmon"))+
  ylab("Monthly Growth Rate of Interstitial")+
  xlab(label = NULL)+theme_bw()