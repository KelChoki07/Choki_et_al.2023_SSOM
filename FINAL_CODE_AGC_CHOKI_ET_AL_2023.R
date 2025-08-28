#  _____  _             _     _          _            _      ___    ___  ___   ____  
# / ____|| |           | |   (_)        | |          | |    |__ \  / _ \|__ \ |___ \ 
#| |     | |__    ___  | | __ _     ___ | |_    __ _ | |       ) || | | |  ) |  __) |
#| |     | '_ \  / _ \ | |/ /| |   / _ \| __|  / _` || |      / / | | | | / /  |__ < 
#| |____ | | | || (_) ||   < | |  |  __/| |_  | (_| || | _   / /_ | |_| |/ /_  ___) |
# \_____||_| |_| \___/ |_|\_\|_|   \___| \__|  \__,_||_|(_) |____| \___/|____||____/ 
                                                                                                                                                                      


#------------------------------------------------------------------------------------------------------------#
# Publication details                                                                                        #
# Title: Conservation potential of non-protected area for sympatric carnivores in Bhutan                     #
# Authors: Karma Choki, Phub Dhendup, Jigme Tenzin, Dago Dorji, Kuenley Tenzin, Tenzin Wangmo, Ugyen Penjor  #
# Journal: Global Ecology and Conservation, Volume 42, 2023                                                  #
# DOI: https://doi.org/10.1016/j.gecco.2023.e02392                                                           #
#------------------------------------------------------------------------------------------------------------#

# Asian golden cat
# Revision: 25 December 2022
# Multi-scale covariates
# Scale optimisation via univariate modeling
# Single-season occupancy model

# Load packages
install.packages("corrplot")
library(unmarked)
library(AICcmodavg)
library(corrplot)
library(MuMIn)
library(ncf)
library(ggplot2)
library(raster)

# Load data 
agc_dh <- read.csv("DetHis_asiaticgoldencat.csv", header=T, row.names=1 , na.strings="NA")
head(agc_dh)

# Load covariates
SiteCovMulti <- read.csv('cov_multiscale_250822_labels.csv')

str(SiteCovMulti) #munt_P- poisson, munt_ZIP- zero inflated poisson 

# Import detection covariate
eFF <- read.csv("Effort2.csv", header=T, row.names=1)
head(eFF)

# Get naive occupancy
(naive_occ <- sum(ifelse(rowSums(agc_dh, na.rm=T)>0,1,0))/nrow(agc_dh))

# Prepare unmarked frame object 
Asiaticgoldencat <- unmarkedFrameOccu(y=agc_dh, siteCovs=SiteCovMulti, obsCovs=list(Effort=eFF))
str(Asiaticgoldencat)
?unmarkedFrameOccu
### ----------------------------------------------------------------------------

# Univariate model selection (scale-optimisation)
e5 <- occu(~Effort ~scale(ele500m), data=Asiaticgoldencat)
e1 <- occu(~Effort ~scale(ele1km), data=Asiaticgoldencat)
e2 <- occu(~Effort ~scale(ele2km), data=Asiaticgoldencat)
e4 <- occu(~Effort ~scale(ele4km), data=Asiaticgoldencat)

e.mods <- list("e5"=e5, "e1"=e1, "e2"=e2, "e4"=e4)
aictab(cand.set = e.mods, second.ord = F) # e5 #second.ord=F bc we are using AIC
#Your sample size here is the site within each elevation scale
#So, it'll be less than total stations bc not all the stations are in the same elev

s5 <- occu(~Effort ~scale(set500m), data=Asiaticgoldencat)
s1 <- occu(~Effort ~scale(set1km), data=Asiaticgoldencat)
s2 <- occu(~Effort ~scale(set2km), data=Asiaticgoldencat)
s4 <- occu(~Effort ~scale(set4km), data=Asiaticgoldencat)

s.mods <- list("s5"=s5, "s1"=s1, "s2"=s2, "s4"=s4)
aictab(cand.set = s.mods, second.ord = F) # s1
#how many pixel had housing point, and average within the given buffer
#focal mean-density- you count all the pixel that has feature and divide it by the total area
#within the buffer
#its density, not distance to the settle

r5 <- occu(~Effort ~scale(riv500m), data=Asiaticgoldencat)
r1 <- occu(~Effort ~scale(riv1km), data=Asiaticgoldencat)
r2 <- occu(~Effort ~scale(riv2km), data=Asiaticgoldencat)
r4 <- occu(~Effort ~scale(riv4km), data=Asiaticgoldencat)

r.mods <- list("r5"=r5, "r1"=r1, "r2"=r2, "r4"=r4)
aictab(cand.set = r.mods, second.ord = F) # r1
#its density, not distance

ro5 <- occu(~Effort ~scale(road500m), data=Asiaticgoldencat)
ro1 <- occu(~Effort ~scale(road1km), data=Asiaticgoldencat)
ro2 <- occu(~Effort ~scale(road2km), data=Asiaticgoldencat)
ro4 <- occu(~Effort ~scale(road4km), data=Asiaticgoldencat)

ro.mods <- list("ro5"=ro5, "ro1"=ro1, "ro2"=ro2, "ro4"=ro4)
aictab(cand.set = ro.mods, second.ord = F) # ro1

t5 <- occu(~Effort ~scale(tree500m), data=Asiaticgoldencat)
t1 <- occu(~Effort ~scale(tree1km), data=Asiaticgoldencat)
t2 <- occu(~Effort ~scale(tree2km), data=Asiaticgoldencat)
t4 <- occu(~Effort ~scale(tree4km), data=Asiaticgoldencat)

t.mods <- list("t5"=t5, "t1"=t1, "t2"=t2, "t4"=t4)
aictab(cand.set = t.mods, second.ord = F) # t2
#same concept like ele within each buffer scale
#high buffer area not necessary mean you will have high forest coverage due to 
#other built up area
#this is unclassified tree/forest cover

tn5 <- occu(~Effort ~scale(tree1500m), data=Asiaticgoldencat)
tn1 <- occu(~Effort ~scale(tree11km), data=Asiaticgoldencat)
tn2 <- occu(~Effort ~scale(tree12km), data=Asiaticgoldencat)
tn4 <- occu(~Effort ~scale(tree14km), data=Asiaticgoldencat)

tn.mods <- list("tn5"=tn5, "tn1"=tn1, "tn2"=tn2, "tn4"=tn4)
aictab(cand.set = tn.mods, second.ord = F) # tn2
#' tn= tree non forest
#' Hansen tree cover layer- use RECLASSIFY tool to classify tree cover into 3 classes-its continues value)
#' 0-20 % (if the pixel has any value b/n 0-20, you count those pixrls within given buffer) consisting of only 20% tree cover within

to5 <- occu(~Effort ~scale(tree2500m), data=Asiaticgoldencat)
to1 <- occu(~Effort ~scale(tree21km), data=Asiaticgoldencat)
to2 <- occu(~Effort ~scale(tree22km), data=Asiaticgoldencat)
to4 <- occu(~Effort ~scale(tree24km), data=Asiaticgoldencat)

to.mods <- list("to5"=to5, "to1"=to1, "to2"=to2, "to4"=to4)
aictab(cand.set = to.mods, second.ord = F) # to5
# class 2 (21-40%)

tc5 <- occu(~Effort ~scale(tree3500m), data=Asiaticgoldencat)
tc1 <- occu(~Effort ~scale(tree31km), data=Asiaticgoldencat)
tc2 <- occu(~Effort ~scale(tree32km), data=Asiaticgoldencat)
tc4 <- occu(~Effort ~scale(tree34km), data=Asiaticgoldencat)

tc.mods <- list("tc5"=tc5, "tc1"=tc1, "tc2"=tc2, "tc4"=tc4)
aictab(cand.set = tc.mods, second.ord = F) # tc5
#Class 3 (above 40%)

# Covariates to retain: ele500m, set1km, riv1km, road1km, tree2km, tree12km, tree2500m, tree3500m
#scale optimization- we are finding the optimum scale for each variable!

# Select prey covariate
pmun <- occu(~Effort ~scale(munt_P), data=Asiaticgoldencat) #your sample size here 
# is not the overall sites but the sites that have AGC detection 
#for muntjac and sambar are solitary animal, its more more appro to use Royle-Nichols (RN model)
psam <- occu(~Effort ~scale(sam_P), data=Asiaticgoldencat)
pwp <- occu(~Effort ~scale(wp_P), data=Asiaticgoldencat) #We used N-mix model bc they are group living group
prey <- list("pmun"=pmun, "psam"=psam, "pwp"=pwp) 
aictab(cand.set=prey, second.ord=F) #

### ----------------------------------------------------------------------------
head(SiteCovMulti)
# Check for correlation among the scale optimised covariates
# Final covariates for multivariate model: ele500m, set1km, riv1km, road1km, tree2km, tree2500m, munt_P

### ----------------------------------------------------------------------------

# Detection probability (we are just taking the constant det prob)
(d1 <- occu(~1 ~(ele500m)+scale(set1km)+scale(riv1km)+scale(tree2500m), data=Asiaticgoldencat))
(d2 <- occu(~scale(Effort) ~scale(ele500m)+scale(set1km)+scale(riv1km)+scale(tree2500m), data=Asiaticgoldencat, starts=c(1,0,0,0,0,0,0)))
det_mod <- list("d1"=d1, "d2"=d2)
aictab(cand.set = det_mod, second.ord = F) # to5

### ----------------------------------------------------------------------------

# Multivariate modelling 
# (mod1 <- occu(~scale(Effort) ~scale(ele500m)+scale(tree2500m)+scale(munRN_top), data=Asiaticgoldencat, starts=c(0, -1, 1, -3, -2, 2)))#, starts=c(0, 0, 0, 0, 0, 0)))
# (mod4 <- occu(~scale(Effort) ~scale(ele500m)+scale(set1km)+scale(tree2500m), data=Asiaticgoldencat))#, starts=c(0, 0, 0, 0, 0, 0)))
# (mod5 <- occu(~scale(Effort) ~scale(tree2500m)+scale(munRN_top), data=Asiaticgoldencat))#, starts=c(0, 1, -1, -1, 1)))
# (mod6 <- occu(~scale(Effort) ~scale(set1km)+scale(tree2500m), starts=c(0, 0, 0, -1, 1), data=Asiaticgoldencat))#, starts=c(0, 0, 0, 0, 0, 0)))
# (mod7 <- occu(~scale(Effort) ~scale(ele500m)+scale(tree2500m), starts=c(0, 0, 0, -2, 1), data=Asiaticgoldencat))#, starts=c(0, 0, 0, 0, 0, 0, 0)))
# (mod8 <- occu(~scale(Effort) ~scale(ele500m)+I(scale(ele500m)^2)+scale(tree2500m), data=Asiaticgoldencat))#, starts=c(0, 0, 0, 0, 0, 0)))
# (mod13 <- occu(~scale(Effort) ~scale(ele500m)+I(scale(ele500m)^2)+scale(munRN_top), data=Asiaticgoldencat))
# (nullmod <- occu(~1 ~1, data=Asiaticgoldencat))
# 
# mod_list <- list("mod4"=mod4, "mod5"=mod5, "mod6"=mod6, "mod7"=mod7, 
#                  "mod8"=mod8, "mod13"=mod13, 
#                  "null"=nullmod)
# aictab(cand.set=mod_list, second.ord=F) 

####################################################################################################

# Constant detection probability

(mod1 <- occu(~1 ~scale(ele500m)+scale(tree2500m)+scale(munRN_top), data=Asiaticgoldencat, starts=c(0, -1, 1, -3, -2)))
(mod4 <- occu(~1 ~scale(ele500m)+scale(set1km)+scale(tree2500m), data=Asiaticgoldencat))#, starts=c(0, 0, 0, 0, 0, 0)))
(mod5 <- occu(~1 ~scale(tree2500m)+scale(munRN_top), data=Asiaticgoldencat))#, starts=c(0, 1, -1, -1, 1)))
(mod6 <- occu(~1 ~scale(set1km)+scale(tree2500m), starts=c(0, 0, 0, -1), data=Asiaticgoldencat))
(mod7 <- occu(~1 ~scale(ele500m)+scale(tree2500m), starts=c(0, 0, 0, -2), data=Asiaticgoldencat))
(mod8 <- occu(~1 ~scale(ele500m)+I(scale(ele500m)^2)+scale(tree2500m), data=Asiaticgoldencat))#, starts=c(0, 0, 0, 0, 0, 0)))
(mod13 <- occu(~1 ~scale(ele500m)+I(scale(ele500m)^2)+scale(munRN_top), data=Asiaticgoldencat))
(nullmod <- occu(~1 ~1, data=Asiaticgoldencat))

mod_list <- list("mod4"=mod4, "mod5"=mod5, "mod6"=mod6, "mod7"=mod7, 
                 "mod8"=mod8, "mod13"=mod13, 
                 "null"=nullmod)
aictab(cand.set=mod_list, second.ord=F) 
#starting values are in seq of psi intercept, cov, det intercept & det cov.
# AGC has non-linear r/n w/ ele- ecological hypotheses- 
# 
####################################################################################################

# Goodness of fit #This is mackenzie-bailey
(occ_gof <- mb.gof.test(mod1, nsim=1000, plot.hist=T)) # c-hat=0.06; p=0.84
(occ_gof_top <- mb.gof.test(mod8, nsim=100, plot.hist=T)) # c-hat=0.07; p=0.83
# here we are checking the GOF of top model
# using our top, we are generating new dataset for 1000 times, so that we can check
#if the expected value are similar to that of our observed data

### ----------------------------------------------------------------------------

# Beta estimates and confidence intervals of the top model
# Change 'level=0.95' for 95% confidence interval
coef(mod8)
confint(mod8, type = "state", level=0.85) #bc we have highly uncertain estimates, we are bringing down the level of CI
confint(mod8, type = "det", level=0.85)

### ----------------------------------------------------------------------------

# Occupancy estimates

# Proportion of area occupied (PAO)
best_mod <- occu(~1 ~scale(ele500m)+I(scale(ele500m)^2)+scale(tree2500m), data=Asiaticgoldencat)
re <- ranef(best_mod)
EBUP <- bup(re, stat="mean")
CI <- confint(re, level=0.95)
rbind(PAO=c(Estimate=sum(EBUP), colSums(CI)) / 31)
#PAO- your PAO estimate cannot be less than your naive estimate (see the lower CI, that should be equal naive estimate)
# the proportion can be less than the sites where you have detected your species
#Whereas occupancy prob- the lower CI can be lower than naive estimate

# Occupancy probability
occ_pred <- predict(best_mod, 
                    newdata=data.frame(cbind("ele500m"=SiteCovMulti$ele500m, 
                                             "tree2500m"=SiteCovMulti$tree2500m)), 
                    type="state", 
                    inf.rm=T) #Here we're using the non-standardize variables (real variable)

cbind(Predicted=mean(occ_pred$Predicted), SE=mean(occ_pred$SE), 
      lower=mean(occ_pred$lower), upper=mean(occ_pred$upper))

# Detection probability
p_pred <- predict(best_mod, 
                  newdata=NULL, 
                  type="det", 
                  inf.rm=T)

cbind(Predicted=mean(p_pred$Predicted), SE=mean(p_pred$SE), 
      lower=mean(p_pred$lower), upper=mean(p_pred$upper))

### ----------------------------------------------------------------------------

# Model averaging #this is used for interpretation and prediction mapping
# Models within delta 6
occu_model_list <- list(occ_1=mod8,
                        occ_2=mod5,
                        occ_3=mod7,
                        occ_4=nullmod,
                        occ_5=mod6,
                        occ_6=mod4)

occ_avg <- model.avg(occu_model_list, fit=T)
coef(occ_avg)
confint(occ_avg, level=0.85)
#golden relationship with elevation is non-linear, 
#The occupancy of AGC peaks at 1500m, which is mid-altitude, whereas the muntjac prefers low elevation
#so mid elevation there is no muntjac, so AGC might have peferred other smaller species, previous studies
#has shown different prey


# Prediction

# Forest cover
# mod5 - you can try different models here (the overall effect does not change)
newDataF <- data.frame(
  tree2500m=seq(range(SiteCovMulti$tree2500m)[1], range(SiteCovMulti$tree2500m)[2], 100),
  munRN_top=mean(SiteCovMulti$munRN_top))

occ.probF <- predict(mod5, type="state", newdata=newDataF, appendData=T, level=0.85)



library(beepr)
library(sys)
beep("fanfare")
sys.sleep(2)

### ----------------------------------------------------------------------------

