####This script will replicate all analysis for the "Effect of mating among relatives on seed set" section of the paper. 
####The end goal is the estimate the number of lethal equivalent alleles in this population using a regression which follows Morton et al., 1956 and Nietlisbach et al., 2018. 


library(tidyverse)
library(car)


cdata <- read.csv("crossdata2.csv")
cdata$SeedSet <- cdata$Filled/cdata$Total
cdata$Individual <- paste(cdata$Family, cdata$Plant)


#Count genotypes in each treatment type
cdata %>% group_by(CrossType) %>% count(Individual) %>% count(Individual) %>% summarize(n = sum(n))
#Count observations of each treatment type
cdata %>% group_by(CrossType) %>% count(Individual) %>% summarize(n = sum(n))

###Remove incompatible crosses or errors
cdata2 <- filter(cdata, SeedSet > 0.22)

#Test difference between outcross polycross and open pollination
cdata2 %>% filter(Year == 2016, xor(CrossType == "W", CrossType == "OP")) %>% group_by(CrossType) %>% summarise(mean(SeedSet, na.rm = T))
cdata2 %>% filter(Year == 2016, xor(CrossType == "W", CrossType == "OP")) %>% summarize(wilcox.test(SeedSet~CrossType)$p.value)
cdata2 %>% filter(Year == 2016, xor(CrossType == "W", CrossType == "OP")) %>% summarize(wilcox.test(SeedSet~CrossType)$statistic)


#Test outcrossed difference between FS and HS families
cdata2 %>% filter(xor(CrossType == "W", CrossType == "OP")) %>% group_by(IBC) %>% summarise(mean(SeedSet, na.rm = T))
cdata2 %>% filter( xor(CrossType == "W", CrossType == "OP")) %>% summarize(wilcox.test(SeedSet~IBC)$p.value)
cdata2 %>% filter( xor(CrossType == "W", CrossType == "OP")) %>% summarize(wilcox.test(SeedSet~IBC)$statistic)


#Test difference between within family single cross and within family polycross
cdata2 %>% filter(xor(CrossType == "B", CrossType == "Y")) %>% group_by(IBC, CrossType) %>% summarise(mean(SeedSet, na.rm = T))
cdata2 %>% filter( xor(CrossType == "B", CrossType == "Y")) %>% group_by(IBC) %>% summarize(wilcox.test(SeedSet~CrossType)$p.value)
cdata2 %>% filter( xor(CrossType == "B", CrossType == "Y")) %>% group_by(IBC) %>% summarize(wilcox.test(SeedSet~CrossType)$statistic)

cdata2 %>% filter( xor(CrossType == "B", CrossType == "Y")) %>% group_by(IBC) %>% summarize(fligner.test(SeedSet~CrossType)$p.value)
cdata2 %>% filter( xor(CrossType == "B", CrossType == "Y")) %>% group_by(IBC) %>% summarize(fligner.test(SeedSet~CrossType)$statistic)


#Combine cross types into just inbred or outbred
cdata2$CrossType2 <- 0
for(i in 1: length(row.names(cdata2))){
  if(cdata2[i, 3] == "W"){
    cdata2[i,11] <- "OUT"
  }else {
    cdata2[i,11] <- "IN"
  }
}

#Remove open pollinated data, and outbred data from FS families
cdata2 <- filter(cdata2, CrossType != "OP")
cdata3 <- filter(cdata2, IBC == 0 | CrossType2 != "OUT")


##Check for interaction between treatment and year
Anova(lm(SeedSet~Year+CrossType2+Year*CrossType2, data = cdata3))


#Calculation of lethal equivalents
means <- cdata3 %>% group_by(CCA) %>% summarise(mean = mean(SeedSet, na.rm = T))


n <- cdata3 %>% group_by(CCA) %>% summarise(n = sum(Total))

w <- ((n[2] * means[2])/ (1-means[2]))
w2 <- as.numeric(w[,1])
CCA <- as.data.frame(means[,1])
CCA <- as.numeric(CCA[,1])
JP <- summary(lm(log(mean)~CCA, data = means, weights = w2))$coefficients

a <- JP[1,1]
b <- JP[2,1]
repeat{
predict <- exp(a+b*CCA)
w3 <- ((n[2] * predict[2])/ (1-predict[2]))
w4 <- as.numeric(w3[,1])

JP2 <- summary(lm(log(mean)~CCA, data = means, weights = w4))$coefficients
a2 <- JP2[1,1]
b2 <- JP2[2,1]

dif <- b-b2
if(dif > 0.0001){
  a <- a2
  b <- b2
}else{ print(paste0("A = ", a2*-1, "    B = ", b2*-1))
  break()}}


##Calculate lethal equivalents, not accounting for self-incompatibility

cdatab <- cdata %>% filter(CrossType != "OP") %>% filter(IBC == 0 | CrossType != "W")

meansb <- cdatab %>% group_by(CCA) %>% summarise(mean = mean(SeedSet, na.rm = T))


nb <- cdatab %>% group_by(CCA) %>% summarise(n = sum(Total))

wb <- ((nb[2] * meansb[2])/ (1-meansb[2]))
w2b <- as.numeric(wb[,1])

Mb <- summary(lm(log(mean)~CCA, data = meansb, weights = w2b))$coefficients

ab <- Mb[1,1]
bb <- Mb[2,1]
repeat{
  predictb <- exp(ab+bb*CCA)
  w3b <- ((nb[2] * predictb[2])/ (1-predictb[2]))
  w4b <- as.numeric(w3b[,1])
  
  M2b <- summary(lm(log(mean)~CCA, data = meansb, weights = w4b))$coefficients
  a2b <- M2b[1,1]
  b2b <- M2b[2,1]
  
  difb <- bb-b2b
  if(difb > 0.0001){
    ab <- a2b
    bb <- b2b
  }else{print(paste0("A = ", a2b*-1, "    B = ", b2b*-1))
    break()}}



##Figure

cdata3 %>% ggplot(aes(y= SeedSet, x= CCA, group = CCA)) + geom_boxplot() + geom_abline(slope = -(1-exp(b2)), intercept = exp(a2), size = .5) + xlab("Coefficient of Coancestry") + ylab("Seed Set") +
  stat_summary(fun=mean, geom="point", shape=23, size=2.2, color="#a3a3a3", fill="#a3a3a3") + theme_clean(base_family = "sans") + 
  theme(axis.title.x =  element_blank(),
        axis.title.y =  element_blank(),
        axis.text.x = element_text(size= 15),
        axis.text.y = element_text(size= 15))


