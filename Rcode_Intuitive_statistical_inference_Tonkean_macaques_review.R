########################################################################################################
############## R script for article "Intuitive Statistical inferences in Tonkean macaques" #############
########################################################################################################
######################################################################################################## 
########################################################################################################

##### Packages needed #####

library(tidyverse) # version 2.0.0
library(ggplot2) # version 3.4.4      
library(lme4) # version 1.1-35.1      
library(dplyr) # version 1.1.4
library(tidyr) # version 1.3.0
library(rstatix) # version 0.7.2
library(gmodels) # version 2.18.1.1
library(lmerTest) # version 3.1-3     
library(readr) # version 2.1.4
library(cowplot) # version 1.1.1
library(car) # version 3.1-2          
library(MASS) # version 7.3-60
library(DHARMa) # version 0.4.6
library(corrplot) # version 0.92
library(performance) # version 0.12.0
library(see) # version 0.8.4
library(patchwork) # version 1.2.0
library(Rmisc) # version 1.5.1
library(glmulti) # version 1.0.8
library(rJava) # version 1.0-11 
library(ggpubr) # version 0.6.0 
library(Hmisc) # version 5.1-3

# Write your directory here

setwd("Write your directory here")
getwd()

##### Upload data #####
data <- read.csv2(file = "RawData_Test - Reviewed.csv", header = TRUE, sep = ';', dec = ',', stringsAsFactors = FALSE)

data$Park <- as.factor(data$Park)
data$ID_individual <- as.factor(data$ID_individual)
data$Age <- as.factor(data$Age)
data$Sex <- as.factor(data$Sex)
data$Phase <- as.factor(data$Phase)
data$Test_condition <- as.factor(data$Test_condition)
data$Session <- as.factor(data$Session)
data$Trial <- as.factor(data$Trial)
data$Experimenter_ID <- as.factor(data$Experimenter_ID)
data$Experimenter_position <- as.factor(data$Experimenter_position)
data$Arms_position <- as.factor(data$Arms_position)
data$Fav_Jar_position <- as.factor(data$Fav_Jar_position)
data$Side_Success <- as.factor(data$Side_Success)
data$Selected_side <- as.factor(data$Selected_side)
data$Selected_reward <- as.factor(data$Selected_reward)
data$Selected_side_binary <- as.numeric(data$Selected_side_binary)
data$Success <- as.numeric(data$Success)
data$ROR <- as.numeric(data$ROR)
data$RQP <- as.numeric(data$RQP)
data$RQB <- as.numeric(data$RQB)
data$Raisin <- as.factor(data$Raisin)

data$Test_condition <- fct_relevel(data$Test_condition, c("cond1","cond3","cond4","cond6","cond2a","cond5a","cond2b","cond5b"))
levels(data$Test_condition)
str(data)
summary(data)

###################
##### Table 2 #####
###################
#Weber's law testing models analyses

## Creating new colum : WL.Diff.Peanut = Higher QP - Lower QP // with QP = peanuts quantity in one jar
## Creating new colum : WL.Ratio.Peanut = Lower QP / (QP jar A + QP jar B) // with QP = peanuts quantity in one jar

data$WL.Diff.Peanuts <- 0
data$WL.Ratio.Peanuts <- 0

for (i in 1:nrow(data)){
  if (data$Test_condition[i]=="cond1" | data$Test_condition[i]=="cond4"){
    data$WL.Diff.Peanuts[i] <- data$Total_peanuts_potA[i] - data$Total_peanuts_potB[i]
    data$WL.Ratio.Peanuts[i] <- data$Total_peanuts_potB[i] / (data$Total_peanuts_potA[i]+data$Total_peanuts_potB[i])
  }
  else if (data$Test_condition[i]!="cond1" & data$Test_condition[i]!="cond4"){
    data$WL.Diff.Peanuts[i] <- data$Total_peanuts_potB[i] - data$Total_peanuts_potA[i]
    data$WL.Ratio.Peanuts[i] <- data$Total_peanuts_potA[i] / (data$Total_peanuts_potB[i]+data$Total_peanuts_potA[i])
  }
}

## Creating new colum : MorePeanut.Selected = Choice for the highest QP (coded 1)

data$MorePeanut.Selected <- 0

for (i in 1:nrow(data)){
  if (data$Test_condition[i]=="cond1" | data$Test_condition[i]=="cond4"){
    if (data$Success[i]==1){
      data$MorePeanut.Selected [i] <- 1
    }
    else if (data$Success[i]==0){
      data$MorePeanut.Selected [i] <- 0
    }
  }
  else if (data$Test_condition[i]!="cond1" & data$Test_condition[i]!="cond4"){
    if (data$Success[i]==1){
      data$MorePeanut.Selected [i] <- 0
    }
    else if (data$Success[i]==0){
      data$MorePeanut.Selected [i] <- 1
    }
  }
}

data$MorePeanut.Selected <- as.factor(data$MorePeanut.Selected)

#check data
str(data)
unique(data$WL.Diff.Peanuts)
unique(data$WL.Ratio.Peanuts)

## Selection of data without cond6 (Ctrl-Opaque)
data.WL <- data %>% 
  filter(Test_condition!="cond6")
#check data
unique(data.WL$Test_condition)

#### Results model (1)

Null.WeberLaw <- glm(MorePeanut.Selected ~ 1, family = binomial, data = data.WL)
Model1 <- glm(MorePeanut.Selected  ~ WL.Diff.Peanuts, family = binomial, data = data.WL)
Anova(Model1)
summary(Model1)

#### Results model (2)

Null.WeberLaw <- glm(MorePeanut.Selected ~ 1, family = binomial, data = data.WL)
Model2 <- glm(MorePeanut.Selected  ~ WL.Ratio.Peanuts, family = binomial, data = data.WL)
Anova(Model2)
summary(Model2)

####################
##### TABLE 3 #####
###################
#Individual performances across conditions --> binomial analyses

## Caculation of performance / individual / conditions
NbSuccess.Cond.Ind <- data %>%group_by(ID_individual,Test_condition) %>%
                              dplyr::summarise(NbSucces=sum(Success),N=n()) %>%
                              dplyr::mutate(Perf= NbSucces/N)
                    
## Creating table with p-value of binomial test for each individual / condition

NbSuccess.Cond.Ind$Pvalue_2sided <- "initialisation"
NbSuccess.Cond.Ind$Pvalue_greater <- "initialisation"
NbSuccess.Cond.Ind

for (i in 1:nrow(NbSuccess.Cond.Ind)){
  n <- NbSuccess.Cond.Ind$NbSucces[i]
  Tot <- NbSuccess.Cond.Ind$N[i]
  Test_Binom2sided <- binom.test(n,Tot,p=0.5,alternative="two.sided")
  Test_BinomGreater <- binom.test(n,Tot,p=0.5,alternative="greater")
  NbSuccess.Cond.Ind$Pvalue_2sided[i] <- Test_Binom2sided[["p.value"]]
  NbSuccess.Cond.Ind$Pvalue_greater[i] <- Test_BinomGreater[["p.value"]]
}
#check data
NbSuccess.Cond.Ind <- as.data.frame(NbSuccess.Cond.Ind)
str(NbSuccess.Cond.Ind)

### Results Table 3
write.table(NbSuccess.Cond.Ind,"Table 3.csv",sep=";",row.names=FALSE)

########################
##### FIGURE 3 & 4 #####
########################

## Adding column Status  complete/incomplete at table with performance / condition / individual to be added in graphs
NbSuccess.Cond.Ind$Status <- "status"

for (i in 1:nrow(NbSuccess.Cond.Ind)){
  if (NbSuccess.Cond.Ind$N[i]>=20){
    NbSuccess.Cond.Ind$Status[i] <- "complete"
  }
  else(NbSuccess.Cond.Ind$Status[i]<- "incomplete")
}

## Creating a table containing n= total number of trials per conditions, to be written in graphs
TableTest <- NbSuccess.Cond.Ind %>%
  group_by(Test_condition)%>%
  dplyr::summarise(NbTrial = sum(N))%>%
  dplyr::mutate(NumberTrial = paste("n=", NbTrial))

## Graphic layout
Th <-theme(axis.title=element_text(12),axis.text.x=element_text(size=12),axis.text.y=element_text(size=15),legend.title=element_text(size=15),legend.text=element_text(size=15),plot.title = element_text(hjust=0.5,size=20))

########################
### Graph Figure 3 #####
########################

##Selecting only conditions that showed statistical inferences
NbSuccess.Cond.Test <- NbSuccess.Cond.Ind %>% filter(Test_condition=="cond1"|Test_condition=="cond3"|Test_condition=="cond4"|Test_condition=="cond6")
Table.Cond.Test <- TableTest %>% filter(Test_condition=="cond1"|Test_condition=="cond3"|Test_condition=="cond4"|Test_condition=="cond6")

x11()
Figure3 <- ggplot(NbSuccess.Cond.Test)+ geom_boxplot(aes(x=Test_condition, y=Perf),width = 0.5)+
  geom_jitter(aes(x=Test_condition, y=Perf, color=ID_individual, shape=Status),size=2.5,width=0.3,height=0)+
  geom_hline(yintercept=0.50, linetype="dashed", color="black", linewidth=0.7)+
  scale_shape_manual(values=c(19,1))+
  scale_color_manual(values = c("purple","blue","chocolate4","deepskyblue4","azure4","darkgreen","green","red","#ffcc00ff","aquamarine","#ff6600ff","chartreuse3","magenta"))+
  geom_label(data=Table.Cond.Test, aes(x=Test_condition, y=1.2,label=NumberTrial))+
  scale_y_continuous(name="Jar A choice rate", breaks=seq(from=0,to=1,by=0.25), limits=c(0,1.3))+
  scale_x_discrete(labels=c("1-Reversed proportions","3-Ctrl-preferred items quantity","4-Ctrl-non preferred items quantity","6-Ctrl-Opaque"))+
  ggtitle("Figure 3")+
  xlab("Conditions")+Th
print(Figure3)

# Saving graph
svg("Figure3.svg", width = 9, height = 9.5) 
print(Figure3)
dev.off()

########################
### Graph Figure 4 #####
########################

##Selecting only conditions that showed inferences based on preferred items quantities
NbSuccess.Cond.Ctrl <- NbSuccess.Cond.Ind %>% filter(Test_condition!="cond1"& Test_condition!="cond3"&Test_condition!="cond4"&Test_condition!="cond6")
Table.Cond.Ctrl <- TableTest %>% filter(Test_condition!="cond1"& Test_condition!="cond3"& Test_condition!="cond4"&Test_condition!="cond6")

x11()
Figure4 <- ggplot(NbSuccess.Cond.Ctrl)+ geom_boxplot(aes(x=Test_condition, y=Perf),width = 0.5)+
  geom_jitter(aes(x=Test_condition, y=Perf, color=ID_individual, shape=Status),size=2.5,width=0.3, height=0)+
  geom_hline(yintercept=0.50, linetype="dashed", color="black", linewidth=0.7)+
  scale_shape_manual(values=c(19,1))+
  scale_color_manual(values = c("purple","blue","chocolate4","deepskyblue4","azure4","darkgreen","green","red","#ffcc00ff","aquamarine","#ff6600ff","chartreuse3","magenta"))+
  geom_label(data=Table.Cond.Ctrl, aes(x=Test_condition, y=1.2,label=NumberTrial))+
  scale_y_continuous(name="Succes rate", breaks=seq(from=0,to=1,by=0.25), limits=c(0,1.3))+
  scale_x_discrete(labels=c("2a-Quantity-proportions decorrelated","5a-Quantity-Ctrl","2b-Quantity-proportions decorrelated","5b-Quantity-Ctrl"))+
  ggtitle("Figure 4")+
  xlab("Conditions")+ Th

print(Figure4)

# Saving graph
svg("Figure4.svg", width = 9, height = 9.5) 
print(Figure4)
dev.off()

#####################################################
##### Group analyses associated to Figure 3 & 4 #####
#####################################################

## Creating a table of performance / condition
Tab.Wilcx.Cond <- NbSuccess.Cond.Ind %>%
  group_by(Test_condition)%>%
  dplyr::summarise(SumSuccess= sum(NbSucces), SumN = sum(N))%>%
  dplyr::mutate(Perf=SumSuccess/SumN)

## two-tailed Wilcoxon signed rank-test comparing median results of each condition chance level mu = 0,5

Condition <- c("cond1", "cond3", "cond4", "cond6", "cond2a", "cond5a", "cond2b", "cond5b")
Tab.Wilcx.Cond$TestWilcoxon <- "p-value"
Tab.Wilcx.Cond$Vpara <- "V"

for (i in 1:length(Condition)){
  NbSuccess.Cond <- NbSuccess.Cond.Ind %>% filter(Test_condition==Condition[i])
  Test_Wilcoxon <- wilcox.test(NbSuccess.Cond$Perf, alternative='two.sided', mu=0.5, paired = FALSE)
  Tab.Wilcx.Cond$TestWilcoxon[i] <- Test_Wilcoxon[["p.value"]]
  Tab.Wilcx.Cond$Vpara[i] <- Test_Wilcoxon[[1]]
}

print(Tab.Wilcx.Cond)

## two-tailed paired Wilcoxon signed rank-test : comparing results between 2 conditions of interest
# switching table
TablePaired <- NbSuccess.Cond.Ind %>% 
  dplyr::select(ID_individual,Test_condition,Perf)%>% 
  pivot_wider(names_from = Test_condition, values_from = Perf)

## Condition 1 vs Condition 4
median(TablePaired$cond1 - TablePaired$cond4, na.rm=TRUE)
wilcox.test(TablePaired$cond1, TablePaired$cond4, alternative='two.sided', paired=TRUE)
wilcox.test(TablePaired$cond1, TablePaired$cond4, alternative='greater', paired=TRUE)

## Condition 1 vs Condition 3
median(TablePaired$cond1 - TablePaired$cond3, na.rm=TRUE)
wilcox.test(TablePaired$cond1, TablePaired$cond3, alternative='two.sided', paired=TRUE)
wilcox.test(TablePaired$cond1, TablePaired$cond3, alternative='greater', paired=TRUE)

## Condition 3 vs Condition 4
median(TablePaired$cond3 - TablePaired$cond4, na.rm=TRUE)
wilcox.test(TablePaired$cond3, TablePaired$cond4, alternative='two.sided', paired=TRUE)
wilcox.test(TablePaired$cond4, TablePaired$cond3, alternative='greater', paired=TRUE)

#########################################################################
##### GLMER to test for effect of experimental factors and learning #####
#########################################################################

Null.ExpeFactor <- lme4::glmer(Success ~ (1|ID_individual), family = binomial, data = data)
Full.ExpeFactor <- lme4::glmer(Success ~ Arms_position + Fav_Jar_position + Experimenter_ID + Session + (1|Experimenter_position) + (1|ID_individual), family = binomial, data = data)

summary(Full.ExpeFactor)
Anova(Full.ExpeFactor)
anova(Null.ExpeFactor,Full.ExpeFactor,test="Chisq")


####################
##### Table 4 #####
###################
#Model selection with group data to compare effect of ROR, RQP and RQB on jar's choice

## Selection of data without cond6 (Ctrl-Opaque)
dataGLM <- data %>% 
  filter(Test_condition!="cond6")

# check data
unique(dataGLM$Test_condition)

## Null model

Null.StrategyResponse.GLMM <- glmer(Selected_side_binary ~ 1|ID_individual, family = binomial, data = dataGLM)

## Model (3)
Model3 <- glmer(Selected_side_binary ~ scale(RQB) + scale(RQP) + scale(ROR) + (1|ID_individual), family = binomial, data = dataGLM)
summary(Model3) # BIC = 1848
vif(Model3) 
check_model(Model3)
check_residuals(Model3)
check_convergence(Model3)

## Model (4)
Model4 <- glmer(Selected_side_binary ~ scale(RQB) + scale(RQP) + (1|ID_individual), family = binomial, data = dataGLM)
summary(Model4) # BIC = 1847.5
vif(Model4) 
check_model(Model4)
check_residuals(Model4)
check_convergence(Model4)


## Model (5)
Model5 <-glmer(Selected_side_binary ~ scale(RQP) + scale(ROR)+ (1|ID_individual), family = binomial, data = dataGLM)
summary(Model5) # AIC = 1849
vif(Model5)
check_model(Model5)
check_residuals(Model5)
check_convergence(Model5)

## Model (6)
Model6 <-glmer(Selected_side_binary ~ scale(RQB) + scale(ROR)+ (1|ID_individual), family = binomial, data = dataGLM)
summary(Model6) # AIC = 1931.8
vif(Model6)
check_model(Model6)
check_residuals(Model6)
check_convergence(Model6)

## Model (3) selected for model selection per individual

Model3 <- glmer(Selected_side_binary ~ scale(RQB) + scale(RQP) + scale(ROR) + (1|ID_individual), family = binomial, data = dataGLM)
summary(Model3) # BIC = 1848


###################
##### Table 5 #####
###################

## Selection of one model for each individual based on BIC criterion

Individual <- c("Abricot","Alaryc","Barnabe","Eric","Ficelle","Horus","Nema","Nereis","Olli","Walt")

for (i in 1:length(Individual)){
  IndGLM <- dataGLM %>% filter(ID_individual==Individual[i])
  Modele.Full <- glmulti(Selected_side_binary ~ ROR + RQB + RQP,
                         data=IndGLM,
                         crit=bic,
                         level=1,
                         method="h",
                         family=binomial,
                         fitfunction = glm,
                         confsetsize = 100)
  SumModele.Full <- summary(Modele.Full)
  print(Individual[i])
  print(SumModele.Full)
}
## Storing and summary of each model

#Abricot
AbricotGLM <- dataGLM %>% filter(ID_individual=="Abricot")
unique(AbricotGLM $ID_individual)
AbricotModel <- glm(Selected_side_binary ~ scale(RQP) + scale(RQB), family = binomial, data = AbricotGLM)
summary(AbricotModel)

#Alaryc
AlarycGLM <- dataGLM %>% filter(ID_individual=="Alaryc")
unique(AlarycGLM$ID_individual)
AlarycModel <- glm(Selected_side_binary ~ scale(RQP) + scale(ROR), family = binomial, data = AlarycGLM)
summary(AlarycModel)

#Barnabe
BarnabeGLM <- dataGLM %>% filter(ID_individual=="Barnabe")
unique(BarnabeGLM$ID_individual)
BarnabeModel <- glm(Selected_side_binary ~ scale(ROR), family = binomial, data = BarnabeGLM)
summary(BarnabeModel)

#Eric
EricGLM <- dataGLM %>% filter(ID_individual=="Eric")
unique(EricGLM$ID_individual)
EricModel <- glm(Selected_side_binary ~ scale(RQP) + scale(RQB), family = binomial, data = EricGLM)
summary(EricModel)

#Ficelle
FicelleGLM <- dataGLM %>% filter(ID_individual=="Ficelle")
unique(FicelleGLM$ID_individual)
FicelleModel <- glm(Selected_side_binary ~ scale(ROR), family = binomial, data = FicelleGLM)
summary(FicelleModel)

#Horus
HorusGLM <- dataGLM %>% filter(ID_individual=="Horus")
unique(HorusGLM$ID_individual)
HorusModel <- glm(Selected_side_binary ~ scale(RQP) + scale(ROR), family = binomial, data = HorusGLM)
summary(HorusModel)

#Nema
NemaGLM <- dataGLM %>% filter(ID_individual=="Nema")
unique(NemaGLM$ID_individual)
NemaModel <- glm(Selected_side_binary ~ scale(RQP), family = binomial, data = NemaGLM)
summary(NemaModel)

#Nereis
NereisGLM <- dataGLM %>% filter(ID_individual=="Nereis")
unique(NereisGLM$ID_individual)
NereisModel <- glm(Selected_side_binary ~ scale(RQP), family = binomial, data = NereisGLM)
summary(NereisModel)

#Olli
OlliGLM <- dataGLM %>% filter(ID_individual=="Olli")
unique(OlliGLM$ID_individual)
OlliModel <- glm(Selected_side_binary ~ scale(RQP) + scale(RQB), family = binomial, data = OlliGLM)
summary(OlliModel)

#Walt
WaltGLM <- dataGLM %>% filter(ID_individual=="Walt")
unique(WaltGLM$ID_individual)
WaltModel <- glm(Selected_side_binary ~ scale(ROR), family = binomial, data = WaltGLM)
summary(WaltModel)

####################
##### Figure S1 ####
####################

GrapeGroupe <- data %>%filter(ID_individual=="Eric"|ID_individual=="Ficelle"|ID_individual=="Horus"|ID_individual=="Olli"|ID_individual=="Walt")%>% 
  filter(Test_condition=="cond1"|Test_condition=="cond2a"|Test_condition=="cond3"|Test_condition=="cond4")%>%
  group_by(Test_condition)
GrapeGroupe$Success <- as.factor(GrapeGroupe$Success)

x11()
GrapheGrape <- ggplot(GrapeGroupe)+ geom_bar(aes(x=Raisin, fill=Success), position = "dodge")+facet_wrap(~Test_condition)+
  ggtitle("Performance in function of obtaining a grape per condition")+
  scale_x_discrete(labels=c("NO", "YES"))+
  theme_minimal()

svg("GrapheGrape.svg", width = 10, height = 10) 
print(GrapheGrape)
dev.off()

###################
##### Table S2 ####
###################

## Performance depending on whether a dried grape was given or not in this session in conditions of interest

GrapeProp <- data%>% filter(ID_individual=="Eric"|ID_individual=="Ficelle"|ID_individual=="Horus"|ID_individual=="Olli"|ID_individual=="Walt")%>%
  filter(Test_condition!="cond2b"&Test_condition!="cond6") %>%
  group_by(Test_condition,Raisin)%>%
  dplyr::summarise(TotSucces=sum(Success),N=n())%>%
  mutate(TotFail = N-TotSucces)%>%
  dplyr::select(Test_condition,Raisin,TotSucces, TotFail)

GrapeProp <- as.data.frame(GrapeProp)

# Chi-square test Condition 1
GrapeProp.Groupe.1 <- GrapeProp%>% filter(Test_condition=="cond1")%>% dplyr::select(TotSucces, TotFail)
chisq.test(GrapeProp.Groupe.1)

# Chi-square test Condition 2a
GrapeProp.Groupe.2a <- GrapeProp%>% filter(Test_condition=="cond2a")%>% dplyr::select(TotSucces, TotFail)
chisq.test(GrapeProp.Groupe.2a)

# Chi-square test Condition 3
GrapeProp.Groupe.3 <- GrapeProp%>% filter(Test_condition=="cond3")%>% dplyr::select(TotSucces, TotFail)
chisq.test(GrapeProp.Groupe.3)

# Chi-square test Condition 4
GrapeProp.Groupe.4 <- GrapeProp%>% filter(Test_condition=="cond4")%>% dplyr::select(TotSucces, TotFail)
chisq.test(GrapeProp.Groupe.4)

# Chi-square test Condition 5a
GrapeProp.Groupe.5a <- GrapeProp%>% filter(Test_condition=="cond5a")%>% dplyr::select(TotSucces, TotFail)
chisq.test(GrapeProp.Groupe.5a)

# Chi-square test Condition 5b
GrapeProp.Groupe.5b <- GrapeProp%>% filter(Test_condition=="cond5b")%>% dplyr::select(TotSucces, TotFail)
chisq.test(GrapeProp.Groupe.5b)

#####################
##### Figure S5 #####
#####################
#Exclusion of control opaque condition 6
dataGLM <- data %>% filter(Test_condition!="cond6")
unique(dataGLM$Test_condition)

# Selection of column of interest ROR, RQP and RQB
SubTableR <- dataGLM %>% dplyr::select(ROR, RQP, RQB)
SubTableR <- as.data.frame(SubTableR)

#Correlation Matrix
mcorR <-cor(SubTableR)
write.table(mcorR,"mcorR.csv",sep=";",row.names=FALSE)
rcorr(as.matrix(SubTableR), type=c("pearson"))

mcorR
FigureS5 <- corrplot(mcorR, type="upper", order="hclust", tl.col="black", tl.srt=45)

# Saving graph
svg("FigureS5.svg", width = 9, height = 9.5) 
print(FigureS5)
dev.off()

#####################
##### Figure S6 #####
#####################

## Mean performance per trial per condition of the group

Learning <- data %>% filter(Phase=="Test") %>%
  group_by(Test_condition,Trial) %>% 
  dplyr::summarise(NbSucces.Trial=sum(Success),N=n())%>%
  dplyr::mutate(Perf.Trial= NbSucces.Trial/N)%>%
  dplyr::mutate(RoundedPerf = round(Perf.Trial, digits =2))
Learning <- as.data.frame(Learning)
Learning$Trial <- as.integer(Learning$Trial)

str(Learning)

## Testing correlation between performance and trial number
# Test conditions 1, 3, 4 and 6
x11()
Learning.1.3.4.6 <- Learning %>% filter(Test_condition=="cond1"|Test_condition=="cond3"|Test_condition=="cond4"|Test_condition=="cond6")

GrapheLearning.1.3.4.6 <- ggplot(Learning.1.3.4.6, aes(x=Trial, y=Perf.Trial)) +
  facet_grid(~ Test_condition) +
  geom_point(shape=1) +    
  geom_smooth(method=lm,   
              se=FALSE, level = .95, colour = "#045a8d")+
  geom_smooth(colour="#bdc9e1", se=FALSE)+
  theme_bw() +
  labs(x = "Trial", y = "Performance") +
  expand_limits(y=c(0,1)) +
  scale_x_continuous(breaks=seq(1, 20, 1)) +
  scale_y_continuous(breaks=seq(0, 1, 0.1))+
  stat_cor(method = "pearson", label.x = 5, label.y = 1.1)

print(GrapheLearning.1.3.4.6)

# Saving graph
svg("GrapheLearning.1.3.4.6.svg", width = 15, height = 5) 
print(GrapheLearning.1.3.4.6)
dev.off()

# Control conditions 2a, 2b, 5a, 5b
Learning.2.5 <- Learning %>% filter(Test_condition=="cond2a"|Test_condition=="cond2b"|Test_condition=="cond5a"|Test_condition=="cond5b")
x11()
GrapheLearning.2.5 <- ggplot(Learning.2.5, aes(x=Trial, y=Perf.Trial)) +
  facet_grid(~ Test_condition) +
  geom_point(shape=1) +    
  geom_smooth(method=lm,   
              se=FALSE, level = .95, colour = "#045a8d")+
  geom_smooth(colour="#bdc9e1", se=FALSE)+
  theme_bw() +
  labs(x = "Trial", y = "Performance") +
  expand_limits(y=c(0,1)) +
  scale_x_continuous(breaks=seq(1, 20, 1)) +
  scale_y_continuous(breaks=seq(0, 1, 0.1))+
  stat_cor(method = "pearson", label.x = 5, label.y = 1.1)

print(GrapheLearning.2.5)

# Saving graph
svg("GrapheLearning.2.5.svg", width = 15, height = 5) 
print(GrapheLearning.2.5)
dev.off()

finalPlot <- ggarrange(GrapheLearning.1.3.4.6,  GrapheLearning.2.5,
                       ncol=4, nrow=2, common.legend = FALSE)
print(finalePlot)


#####################
##### Figure S7 #####
#####################

### Graph S7A

## Storage of models according to each strategies' group
StockMod.A <- list(AbricotModel,EricModel,OlliModel) #Strategy A: RQP + RQB
StockMod.B <- list(AlarycModel, HorusModel) #Strategy B: RQP + ROR
StockMod.C <- list(BarnabeModel, FicelleModel, WaltModel) #Strategy C: ROR
StockMod.D <- list(NereisModel, NemaModel) #Strategy D: RQP

##Storage of identity of individuals according to each strategies' group
Ind.A <- c("Abricot","Eric","Olli")
Ind.B<- c("Alaryc","Horus")
Ind.C <- c("Barnabe","Ficelle","Walt")
Ind.D <- c("Nema","Nereis")

## Storage Predictions of Model Strategy C : ROR
StockPredictionROR.C <- list()

for (i in 1:length(Ind.C)){
  dataInd <- dataGLM %>% filter(ID_individual==Ind.C[i])
  pred_data <- expand.grid(
    ROR = seq(min(dataInd$ROR), max(dataInd$ROR), length.out = 100)
  )
  pred_data <- as.data.frame(pred_data)
  predictions <- predict(StockMod.C[[i]], newdata = pred_data, type = "link", se.fit = TRUE, allow.new.levels = TRUE)
  pred_data$fit <- predictions$fit
  pred_data$se.fit <- predictions$se.fit
  pred_data <- pred_data %>%
    mutate(
      lower = plogis(fit - 1.96 * se.fit),
      upper = plogis(fit + 1.96 * se.fit),
      pred = plogis(fit)
    )
  StockPredictionROR.C <- c(StockPredictionROR.C, list(pred_data))
}
# Rename according to individuals' names
names(StockPredictionROR.C) <- Ind.C

## Storage Predictions of Model Strategy B : RQP + ROR
StockPredictionROR.B <- list()

for (i in 1:length(Ind.B)){
  dataInd <- dataGLM %>% filter(ID_individual==Ind.B[i])
  pred_data <- expand.grid(
    ROR = seq(min(dataInd$ROR), max(dataInd$ROR), length.out = 100),
    RQP= mean(dataInd$RQP)# Fixer z à sa moyenne pour la prédiction
  )
  pred_data <- as.data.frame(pred_data)
  predictions <- predict(StockMod.B[[i]], newdata = pred_data, type = "link", se.fit = TRUE, allow.new.levels = TRUE)
  pred_data$fit <- predictions$fit
  pred_data$se.fit <- predictions$se.fit
  pred_data <- pred_data %>%
    mutate(
      lower = plogis(fit - 1.96 * se.fit),
      upper = plogis(fit + 1.96 * se.fit),
      pred = plogis(fit)
    )
  StockPredictionROR.B <- c(StockPredictionROR.B, list(pred_data))
}
# Rename according to individuals' names
names(StockPredictionROR.B) <- Ind.B

## Table storage of performance depending on ROR per session per individual 
dataGLM$Selected_side_binary <- as.numeric(dataGLM$Selected_side_binary)
unique(dataGLM$Selected_side_binary)
Table.ROR <- dataGLM %>%
  filter(ID_individual!="Dory"&ID_individual!="Jeanne"&ID_individual!="Wallace")%>%
  group_by(ID_individual, Session, ROR)%>%
  dplyr::summarise(Somme_Success=sum(Selected_side_binary), N=n())%>%
  mutate(Perf_Session = Somme_Success/N)

## Graph Figure S7A
x11()
Graphe_ROR <- ggplot()+ 
  geom_jitter(data=Table.ROR,aes(x=ROR,y=Perf_Session, color=ID_individual),width=0.2, height=0, size=2 ) + #plot of performance per session per individual
  scale_color_manual(values = c("purple","blue","chocolate4","azure4","darkgreen","green","#ffcc00ff","aquamarine","#ff6600ff","magenta"))+
  #Strategy C 
  geom_line(data = StockPredictionROR.C[[1]], aes(x = ROR, y = pred), color = "chocolate4", linewidth = 1)+ #fitted model Barnabe
  geom_line(data = StockPredictionROR.C[[2]], aes(x = ROR, y = pred), color = "darkgreen", linewidth = 1)+ #fitted model Ficelle
  geom_line(data = StockPredictionROR.C[[3]], aes(x = ROR, y = pred), color = "magenta", linewidth = 1)+ # fitted model Walt
  #Strategy B
  geom_line(data = StockPredictionROR.B[[1]], aes(x = ROR, y = pred), color = "blue", linewidth = 1)+ # fitted model Alaryc
  geom_line(data = StockPredictionROR.B[[2]], aes(x = ROR, y = pred), color = "green", linewidth = 1)+ # fitted model Horus
  theme_minimal()+ theme(aspect.ratio = 1)+
  scale_y_continuous(name="Choice rate per session", breaks=seq(from=0,to=1,by=0.25), limits=c(0,1))

print(Graphe_ROR)

# Saving graph
svg("Graphe_ROR.svg", width = 10, height = 10) 
print(Graphe_ROR)
dev.off()

### Graph S7B

## Storage of models according to each strategies' group
StockMod.A <- list(AbricotModel,EricModel,OlliModel) #Strategy A: RQP + RQB
StockMod.B <- list(AlarycModel, HorusModel) #Strategy B: RQP + ROR
StockMod.C <- list(BarnabeModel, FicelleModel, WaltModel) #Strategy C: ROR
StockMod.D <- list(NereisModel, NemaModel) #Strategy D: RQP

##Storage of identity of individuals according to each strategies' group
Ind.A <- c("Abricot","Eric","Olli")
Ind.B<- c("Alaryc","Horus")
Ind.C <- c("Barnabe","Ficelle","Walt")
Ind.D <- c("Nema","Nereis")

## Storage Predictions of Model Strategy A : RQP + RQB
StockPredictionRQP.A <- list()

for (i in 1:length(Ind.A)){
  dataInd <- dataGLM %>% filter(ID_individual==Ind.A[i])
  pred_data <- expand.grid(
    RQP = seq(min(dataInd$RQP), max(dataInd$RQP), length.out = 100),
    RQB= mean(dataInd$RQB)# Fixer z à sa moyenne pour la prédiction
  )
  pred_data <- as.data.frame(pred_data)
  predictions <- predict(StockMod.A[[i]], newdata = pred_data, type = "link", se.fit = TRUE, allow.new.levels = TRUE)
  pred_data$fit <- predictions$fit
  pred_data$se.fit <- predictions$se.fit
  pred_data <- pred_data %>%
    mutate(
      lower = plogis(fit - 1.96 * se.fit),
      upper = plogis(fit + 1.96 * se.fit),
      pred = plogis(fit)
    )
  StockPredictionRQP.A <- c(StockPredictionRQP.A, list(pred_data))
}
# Rename according to individuals' names
names(StockPredictionRQP.A) <- Ind.A

## Storage Predictions of Model Strategy B : RQP + ROR
StockPredictionRQP.B <- list()

for (i in 1:length(Ind.B)){
  dataInd <- dataGLM %>% filter(ID_individual==Ind.B[i])
  pred_data <- expand.grid(
    RQP = seq(min(dataInd$RQP), max(dataInd$RQP), length.out = 100),
    ROR= mean(dataInd$ROR)
  )
  pred_data <- as.data.frame(pred_data)
  predictions <- predict(StockMod.B[[i]], newdata = pred_data, type = "link", se.fit = TRUE, allow.new.levels = TRUE)
  pred_data$fit <- predictions$fit
  pred_data$se.fit <- predictions$se.fit
  pred_data <- pred_data %>%
    mutate(
      lower = plogis(fit - 1.96 * se.fit),
      upper = plogis(fit + 1.96 * se.fit),
      pred = plogis(fit)
    )
  StockPredictionRQP.B <- c(StockPredictionRQP.B, list(pred_data))
}
# Rename according to individuals' names
names(StockPredictionRQP.B) <- Ind.B

## Storage Predictions of Model Strategy D : RQP
StockPredictionRQP.D <- list()

for (i in 1:length(Ind.D)){
  dataInd <- dataGLM %>% filter(ID_individual==Ind.D[i])
  pred_data <- expand.grid(
    RQP = seq(min(dataInd$RQP), max(dataInd$RQP), length.out = 100)
  )
  pred_data <- as.data.frame(pred_data)
  predictions <- predict(StockMod.D[[i]], newdata = pred_data, type = "link", se.fit = TRUE, allow.new.levels = TRUE)
  pred_data$fit <- predictions$fit
  pred_data$se.fit <- predictions$se.fit
  pred_data <- pred_data %>%
    mutate(
      lower = plogis(fit - 1.96 * se.fit),
      upper = plogis(fit + 1.96 * se.fit),
      pred = plogis(fit)
    )
  StockPredictionRQP.D <- c(StockPredictionRQP.D, list(pred_data))
}
# Rename according to individuals' names
names(StockPredictionRQP.D) <- Ind.D

## Table storage of performance depending on RQP per session per individual 
Table.RQP <- dataGLM %>%
  filter(ID_individual!="Dory"&ID_individual!="Jeanne"&ID_individual!="Wallace")%>%
  group_by(ID_individual, Session, RQP)%>%
  dplyr::summarise(Somme_Success=sum(Selected_side_binary), N=n())%>%
  mutate(Perf_Session = Somme_Success/N)

## Graph Figure S7B
x11()
Graphe_RQP <-  ggplot()+
  geom_jitter(data=Table.RQP,aes(x=RQP,y=Perf_Session, color=ID_individual),width=0.03, height=0, size=2 ) + #plot of performance per session per individual
  scale_color_manual(values = c("purple","blue","chocolate4","azure4","darkgreen","green","#ffcc00ff","aquamarine","#ff6600ff","magenta"))+
  #Strategy A 
  geom_line(data = StockPredictionRQP.A[[1]], aes(x = RQP, y = pred), color = "purple", linewidth = 1)+ #fitted model Abricot
  geom_line(data = StockPredictionRQP.A[[2]], aes(x = RQP, y = pred), color = "azure4", linewidth = 1)+ #fitted model Eric
  geom_line(data = StockPredictionRQP.A[[3]], aes(x = RQP, y = pred), color = "#ff6600f0", linewidth = 1)+ # fitted model Olli
  #Strategy B
  geom_line(data = StockPredictionRQP.B[[1]], aes(x = RQP, y = pred), color = "blue", linewidth = 1)+ # fitted model Alaryc
  geom_line(data = StockPredictionRQP.B[[2]], aes(x = RQP, y = pred), color = "green", linewidth = 1)+ # fitted model Horus
  #Strategy D
  geom_line(data = StockPredictionRQP.D[[1]], aes(x = RQP, y = pred), color = "#ffcc00f0", linewidth = 1)+ # fitted model Nema
  geom_line(data = StockPredictionRQP.D[[2]], aes(x = RQP, y = pred), color = "aquamarine", linewidth = 1)+ # fitted model Nereis
  theme_minimal()+coord_fixed(ratio=1)
print(Graphe_RQP)

# Saving graph
svg("Graphe_RQP.svg", width = 10, height = 10) 
print(Graphe_RQP)
dev.off()

### Graph S7C

## Storage of models according to each strategies' group
StockMod.A <- list(AbricotModel,EricModel,OlliModel) #Strategy A: RQP + RQB
StockMod.B <- list(AlarycModel, HorusModel) #Strategy B: RQP + ROR
StockMod.C <- list(BarnabeModel, FicelleModel, WaltModel) #Strategy C: ROR
StockMod.D <- list(NereisModel, NemaModel) #Strategy D: RQP

##Storage of identity of individuals according to each strategies' group
Ind.A <- c("Abricot","Eric","Olli")
Ind.B<- c("Alaryc","Horus")
Ind.C <- c("Barnabe","Ficelle","Walt")
Ind.D <- c("Nema","Nereis")

## Storage Predictions of Model Strategy A : RQP + RQB
StockPredictionRQB.A <- list()

for (i in 1:length(Ind.A)){
  dataInd <- dataGLM %>% filter(ID_individual==Ind.A[i])
  pred_data <- expand.grid(
    RQB = seq(min(dataInd$RQB), max(dataInd$RQB), length.out = 100),
    RQP= mean(dataInd$RQP)
  )
  pred_data <- as.data.frame(pred_data)
  predictions <- predict(StockMod.A[[i]], newdata = pred_data, type = "link", se.fit = TRUE, allow.new.levels = TRUE)
  pred_data$fit <- predictions$fit
  pred_data$se.fit <- predictions$se.fit
  pred_data <- pred_data %>%
    mutate(
      lower = plogis(fit - 1.96 * se.fit),
      upper = plogis(fit + 1.96 * se.fit),
      pred = plogis(fit)
    )
  StockPredictionRQB.A <- c(StockPredictionRQB.A, list(pred_data))
}

# Rename according to individuals' names
names(StockPredictionRQB.A) <- Ind.A

## Table storage of performance depending on RQB per session per individual 
Table.RQB <- dataGLM %>%
  filter(ID_individual!="Dory"&ID_individual!="Jeanne"&ID_individual!="Wallace")%>%
  group_by(ID_individual, Session, RQB)%>%
  dplyr::summarise(Somme_Success=sum(Selected_side_binary), N=n())%>%
  mutate(Perf_Session = Somme_Success/N)

## Graph Figure S7C
x11()
Graphe_RQB <- ggplot() +
  geom_jitter(data=Table.RQB,aes(x=RQB,y=Perf_Session, color=ID_individual),width=0.03, height=0, size=2 ) + #plot des perf par session des ind
  scale_color_manual(values = c("purple","blue","chocolate4","azure4","darkgreen","green","#ffcc00ff","aquamarine","#ff6600ff","magenta"))+
  #Strategie A 
  geom_line(data = StockPredictionRQB.A[[1]], aes(x = RQB, y = pred), color = "purple", linewidth = 1)+ #courbe du modele Abricot
  geom_line(data = StockPredictionRQB.A[[2]], aes(x = RQB, y = pred), color = "azure4", linewidth = 1)+ #courbe du modele Eric
  geom_line(data = StockPredictionRQB.A[[3]], aes(x = RQB, y = pred), color = "#ff6600f0", linewidth = 1)+ # courbe du modele Olli
  theme_minimal()+coord_fixed(ratio=1)
print(Graphe_RQB)

# Saving graph
svg("Graphe_RQB.svg", width = 10, height = 10) 
print(Graphe_RQB)
dev.off()






