library(pacman)
pacman::p_load(readr, dplyr, lubridate, ggplot2, geepack, lme4, lmtest, gee, MuMIn, car)
GLU=read_csv("C:/Users/andre/Desktop/Biostat 699/glucose.csv")
DEMO=read_csv("C:/Users/andre/Desktop/Biostat 699/demographics.csv")
#Data cleaning
DEMO=rename(DEMO,ID=id,Treatment=`Treatment group`,DiabDiag=`Date Diabetes Diagnosis`,
            heightIN=`Baseline Height in Inches`,heightM=`Baseline Height in Meters`,BP=`Baseline Blood Pressure`,
            SBP=`Baseline Systolic blood pressure`,DBP=`Baseline diastolic blood pressure`,weight=`Weight in kg at baseline`,
            BMI=`BMI at baseline`,HBA1C=`Baseline HbA1c`, Base1=`Date of Baseline for Drug 1`, Base2=`Date of Baseline for Drug 2`)
DEMO$Treatment=as.factor(DEMO$Treatment)
DEMO$Sex=as.factor(DEMO$Sex)
DEMO$Race=as.factor(DEMO$Race)
DEMO$Ethnicity=as.factor(DEMO$Ethnicity)

levels(DEMO$Treatment)=c("A","B")
levels(DEMO$Sex)=c("F","M")
levels(DEMO$Race)=c("Asian","Black","Unknown","White")
levels(DEMO$Ethnicity)=c("HSP","NonHSP")

GLU=rename(GLU, glucose=`Historic Glucose(mg/dL)`)

#Merge datasets
ALL=full_join(GLU,DEMO,by="ID")

full=unique(GLU$ID)
glucosedata=c(rep(0,length(DEMO$ID)))
DEMO=cbind(DEMO,glucosedata)
for (i in 1:length(full)){
  DEMO$glucosedata[which(DEMO$ID==full[i])]=1
}

CompDemo=filter(DEMO,glucosedata=='1')

#Summary statistics

unique(GLU$ID)
length(unique(GLU$ID))#25 patients with glucose tracked
barplot(table(DEMO$Treatment)) #19 in group A, 26 in group B
barplot(table(CompDemo$Treatment))#11 in group A, 14 in group B
barplot(table(DEMO$Sex))#19 females, 26 males
barplot(table(CompDemo$Sex))#11 females, 14 males
barplot(table(DEMO$Race))#30 white, 12 black, 2 Asian, 1 unknown
barplot(table(CompDemo$Race))#17 white, 6 black, 1 Asian, 1 unknown
barplot(table(DEMO$Ethnicity))#42 nonHispanic, 3 Hispanic
barplot(table(CompDemo$Ethnicity))#23 nonHispanic, 2 Hispanic

TIR=c(rep(0,length(GLU$ID)))
TAR=c(rep(0,length(GLU$ID)))
TBR=c(rep(0,length(GLU$ID)))
GLU=cbind(GLU,TIR,TAR,TBR)
for (i in 1:length(GLU$ID)){
  if (GLU$glucose[i]<70){
    GLU$TBR[i]=1
  } else if (GLU$glucose[i]>180){
    GLU$TAR[i]=1
  } else {
    GLU$TIR[i]=1
  }
}

GLUfull=filter(GLU, ID!=35 & ID!=41 & ID!=42 & ID!=44)
measurement=c(rep(0,length(GLUfull$Date)))
stage=c(rep(1,length(GLUfull$Date)))
GLUfull=cbind(GLUfull,measurement,stage)

initial=0
cutpoint=1
for (j in 1:length(GLUfull$Date)){
  if (j!=1 && GLUfull$ID[j]!=GLUfull$ID[j-1]){
    initial=0
    GLUfull$measurement[j]=0
    initial=initial+1
    cutpoint=j
    GLUfull$stage[j:length(GLUfull$Date)]=1
  } else {
    if (mdy(GLUfull$Date[j])-mdy(GLUfull$Date[cutpoint])<30){
      GLUfull$measurement[j]=initial
      initial=initial+1
    } else {
      initial=0
      GLUfull$measurement[j]=0
      initial=initial+1
      cutpoint=j
      GLUfull$stage[j:length(GLUfull$Date)]=2
    }
  }
}

GLUfull$stage=as.factor(GLUfull$stage)
COMP=left_join(GLUfull,DEMO,by="ID")

drug = c(rep(NULL,length(GLUfull$ID)))
for (i in 1:length(COMP$ID)){
  if (COMP$Treatment[i]=='A' & COMP$stage[i]=='1'){
    drug[i]='X'
  } else if (COMP$Treatment[i]=='A' & COMP$stage[i]=='2'){
    drug[i]='Y'
  } else if (COMP$Treatment[i]=='B' & COMP$stage[i]=='1'){
    drug[i]='Y'
  } else {
    drug[i]='X'
  }
}
COMP=cbind(COMP,drug)
COMP$drug=as.factor(COMP$drug)

ids=unique(COMP$ID)
for (i in 1:length(ids)){
  plot=ggplot(data=filter(COMP, ID==ids[i]), aes(x=measurement, y=glucose, group=drug))+geom_line(aes(color=drug), size=1) +ggtitle(
    paste0("Time Series of Glucose Levels for Patient ID ", toString(ids[i]))) + xlab("Measurements Since Baseline") +ylab(
      "glucose (mg/dL)") + theme_bw() + theme(plot.background=element_blank(), panel.grid.major= element_blank(), panel.grid.minor= element_blank()) + 
    theme(plot.title = element_text(hjust=0.5, size=18))+theme(axis.text=element_text(size=12),axis.title=element_text(size=16,))
  print(plot)
}


#Collapsing over day
Day=c(rep(0, length(COMP$ID)))
for (i in 1:length(COMP$ID)){
  #Day[i]=COMP$measurement[i]%/%96 + 1
  if (COMP$stage[i]=='1'){
    Day[i]=as.numeric(mdy(COMP$Date[i])-mdy(COMP$Base1[i]))
  } else {
    Day[i]=as.numeric(mdy(COMP$Date[i])-mdy(COMP$Base2[i]))
  }
}
COMP=cbind(COMP, Day)

avg=c()
std=c()
cv=c()
ptir=c()
ptar=c()
ptbr=c()
ID=c()
Treatment=c()
drug=c()
day=c()
stage=c()
intervals=c()
k=1
for (i in 1:length(ids)){
  temp=filter(COMP, ID==ids[i] & stage=='1')
  dayvec=unique(temp$Day)
  for (j in 1:length(dayvec)){
    newtemp=filter(temp, Day==dayvec[j])
    avg[k]=mean(newtemp$glucose)
    std[k]=sd(newtemp$glucose)
    cv[k]=std[k]/avg[k]*100
    ptir[k]=sum(newtemp$TIR)/length(newtemp$ID)*100
    ptar[k]=sum(newtemp$TAR)/length(newtemp$ID)*100
    ptbr[k]=sum(newtemp$TBR)/length(newtemp$ID)*100
    ID[k]=ids[i]
    Treatment[k]=as.character(newtemp$Treatment[1])
    drug[k]=as.character(newtemp$drug[1])
    day[k]=dayvec[j]
    stage=1
    intervals[k]=length(newtemp$ID)/96
    k = k + 1
  }
  temp=filter(COMP, ID==ids[i] & stage=='2')
  dayvec=unique(temp$Day)
  for (j in 1:length(dayvec)){
    newtemp=filter(temp, Day==dayvec[j])
    avg[k]=mean(newtemp$glucose)
    std[k]=sd(newtemp$glucose)
    cv[k]=std[k]/avg[k]*100
    ptir[k]=sum(newtemp$TIR)/length(newtemp$ID)*100
    ptar[k]=sum(newtemp$TAR)/length(newtemp$ID)*100
    ptbr[k]=sum(newtemp$TBR)/length(newtemp$ID)*100
    ID[k]=ids[i]
    Treatment[k]=as.character(newtemp$Treatment[1])
    drug[k]=as.character(newtemp$drug[1])
    day[k]=dayvec[j]
    stage=2
    intervals[k]=length(newtemp$ID)/96
    k = k + 1
  }
}

DAYS=data_frame(ID, day, drug, avg, std, cv, ptir, ptar, ptbr, intervals)
DAYS$drug=as.factor(DAYS$drug)
DAYS=left_join(DAYS, CompDemo, by='ID')

ggplot(data=DAYS, aes(x=day, y=avg, group=ID)) + geom_line() + facet_grid(.~drug)
ggplot(data=DAYS, aes(x=day, y=std, group=ID)) + geom_line() + facet_grid(.~drug)
ggplot(data=DAYS, aes(x=day, y=cv, group=ID)) + geom_line() + facet_grid(.~drug)


#check for collinearity
cor(DAYS$avg, DAYS$Age)
cor(DAYS$avg, DAYS$BMI)
cor(DAYS$avg, DAYS$day)
cor(DAYS$std, DAYS$Age)
cor(DAYS$std, DAYS$BMI)
cor(DAYS$std, DAYS$day)
cor(DAYS$cv, DAYS$Age)
cor(DAYS$cv, DAYS$BMI)
cor(DAYS$cv, DAYS$day)

#GEE--compound symmetry var-cov
model1A=geeglm(avg~drug+Treatment+Age+Sex+BMI+SBP+DBP+day, id=ID, family='gaussian', data=DAYS, corstr='exchangeable',  weights=intervals)
summary(model1A)$coefficients
model2A=geeglm(std~drug+Treatment+Age+Sex+BMI+SBP+DBP+day, id=ID, family='gaussian', data=DAYS, corstr='exchangeable',  weights=intervals)
summary(model2A)$coefficients
model3A=geeglm(cv~drug+Treatment+Age+Sex+BMI+SBP+DBP+day, id=ID, family='gaussian', data=DAYS, corstr='exchangeable',  weights=intervals)
summary(model3A)$coefficients

#GEE--autoregressive var-cov
model1B=geeglm(avg~drug+Treatment+Age+Sex+BMI+SBP+DBP+day, id=ID, family='gaussian', data=DAYS, corstr='ar1',  weights=intervals)
summary(model1B)$coefficients
model2B=geeglm(std~drug+Treatment+Age+Sex+BMI+SBP+DBP+day, id=ID, family='gaussian', data=DAYS, corstr='ar1',  weights=intervals)
summary(model2B)$coefficients
model3B=geeglm(cv~drug+Treatment+Age+Sex+BMI+SBP+DBP+day, id=ID, family='gaussian', data=DAYS, corstr='ar1',  weights=intervals)
summary(model3B)$coefficients

#GEE--independent var-cov
model1C=geeglm(avg~drug+Treatment+Age+Sex+BMI+SBP+DBP+day, id=ID, family='gaussian', data=DAYS, corstr='independence',  weights=intervals)
summary(model1C)$coefficients
model2C=geeglm(std~drug+Treatment+Age+Sex+BMI+SBP+DBP+day, id=ID, family='gaussian', data=DAYS, corstr='independence',  weights=intervals)
summary(model2C)$coefficients
model3C=geeglm(cv~drug+Treatment+Age+Sex+BMI+SBP+DBP+day, id=ID, family='gaussian', data=DAYS, corstr='independence',  weights=intervals)
summary(model3C)$coefficients

model.sel(model1A, model1B, model1C, rank=QIC)
model.sel(model2A, model2B, model2C, rank=QIC)
model.sel(model3A, model3B, model3C, rank=QIC)
#Use ar1 var-cov for models 2&3, exchangeable for 1
#GEE models with variable p-values
modelA=geeglm(avg~drug+Treatment+Age+Sex+BMI+SBP+DBP+day, id=ID, family='gaussian', data=DAYS, corstr='independence', weights=intervals)
modelB=geeglm(std~drug+Treatment+Age+Sex+BMI+SBP+DBP+day, id=ID, family='gaussian', data=DAYS, corstr='ar1', weights=intervals)
modelC=geeglm(cv~drug+Treatment+Age+Sex+BMI+SBP+DBP+day, id=ID, family='gaussian', data=DAYS, corstr='ar1', weights=intervals)
summary(modelA)$coefficients
summary(modelB)$coefficients
summary(modelC)$coefficients

#check for multicollinearity
vif(modelA)
vif(modelB)
vif(modelC)

#check residual diagnostics
ResidDataA=cbind.data.frame(DAYS,"Residuals"=residuals(modelA), "Fitted" =fitted(modelA))
ggplot(data=ResidDataA) + ggtitle("Model Residuals vs. Fitted Average Daily Glucose Levels, GEE") + geom_point(aes(x=Fitted, y=Residuals)) +
  theme_bw() + theme(plot.background=element_blank(), panel.grid.major= element_blank(), panel.grid.minor= element_blank()) + 
  theme(plot.title = element_text(hjust=0.5, size=14))
cor.test(residuals(modelA),fitted(modelA))

ResidDataB=cbind.data.frame(DAYS,"Residuals"=residuals(modelB), "Fitted" =fitted(modelB))
ggplot(data=ResidDataB) + ggtitle("Model Residuals vs. Fitted Daily Glucose Level Standard Deviation, GEE") + geom_point(aes(x=Fitted, y=Residuals)) +
  theme_bw() + theme(plot.background=element_blank(), panel.grid.major= element_blank(), panel.grid.minor= element_blank()) + 
  theme(plot.title = element_text(hjust=0.5, size=14))
cor.test(residuals(modelB),fitted(modelB))

ResidDataC=cbind.data.frame(DAYS,"Residuals"=residuals(modelC), "Fitted" =fitted(modelC))
ggplot(data=ResidDataC) + ggtitle("Model Residuals vs. Fitted Daily Glucose Coefficient of Variation, GEE") + geom_point(aes(x=Fitted, y=Residuals)) +
  theme_bw() + theme(plot.background=element_blank(), panel.grid.major= element_blank(), panel.grid.minor= element_blank()) + 
  theme(plot.title = element_text(hjust=0.5, size=14))
cor.test(residuals(modelC),fitted(modelC))