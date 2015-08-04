####################################################
####Cumulative (proportional odds) logit analysis
####################################################
rm(list=ls())
gc()
graphics.off()
op<-par()
if(!require("MASS"))
  (install.packages("MASS"))
if(!require("ggplot2"))
  (install.packages("ggplot2"))
if(!require("kernlab"))
  (install.packages("kernlab"))
if(!require("e1071"))
  (install.packages("e1071"))
if(!require("gbm"))
  (install.packages("gbm"))
if(!require("plyr"))
  (install.packages("plyr"))
if(!require("dplyr"))
  (install.packages("dplyr"))
if(!require("reshape2"))
  (install.packages("reshape2"))
if(!require("nnet"))
  (install.packages("nnet"))
if(!require("RCurl"))
  (install.packages("RCurl"))
if(!require("devtools"))
  (install.packages("devtools"))
if(!require("reshape"))
  (install.packages("reshape"))
if(!require("RSNNS"))
  (install.packages("RSNNS"))
if(!require("Rcpp"))
  (install.packages("Rcpp"))
if(!require("neuralnet"))
  (install.packages("neuralnet"))


###Import Data
RAW_DATA<-read.csv("Raw_Data.csv")
RAW_DATA$Productivity.Group<-as.character(RAW_DATA$Productivity.Group)
RAW_DATA$Productivity.Group[ RAW_DATA$Productivity.Group=="B"]="A"
RAW_DATA$Productivity.Group[ RAW_DATA$Productivity.Group=="C"]="B"
RAW_DATA$Productivity.Group[ RAW_DATA$Productivity.Group=="D"]="C"
RAW_DATA$Race<-as.character(RAW_DATA$Race)
RAW_DATA$Race[ RAW_DATA$Race=="CAU"]="WH"


min_thresh<-min(RAW_DATA$Total.Genotype.Frequency[ RAW_DATA$Total.Genotype.Frequency!=0],na.rm=T)
RAW_DATA$Total.Genotype.Frequency[ RAW_DATA$Total.Genotype.Frequency==0]<-min_thresh
RAW_DATA$Total.Genotype.Frequency[ is.na(RAW_DATA$Total.Genotype.Frequency)]=min_thresh

min_thresh<-min(RAW_DATA$Frequency.1[ RAW_DATA$Frequency.1!=0],na.rm=T)
RAW_DATA$Frequency.1[ RAW_DATA$Frequency.1==0]<-min_thresh
RAW_DATA$Frequency.1[ is.na(RAW_DATA$Frequency.1)]=min_thresh

min_thresh<-min(RAW_DATA$Frequency.2[ RAW_DATA$Frequency.2!=0],na.rm=T)
RAW_DATA$Frequency.2[ RAW_DATA$Frequency.2==0]<-min_thresh
RAW_DATA$Frequency.2[ is.na(RAW_DATA$Frequency.2)]=min_thresh

IDs<-RAW_DATA$RID[1:1306]

permute<-rbinom(nrow(RAW_DATA),1,prob=0.5)==1
RAW_DATA[permute,]<-RAW_DATA[permute,c(1:6,8,7,10,9,11:13)]
DATA<-RAW_DATA[c("RID","Race","Frequency.1",
                 "Frequency.2","Rank.1","Rank.2",
                 "Total.Genotype.Frequency","Productivity.Group")]

DATA$Productivity.Group<-ordered(DATA$Productivity.Group,levels=c("C","B","A"))
colnames(DATA)<-c("RID","Race","H1","H2","Rank_H1","Rank_H2","GF","Productivity")
DATA<-DATA[complete.cases(DATA),]###remove NA's
DATA<-DATA[DATA$GF!=0,]
rm(RAW_DATA)

###Format Predictors to be ln-scale
DATA$GF<-log(DATA$GF)
DATA$H1<-log(DATA$H1)
DATA$H2<-log(DATA$H2)
#DATA[c("H1","H2","GF")]<-normalizeData(DATA[c("H1","H2","GF")],type="0_1")
summary(DATA)

####Split Train and Test
set.seed(1103)
train_idx<-sample(1:nrow(DATA),floor(nrow(DATA)*0.5))
#train_idx<-which(DATA$RID %in% IDs)
logical<-rep(FALSE,nrow(DATA))
logical[train_idx]<-TRUE
train_idx<-logical
test_idx<-!logical
rm(logical)
TRAIN<-DATA[train_idx,]
TEST<-DATA[test_idx,]



####summarize deomgraphics data
table(TRAIN$Race)
table(TEST$Race)

ftable(TRAIN$Race,TRAIN$Productivity)
ftable(TEST$Race,TEST$Productivity)


###Try a basic Proportional Odds Model on the genotype Frequencies
###first fit weights to the different classes
ww<-table(TRAIN$Productivity,TRAIN$Race)
ww<-as.data.frame(ww)
for(race in unique(ww$Var2)){
idx<-ww$Var2==race
ww$Freq[idx]<-ww$Freq[idx]/sum(ww$Freq[idx])
ww$Freq[idx]<-1/ww$Freq[idx]
}

###modify the weights
#ww$Freq[ ww$Var1=="B"]=ww$Freq[ ww$Var1=="B"]*1


TRAIN$weights_fit<-0
for(i in 1:nrow(ww)){
  idx<-TRAIN$Race==ww$Var2[i] & as.character(TRAIN$Productivity)==ww$Var1[i]
  TRAIN$weights_fit[idx]<-ww$Freq[i]
}

fit<-polr(Productivity~GF*Race,weights=weights_fit,data=TRAIN)
#fit<-polr(Productivity~GF*Race,data=TRAIN)
summary(fit)
probs<-fitted(fit)
display<-round(probs,2)
geno_values<-seq(min(DATA$GF),max(DATA$GF),length.out=10000)

###################################
####TRAIN Data Fit
####make charts for different races
###################################
par_disp<-unique(TRAIN$Race)
#par(op)
#par(mfrow=c(3,2))
race=par_disp[1]


correct<-data.frame()
for(race in par_disp){ 
  CUT_TRAIN<-TRAIN[TRAIN$Race==race, ]
  class_pred<-predict(fit,newdata=CUT_TRAIN)
  out<-data.frame(table(class_pred,CUT_TRAIN$Productivity))
  for(id in unique(out$Var2)){
    out$Freq[out$Var2==id]<-out$Freq[out$Var2==id]/sum(out$Freq[out$Var2==id])
  }
  
  out<-out[ out$class_pred==out$Var2,]
  c<-paste0(out$Var2,"=",round(out$Freq,2)*100,"%",collapse="; ")
  correct<-rbind(correct,data.frame(pct_correct=c,Race=race))
  
}


####set up race grids for prediction values
grid_values<-data.frame()

for(race in par_disp){ 
  grid<-data.frame("GF"=geno_values,"Race"=race)  
  preds<-as.data.frame(predict(fit,newdata=grid,type="probs"))
  grid$Race=paste0(race,": (",correct$pct_correct[correct$Race==race],")")
  grid_values<-rbind(grid_values,cbind(grid,preds)) 
}

grid<-melt(grid_values,id=c("GF","Race"))

colnames(grid)[3]<-"Prognosis"


####calculate cutoff values from grid
cutoff_values<-data.frame()
par_disp=unique(grid_values$Race)
for(race in par_disp){ 
 dat<-subset(grid_values,Race==race)
 dat<-dat[c("C","B","A")]
 Cutoff<-data.frame(ncol(dat)-1)
 for(i in 1:2){
   CUT<-array(1:10000)
   for(j in 1:10000){
     CUT[j]=dat[j,i+1]>dat[j,i]}
   intersect=10000-sum(CUT)
   Cutoff[i]<-(geno_values[intersect])
 }

cutoff_values<-rbind(cutoff_values,data.frame(cut=t(Cutoff),Race=race))   
}

cutoff_natural_scale<-cutoff_values
cutoff_natural_scale$cut<-exp(cutoff_natural_scale$cut)
print(cutoff_natural_scale)

cutoff_natural_scale$cut<-signif(cutoff_natural_scale$cut,2)
print(cutoff_natural_scale)


class(grid$Race)

TRAIN$Race<-as.character(TRAIN$Race)
mapping<-unique(grid$Race)
for ( b in unique(TRAIN$Race)){
  TRAIN$Race[TRAIN$Race==b]=mapping[grep(b,mapping)]
    
}

TRAIN$value<-0
TRAIN$value[ TRAIN$Productivity=="C"]<-1/6
TRAIN$value[ TRAIN$Productivity=="B"]<-1/6+1/3
TRAIN$value[ TRAIN$Productivity=="A"]<-1/6+2/3



plot_boundaries_train<-ggplot(data=grid)+geom_line(aes(x=GF,y=value,lty=Prognosis),lwd=1.25)+
  facet_wrap(~Race)+geom_vline(data=cutoff_values,aes(xintercept=cut),lwd=1,lty=1)+
  geom_text(data=TRAIN,aes(x=GF,y=value,label=Productivity),size=5,colour="black")+
  ylab("Prognosis Probability")+xlab("Log(Genotype Frequency)")+
  theme(strip.text.x = element_text(size = 15),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15))+
  scale_colour_grey(start=0,end=1)

print(plot_boundaries_train)

print(cutoff_values)



###################################
####TEST Data Fit
####make charts for different races
###################################
par_disp<-unique(TEST$Race)
#par(op)
#par(mfrow=c(3,2))
race=par_disp[1]


####get the confusion matrix for all combinations
all<-data.frame(prediction=predict(fit,newdata=TEST))
all$prediction<-all$prediction
all$truth<-TEST$Productivity
all$race<-TEST$Race

result<-ftable(all$race,all$truth,all$prediction)
print(result)
result<-data.frame(result)
colnames(result)<-c("Race","Truth","Prediction","Freq")

for(u in unique(result$Race)){
  for(b in unique(result$Truth)){
  idx<-result$Race==u & result$Truth==b
  result$Freq[idx]<-result$Freq[idx]/sum(result$Freq[idx])
  }
}

plot_dat<-arrange(result,Race,Truth,Prediction)
plot_dat<-ddply(plot_dat,c("Race","Truth"),transform,label_y=cumsum(Freq)-0.5*Freq)


plot_bars<-ggplot(data=plot_dat,aes(x=Truth,y=Freq,fill=Prediction))+
  geom_bar(stat="identity")+facet_wrap(~Race)+
  geom_text(aes(y=label_y,label=round(Freq*100,1)))+
  ggtitle("Percent Concordance on Validation Data by Class")+
  xlab("True Search Productivity Classification")+ylab("Percent")
print(plot_bars)



correct<-data.frame()
for(race in par_disp){ 
  CUT_TEST<-TEST[TEST$Race==race, ]
  class_pred<-predict(fit,newdata=CUT_TEST)
  out<-data.frame(table(class_pred,CUT_TEST$Productivity))
  for(id in unique(out$Var2)){
    out$Freq[out$Var2==id]<-out$Freq[out$Var2==id]/sum(out$Freq[out$Var2==id])
  }
  
  out<-out[ out$class_pred==out$Var2,]
  c<-paste0(out$Var2,"=",round(out$Freq,2)*100,"%",collapse="; ")
  correct<-rbind(correct,data.frame(pct_correct=c,Race=race))
}


####set up race grids for prediction values
grid_values<-data.frame()

for(race in par_disp){ 
  grid<-data.frame("GF"=geno_values,"Race"=race)  
  preds<-as.data.frame(predict(fit,newdata=grid,type="probs"))
  grid$Race=paste0(race,": (",correct$pct_correct[correct$Race==race],")")
  grid_values<-rbind(grid_values,cbind(grid,preds)) 
}

grid<-melt(grid_values,id=c("GF","Race"))

colnames(grid)[3]<-"Prognosis"


####calculate cutoff values from grid
cutoff_values<-data.frame()
par_disp=unique(grid_values$Race)
for(race in par_disp){ 
  dat<-subset(grid_values,Race==race)
  dat<-dat[c("C","B","A")]
  Cutoff<-data.frame(ncol(dat)-1)
  for(i in 1:2){
    CUT<-array(1:10000)
    for(j in 1:10000){
      CUT[j]=dat[j,i+1]>dat[j,i]}
    intersect=10000-sum(CUT)
    Cutoff[i]<-(geno_values[intersect])
  }
  
  cutoff_values<-rbind(cutoff_values,data.frame(cut=t(Cutoff),Race=race))   
}



class(grid$Race)

TEST$Race<-as.character(TEST$Race)
mapping<-unique(grid$Race)
for ( b in unique(TEST$Race)){
  TEST$Race[TEST$Race==b]=mapping[grep(b,mapping)]
  
}

TEST$value<-0
TEST$value[ TEST$Productivity=="C"]<-1/6
TEST$value[ TEST$Productivity=="B"]<-1/6+1/3
TEST$value[ TEST$Productivity=="A"]<-1/6+2/3




plot_boundaries_test<-ggplot(data=grid)+geom_line(aes(x=GF,y=value,lty=Prognosis),lwd=1.25)+
  facet_wrap(~Race)+geom_vline(data=cutoff_values,aes(xintercept=cut),lwd=1,lty=1)+
  geom_text(data=TEST,aes(x=GF,y=value,label=Productivity),size=5,colour="black")+
  ylab("Prognosis Probability")+xlab("Log(Genotype Frequency)")+
  theme(strip.text.x = element_text(size = 20),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=20),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15))+
  scale_colour_grey(start=0,end=1)

print(plot_boundaries_test)

print(cutoff_values)



##########################################################
########Make Plots of concordant searches for clinical data
########10/10's
##########################################################
CLINICAL<-read.csv("clinical_results_4_29_2015.csv",stringsAsFactors=F)
CLINICAL$Match<-0
CLINICAL$Match[ CLINICAL$Updated.10.10=="Y"]<-1

CLINICAL$GF<-CLINICAL$Genotype.Frequency
CLINICAL$GF[is.na(CLINICAL$GF) | (CLINICAL$GF==0)]<-min(CLINICAL$GF[CLINICAL$GF!=0],na.rm=T)
CLINICAL$GF<-log(CLINICAL$GF)
CLINICAL$Prediction<-predict(fit,newdata=CLINICAL)

ftable(CLINICAL[c("Prediction","Match","Race")])

cutoff_values_2<-cutoff_values

x<-strsplit(as.character(cutoff_values_2$Race),":")
x<-matrix(unlist(x),ncol=2,byrow = T)
cutoff_values_2$Race<-x[,1]
cutoff_values_2<-cutoff_values_2[ cutoff_values_2$Race!="UNK",]

h_lines<-data.frame()
u="API"
b="A"
for(u in unique(CLINICAL$Race)){
  for(b in unique(CLINICAL$Prediction)){
    idx<-CLINICAL$Race==u & CLINICAL$Prediction==b
    match_rate<-mean(CLINICAL$Match[idx])
    vals<-cutoff_values_2[cutoff_values_2$Race==u,]
    if(b=="A"){
      xmin<-max(vals$cut)
      xmax<-max(CLINICAL$GF)
    }
    if(b=="B"){
      xmin<-min(vals$cut)
      xmax<-max(vals$cut)
    }
    if(b=="C"){
      xmin<-min(CLINICAL$GF)
      xmax<-min(vals$cut)
    }
    
    xcoord<-c(xmin,xmax)
    xavg<-mean(c(xmin,xmax))
    h_lines<-rbind(h_lines,data.frame("Race"=u,"Class"=b,"Match_Rate"=match_rate,
                                      "END_Pts"=xcoord,"xcenter"=xavg))
  }
}

h_lines$Y_pos<-h_lines$Match_Rate
h_lines$Y_pos[ h_lines$Class %in% c("A")]<-h_lines$Y_pos[ h_lines$Class %in% c("A")]-0.1
h_lines$Y_pos[ h_lines$Class %in% c("C","B")]<-h_lines$Y_pos[ h_lines$Class %in% c("C","B")]+0.1

annotated_text<-data.frame()
for(i in unique(h_lines$Race)){
  for(b in unique(h_lines$Class)){
    idx<-which(h_lines$Race==i & h_lines$Class==b)[1]
    annotated_text<-rbind(annotated_text,h_lines[idx,])
  }
}

annotated_text$label<-round(annotated_text$Match_Rate*100,1)
annotated_text$label<-paste0(annotated_text$label,"%")
annotated_text$Class_desc<-as.character(annotated_text$Class)
annotated_text$Class_desc[annotated_text$Class_desc=="A"]="Good"
annotated_text$Class_desc[annotated_text$Class_desc=="B"]="Fair"
annotated_text$Class_desc[annotated_text$Class_desc=="C"]="Poor"



plot_clinical<-ggplot(data=CLINICAL)+geom_point(aes(y=Match,x=GF),position=position_jitter(w = 0, h = 0.01))+
  facet_wrap(~Race)+geom_vline(data=cutoff_values_2,aes(xintercept=cut),lwd=1,lty=1)+
  geom_line(data=h_lines,aes(x=END_Pts,y=Match_Rate,lty=Class),lwd=1)+
  geom_text(data=annotated_text,aes(x=xcenter,y=Y_pos,label=label),vjust=-0.4,size=7)+
  geom_text(data=annotated_text,aes(x=xcenter,y=Y_pos,label=Class_desc),vjust=1,size=7)+
  theme(legend.position="none")+
  scale_y_continuous(breaks=c(0,1),labels=c("No 10/10 Donor","10/10 Donor"))+
  ylab(NULL)+xlab("Log(Genotype Frequency)")+
  theme(strip.text.x = element_text(size = 20),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15))+
  scale_colour_grey(start=0,end=1)
  

print(plot_clinical)



##########################################################
########Make Plots of concordant searches for clinical data
########at least 9/10
##########################################################
CLINICAL<-read.csv("clinical_results_4_29_2015.csv",stringsAsFactors=F)
CLINICAL$Match<-0
CLINICAL$Match[ CLINICAL$At.least.a.9.10=="Y"]<-1

CLINICAL$GF<-CLINICAL$Genotype.Frequency
CLINICAL$GF[is.na(CLINICAL$GF) | (CLINICAL$GF==0)]<-min(CLINICAL$GF[CLINICAL$GF!=0],na.rm=T)
CLINICAL$GF<-log(CLINICAL$GF)
CLINICAL$Prediction<-predict(fit,newdata=CLINICAL)

cutoff_values_2<-cutoff_values

x<-strsplit(as.character(cutoff_values_2$Race),":")
x<-matrix(unlist(x),ncol=2,byrow = T)
cutoff_values_2$Race<-x[,1]
cutoff_values_2<-cutoff_values_2[ cutoff_values_2$Race!="UNK",]

h_lines<-data.frame()
for(u in unique(CLINICAL$Race)){
  for(b in unique(CLINICAL$Prediction)){
    idx<-CLINICAL$Race==u & CLINICAL$Prediction==b
    match_rate<-mean(CLINICAL$Match[idx])
    vals<-cutoff_values_2[cutoff_values_2$Race==u,]
    if(b=="A"){
      xmin<-max(vals$cut)
      xmax<-max(CLINICAL$GF)
    }
    if(b=="B"){
      xmin<-min(vals$cut)
      xmax<-max(vals$cut)
    }
    if(b=="C"){
      xmin<-min(CLINICAL$GF)
      xmax<-min(vals$cut)
    }
    
    xcoord<-c(xmin,xmax)
    xavg<-mean(c(xmin,xmax))
    h_lines<-rbind(h_lines,data.frame("Race"=u,"Class"=b,"Match_Rate"=match_rate,
                                      "END_Pts"=xcoord,"xcenter"=xavg))
  }
}

h_lines$Y_pos<-h_lines$Match_Rate
h_lines$Y_pos[ h_lines$Class %in% c("A","B")]<-h_lines$Y_pos[ h_lines$Class %in% c("A","B")]-0.1
h_lines$Y_pos[ h_lines$Class %in% c("C")]<-h_lines$Y_pos[ h_lines$Class %in% c("C")]+0.1

annotated_text<-data.frame()
for(i in unique(h_lines$Race)){
  for(b in unique(h_lines$Class)){
    idx<-which(h_lines$Race==i & h_lines$Class==b)[1]
    annotated_text<-rbind(annotated_text,h_lines[idx,])
  }
}

annotated_text$label<-round(annotated_text$Match_Rate*100,1)
annotated_text$label<-paste0(annotated_text$label,"%")
annotated_text$Class_desc<-as.character(annotated_text$Class)
annotated_text$Class_desc[annotated_text$Class_desc=="A"]="Good"
annotated_text$Class_desc[annotated_text$Class_desc=="B"]="Fair"
annotated_text$Class_desc[annotated_text$Class_desc=="C"]="Poor"



plot_clinical<-ggplot(data=CLINICAL)+geom_point(aes(y=Match,x=GF),position=position_jitter(w = 0, h = 0.01))+
  facet_wrap(~Race)+geom_vline(data=cutoff_values_2,aes(xintercept=cut),lwd=1,lty=1)+
  geom_line(data=h_lines,aes(x=END_Pts,y=Match_Rate,lty=Class),lwd=1)+
  geom_text(data=annotated_text,aes(x=xcenter,y=Y_pos,label=label),vjust=-0.4,size=7)+
  geom_text(data=annotated_text,aes(x=xcenter,y=Y_pos,label=Class_desc),vjust=1,size=7)+
  theme(legend.position="none")+
  scale_y_continuous(breaks=c(0,1),labels=c("No 9/10 Donor","at least 9/10 Donor"))+
  ylab(NULL)+xlab("Log(Genotype Frequency)")+
  theme(strip.text.x = element_text(size = 20),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size=15),
        axis.title.x = element_text(size=20),
        axis.title.y = element_text(size=15),
        legend.text=element_text(size=15),
        legend.title=element_text(size=15))+
  scale_colour_grey(start=0,end=1)


print(plot_clinical)

