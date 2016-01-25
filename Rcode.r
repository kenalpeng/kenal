# =======================================================================================================
# 1.Hazard Model with log(time)
# Estimation sample: 2000Q1~ 2012Q4;   
# Estimation Outsample forecast evalautaion: 2008Q1~2012Q4

# Auto.arima forecast: 
#    1). Company variables: past two years 
#    2). Macro variables: (1)None v.s. 
#                         (2)Auto.arima forecast of past two years (period=8)  

# 2. Two step : First step is created dataspace for variable then use Linear interpolation and Auto.arima(8 periods)
#               to create missing value and forecast.  
 
 
rm(list=ls());
install.packages("INLA", repos="http://www.math.ntnu.no/inla/R/stable")
library(survival)
library(ROCR)
library(forecast)
library(moments)
library(INLA)
library(Cairo)
#inla.upgrade(testing=TRUE)
# =================================================================================
# ================================== Step 1 =======================================
# ================================   Setup    =====================================
deventcode=1 # 1: ues news default event   ; 0: use last data be default time

# =================================== Input informatin ======================================================
vars=c(1:24)
# Code map for def and nondef firm data 
# -------------------------------------------------------------------------------------------------------------------------------------------
#      [1]    [2]       [3]     [4]        [5]   [6]    [7]    [8]    [9]      [10]
# [1]  "firm" "yyyymm"  "code"  "datatime" "TA"  "TLTA" "TLTE" "CLTL" "OCF2TA" "CLTA" 
# [11] "LR"   "LiquidR" "RE2TA" "ROA"      "EPS" "OCR"  "RT"   "DAT"  "IT"     "TA.1" 
# [21] "CL"   "TL" "CA" "TE"    "OCF"      "NS"  "OC"   "GP"   "inv"  "lend"      
#  
# =================================================================================================================================
yyyymm=201303 ; # The End time of dataspace
newsyear=c(2008:2012) # The years of News data
nondef=read.table('F:/0_Study/6_Theses/code/main_code/firmdata/20131220/nondef20131220.txt',header=TRUE)
def=read.table('F:/0_Study/6_Theses/code/main_code/firmdata/20131220/def20131220.txt',header=TRUE) 
Devent=read.table('F:/0_Study/6_Theses/code/Cdata/defaultevent.txt',header=TRUE) #違約公司時間點(TCRI違約新聞事件)
Devent$code=1
marc_mark=read.table('F:/0_Study/6_Theses/code/Cdata/marc_mark.txt',header=TRUE) 
moneyS=read.table('F:/0_Study/6_Theses/code/Cdata/money.txt',header=TRUE)

# Number of default firm at TCRI rating
yyyyNdef0=read.table('F:/0_Study/6_Theses/code/main_code/firmdata/20131220/yyyyNdef.txt',header=TRUE)
TCRI=read.table('F:/0_Study/6_Theses/code/main_code/TCRI.txt',header=TRUE) 
TEJ10=read.table('F:/0_Study/6_Theses/code/Cdata/TEJ10.txt',header=TRUE) # 產業分類碼(未違約) 
TEJ11=read.table('F:/0_Study/6_Theses/code/Cdata/TEJ11.txt',header=TRUE) # 產業分類碼(違約)
TimeMapping=read.table("F:/0_Study/6_Theses/code/Cdata/TimeMapping.txt",header=TRUE) ;# TimeMapping.txt

firm1_keyword=read.table("F:/0_Study/6_Theses/code/main_code/firm_code1.txt");
firm2_keyword=read.table("F:/0_Study/6_Theses/code/main_code/firm_code2.txt");
firm2_keyword=c(as.matrix(firm2_keyword))
keywordpath=sprintf("F:/0_Study/6_Theses/code/main_code/GBN/GBNkeyword36_%s.txt",newsyear)
period=8  # forcasting period
# ==========================================================================
outputname=c('ParametorEstimateT','ParametorEstimate_WithoutMarcT','logisticIntensityT',
               'summaryT','RMSE','Aver_intensity','News')
temppath=sprintf("F:/0_Study/6_Theses/code/temp.txt")
outpath=sprintf('F:/0_Study/6_Theses/code/output/%s_%s_1211.txt',outputname,(deventcode+1))
pdfoutput=sprintf('F:/0_Study/6_Theses/code/output/%s_%s_1211.pdf',outputname,(deventcode+1))

# ==========================================================================
# ===========================  Model (DHM & MCMC) =========================================
myfun10=code ~ I(log(TimeMapping)) + TLTA_1 + LiquidR + RE2TA + EPS + OCR 
myfun11=code ~ I(log(TimeMapping)) + TLTA_1 + LiquidR + RE2TA + EPS + OCR + unemployment+ stock

# ==================================  End setup ===========================================
# =========================================================================================
nondef=nondef[nondef$yyyymm<yyyymm,] ; def=def[def$yyyymm<yyyymm,]
nondef_NA=nondef[,vars] ; def_NA=def[,vars]
nondef=na.omit(nondef) ; def=na.omit(def)
 
# ========================================== Function =====================================================

## ================= data summary ===============
# All data
NAsumary=function(data){
nobNA=matrix(numeric(ncol(data)*3),nrow=3)
dataspace=data.frame(data)
for(i in 1:ncol(data)){
tmp=length(data[,i])-length(na.omit(data[,i]))
nobNA[1,i]=length(data[,i]);
nobNA[2,i]=tmp;
}
nobNA[3,]=nobNA[1,]-nobNA[2,]
colnames(nobNA)=colnames(data)
rownames(nobNA)=c('All_nobs','NA_nobs','Diff._nobs')
return(nobNA)
}

## Data summary
# =============================================================================
var_summary=function(d){
var_summary=matrix(numeric(ncol(d)*8),nrow=8)
for(i in 1:ncol(d)){
var_summary[1,i]= length(na.omit(d[,i]))
var_summary[2,i]= mean(na.omit(d[,i]))
var_summary[3,i]= sd(na.omit(d[,i]))
var_summary[4,i]= median(na.omit(d[,i]))
var_summary[5,i]= max(na.omit(d[,i]))
var_summary[6,i]= min(na.omit(d[,i]))
var_summary[7,i]= skewness(na.omit(d[,i]))
var_summary[8,i]= kurtosis(na.omit(d[,i]))
}
colnames(var_summary)=colnames(d) ;
rownames(var_summary)=c('Number','mean','sd','median','max','min','skewness','kurtosis')
return(var_summary)
}
# =============================================================================
# ======= Create new forcasting data ; use auto.arima and Order of first-differencing is 1
forcastVar=function(firm,Cdata,period){
 firm=firm[firm$V2>7,]
 newvars={}
 for(m in 1:nrow(firm)){
  newvars0=matrix(numeric(period*ncol(Cdata)),ncol=ncol(Cdata))
  tmp=as.matrix(Cdata[Cdata$firm==firm$V1[m],])
  tmp0=rbind(tmp,newvars0)
    for(j in 1:period){
      for(i in 5:(ncol(tmp0))){
	   p=tmp0[2:(nrow(tmp)+j-1),i]-tmp0[1:(nrow(tmp)+j-2),i]
	   if(sum(p)==0){tmp0[(nrow(tmp)+j),i]=mean(tmp0[2:(nrow(tmp)+j-1),i])}
       if(sum(p)!=0){ 
	   vari0=auto.arima(as.numeric(tmp0[1:(nrow(tmp)+j-1),i]))
       vari=length(vari0$model$phi)+length(vari0$model$theta)+length(vari0$model$Delta)
       #vari=forecast(auto.arima(tmp0[1:(nrow(tmp)+j-1),i]),h=1)$mean
          if(vari!=0){tmp0[(nrow(tmp)+j),i]=forecast(auto.arima(as.numeric(tmp0[1:(nrow(tmp)+j-1),i])),h=1)$mean}
          if(vari==0){tmp0[(nrow(tmp)+j),i]=forecast(auto.arima(as.numeric(tmp0[1:(nrow(tmp)+j-1),i]),d=1),h=1)$mean}
          }
	   } 
    tmp0[(nrow(tmp)+1):(nrow(tmp)+period),1]=rep(tmp[1,1],period)
    tmp0[(nrow(tmp)+1):(nrow(tmp)+period),c(2,4)]=as.matrix(TimeMapping0[TimeMapping0$mapping > tmp[nrow(tmp),4] & 
                                                             TimeMapping0$mapping < tmp[nrow(tmp),4]+period+1,])
   }  
newvars=rbind(newvars,tmp0[(nrow(tmp)+1):(nrow(tmp)+period),])
}
forcastVar=newvars
}
# ============================================================================================================
#  combine data files (s1 & s2) as one. (yyyymm  vars.)
# ===============================================
CombidD <- function(s1,s2){
CombidD1=cbind(s1,matrix(rep(NA,nrow(s1)*(ncol(s2)-1)),ncol=(ncol(s2)-1)))
n=ncol(s1)+1 ;
for(i in 1:nrow(CombidD1)){
tmp=as.matrix(s2[s2[,1]==CombidD1[i,1],2:ncol(s2)])
  if(nrow(tmp)==1){
    CombidD1[i,n:(n+ncol(tmp)-1)]=as.matrix(tmp)
   }
 }
 colnames(CombidD1)=c(colnames(s1),colnames(s2)[2:ncol(s2)])
 return(data.frame(CombidD1))
}

CombidD0 <- function(s1,s2,chara){
CombidD1=cbind(s1,matrix(rep(chara,nrow(s1)*(ncol(s2)-1)),ncol=(ncol(s2)-1)))
n=ncol(s1)+1 ;
for(i in 1:nrow(CombidD1)){
tmp=as.matrix(s2[s2[,1]==CombidD1[i,1],2:ncol(s2)])
  if(nrow(tmp)==1){
    CombidD1[i,n:(n+ncol(tmp)-1)]=as.matrix(tmp)
   }
 }
 colnames(CombidD1)=c(colnames(s1),colnames(s2)[2:ncol(s2)])
 return(data.frame(CombidD1))
}

# s1 : firm_code yyyymm vars. ; s2 : firm_code yyyymm vars.
CombidD1 <- function(s1,s2,chara){
CombidD1=cbind(s1,matrix(rep(chara,nrow(s1)*(ncol(s2)-2)),ncol=(ncol(s2)-2)))
n=ncol(s1)+1 ;
for(i in 1:nrow(CombidD1)){
tmp=as.matrix(s2[s2[,1]==CombidD1[i,1] & s2[,2]==CombidD1[i,2],3:ncol(s2)])
  if(nrow(tmp)==1){
    CombidD1[i,n:(n+ncol(tmp)-1)]=as.matrix(tmp)
   }
 }
 colnames(CombidD1)=c(colnames(s1),colnames(s2)[3:ncol(s2)])
 return(data.frame(CombidD1))
}
#write.table(CombidD1,'F:/1_Lin/2013/data/alldataCb.txt',row.names = F)

# =======================================================================================================
# Newsdata : data space with three part. The first is date (yyyymmdd). The second is keyword. 
#            The three is times of 'Good','Bad' and 'Neutral' in keyword.txt from keyword.r program.
# word : The word of want to calculate the prob. which must contain Big5 column "matrix".
# Ex. for marc word : Inf=Information(newsdata=news_tmpQ1,word=matrix(Marc_keyword,ncol=1)) 
# Ex. for firm word : firmname=cbind(firm1_keyword,firm2_keyword)
#                     Inf=Information(newsdata=news_tmpQ1,word=firmname)

Information = function(newsdata,word){
firm_prob={} ;prob_tmp={};
if(ncol(word)>1){
   word=cbind(1:nrow(word),word)
   colnames(word)=c('map','firm','Big5')
   word=data.frame(word)
 for(i in 1:nrow(word)){
  if (sum(colnames(newsdata)==word$Big5[i])!=0){
      tmp=cbind(c(1:nrow(newsdata)),newsdata[,colnames(newsdata)==word$Big5[i]])
	  tmp0=tmp[tmp[,2]!=0,1] ; 
	  if(length(tmp0)!=0){
	  tmp1=sum(newsdata$GNews[tmp0]) ;  tmp2=sum(newsdata$BNews[tmp0]) ; tmp3=sum(newsdata$NNews[tmp0]) ; tmp4=length(tmp0)
	 #tmp1=sum(newsdata$GNews[tmp0])/length(tmp0)  ;  tmp2=sum(newsdata$Bews[tmp0])/length(tmp0) ; tmp3=sum(newsdata$Uews[tmp0])/length(tmp0)
	 #tmp1=sum(newsdata$GNews[tmp0])/length(tmp0)  ;  tmp2=sum(newsdata$BNews[tmp0])/length(tmp0) ; tmp3=sum(newsdata$UNews[tmp0])/length(tmp0)
	 prob_tmp=cbind(word$firm[i],tmp1,tmp2,tmp3,tmp4)
	  } else {
          prob_tmp=cbind(word$firm[i],0,0,0,0)
		  }
     } else {
      prob_tmp=cbind(word$firm[i],0,0,0,0)
   }	 
   firm_prob=rbind(firm_prob,prob_tmp)
}
firm_prob=data.frame(cbind(firm_prob[,1],round(newsdata$yyyymmdd[nrow(newsdata)]/100),firm_prob[,2:ncol(firm_prob)]))
colnames(firm_prob)=c('firm','yyyymm','Good_News','Bad_News','Neutral_News','All_News')
} else {
   word=cbind(1:nrow(word),word)
   colnames(word)=c('map','Big5')
   word=data.frame(word)
 for(i in 1:nrow(word)){
  if (sum(colnames(newsdata)==word$Big5[i])!=0){
      tmp=cbind(c(1:nrow(newsdata)),newsdata[,colnames(newsdata)==word$Big5[i]])
	  tmp0=tmp[tmp[,2]!=0,1] ; 
        if(length(tmp0)!=0){ 	  
	      tmp1=sum(newsdata$GNews[tmp0])   ;  tmp2=sum(newsdata$BNews[tmp0])  ; tmp3=sum(newsdata$NNews[tmp0]) ; tmp4=length(tmp0)
	      prob_tmp=c(sprintf('%s',word$Big5[i]),c(tmp1,tmp2,tmp3,tmp4))
	     } else {
          prob_tmp=c(sprintf('%s',word$Big5[i]),0,0,0,0)
		  }
     } else {
      prob_tmp=c(sprintf('%s',word$Big5[i]),0,0,0,0)
   }	 
   firm_prob=rbind(firm_prob,prob_tmp)
}
yyyymm=round(newsdata$yyyymmdd[nrow(newsdata)]/100)
firm_prob=data.frame(cbind(matrix(firm_prob[,1],ncol=1),yyyymm,matrix(as.numeric(firm_prob[,c(2:ncol(firm_prob))]),ncol=(ncol(firm_prob)-1))))
colnames(firm_prob)=c('Big5','yyyymm','Good_News','Bad_News','Neutral_News','All_News');
firm_prob$Good_News=as.numeric(matrix(firm_prob$Good_News,ncol=1)) 
firm_prob$Bad_News=as.numeric(matrix(firm_prob$Bad_News,ncol=1)) 
firm_prob$Neutral_News=as.numeric(matrix(firm_prob$Neutral_News,ncol=1)) 
}
return(firm_prob)
}

# ===============================================================================================================
FontsCE<-function(x){
        FC<-c( "新細明体","細明体","標楷体","黑体","宋体","新宋体","仿宋","楷体","仿宋_GB2312","楷体_GB2312","微軟正黑体","微軟雅黑",
    "隸書","幼圓","華文細黑","華文楷体","華文宋体","華文中宋", "華文仿宋","方正舒体" ,"方正姚体","華文彩云","華文琥珀","華文??", "華文行楷","華文新魏" )
        FE<-c( "PMingLiU" ,"MingLiU","DFKai-SB","SimHei","SimSun" ,"NSimSun","FangSong","KaiTi",
    "FangSong_GB2312","KaiTi_GB2312","Microsoft JhengHei","Microsoft YaHei","LiSu","YouYuan","STXihei","STKaiti",
    "STSong","STZhongsong","STFangsong","FZShuTi","FZYaoti","STCaiyun","STHupo","STLiti","STXingkai","STXinwei")
        CFonts=data.frame(FC,FE,stringsAsFactors=F)
        n<-which(CFonts==x)
            result<-CFonts[n,2]
            return(result)
            }
#FontsCE("仿宋")

roc_curve=function(modelfit,code,main,windows){
if(windows==1){windows()}
a=round(performance(prediction(modelfit, code),'auc')@ y.values[[1]],3)
plot(performance(prediction(modelfit, code),"tpr","fpr"),ylab='sensitivity:true positive rate',
     xlab='(1-specificity):false positive rate',main=main,col=4) 
title(sub=paste("Area under ROC:",round(a,5)),cex.sub = 1.2, font.sub = 3, col.sub = "red")
}
# =================================================================================================================

# ============================================= End function ======================================================
# =================================================================================================================
marc_mark=CombidD(marc_mark,moneyS)
# =========== combind with marc_mark and timemapping =================================
tmp=cbind(nondef[,1:2],matrix(numeric(nrow(nondef)*ncol(marc_mark)),ncol=ncol(marc_mark)))
write.table(table(tmp[,2]),temppath,row.names = F,col.names = F)
yyyymm=read.table(temppath,header=F)
for(i in 1:nrow(yyyymm)){
tmp[tmp[,2]==yyyymm$V1[i],3]=TimeMapping[TimeMapping[,1]==yyyymm$V1[i],2]
tmp[tmp[,2]==yyyymm$V1[i],(4:(ncol(tmp)))]=marc_mark[marc_mark[,1]==yyyymm$V1[i],(2:ncol(marc_mark))]
}
colnames(tmp)=c("firm","yyyymm","TimeMapping",colnames(marc_mark[,2:ncol(marc_mark)]))
nondef_C=cbind(nondef[,1:3],tmp[,3:ncol(tmp)],nondef[,4:ncol(nondef)])
tmp=cbind(def[,1:2],matrix(numeric(nrow(def)*ncol(marc_mark)),ncol=ncol(marc_mark)))
write.table(table(tmp[,2]),temppath,row.names = F,col.names = F)
yyyymm=read.table(temppath,header=F)
for(i in 1:nrow(yyyymm)){
tmp[tmp[,2]==yyyymm$V1[i],3]=TimeMapping[TimeMapping[,1]==yyyymm$V1[i],2]
tmp[tmp[,2]==yyyymm$V1[i],(4:(ncol(tmp)))]=marc_mark[marc_mark[,1]==yyyymm$V1[i],(2:ncol(marc_mark))]
}
colnames(tmp)=c("firm","yyyymm","TimeMapping",colnames(marc_mark[,2:ncol(marc_mark)]))
def_C=cbind(def[,1:3],tmp[,3:ncol(tmp)],def[,4:ncol(def)])
write.table(table(def_C[,1]),temppath,row.names = F,col.names = F)
firm=read.table(temppath,header=F)
def_C$code=0
tmp0={}
for(i in 1:nrow(firm)){
tmp=def_C[def_C$firm==firm$V1[i],]
tmp$code[nrow(tmp)]=1
tmp0=rbind(tmp0,tmp)
}
def_C0=tmp0
alldataspace=rbind(nondef_C,def_C0)
write.table(table(TCRI[,1]),temppath,row.names = F,col.names = F)
firmTCRI=read.table(temppath,header=F)

#======================================= firm summary    ======================================
allfirm=length(table(nondef_NA$firm))+length(table(def_NA$firm)) # All firms
allfirm0=length(table(nondef_NA$firm)) # Non-default firms
allfirm1=length(table(def_NA$firm)) # Default firms

# length(table(alldataspace$firm))  
# k=alldataspace[alldataspace$code==1,]
# length(table(k$firm)) 
# length(table(alldataspace$firm))-length(table(k$firm))

# =================================================================================
# ================================== Step 2 =======================================
# =================================================================================
if(deventcode==1){
alldataspace0=na.omit(alldataspace) ;
alldataspace0=cbind(CombidD1(alldataspace0[,1:2],Devent,chara=0),alldataspace[,4:ncol(alldataspace)]);
}else{ alldataspace0=na.omit(alldataspace) }
alldataspace0$LR_1=1/alldataspace0$LR
alldataspace0$TLTA_1=1/alldataspace0$TLTA
alldataspace0$CLTA_1=1/alldataspace0$CLTA

# ================================================================================
rocplot=alldataspace0[,1:3]

dmyfun10=glm(formula = myfun10 , family = quasibinomial, data = alldataspace0)
dmyfun11=glm(formula = myfun11 , family = quasibinomial, data = alldataspace0)

# Created timemapping0 for forcasting
k=matrix(numeric(period*2),ncol=2)
k0=rep(c(3,6,9,12),period)
yyyy0=round(TimeMapping[nrow(TimeMapping),1]/100)
mm0=TimeMapping[nrow(TimeMapping),1]-yyyy0*100
k00=rep(c(yyyy0)+c(0:period),each = 4)
if(mm0==3){
k[1:nrow(k),1]=k00[2:(2+period-1)]*100+k0[2:(2+period-1)]
}
if(mm0==6){
k[1:nrow(k),1]=k00[3:(3+period-1)]*100+k0[3:(3+period-1)]
}
if(mm0==9){
k[1:nrow(k),1]=k00[4:(4+period-1)]*100+k0[4:(4+period-1)]
}
if(mm0==12){
k[1:nrow(k),1]=k00[5:(5+period-1)]*100+k0[5:(5+period-1)]
}
k[,2]=c((TimeMapping[nrow(TimeMapping),2]+1):(TimeMapping[nrow(TimeMapping),2]+period))
TimeMapping0=data.frame(rbind(as.matrix(TimeMapping),k));
# =============================================================================================
## plot Cumulative intensity for model 
fitvalue={} ; 
fitvalue=cbind(alldataspace0$firm,alldataspace0$yyyymm)
colnames(fitvalue)=c("firm","yyyymm")
fitvalue=data.frame(fitvalue)
fitvalue$Model10=fitted(dmyfun10) # user-define ; without Marc
fitvalue$Model11=fitted(dmyfun11) # user-define ; with Marc
fitvalue=data.frame(fitvalue)
write.table(table(fitvalue$firm),temppath,row.names = F ,col.names = F)
firm=read.table(temppath,header=F)
# =======================================================================================================
# =====================================    EB       =====================================================
# =======================================================================================================
fitvalueINLA={} ; coef_INLA={}
fitvalueINLA=cbind(alldataspace0$firm,alldataspace0$yyyymm)
colnames(fitvalueINLA)=c("firm","yyyymm")
fitvalueINLA=data.frame(fitvalueINLA)
fitvalueINLA$Model10=fitted(dmyfun10) # user-define ; without Marc
fitvalueINLA$Model11=fitted(dmyfun11) # user-define ; with Marc
fitvalueINLA=data.frame(fitvalueINLA)
aaa=inla(myfun10,data=alldataspace0,family="binomial",control.predictor=list(compute=T))
fitvalueINLA$INLA_Model10=aaa$summary.fitted.values[,1]
coef_INLA=rbind(coef_INLA,aaa$summary.fixed) # coef of model10
aaa=inla(myfun11,data=alldataspace0,family="binomial",control.predictor=list(compute=T))
fitvalueINLA$INLA_Model11=aaa$summary.fitted.values[,1]
coef_INLA=rbind(coef_INLA,aaa$summary.fixed) # coef of model11
# Roc for EB

windows()
par(mfrow=c(2,2))
roc_curve(modelfit=fitted(dmyfun10),code=alldataspace0$code,main="Model 1 ",windows=0) # ROC curve for DHM without Market vars.
roc_curve(modelfit=fitted(dmyfun11),code=alldataspace0$code,main="Model 2 ",windows=0) # ROC curve for DHM with Market vars.
roc_curve(modelfit=fitvalueINLA$INLA_Model10,code=alldataspace0$code,main="Model 3 ",windows=0)# ROC curve for INLA without Market vars.
roc_curve(modelfit=fitvalueINLA$INLA_Model11,code=alldataspace0$code,main="Model 4 ",windows=0)# ROC curve for INLA with Market vars.


# =======================================================================================================
# ================================ Auto.arima firm and Marc vars. ============================================
# Note: If the first time use this program, you must use the following two lines of command to forcasting
#       finance data else you can direct read the result from the first.
# ============================================================================================================ 
#Cdata1new=data.frame(forcastVar(firm=firm,Cdata=alldataspace0[,c(1:4,(4+ncol(marc_mark)+1):ncol(alldataspace0))],period=period)) # forcasting : newdata for all firm
#write.table(Cdata1new,"F:/0_Study/6_Theses/code/main_code/Cdata1new.txt",row.names = F)
Cdata1new=read.table("F:/0_Study/6_Theses/code/main_code/Cdata1new.txt",header=T)
Cdata1new$TLTA_1=1/Cdata1new$TLTA;
# forcasting Marc variables
  tmp=as.matrix(marc_mark)
  newvars0=matrix(numeric(period*ncol(marc_mark)),ncol=ncol(marc_mark))
  FstartT=TimeMapping0[TimeMapping0$time==marc_mark$yyyymm[nrow(marc_mark)],2]+1
  FendT=FstartT+period-1
  newvars0[,1]=TimeMapping0[FstartT:FendT,1]
  NewMarc=rbind(tmp,newvars0)
  for(i in 2:ncol(marc_mark)){
    for(j in 1:period){
	  vari=forecast(auto.arima(NewMarc[1:(nrow(tmp)+j-1),i]),h=1)$mean
      NewMarc[(nrow(tmp)+j),i]=vari
       }  
   }
 NewMarc=data.frame(NewMarc)
 Cdata1new_tmp=cbind(Cdata1new,matrix(rep(NA,nrow(Cdata1new)*(ncol(NewMarc)-1)),ncol=(ncol(NewMarc)-1)))
 write.table(table(Cdata1new_tmp$yyyymm),temppath,row.names = F ,col.names = F)
 yyyymm=read.table(temppath,header=F)
for(i in 1:nrow(yyyymm)){
Cdata1new_tmp[Cdata1new_tmp$yyyymm==yyyymm$V1[i],(ncol(Cdata1new)+1):ncol(Cdata1new_tmp)]=NewMarc[NewMarc$yyyymm==yyyymm$V1[i],2:ncol(NewMarc)] 
}
colnames(Cdata1new_tmp)=c(colnames(Cdata1new),colnames(NewMarc[,2:ncol(NewMarc)]))  
Cdata1new=data.frame(Cdata1new_tmp)   
# =====================================================================================
# predict newdata form 201212~201412
allnew10=predict(dmyfun10,Cdata1new,type = "response") ; allnew11=predict(dmyfun11,Cdata1new,type = "response")
allnew=cbind(Cdata1new[,c(1:2)],allnew10,allnew11)
colnames(allnew)=colnames(fitvalue) ; fitvalue=rbind(fitvalue,allnew)

write.table(table(nondef$firm),temppath,row.names = F ,col.names = F)
nondef_firm=read.table(temppath,header=F)

fitvalue0=fitvalue[fitvalue$firm%in%nondef_firm$V1,] # fit value of nondefault firm
fitvalueINLA0=fitvalueINLA[fitvalueINLA$firm%in%nondef_firm$V1,] # fit value of nondefault firm
# ======= plot all intensity by yyyymm
write.table(table(fitvalue$yyyymm),temppath,row.names = F ,col.names = F)
yyyymm=read.table(temppath,header=F)
y=matrix(numeric((ncol(fitvalue)-2)*nrow(yyyymm)),nrow=nrow(yyyymm))
yINLA=matrix(numeric((ncol(fitvalueINLA)-2)*nrow(yyyymm)),nrow=nrow(yyyymm))

for(i in 1:nrow(yyyymm)){
y[i,]=apply(fitvalue[fitvalue$yyyymm==yyyymm$V1[i],3:ncol(fitvalue)],2,mean)
yINLA[i,]=apply(fitvalueINLA[fitvalueINLA$yyyymm==yyyymm$V1[i],3:ncol(fitvalueINLA)],2,mean)
}
yyyymmINLA=cbind(yyyymm$V1,yINLA) ;
yyyymm=cbind(yyyymm$V1,y) ; 
yyyymmINLA=data.frame(yyyymmINLA) ; 
yyyymm=data.frame(yyyymm) ; 
colnames(yyyymm)=c("V1",colnames(fitvalue)[3:ncol(fitvalue)]) ; 
colnames(yyyymmINLA)=c("V1",colnames(fitvalueINLA)[3:ncol(fitvalueINLA)])

yyyy=round(yyyymm$V1[1]/100)
if ((yyyymm$V1[1]-yyyy*100)==3){QQ=1}
if ((yyyymm$V1[1]-yyyy*100)==6){QQ=2}
if ((yyyymm$V1[1]-yyyy*100)==9){QQ=3}
if ((yyyymm$V1[1]-yyyy*100)==12){QQ=4}
starttime=TimeMapping0[TimeMapping0$time==yyyymm$V1[1],2]
endtime=TimeMapping0[TimeMapping0$time==yyyymm$V1[nrow(yyyymm)],2]
yy=data.frame(cbind(TimeMapping0[starttime:endtime,1],matrix(rep(NA,(endtime-starttime+1)*(ncol(yyyymm)-1)),ncol=(ncol(yyyymm)-1))) )
yyINLA=data.frame(cbind(TimeMapping0[starttime:endtime,1],matrix(rep(NA,(endtime-starttime+1)*(ncol(yyyymmINLA)-1)),ncol=(ncol(yyyymmINLA)-1))) )
for(i in 1:nrow(yy)){
if (length(yyyymm[yyyymm[,1]==yy[i,1],2])!=0){yy[i,2:ncol(yy)]=yyyymm[yyyymm[,1]==yy[i,1],2:ncol(yyyymm)]}
if (length(yyyymmINLA[yyyymmINLA[,1]==yy[i,1],2])!=0){yyINLA[i,2:ncol(yyINLA)]=yyyymmINLA[yyyymmINLA[,1]==yy[i,1],2:ncol(yyINLA)]}
}

write.table(table(fitvalue0$yyyymm),temppath,row.names = F ,col.names = F)
yyyymm0=read.table(temppath,header=F)
y0=matrix(numeric((ncol(fitvalue0)-2)*nrow(yyyymm0)),nrow=nrow(yyyymm0))
yINLA0=matrix(numeric((ncol(fitvalueINLA0)-2)*nrow(yyyymm0)),nrow=nrow(yyyymm0))
for(i in 1:nrow(yyyymm0)){
y0[i,]=apply(fitvalue0[fitvalue0$yyyymm==yyyymm0$V1[i],3:ncol(fitvalue0)],2,mean)
yINLA0[i,]=apply(fitvalueINLA0[fitvalueINLA0$yyyymm==yyyymm0$V1[i],3:ncol(fitvalueINLA0)],2,mean)
}
yyyymmINLA0=cbind(yyyymm$V1,yINLA0) ;
yyyymm0=cbind(yyyymm$V1,y0) ; 
yyyymmINLA0=data.frame(yyyymmINLA0) ; 
yyyymm0=data.frame(yyyymm0) ; 
colnames(yyyymm)=c("V1",colnames(fitvalue0)[3:ncol(fitvalue0)]) ; 
colnames(yyyymmINLA0)=c("V1",colnames(fitvalueINLA0)[3:ncol(fitvalueINLA0)])
yy0=data.frame(cbind(TimeMapping0[starttime:endtime,1],matrix(rep(NA,(endtime-starttime+1)*(ncol(yyyymm)-1)),ncol=(ncol(yyyymm)-1))) )
yyINLA0=data.frame(cbind(TimeMapping0[starttime:endtime,1],matrix(rep(NA,(endtime-starttime+1)*(ncol(yyyymmINLA0)-1)),ncol=(ncol(yyyymmINLA0)-1))) )
for(i in 1:nrow(yy0)){
if (length(yyyymm0[yyyymm0[,1]==yy0[i,1],2])!=0){yy0[i,2:ncol(yy0)]=yyyymm0[yyyymm0[,1]==yy0[i,1],2:ncol(yyyymm0)]}
if (length(yyyymmINLA0[yyyymmINLA0[,1]==yy0[i,1],2])!=0){yyINLA0[i,2:ncol(yyINLA0)]=yyyymmINLA0[yyyymmINLA0[,1]==yy0[i,1],2:ncol(yyINLA0)]}
}

jpeg(file=sprintf('%s_withINLA.jpeg',pdfoutput[6]),quality = 100,width = 1000, height = 800,antialias = "cleartype") ;
par(mfrow=c(2,2))
ts.plot(ts(yy[,2],start=c(yyyy,QQ),freq=4),ts(yyINLA[,4],start=c(yyyy,QQ),freq=4),
                                  ylab='Average Intensity',main='Without Marc vars.',lty=c(1,1343),col=c(2,4), lwd=c(2,2))
abline(v = c(2008.75,2012), col = c(3,2), lty=2 )
legend("topright",c("Model 1",'Model 3',"Predict value"),lty=c(1,1343,1),col=c(2,4,6),lwd=c(2,2,2),cex =0.9);
x0=c(2013,2013.25,2013.5,2013.75,2014,2014.25,2014.5,2014.75)
y0=yy[(nrow(yy)-period+1):nrow(yy),2]
matlines(x=x0 ,y=y0,col=6,lty=1,lwd=2)
# ===========================================
ts.plot(ts(yy[,3],start=c(yyyy,QQ),freq=4),ts(yyINLA[,5],start=c(yyyy,QQ),freq=4),
                                  ylab='Average Intensity',main='With Marc vars.',lty=c(1,1343),col=c(2,4), lwd=c(2,2))
abline(v = c(2008.75,2012), col = c(3,2), lty=2 )
legend("topright",c("Model 2",'Model 4',"Predict value"),lty=c(1,1343,1),col=c(2,4,6),lwd=c(2,2,2),cex =0.9);
x0=c(2013,2013.25,2013.5,2013.75,2014,2014.25,2014.5,2014.75)
y0=yy[(nrow(yy)-period+1):nrow(yy),3]
matlines(x=x0 ,y=y0,col=6,lty=1,lwd=2)
dev.off();

# ====== For each firm  incouding default and nondefault firm
firm_mapping=cbind(firm1_keyword,firm2_keyword)
colnames(firm_mapping)=c('firm','firm_Big5')

CairoPDF(file=pdfoutput[3],width = 11, height = 11) ; 
for(i in 1:nrow(firm)){
p=nrow(Cdata1new[Cdata1new$firm==firm$V1[i],]) ; 
if (p!=0){
newdat=data.frame(Cdata1new[Cdata1new$firm==firm$V1[i],])
firmnew10=predict(dmyfun10,newdat,type = "response")
firmnew11=predict(dmyfun11,newdat,type = "response")
firmnew=cbind(Cdata1new[Cdata1new$firm==firm$V1[i],c(1:2)],cbind(cbind(firmnew10,firmnew11)))
tmp=fitvalue[fitvalue$firm==firm$V1[i],] ; 
tmp0=fitvalueINLA[fitvalueINLA$firm==firm$V1[i],c(1,2,5,6)] ;
tmp[,3]=cumsum(tmp[,3]) ;tmp[,4]=cumsum(tmp[,4]) ;
tmp0[,3]=cumsum(tmp0[,3]) ;tmp0[,4]=cumsum(tmp0[,4])
x0=round(c(tmp$yyyymm)/100)+(c(tmp$yyyymm)-round(c(tmp$yyyymm)/100)*100)/12
matplot(x=x0,tmp[,c(3,4)],ylab='Cumulative Intensity',type='l',xlab='Year',col=c(2,4),lty=1,lwd=1,
        main=sprintf('For firm %s(%i)',as.character(firm_mapping[firm_mapping[,1]==firm$V1[i],2]),firm$V1[i]),family='SimSun')
matlines(x=x0[1:nrow(tmp0)],tmp0[,c(3,4)],col=c(2,4),lty=2,lwd=1)
legend("topleft",c("with Marc vars.","without Marc vars.","with Marc vars.(INLA)",
       "without Marc vars.(INLA)"),lwd=1,lty=c(1,1,2,2),col=c(2,4,2,4));
y0=tmp[(nrow(tmp)-period):nrow(tmp),c(3,4)]
xtmp=(c(tmp$yyyymm[(nrow(tmp)-period):nrow(tmp)])-round(c(tmp$yyyymm[(nrow(tmp)-period):nrow(tmp)])/100)*100)/12
x0=round(c(tmp$yyyymm[(nrow(tmp)-period):nrow(tmp)])/100)+xtmp
matlines(x=x0 ,y=y0,col=6,lty=1,lwd =1)
  }	  
}
 dev.off();
 
# ================================ RMSE : model performance ============================================= 
insampleDate=200812      
LRName10=c('Intercept','I(log(TimeMapping))','TLTA_1','LiquidR','RE2TA','EPS','OCR')
LRName11=c('Intercept','I(log(TimeMapping))','TLTA_1','LiquidR','RE2TA','EPS','OCR','unemployment','stock')

l=length(LRName10) ; h=length(LRName11)
Cdata1= alldataspace0   # No NA data space

write.table(table(Cdata1$yyyymm),temppath,row.names = F ,col.names = F)
yyyymm=read.table(temppath,header=F)
k=sum((yyyymm$V1==insampleDate)*c(1:nrow(yyyymm)))
tmpyyyymm=yyyymm$V1[k:(nrow(yyyymm)-1)]
RMSE1=numeric(length(tmpyyyymm)) ; RMSE2=numeric(length(tmpyyyymm)) # RMSE1 : model without marc ; RMSE2 : model with marc
logisticP10=matrix(numeric(l*2*length(tmpyyyymm)),nrow=l) ;logisticP11=matrix(numeric(h*2*length(tmpyyyymm)),nrow=h)

for(i in 1:length(tmpyyyymm)){
insampleDate0=tmpyyyymm[i]
insample1=Cdata1[Cdata1[,2]<=insampleDate0,]; 
outsample1=Cdata1[Cdata1[,2]>insampleDate0,]; 
in_fit10=glm(myfun10 ,data=insample1,family= quasibinomial)
in_fit11=glm(myfun11 ,data=insample1,family= quasibinomial)
logisticP10[,((i-1)*2+1):(2*i)]=summary(in_fit10)$coeff[,c(1,4)];
logisticP11[,((i-1)*2+1):(2*i)]=summary(in_fit11)$coeff[,c(1,4)];
outsample1$predictValue10 = predict(in_fit10,newdata=outsample1,type = "response")
outsample1$predictValue11 = predict(in_fit11,newdata=outsample1,type = "response")
write.table(table(outsample1$firm),temppath,row.names = F ,col.names = F)
n=read.table(temppath,header=F)
tmp={}
for(j in 1:nrow(n)){
tmp0=fitvalue[fitvalue$firm==n$V1[j],] 
tmp1=outsample1[outsample1$firm==n$V1[j],]
tmp1$fitValue10=tmp0[tmp0$yyyymm%in%tmp1$yyyymm,3] ;tmp1$fitValue11=tmp0[tmp0$yyyymm%in%tmp1$yyyymm,4]
tmp=rbind(tmp,tmp1)
}
outsample1=tmp
outsample1$diff10=outsample1$fitValue10-outsample1$predictValue10
outsample1$diff11=outsample1$fitValue11-outsample1$predictValue11
RMSE1[i]=(median(outsample1$diff10^2))^0.5 ; RMSE2[i]=(median(outsample1$diff11^2))^0.5
}
colnames(logisticP10)=rep(tmpyyyymm,each=2)
rownames(logisticP10)=LRName10
colnames(logisticP11)=rep(tmpyyyymm,each=2)
rownames(logisticP11)=LRName11

pdf(file=pdfoutput[5]) 
T1=ts(RMSE1,start=c(round(tmpyyyymm[1]/100),(tmpyyyymm[1]-round(tmpyyyymm[1]/100)*100)/3),frequency = 4)
T2=ts(RMSE2,start=c(round(tmpyyyymm[1]/100),(tmpyyyymm[1]-round(tmpyyyymm[1]/100)*100)/3),frequency = 4)
ts.plot(T1,T2,main=sprintf("RMSE from %d to %d",insampleDate,tmpyyyymm[length(tmpyyyymm)])
		,ylab='Value',lty=c(1,1343),col=c(2,4), lwd=c(2,2))
legend("topright",c("Model 1",'Model 2'),lty=c(1,1343),col=c(2,4),lwd=c(2,2));
dev.off();

# ========== Summary output =================================
sink(file=outpath[4]);
      print("Dataspace summary all Data");
      print( var_summary(alldataspace0[,4:(ncol(alldataspace0)-5)]) );
	  print("===============================================================================================================================================================================================================================")
	  print("Dataspace summary without na.omit ");
      print( var_summary(Cdata1[,4:(ncol(Cdata1)-5)]) );
	  print("===============================================================================================================================================================================================================================")
	  print("DHR without marc ");
	  summary(dmyfun10)
	  print("===============================================================================================================================================================================================================================")
	  print("DHR with marc ");
	  summary(dmyfun11)
	  print("===============================================================================================================================================================================================================================")
	  sprintf("Parameter estimation for Logistic Model from %d to %d(Without Marc vars.)",insampleDate,tmpyyyymm[length(tmpyyyymm)]);
	  logisticP10
	  print("===============================================================================================================================================================================================================================")
	  sprintf("Parameter estimation for Logistic Model from %d to %d",insampleDate,tmpyyyymm[length(tmpyyyymm)]);
	  logisticP11
	  print("=====================================================================================================================")
	  sprintf('Parameter estimation for Logistic Model' );
	  coef_INLA[1:length(LRName10),]
	  print("=====================================================================================================================")
	  coef_INLA[(length(LRName10)+1):(length(LRName10)+length(LRName11)),]
	  sink(type="message");
sink(); 

# ======================================================================================================================
# ===================================== News infromation ===============================================================
# =======================================================================================================
# Newsdata : data space with three part. The first is date (yyyymmdd). The second is keyword. 
#            The three is times of 'Good','Bad' and 'Neutral' in keyword.txt from keyword.r program.
# word : The word of want to calculate the prob. which must contain Big5 column "matrix".
# Ex. for marc word : Inf=Information(newsdata=news_tmpQ1,word=matrix(Marc_keyword,ncol=1)) 
# Ex. for firm word : firmname=cbind(firm1_keyword,firm2_keyword)
#                     Inf=Information(newsdata=news_tmpQ1,word=firmname)
# =======================================================================================================

# without daily weight
# keywordpath[i] : i = 2008,2009,2010,...,2012
month=c(3,6,9,12)

## ======================  Information for firm   ============================
firm_mapping$firm_Big5=as.character(firm_mapping[,2])
firm_Inf={} ;
for(k in 1:length(keywordpath)){
Hpyeartmp={}
keyword=read.table(keywordpath[k],header=TRUE) ; 
news_tmpQ1=keyword[keyword$yyyymmdd>(newsyear[k]-1)*10000+month[4]*100 & keyword$yyyymmdd<newsyear[k]*10000+month[1]*100+32,]
news_tmpQ2=keyword[keyword$yyyymmdd>newsyear[k]*10000+month[1]*100+32 & keyword$yyyymmdd<newsyear[k]*10000+month[2]*100+32,]
news_tmpQ3=keyword[keyword$yyyymmdd>newsyear[k]*10000+month[2]*100+32 & keyword$yyyymmdd<newsyear[k]*10000+month[3]*100+32,]
news_tmpQ4=keyword[keyword$yyyymmdd>newsyear[k]*10000+month[3]*100+32 & keyword$yyyymmdd<newsyear[k]*10000+month[4]*100+32,]
HpQ1=Information(news=news_tmpQ1,word=firm_mapping) 
HpQ2=Information(news=news_tmpQ2,word=firm_mapping)
HpQ3=Information(news=news_tmpQ3,word=firm_mapping)
HpQ4=Information(news=news_tmpQ4,word=firm_mapping)
Hpyeartmp=rbind(HpQ1,HpQ2,HpQ3,HpQ4)
firm_Inf=rbind(firm_Inf,Hpyeartmp)
 }

write.table(table(firm_Inf$yyyymm),temppath,row.names = F ,col.names = F)
yyyymm=read.table(temppath,header=F)
tmp={}
for(i in 1:nrow(firm_mapping)){
tmp0=firm_Inf[firm_Inf[,1]%in%firm_mapping[i,1],]
tmp=rbind(tmp,tmp0)
}
 firm_Inf=tmp

## ======================  Information for marc   ============================
marc=c(1:53)+2 ; # GBNkeyword33

keyword=read.table(keywordpath[1],header=TRUE)
marc_mapping=matrix(colnames(keyword)[marc],ncol=1) # The names of marc. variables had changed.
marc_Inf={} ; 
for(k in 1:length(keywordpath)){
Hpyeartmp={}
keyword=read.table(keywordpath[k],header=TRUE) ; 
news_tmpQ1=keyword[keyword$yyyymmdd>(newsyear[k]-1)*10000+month[4]*100 & keyword$yyyymmdd<newsyear[k]*10000+month[1]*100+32,]
news_tmpQ2=keyword[keyword$yyyymmdd>newsyear[k]*10000+month[1]*100+32 & keyword$yyyymmdd<newsyear[k]*10000+month[2]*100+32,]
news_tmpQ3=keyword[keyword$yyyymmdd>newsyear[k]*10000+month[2]*100+32 & keyword$yyyymmdd<newsyear[k]*10000+month[3]*100+32,]
news_tmpQ4=keyword[keyword$yyyymmdd>newsyear[k]*10000+month[3]*100+32 & keyword$yyyymmdd<newsyear[k]*10000+month[4]*100+32,]
HpQ1=Information(news=news_tmpQ1,word=marc_mapping) 
HpQ2=Information(news=news_tmpQ2,word=marc_mapping)
HpQ3=Information(news=news_tmpQ3,word=marc_mapping)
HpQ4=Information(news=news_tmpQ4,word=marc_mapping)
Hpyeartmp=rbind(HpQ1,HpQ2,HpQ3,HpQ4)
marc_Inf=rbind(marc_Inf,Hpyeartmp)
 }
 
write.table(table(marc_Inf$yyyymm),temppath,row.names = F ,col.names = F)
yyyymm=read.table(temppath,header=F)
tmp={}
for(i in 1:nrow(marc_mapping)){
tmp0=marc_Inf[marc_Inf[,1]%in%marc_mapping[i,1],]
tmp=rbind(tmp,tmp0)
}
marc_Inf=tmp

# ============= firm effect from News ======================= 
Cdata2008=Cdata1[Cdata1$yyyymm>200712,]; 
Cdata20081=CombidD1(Cdata2008,firm_Inf,chara=0);
Cdata2008=Cdata20081
Cdata2008$Good_NewsR=Cdata2008$Good_News/Cdata2008$All_News ;#Cdata2008$Good_NewsR_1=Cdata2008$Good_News_1/Cdata2008$All_News_1 
Cdata2008$Bad_NewsR=Cdata2008$Bad_News/Cdata2008$All_News ; #Cdata2008$Bad_NewsR_1=Cdata2008$Bad_News_1/Cdata2008$All_News_1
Cdata2008[Cdata2008$All_News==0,(ncol(Cdata2008)-1):ncol(Cdata2008)]=0
mm1=(glm(formula = code ~ I(log(TimeMapping))+ Good_NewsR + Bad_NewsR  ,family = quasibinomial, data = Cdata2008))
mm2=(glm(formula = code ~I(log(TimeMapping)) + All_News ,family = quasibinomial, data = Cdata2008))

# ============= firm & marc effect from News =======================
marc_fitVar=c(1:nrow(marc_mapping)) 
coeftmp1={} ; coeftmp2={} ; marc_fit0={} ;  
for(j in 1:nrow(marc_mapping)) {
Cdata20082=Cdata20081
marctemp0=data.frame(marc_Inf[marc_Inf$Big5==marc_mapping[j,1],2:ncol(marc_Inf)])
marctemp=data.frame(matrix(as.numeric(matrix(marctemp0[,1],ncol=1)),ncol=1))
marctemp$Good_News=marctemp0$Good_News/matrix(as.numeric(matrix(marctemp0$All_News,ncol=1)),ncol=1)
marctemp$Bad_News=marctemp0$Bad_News/matrix(as.numeric(matrix(marctemp0$All_News,ncol=1)),ncol=1)
marctemp[marctemp0$All_News==0,2:3]=0 ;
colnames(marctemp)=c('yyyymm','Good_NewsR','Bad_NewsR')
Cdata20082=data.frame(cbind(Cdata20082,matrix(rep(0,nrow(Cdata20082)*(ncol(marctemp)-1)),ncol=(ncol(marctemp)-1)) )  )
for(i in 1:nrow(yyyymm)){
Cdata20082[Cdata20082[,2]==yyyymm[i,1],(ncol(Cdata20082)-ncol(marctemp)+2):ncol(Cdata20082)] = 
matrix(rep(as.numeric(marctemp[marctemp$yyyymm==yyyymm[i,1],2:ncol(marctemp)]),nrow(Cdata20082[Cdata20082[,2]==yyyymm[i,1],(ncol(Cdata20082)-ncol(marctemp)+2):ncol(Cdata20082)])),ncol=2,byrow = T)
}
m10=glm(formula = code ~ I(log(TimeMapping)) + X1 + X2 , family = quasibinomial, data = Cdata20082) ; m1=summary(m10)
m20=glm(formula = code ~ X1+ X2  , family = quasibinomial, data = Cdata20082)  ; m2=summary(m20)
coeftmp1=rbind(coeftmp1,cbind(j,round(m1$coeff[,c(1,4)],4)) )
coeftmp2=rbind(coeftmp2,cbind(j,round(m2$coeff[,c(1,4)],4)) )
  if(sum(marc_fitVar%in%j)!=0){
    marc_fit0=cbind(marc_fit0,matrix(fitted(m10),ncol=1)) 
  }
 }
colnames(marc_fit0) = marc_mapping[marc_fitVar,1]
write.table(as.matrix(round(apply(marc_fit0,2,sd),5)),temppath,row.names = T ,col.names = F)
marc_fitSD=read.table(temppath,header=F)
marc_firm_fit=cbind(marc_fit0[,order(marc_fitSD$V2,decreasing = T)[1:5]],fitted(mm1))  #(cbinde firm and MARC keyword )large -> small sd
marc_firm_fit1=cbind(marc_fit0[,order(marc_fitSD$V2,decreasing = F)[1:5]],fitted(mm1)) #(cbinde firm and MARC keyword ) small-> large sd
# ==========================================================================================================================
tmp=t(apply(marc_firm_fit,1,cumsum)) ; tmp1=t(apply(marc_firm_fit1,1,cumsum))
h=alldataspace0[alldataspace0$yyyymm>200712,]
aaa=inla(myfun10,data=h,family="binomial",offset=tmp[,ncol(tmp)],control.predictor=list(compute=T))
coef_INLA_news10=aaa$summary.fixed # coef of model11
firm_marc_fitPD=cbind(Cdata2008[,1:2],aaa$summary.fitted.values[,1] )  ;
aaa=inla(myfun11,data=h,family="binomial",offset=tmp1[,ncol(tmp1)],control.predictor=list(compute=T),control.compute=list(dic=TRUE, cpo=TRUE))
coef_INLA_news11=aaa$summary.fixed # coef of model11
firm_marc_fitPD1=cbind(Cdata2008[,1:2],aaa$summary.fitted.values[,1] )
colnames(firm_marc_fitPD)=c('firm','yyyymm','fitPD') ; colnames(firm_marc_fitPD1)=c('firm','yyyymm','fitPD')
firm_marc_fitPD=data.frame(firm_marc_fitPD ) ; firm_marc_fitPD1=data.frame(firm_marc_fitPD1)

write.table(table(firm_marc_fitPD$yyyymm),temppath,row.names = F ,col.names = F)
yyyymm=read.table(temppath,header=F)
for(i in 1:nrow(yyyymm)){
yyyymm$V2[i]=mean(firm_marc_fitPD[firm_marc_fitPD$yyyymm==yyyymm$V1[i],ncol(firm_marc_fitPD)])
yyyymm$V3[i]=mean(firm_marc_fitPD1[firm_marc_fitPD1$yyyymm==yyyymm$V1[i],ncol(firm_marc_fitPD1)])
}
yyyymm=data.frame(CombidD(yyyymmINLA,yyyymm))
yyyymm=yyyymm[yyyymm$V1<=201212,] ; yyyymm[is.na(yyyymm$V2),c((ncol(yyyymm)-1):ncol(yyyymm))]=0
k=yyyymm[yyyymm$V2==0,]$Model10 ; k1=yyyymm[yyyymm$V2==0,]$Model11
Model10PD=c(k,yyyymm[yyyymm$V2!=0,]$V2) ; Model10PD1=c(k,yyyymm[yyyymm$V3!=0,]$V3)
Model11PD=c(k1,yyyymm[yyyymm$V2!=0,]$V2) ; Model11PD1=c(k1,yyyymm[yyyymm$V3!=0,]$V3)
yyyymm=cbind(yyyymm$V1,yyyymm$Model10,yyyymm$Model11,Model10PD,Model10PD1,Model11PD,Model11PD1)
colnames(yyyymm)=c('V1','Model10','Model11','Model10PD','Model10PD1','Model11PD','Model11PD1') ; yyyymm=data.frame(yyyymm)
yyyy=round(yyyymm$V1[1]/100)
if ((yyyymm$V1[1]-yyyy*100)==3){QQ=1}
if ((yyyymm$V1[1]-yyyy*100)==6){QQ=2}
if ((yyyymm$V1[1]-yyyy*100)==9){QQ=3}
if ((yyyymm$V1[1]-yyyy*100)==12){QQ=4}

starttime=TimeMapping0[TimeMapping0$time==yyyymm$V1[1],2]
endtime=TimeMapping0[TimeMapping0$time==yyyymm$V1[nrow(yyyymm)],2]
yy=cbind(TimeMapping0[starttime:endtime,1],matrix(rep(NA,(endtime-starttime+1)*(ncol(yyyymm)-1)),ncol=(ncol(yyyymm)-1)) )
for(j in 2:ncol(yyyymm)){
for(i in 1:nrow(yy)){
if (length(yyyymm[yyyymm$V1==yy[i,1],j])!=0){yy[i,j]=yyyymm[yyyymm$V1==yy[i,1],j]}
 }
}
colnames(yy)=colnames(yyyymm)

jpeg(file=sprintf('%s_withNews1.jpeg',pdfoutput[6]),quality = 100,width = 1000, height = 800,antialias = "cleartype") ;
par(mfcol=c(2,2))
ts.plot(ts(yy[,2],start=c(yyyy,QQ),freq=4),ts(yy[,4],start=c(yyyy,QQ),freq=4),
                                  ylab='Average Intensity',main='All Firm without Marc Vars.(L)',lty=c(1,1343),col=c(2,4), lwd=c(2,2))
abline(v = c(2008.75,2012),h=0, col = c(2,3), lty=2 )
legend("topright",c("Model 1",'Model 5'),lty=c(1,1343),col=c(2,4),lwd=c(2,2),cex =1.2);
ts.plot(ts(yy[,2],start=c(yyyy,QQ),freq=4),ts(yy[,5],start=c(yyyy,QQ),freq=4),
                                  ylab='Average Intensity',main='All Firm without Marc Vars.(S)',lty=c(1,1343),col=c(2,4), lwd=c(2,2))
abline(v = c(2008.75,2012),h=0, col = c(2,3), lty=2 )
legend("topright",c("Model 1",'Model 5'),lty=c(1,1343),col=c(2,4),lwd=c(2,2),cex =1.2);
ts.plot(ts(yy[,3],start=c(yyyy,QQ),freq=4),ts(yy[,6],start=c(yyyy,QQ),freq=4),
                                  ylab='Average Intensity',main='All Firm with Marc Vars.(L)',lty=c(1,1343),col=c(2,4), lwd=c(2,2))
abline(v = c(2008.75,2012),h=0, col = c(2,3), lty=2 )
legend("topright",c("Model 2",'Model 6'),lty=c(1,1343),col=c(2,4),lwd=c(2,2),cex =1.2);

ts.plot(ts(yy[,3],start=c(yyyy,QQ),freq=4),ts(yy[,7],start=c(yyyy,QQ),freq=4),
                                  ylab='Average Intensity',main='All Firm with Marc Vars.(S)',lty=c(1,1343),col=c(2,4), lwd=c(2,2))
abline(v = c(2008.75,2012),h=0 ,col = c(2,3), lty=2 )
legend("topright",c("Model 2",'Model 6'),lty=c(1,1343),col=c(2,4),lwd=c(2,2),cex =1.2);
dev.off();

# ====== For each firm  incouding default and nondefault firm
# ====== default intensity with News =========================
# PD : large -> small sd for marc keyword (L) ; PD1 : small -> large sd for marc keyword (S)
marc_fitPD=data.frame(cbind(firm_marc_fitPD,firm_marc_fitPD1[,3])) ; colnames(marc_fitPD)=c('firm','yyyymm','fitPD','fitPD1')
fitvalue_news=CombidD1(fitvalue,marc_fitPD,chara=0) ; fitvalue_news=fitvalue_news[fitvalue_news$yyyymm<201303,]
fitvalue_news$M10_PD=0 ; fitvalue_news$M10_PD1=0 ; fitvalue_news$M11_PD=0 ; fitvalue_news$M11_PD1=0 ;
fitvalue_news[fitvalue_news$yyyymm<200803,]$M10_PD= fitvalue_news[fitvalue_news$yyyymm<200803,]$Model10; 
fitvalue_news[fitvalue_news$yyyymm>=200803,]$M10_PD= fitvalue_news[fitvalue_news$yyyymm>=200803,]$fitPD; 
fitvalue_news[fitvalue_news$yyyymm<200803,]$M10_PD1= fitvalue_news[fitvalue_news$yyyymm<200803,]$Model10
fitvalue_news[fitvalue_news$yyyymm>=200803,]$M10_PD1= fitvalue_news[fitvalue_news$yyyymm>=200803,]$fitPD1;
fitvalue_news[fitvalue_news$yyyymm<200803,]$M11_PD= fitvalue_news[fitvalue_news$yyyymm<200803,]$Model11
fitvalue_news[fitvalue_news$yyyymm>=200803,]$M11_PD= fitvalue_news[fitvalue_news$yyyymm>=200803,]$fitPD;
fitvalue_news[fitvalue_news$yyyymm<200803,]$M11_PD1= fitvalue_news[fitvalue_news$yyyymm<200803,]$Model11
fitvalue_news[fitvalue_news$yyyymm>=200803,]$M11_PD1= fitvalue_news[fitvalue_news$yyyymm>=200803,]$fitPD1;

CairoPDF(file=sprintf('%s_withNews1.pdf',pdfoutput[3]),width = 11, height = 11) ; 
for(k in 1:nrow(firm)){
p=nrow(fitvalue_news[fitvalue_news$firm==firm$V1[k],]) ; 
 if (p>=4){
  tmp=fitvalue_news[fitvalue_news$firm==firm$V1[k],2:ncol(fitvalue_news)] 
  starttime=TimeMapping0[TimeMapping0$time==tmp$yyyymm[1],2]
  endtime=TimeMapping0[TimeMapping0$time==tmp$yyyymm[nrow(tmp)],2]
  yy=cbind(TimeMapping0[starttime:endtime,1],matrix(rep(NA,(endtime-starttime+1)*(ncol(tmp)-1)),ncol=(ncol(tmp)-1)) )
 for(j in 2:ncol(tmp)){
  for(i in 1:nrow(yy)){
   if (length(tmp[tmp$yyyymm==yy[i,1],j])!=0){yy[i,j]=tmp[tmp$yyyymm==yy[i,1],j]}
   }
  }              
 for(j in 2:ncol(tmp)){
  for(i in 1:nrow(yy)){
   if (is.na(yy[i,j])==1){
     if (yy[i,1]==200012){
	    yy[i,j]=0
	   } else {
	    ii=0 
        for(ij in 1:nrow(yy)) {
		 ttmp=yy[(i+ij),j]
		 if (is.na(ttmp)==1){ii=ii+ij
		   }else{break} 
		 }
	   }
     yy[i,j]=mean(yy[(i-1),j],ttmp)
	 }
   }
 }
 tmp=data.frame(cbind(yy[,1],apply(yy[,2:ncol(yy)],2,cumsum)))
colnames(tmp)=colnames(fitvalue_news)[2:ncol(fitvalue_news)]
 x0=round(c(tmp$yyyymm)/100)+(c(tmp$yyyymm)-round(c(tmp$yyyymm)/100)*100)/12
par(mfcol=c(2,1))
matplot(x=x0,tmp[,c(2,6,7)],ylab='Cumulative Intensity',type='l',xlab='Year',col=c(2:4),lty=c(1:3),lwd=1,
        main=sprintf('For firm %s(%i)',as.character(firm_mapping[firm_mapping[,1]==firm$V1[k],2]),firm$V1[k]),family='KaiTi_GB2312')
abline(v = c(2008.75), col = c(2), lty=1 )
legend("topleft",c("without Marc vars.","with News (L).","with News (S)"),col=c(2:4),lty=c(1:3),lwd=1);
matplot(x=x0,tmp[,c(3,8,9)],ylab='Cumulative Intensity',type='l',xlab='Year',col=c(2:4),lty=c(1:3),lwd=1,
        main=sprintf('For firm %s(%i)',as.character(firm_mapping[firm_mapping[,1]==firm$V1[k],2]),firm$V1[k]),family='KaiTi_GB2312')
abline(v = c(2008.75), col = c(2), lty=1 )
legend("topleft",c("with Marc vars.","with News (L).","with News (S)"),col=c(2:4),lty=c(1:3),lwd=1);
  }	  
}
 dev.off();

# ==========News Summary output =================================
sink(file=outpath[7]);
      print("firm effect from News");
      print( summary(mm1) );
	  print("firm effect from News(All News)");
      print( summary(mm2) );
	  print("=====================================================")
	  print("news and model 1 ");
	  coef_INLA_news10
	  print("=====================================================")
	  print("News amd model 2 ");
	  coef_INLA_news11
	  print("=====================================================")
	  print("firm & marc effect from News ");
      print( coeftmp1 );
	  print("=====================================================")
	  sprintf("firm & marc effect from News without TimeMapping ");
	  coeftmp2
sink(type="message");
sink(); 

jpeg(file='F:/0_Study/6_Theses/code/output/TCRI_rating.jpeg' ,quality = 100,width = 700, height = 700,antialias = "cleartype") ;
par(mfcol=c(3,3))
for(j in 1:9){
yyyyNdef=yyyyNdef0[yyyyNdef0$TCRI==j,]
write.table(table(yyyyNdef$yyyy),temppath,row.names = F ,col.names = F)
yyyyNdef00=read.table(temppath,header=F)
for(i in 1:nrow(yyyyNdef00)){
yyyyNdef00$V2[i]=sum(yyyyNdef[yyyyNdef$yyyy==yyyyNdef00$V1[i],3])/sum(yyyyNdef[yyyyNdef$yyyy==yyyyNdef00$V1[i],2])
}
ts.plot(ts(yyyyNdef00$V2,start=1999),main=sprintf( 'Rating of %d',j),ylab='default Rate ')
abline(v = c(2005,2008), col = c(3,2), lty=2 )
legend("topright",c("2005",'2008'),lty=2,col=c(3,2),lwd=1,cex =1);
}
 dev.off();

# ============================================================================================= 

fitvalue_news0=CombidD1(alldataspace0[,1:3],fitvalue_news,chara=0)
windows()
par(mfcol=c(2,2))
roc_curve(modelfit=fitvalue_news0$M10_PD,code=fitvalue_news0$code,main="Model 5 (L)",windows=0) # Without Marc Vars and (L)
roc_curve(modelfit=fitvalue_news0$M11_PD,code=fitvalue_news0$code,main="Model 6 (L)",windows=0) # With Marc Vars and (L)
roc_curve(modelfit=fitvalue_news0$M10_PD1,code=fitvalue_news0$code,main="Model 5 (S)",windows=0)# Without Marc Vars and (S)
roc_curve(modelfit=fitvalue_news0$M11_PD1,code=fitvalue_news0$code,main="Model 6 (S)",windows=0)# With Marc Vars and (S)

fitvalue_news2008=fitvalue_news0[fitvalue_news0$yyyymm>=200803,]
windows()
par(mfrow=c(3,2))
roc_curve(modelfit=fitvalue_news2008$Model10,code=fitvalue_news2008$code,main='Without Marc Vars(2008Q1~2012Q4)',windows=0)
roc_curve(modelfit=fitvalue_news2008$Model11,code=fitvalue_news2008$code,main='With Marc Vars(2008Q1~2012Q4)',windows=0)
roc_curve(modelfit=fitvalue_news2008$M10_PD,code=fitvalue_news2008$code,main='Without Marc Vars and (L)(2008Q1~2012Q4)',windows=0)
roc_curve(modelfit=fitvalue_news2008$M11_PD,code=fitvalue_news2008$code,main='With Marc Vars and (L)(2008Q1~2012Q4)',windows=0)
roc_curve(modelfit=fitvalue_news2008$M10_PD1,code=fitvalue_news2008$code,main='Without Marc Vars and (S)(2008Q1~2012Q4)',windows=0)
roc_curve(modelfit=fitvalue_news2008$M11_PD1,code=fitvalue_news2008$code,main='With Marc Vars and (S)(2008Q1~2012Q4)',windows=0)


 # =========================================== The TS plot of default  ================================================
def_time=read.table('F:/0_Study/6_Theses/code/main_code/firmdata/20130916/deffirmtime.txt',header=TRUE) # Time to delist by firm
def_time$yyyy=round(def_time$yyyymm/100) ; def_time$mm=def_time$yyyymm-(def_time$yyyy)*100
def_time$QQ=0 ; 
def_time[def_time$mm>=1 & def_time$mm<=3,5]=3
def_time[def_time$mm>=4 & def_time$mm<=6,5]=6
def_time[def_time$mm>=7 & def_time$mm<=9,5]=9
def_time[def_time$mm>=10 & def_time$mm<=12,5]=12
def_time$yyyyQQ=def_time$yyyy*100+def_time$QQ ;  def_time=cbind(def_time$firm,def_time$yyyyQQ,def_time$yyyy) 
colnames(def_time)=c('firm','deftime','defyear') ; def_time=data.frame(def_time)
def_time=def_time[def_time$deftime>=200003,] ; def_time=def_time[def_time$deftime<=201212,]

write.table(table(def_time$defyear),temppath,row.names = F ,col.names = F)
defyear=read.table(temppath,header=F)
allfirmT=rbind(nondef[,1:2],def[,1:2]) ; allfirmT$year=round(allfirmT$yyyymm/100)

# ============  Delisting   ==================
write.table(table(allfirmT$year),temppath,row.names = F ,col.names = F)
allfirmTY=read.table(temppath,header=F)
defyear$rate=(defyear$V2/allfirmTY$V2)*100
windows()
par(mfrow=c(2,2))
ts.plot(ts(defyear$rate,start=(2000)),ylab='Default rate (%)',main='Delisting (Year)')

write.table(table(def_time$deftime),temppath,row.names = F ,col.names = F)
defQ=read.table(temppath,header=F)
write.table(table(allfirmT$yyyymm),temppath,row.names = F ,col.names = F)
allfirmTQ=read.table(temppath,header=F)
defQ=CombidD0(allfirmTQ,defQ,chara=0) ; colnames(defQ)=c('yyyymm','allfirm','deffirm') ; defQ=data.frame(defQ)
defQ$rate=(defQ$deffirm/defQ$allfirm)*100
windows()
ts.plot(ts(defQ$rate,start=c(2000,1),freq=4),ylab='Default rate (%)',main='Delisting (Quarterly)')

# ============== TCRI default =================
Devent$yyyy=round(Devent$yyyymm/100)
write.table(table(Devent$yyyy),temppath,row.names = F ,col.names = F)
TCRIdefY=read.table(temppath,header=F)
TCRIdefY=CombidD0(allfirmTY,TCRIdefY,chara=0) ; colnames(TCRIdefY)=c('yyyy','allfirm','deffirm') ; TCRIdefY=data.frame(TCRIdefY)
TCRIdefY$rate=(TCRIdefY$deffirm/TCRIdefY$allfirm)*100
windows()
ts.plot(ts(TCRIdefY$rate,start=c(2000)),ylab='Default rate (%)',main='TCRI (Year)')

write.table(table(Devent$yyyymm),temppath,row.names = F ,col.names = F)
TCRIdefQ=read.table(temppath,header=F)
write.table(table(allfirmT$yyyymm),temppath,row.names = F ,col.names = F)
allfirmTQ=read.table(temppath,header=F)
TCRIdefQ=CombidD0(allfirmTQ,TCRIdefQ,chara=0) ; colnames(TCRIdefQ)=c('yyyymm','allfirm','deffirm') ; TCRIdefQ=data.frame(TCRIdefQ)
TCRIdefQ$rate=(TCRIdefQ$deffirm/TCRIdefQ$allfirm)*100
windows()
ts.plot(ts(TCRIdefQ$rate,start=c(2000,1),freq=4),ylab='Default rate (%)',main='TCRI (Quarterly)')

# ==================================================================================================
fitvalue_news$code=0 ; 
for(i in 1:nrow(def_time)){
if(sum(fitvalue_news$firm%in%def_time$firm[i])!=0){
fitvalue_news[fitvalue_news$firm==def_time$firm[i] & fitvalue_news$yyyymm==def_time$deftime[i],11]=1
 } 
}
fitvalue_newsE=fitvalue_news[fitvalue_news$yyyymm>200712,]

# ====================================     ROC　　===================================
windows()
par(mfrow=c(3,2))
   a=round(performance(prediction(fitvalue_news$Model10, fitvalue_news$code),'auc')@ y.values[[1]],3)
plot(performance(prediction(fitvalue_news$Model10, fitvalue_news$code),"tpr","fpr"),main=sprintf('Without Marc Vars ROC=%s',a)) 
   a=round(performance(prediction(fitvalue_news$Model11, fitvalue_news$code),'auc')@ y.values[[1]],3)
plot(performance(prediction(fitvalue_news$Model11, fitvalue_news$code),"tpr","fpr"),main=sprintf('With Marc Vars ROC=%s',a)) 
   a=round(performance(prediction(fitvalue_news$M10_PD1, fitvalue_news$code),'auc')@ y.values[[1]],3)
plot(performance(prediction(fitvalue_news$M10_PD1, fitvalue_news$code),"tpr","fpr"),main=sprintf('Without Marc Vars and (S)ROC=%s',a) ) 
   a=round(performance(prediction(fitvalue_news$M11_PD1, fitvalue_news$code),'auc')@ y.values[[1]],3)
plot(performance(prediction(fitvalue_news$M11_PD1, fitvalue_news$code),"tpr","fpr"),main=sprintf('With Marc Vars and (S)ROC=%s',a)) 
   a=round(performance(prediction(fitvalue_news$M10_PD, fitvalue_news$code),'auc')@ y.values[[1]],3)
plot(performance(prediction(fitvalue_news$M10_PD, fitvalue_news$code),"tpr","fpr"),main=sprintf('Without Marc Vars and (L) ROC=%s',a) ) 
   a=round(performance(prediction(fitvalue_news$M11_PD, fitvalue_news$code),'auc')@ y.values[[1]],3)
plot(performance(prediction(fitvalue_news$M11_PD, fitvalue_news$code),"tpr","fpr"),main=sprintf('With Marc Vars and (L)ROC=%s',a) ) 

 



