setwd("D:\\study\\sub2\\brain_age\\variables_for_normative_modeling\\")
source('Dependencies\\helpers.R') 
library(ggpointdensity)
library(cowplot)
library(patchwork)
library(qqplotr)
library(R.matlab)
library(Hmisc)
library(rstatix)
library(ggpubr)
library(scoring)
library(gamlss)
library(Hmisc)
library(pracma)
#modality<-c( "t1psighc","t1psigsub1","t1npsighc","t1npsigsub1", "t2psigsub2","t2psighc","t2npsighc","t2npsigsub2")#main
modality<-c("t2psighcnew","t2psigsub2new","t1psighcnew","t1psigsub1new","t2npsighcnew","t2npsigsub2new","t1npsighcnew","t1npsigsub1new")#srpbs 25-75
#modality<-c("t1psighc","t1npsighc","t1psighc2","t1npsighc2") #permut(11-65)/boot(10-65)
#modality<-c("osext1psighc","osext1npsigsub1","osext1npsighc","osext2psighc","osext2npsighc","osext2npsigsub2","osext2psigsub2","osext2npsigsub2") #validation
#modality<-c("ot1psighc","ot1npsighc","ot1npsigsub1","ot2psighc","ot2npsighc","ot1psigsub1","ot2psigsub2","ot2npsigsub2") #validation
#modality<-c("omfdt1psighc","omfdt1npsighc","omfdt2psighc","omfdt2npsighc","omfdt1psigsub1","omfdt1npsigsub1","omfdt2psigsub2","omfdt2npsigsub2") #validation
modality<-c("t1asigsub2")#abs (10-65)
type<-c("2")
type<-seq(953,1000)
#type<-seq(1)

# run analysis for each modality and each type (i.e., cross-sectional and longitudinal)
for (idx_m in 1:length(modality)){ 
  for (idx_type in 1:length(type)){
    setwd("D:\\study\\sub2\\brain_age\\variables_for_normative_modeling\\permut_abs\\")
    #setwd("D:\\study\\sub2\\brain_age\\variables_for_normative_modeling\\permut\\")
    #setwd("D:\\study\\sub2\\brain_age\\variables_for_normative_modeling\\vali\\")
    fname<-paste(modality[idx_m],"_",type[idx_type],"_lo_vars_for_R_","cross.mat",sep="")
    M<-readMat(fname)
    # read data
    mydata <- as.data.frame(cbind(M$age,M$sex,M$covars,M$y))
    names(mydata)<-c("age","sex","covar","phenotype");
    # Model training
    mdl <- gamlss(phenotype ~ fp(age,npoly=1)+sex+pb(covar) , # formula for mu
                  sigma.fo=~fp(age,npoly=1),
                  nu.fo=~1,
                  nu.tau=~1,
                  family=BCT(),
                  data=mydata,
                  control = gamlss.control(n.cyc = 500))
    #setwd("D:\\study\\sub2\\brain_age\\variables_for_normative_modeling\\output_normative_modeling\\all\\")
    #outname<-paste("GAMLSS_",modality[idx_m],"_",type[idx_type],"_predicted.png",sep="")
    #png(outname)             # 分辨率
    #plot(mdl)
    #dev.off()

if (mdl$converged) {
  print("Model Converged : YES")
  converged<-1
  setwd("D:\\study\\sub2\\brain_age\\variables_for_normative_modeling\\")
      for (gender in 0:1) {
        print(paste("gender",gender,sep=" "))
        min_age<-10
        max_age<-65
        age_test<-linspace(min_age,max_age,1000)
        sex_test<-matrix(data=gender,nrow=1000,ncol=1)
        covar_test<-matrix(data=mean(mydata$covar),nrow=1000,ncol=1)
        y_test<-matrix(data=mean(mydata$phenotype),nrow=1000,ncol=1) # for the sake of the data structure required for the function
        # but not contribute to the prediction
        data_test<-as.data.frame(cbind(age_test,sex_test,covar_test,y_test))
        names(data_test)<-c("age","sex","covar","phenotype")
        params<-predictAll(mdl,data=mydata,
                           newdata=data_test,output='matrix',type="response",
                           y.value="median",what=c("mu", "sigma", "nu","tau"))


        #main percentiles
        plot(phenotype~age, data=mydata,  col="lightgray", pch=10)
        predictions_quantiles<-matrix(data=0,ncol=5,nrow=1000)
        quantiles <- pnorm(c(-2:2))
        for (i in 1:length(quantiles)){
          Qua <- getQuantile(mdl, quantile=quantiles[i],term="age",fixed.at=list(sex=gender))
          out<-curve(Qua, min_age, max_age,  lwd=2, lty=1, add=T,col="red",n = 1000)
          predictions_quantiles[,i]=as.vector(out$y)
        }


        # save to mat files
        outname<-paste("output_normative_modeling/permut_abs/GAMLSS_",modality[idx_m],"_",type[idx_type],"_predicted_sex",gender,"_WHOLE_SAMPLE.mat",sep="")
        #outname<-paste("output_normative_modeling/permut/GAMLSS_",modality[idx_m],"_",type[idx_type],"_predicted_sex",gender,"_WHOLE_SAMPLE.mat",sep="")
        #outname<-paste("output_normative_modeling/vali/GAMLSS_",modality[idx_m],"_",type[idx_type],"_predicted_sex",gender,"_WHOLE_SAMPLE.mat",sep="")
        writeMat(outname, predictions_quantiles=as.matrix(predictions_quantiles), age=age_test)
      }
      # #for sex validatoin
      #   min_age<-10
      #   max_age<-65
      #   age_test<-linspace(min_age,max_age,1000)
      #   sex_test<-matrix(data=idx_type-1,nrow=1000,ncol=1)
      #   covar_test<-matrix(data=mean(mydata$covar),nrow=1000,ncol=1)
      #   y_test<-matrix(data=mean(mydata$phenotype),nrow=1000,ncol=1) # for the sake of the data structure required for the function
      #   # but not contribute to the prediction
      #   data_test<-as.data.frame(cbind(age_test,sex_test,covar_test,y_test))
      #   names(data_test)<-c("age","sex","covar","phenotype")
      #   params<-predictAll(mdl,data=mydata,
      #                      newdata=data_test,output='matrix',type="response",
      #                      y.value="median",what=c("mu", "sigma", "nu","tau"))
      # 
      # 
      #   #main percentiles
      #   plot(phenotype~age, data=mydata,  col="lightgray", pch=10)
      #   predictions_quantiles<-matrix(data=0,ncol=5,nrow=1000)
      #   quantiles <- pnorm(c(-2:2))
      #   for (i in 1:length(quantiles)){
      #     Qua <- getQuantile(mdl, quantile=quantiles[i],term="age",fixed.at=list(sex=gender))
      #     out<-curve(Qua, min_age, max_age,  lwd=2, lty=1, add=T,col="red",n = 1000)
      #     predictions_quantiles[,i]=as.vector(out$y)
      #   }
      # 
      #   outname<-paste("output_normative_modeling/vali/GAMLSS_",modality[idx_m],"_",type[idx_type],"_predicted_WHOLE_SAMPLE.mat",sep="")
      #   writeMat(outname, predictions_quantiles=as.matrix(predictions_quantiles), age=age_test)
      # # #for sex validatoin
     } else {
       print("Model Cnverged: NO")
       converged<-0
    }
    #gvd1<-deviance(mdl,what="G") # validation global deviance
    #iterations<-mdl$iter # number of iterations taken to converge
    #params<-predictAll(mdl,data=mydata, newdata=mydata, 
                      # output='matrix',type="response",
                       #y.value="median",what=c("mu", "sigma", "nu", "tau"))
    
    #quantiles <- pBCT(params[,1],mu=params[,2],sigma=params[,3],nu=params[,4],tau=params[,5]); #for BCT  
    #z_randomized<-qnorm(quantiles, mean = 0, sd =1);
    #ks_out<-ks.test(z_randomized,"pnorm")
    #path<-"output_normative_modeling";
    #dir.create(path)
    #outname2<-paste("output_normative_modeling/fitstats_GAMLSS_",modality[1],"_",type[1],"_WHOLE_SAMPLE.mat",sep="")
    #print(paste("KS test, p=",ks_out$p.value,sep=""))
    #writeMat(outname2,  AIC1=mdl$aic, gvd1=gvd1,  BIC=mdl$sbc, zscores=z_randomized, ks_p=ks_out$p.value,ks_d=ks_out$statistic, converged=converged, P_deviance=mdl$P.deviance)
    ### predict hypothetical data ###

  }
  
}
