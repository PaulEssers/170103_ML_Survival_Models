set.seed(1)
source("SurvivalCVfunctions.R")
quietLibrary("knitr",quietly=T)
options(stringsAsFactors = F)
quietLibrary("survival",quietly=T)

######################################################
###             Load RadPlat Data                  ###
######################################################


load("../OriginalData/161220_RadPlat.Rdata")
row.names(gene_info)<-row.names(expr.RP.rpkm)
clinical.RP<-radplat.clinical
rm(expr.RP.qn,radplat.clinical)
row.names(clinical.RP)<-clinical.RP$RNA_lab

colnames(expr.RP.rpkm)<-alphaNumeric(colnames(expr.RP.rpkm))
colnames(expr.RP)<-alphaNumeric(colnames(expr.RP))
clinical.RP$RNA_lab<-alphaNumeric(clinical.RP$RNA_lab)

# Remove HPV+ patients
clinical.RP<-clinical.RP[clinical.RP$Reason_exclusion!="Other",]
hpv.neg<-which(clinical.RP$HPV.status=="Negative")
clinical.RP<-clinical.RP[hpv.neg,]
expr.RP<-expr.RP[,clinical.RP$RNA_lab]
expr.RP.rpkm<-expr.RP.rpkm[,clinical.RP$RNA_lab]

# remove non-varying genes, these may cause trouble (?)
varix <- rownames(expr.RP.rpkm)[(apply(expr.RP.rpkm, 1, var, na.rm=TRUE))!=0]
expr.RP.rpkm <- expr.RP.rpkm[varix, ]

#### Rename Data
geneExpr.raw = expr.RP; geneExpr.rpkm=expr.RP.rpkm; 
response<-Surv(clinical.RP$DM.cont,clinical.RP$DM.event)
k = 5; ncore = repeats=10; 

## Calculate expression of signatures
out.file="../OriginalData/170105_SignExpr.Rdata"
if(file.exists(out.file)){load(out.file)}else{
  # load the file with signatures:
  load("~/Databases/Enrichr/Enrichr_Interpretable.Rdata")
  # for some reason there were duplicated signatures in here, remove them because it will cause renaming errors later
  signatures<-signatures.ensembl[unique(names(signatures.ensembl))]
  signExpr<-calculateSignatureExpression(geneExpr=geneExpr.rpkm,signatures)
  save(signExpr,file=out.file)
}


# make sure all column names are correct and the same
colnames(signExpr)<-alphaNumeric(colnames(signExpr))
colnames(geneExpr.rpkm)<-alphaNumeric(colnames(geneExpr.rpkm))
colnames(geneExpr.raw)<-alphaNumeric(colnames(geneExpr.raw))
all.equal(colnames(geneExpr.rpkm),colnames(signExpr))
all.equal(colnames(geneExpr.raw),colnames(signExpr))

# Clinical Variables to use in multivariate Cox model variable selection
clinical.RP$over65<-clinical.RP$AgeDiag>65
clinVar=clinical.RP[,c("Sex","over65","TsiteGen")]

# Time to use as cut-off for variable selection methods that depend on classification
time=2.5*365.25

# Variable Selection methods and machine learning methods used:
# varsel_methods <- c("coeffVar", "uniCox", "multiCox", "DESeq2", "Boruta", "VSURF","mRMR","sSVM_IFR","coeffVar_enr", "uniCox_enr", "multiCox_enr", "Boruta_enr", "VSURF_enr","mRMR_enr","sSVM_IFR_enr")
ml_methods     <- c("rfsrc","superpc","coxnet")

## CV loop
out.file="170217_CV50_RP_DM.Rdata"
if(file.exists(out.file)){load(out.file)}else{
  # Create the data partitions
  require(caret)
  folds<-caret::createMultiFolds(response,k=k,times=repeats)
  
  # setup doSNOW Backend
  quietLibrary("foreach",quietly=T)
  quietLibrary("doSNOW",quietly=T)
  cl <- makeCluster(ncore, outfile="")
  registerDoSNOW(cl)
  
  # loop
  CV.result<-list()
  CV.result<-foreach(fold=1:(k*repeats)) %dopar% {
    tryCatch({
      x<-folds[[fold]]
      # selected_variables<-varselect(x=x, geneExpr.raw=geneExpr.raw, geneExpr.rpkm=geneExpr.rpkm, signExpr=signExpr,response=response)
      # fold.result<-testThresholdCV_Models(x=x,vars=selected_variables,geneExpr.rpkm=geneExpr.rpkm,signExpr=signExpr,response=response)
      selected.vars<-varselect(x=x, geneExpr.raw=geneExpr.raw, geneExpr.rpkm=geneExpr.rpkm, signExpr=signExpr,response=response, clinVar=clinVar, time=time)
      fold.result<-testModels(x=x,vars=selected.vars,geneExpr.rpkm=geneExpr.rpkm,signExpr=signExpr,response=response,models=ml_methods)
      
      return(fold.result)
    }, error=function(e){message("ERROR :",conditionMessage(e), "\n")})
  }
  stopCluster(cl)
  unregister()
  
  save(CV.result,file=out.file)
  result<-matrix(unlist(plyr::ldply(CV.result,cbind)),ncol=5,byrow = F)
  colnames(result)<-c("CI","Dxy","S.D.","model","selection")
  result<-data.frame(result)
  result$CI  <-as.numeric(as.character(result$CI))
  result$Dxy <-as.numeric(as.character(result$Dxy))
  result$S.D.<-as.numeric(as.character(result$S.D.))
  
  save(CV.result,result,file=out.file)
}
# Made a mistake in the code and reversed predictions
# result$CI= 1-result$CI
# result$Dxy= -result$Dxy

pdf("170217_RadPlat_CV_DM.pdf",height=6,width=8)
result<-result[!is.na(result$CI),]
plotCVresults2(result)
plotTopCombinations2(result)
# for some models, no attempt worked, this gives some trouble
#success<-unique(result$model)[c(1,3,4)]
#selected.models<-plotTopCombinations(result[result$model %in% success,],nr=20)
dev.off()

########################################################
####                 Make models                     ###
########################################################

# rename to new-style variable selection
# selected.models$selection<-c("mRMR_gene","mRMR_both","coeffVar_both","coeffVar_gene","sSVMifr_gene")
# 
# sz=20 # number of models in each ensemble
# 
# # Multi-core loop to make the models:
# ncore=5
# library(foreach)
# library(doSNOW)
# cl <- makeCluster(ncore, outfile="")
# registerDoSNOW(cl)
# 
# 
# # loop
# foreach(i=1:nrow(selected.models)) %dopar% {
#   tryCatch({
#     out.file=paste0("../IntermediateData/170130_Surv_",selected.models$selection[i],"_",selected.models$model[i],sz,".Rdata")
#     if(!file.exists(out.file)){
#       mm.ens<-buildEnsemble(geneExpr=geneExpr.rpkm, geneExpr.raw=geneExpr.raw,signExpr=signExpr,response, vselmethod = selected.models$selection[i],model=selected.models$model[i],size=sz)
#       mm.ens<-mm.ens[!sapply(mm.ens,is.null)] # remove failed attempts at building the model
#       class(mm.ens)<-"SurvModelEnsemble"
#       save(mm.ens,file=out.file)
#     }
#     message("finished building model ensemble")
#   }, error=function(e){message("ERROR :",conditionMessage(e), "\n")})
# }
# 
# 
# stopCluster(cl)


#############################################
###    Use on new patient data            ###
#############################################


