## Functions
quietLibrary<-function(x,quietly=T){
  suppressMessages(suppressWarnings(require(x,character.only = T)))
}

# needed to rename signature names to R dimnames standards
is.letter <- function(x) grepl("[[:alpha:]]", x)
alphaNumeric<-function(s){
  first<-sapply(s,function(x){strsplit(x,"")[[1]][1]})
  new<-sapply(1:length(s),function(i){ ifelse(is.letter(first[i]),s[i],paste0("X",s[i]))})
  return(gsub("[^[:alnum:]]", "_", new))
}

calculateSignatureExpression<-function(geneExpr,signatures){
  quietLibrary("GSVA", quietly=T)
  # RPKM values are the input
  sigExpression<-gsva(as.matrix(geneExpr), signatures, min.sz=5, max.sz=500, rnaseq=TRUE, method="gsva",parallel.sz=4)
  # change between algorithms with method parameter: method=c("gsva", "ssgsea", "zscore", "plage")
  signExpr<-sigExpression$es.obs
  row.names(signExpr)<-alphaNumeric(row.names(signExpr))
  return(data.frame(signExpr))
}

## functions needed for variable selection

# Coefficient of Variation
coeffVar<-function(x){sd(x,na.rm=T)/abs(mean(x,na.rm=T))}

# find genes which have measurements in a small subset of samples
proportionZeros<-function(x,cutoff=0.5){mean(x==0)<cutoff}

### Variable Selection

coeffVar_varselect<-function(expr, topnvar=4000){
  message("Calculating coefficient of variance")
  # Select genes and signatures separately, because their variances are calculated differently
  genes.id<-grep("ENS",colnames(expr))
  
  # If there are genes in expr, get the top genes
  if(length(genes.id)>0){
    genes<-expr[,genes.id]
    # remove variables with too many 0-measurements, then report back the most variable genes
    genes.vars <- colnames(genes)[apply(genes, 2, proportionZeros)]
    genes.vars <- colnames(genes[,genes.vars])[order(apply(genes[,genes.vars], 2, coeffVar), decreasing=TRUE)][1:(topnvar/2)] 
  }else{
    genes.vars<-NULL
  }
  # if the number of genes is the total number of variables, skip the signature step.
  if(length(genes.id)==ncol(expr)){return(genes.vars)}
    
  if(length(genes.id)==0){
    signatures<-expr
  }else{
    signatures<-expr[,-genes.id]
  }
  sig.vars <- colnames(signatures)[apply(signatures, 2, proportionZeros)]
  sig.vars <- colnames(signatures[,sig.vars])[order(apply(signatures[,sig.vars], 2, coeffVar), decreasing=TRUE)][1:(topnvar/2)]
  return(c(genes.vars,sig.vars))
}
# coeffVar_varselect(t(rbind(expr.UTSCC.rpkm,enr.UTSCC)))

cox_varselect<-function(expr, response, clinical.vars=NULL,fdr=0.1){
  # make sure names of the survival object and the clinical variables data.frame are identical
  message("Calculating performance in cox models")
  quietLibrary("survival", quietly=T)
  out.df<-data.frame(matrix(ncol=5,nrow=ncol(expr)),row.names=colnames(expr))
  
  if(!is.null(clinical.vars)){
    temp.df<-cbind(clinical.vars,response)
  }else{
    temp.df<-data.frame(response=response)
  }
  
  for(gene in colnames(expr)){
    tryCatch({
      # message(paste0(gene,"\r"))
      temp.df$gene=expr[,gene]
      cs<-summary(coxph(response~.,data=temp.df))
      out.df[gene,]<-cs$coefficients["gene",]
    #}, error=function(e){message("Error in single gene Cox model\n")})
    }, error=function(e){})
  }
  # backup<-out.df
  # out.df<-backup
  out.df<-out.df[is.finite(out.df[,5]),]
  
  out.df$p.adj<-p.adjust(out.df[,5],method="fdr")
  vars<-row.names(out.df[out.df$p.adj<fdr,])
  return(vars)
}

survival2category<-function(expr, response, time=2.5*365.25){
  # assumes expr and response are in the same order!
  response<-as.matrix(response)
  message(paste0("Response Dimensions:", dim(response)))
  class<-ifelse(response[,2]==0 & response[,1]>time,"Negative",
         ifelse(response[,2]==0 & response[,1]<time,"Censored",
         ifelse(response[,2]==1 & response[,1]>time,"Negative","Positive")))
  class<-factor(class)
  
  out.obj<-list()
  out.obj[["expr"]]<-expr[class!="Censored",]
  out.obj[["class"]]<-factor(class[class!="Censored"])
  return(out.obj)
}


DESeq2_varselect<-function(expr, response,logfc=1,fdr=0.1, time=2.5*365.25){
  message("Running DESeq2")
  quietLibrary("DESeq2" , quietly=T)
  
  # remove censored patients at timepoint
  nonCensored<-survival2category(t(expr), response=response, time)
  expr<-t(nonCensored[["expr"]])
  response=data.frame(Class=nonCensored[["class"]] )

  DESeq2Table <- DESeqDataSetFromMatrix(countData=round(expr,0), colData=response, design= ~Class)
  DESeq2Table <- estimateSizeFactors(DESeq2Table)
  DESeq2Table <- estimateDispersions(DESeq2Table)
  DESeq2Table <- nbinomWaldTest(DESeq2Table)

  message("finished DESeq2 run")

  numberGroups<-length(levels(response$Class))
  if(numberGroups==2){
    # need to analyze only one group in case of 2 groups total
    resnames<-resultsNames(DESeq2Table)[3]
  }else{
    resnames<-resultsNames(DESeq2Table)[-1]
  }
  signature<-sapply(resnames,function(name){
    res<-results(DESeq2Table,name=name)
    res<-res[!is.na(res$padj),]
    sig<-res[which(res$padj<fdr & abs(res$log2FoldChange)>logfc),]
    return(row.names(sig))
  })
  return(as.character(unique(unlist(signature))))
}
# test.DESeq2<-DESeq2_varselect(expr=expr.RP, response=Surv(clinical.RP$DM.cont,clinical.RP$DM.event))
# 
Boruta_varselect<-function(expr, response,topnvar=8000,verbose=1, time=2.5*365.25){
  message("Selecting Variables with Boruta")
  quietLibrary("Boruta", quietly=T)

  # remove censored patients at timepoint
  nonCensored<-survival2category(expr=expr, response=response, time)
  expr=nonCensored[["expr"]]
  response=nonCensored[["class"]]

  varix<-coeffVar_varselect(expr,topnvar=topnvar)
  expr <- expr[,varix]

  bor<-Boruta(x=expr,y=response,maxRuns = 250, doTrace=verbose,holdHistory = F, num.threads=1)
  return(names(bor$finalDecision[bor$finalDecision!="Rejected"]))
}
# test.bor<-Boruta_varselect(expr=t(expr.RP), response=Surv(clinical.RP$DM.cont,clinical.RP$DM.event))

VSURF_varselect<-function(expr, response, time=2.5*365.25){
  message("Selecting Variables with VSURF")
  quietLibrary("VSURF", quietly=T)
  
  # remove censored patients at timepoint
  message("Converting survival data to categories for VSURF")
  nonCensored<-survival2category(expr=expr, response=response, time=time)
  expr=nonCensored[["expr"]]
  response=nonCensored[["class"]]
  
  message("Finished converting to categories for VSURF")
  # reduce nr of variables, or it will take forever. The function will return only half if there are no signatures.
  varix<-coeffVar_varselect(expr,topnvar=8000)
  expr<-expr[,varix]
  # variable selection using VSURF
  vsurf.out<-VSURF(x=expr,y=response, nfor.thres = 100, parallel = FALSE)
  # unregister() # unregister parallel backend used by vsurf or it will lead to trouble
  return(colnames(expr)[vsurf.out$varselect.interp])
}
# # test<-VSURF_varselect(t(expr.UTSCC.rpkm),response=phenoDat$Class)
 
mRMR_varselect<-function(expr, response, n=20, sl=100, time=2.5*365.25){
  message("Running mRMRe")
  tryCatch({
    quietLibrary("mRMRe", quietly=T)

    nonCensored<-survival2category(expr=expr, response=response, time)
    expr=nonCensored[["expr"]]
    response=nonCensored[["class"]]
    
    # mRMRe requires ordered factors from some reason
    response=factor(response,ordered = T)

    # Force using only one core
    set.thread.count(1)

    # preselect the top nvar variable genes.
    subs<-coeffVar_varselect(expr,topnvar=8000)
    df<- data.frame(expr[,subs])
    
    mdata <- mRMR.data(data=df,strata=response)
    ens1.m <- mRMR.classic(mdata, feature_count=sl, target_indices=1)
    ens1.fs <- apply(solutions(ens1.m)[[1]], 2, function(x, y) { return(y[x]) }, y=mRMRe::featureNames(mdata))
    return(ens1.fs[,1])

  }, error=function(e){message("ERROR :",conditionMessage(e), "\n")})
}
# mRMR_varselect(t(expr.UTSCC.rpkm),response=phenoDat$Class)

sSVMifr_varselect<-function(expr,resp,n=100 ,k=3,threshold=0.9,plot=F,output=1, time=2.5*365.25){
  quietLibrary("penalizedSVM", quietly=T)
  # output =1 display iteration nr, =2 get all the output from lpsvm

  nonCensored<-survival2category(expr=expr, response=resp, time)
  expr=nonCensored[["expr"]]
  resp=nonCensored[["class"]]
  
  resp<-factor(resp)
  if(length(levels(resp))!=2){
    message("This only works for 2 classes!")
    return(NULL)
  }
  resp<-ifelse(resp==levels(resp)[1],-1,1)
  gs=list()
  res<-data.frame(iter=1:n,trainCorr=0,testCorr=0)

  for(i in 1:n){
    if(output>0){cat("\r","iteration no: ",i)}
    mm<-lpsvm(expr[,!(colnames(expr) %in% unlist(gs))],resp,k=k,output=max(0,output-1))
    res$trainCorr[i]<-mm$trainCorr
    res$testCorr[i] <-mm$testCorr
    gs[[i]]<- colnames(expr[,!(colnames(expr) %in% unlist(gs))])[mm$xind]

    if(i>10){
      # break if performance is under the treshold for 3 iterations in a row
      underThreshold= res$testCorr< threshold*quantile(res$testCorr[1:10],0.50)
      stop<-sapply(res$iter,function(j){
        if(j>n-3){return(F)}else{return(sum(underThreshold[j:(j+2)])>2  )}
      })
      thr<-min(res$iter[stop==T])

      if(thr<i){message("Performance is under threshold. Finished!"); break}
    }
  }
  #remove the last two iterations if the loop stopped before reaching the maximum nr of iterations
  if(length(gs)<n){gs=gs[1:(length(gs)-2)]}

  if(plot==T){
    plot(res[,2],type="b",ylim=c(0,100),ylab="correctness",xlim=c(0,length(gs)+2))
    lines(res[,3],col="red",type="b")
    abline(h=threshold*quantile(res$testCorr[1:10],0.50))
    # using quantiles avoids the line being dependend on outliers
  }
  return(unlist(gs))
}


varselect<-function(x, geneExpr.raw=NA, geneExpr.rpkm, signExpr,response, clinVar=NULL, methods=c("coeffVar_gene", "coeffVar_sign","coeffVar_both","uniCox_gene","uniCox_sign","uniCox_both", "multiCox_gene","multiCox_sign","multiCox_both", "DESeq2", "Boruta_gene","Boruta_sign","Boruta_both", "VSURF_gene","VSURF_sign","VSURF_both","mRMR_gene","mRMR_sign","mRMR_both","sSVMifr_gene","sSVMifr_sign","sSVMifr_both"),time=2.5*365.25){
  # x = samples in the current fold
  combExpr<-t(rbind(geneExpr.rpkm,signExpr))

  selected<-list()
  if("DESeq2" %in% methods){selected[["DESeq2"]]<-DESeq2_varselect(geneExpr.raw[,x], response[x,], time=time) }
  
  # all methods with only gene expression
  if("coeffVar_gene" %in% methods){selected[["coeffVar_gene"]]<-coeffVar_varselect(t(geneExpr.rpkm)[x,]) }
  if("uniCox_gene" %in% methods){selected[["uniCox_gene"]]<-cox_varselect(t(geneExpr.rpkm)[x,],response[x,],clinical.vars = NULL)}
  if("multiCox_gene" %in% methods){selected[["multiCox_gene"]]<-cox_varselect(t(geneExpr.rpkm)[x,],response[x,],clinical.vars = clinVar[x,])}
  if("Boruta_gene" %in% methods){selected[["Boruta_gene"]]<-Boruta_varselect(t(geneExpr.rpkm)[x,], response[x,], time=time) }
  if("VSURF_gene" %in% methods){selected[["VSURF_gene"]]<-VSURF_varselect(t(geneExpr.rpkm)[x,], response[x,], time=time) }
  if("mRMR_gene" %in% methods){selected[["mRMR_gene"]]<-mRMR_varselect(t(geneExpr.rpkm)[x,], response[x,], time=time) }
  if("sSVMifr_gene" %in% methods){selected[["sSVMifr_gene"]]<-sSVMifr_varselect(t(geneExpr.rpkm)[x,], response[x,], time=time) }

  # all methods with gvsa signature expression (except DESeq2, there it is not possible to do this)
  if("coeffVar_sign" %in% methods){selected[["coeffVar_sign"]]<-coeffVar_varselect(t(signExpr)[x,]) }
  if("uniCox_sign" %in% methods){selected[["uniCox_sign"]]<-cox_varselect(t(signExpr)[x,],response[x,],clinical.vars = NULL)}
  if("multiCox_sign" %in% methods){selected[["multiCox_sign"]]<-cox_varselect(t(signExpr)[x,],response[x,],clinical.vars = clinVar[x,])}
  if("Boruta_sign" %in% methods){selected[["Boruta_sign"]]<-Boruta_varselect(t(signExpr)[x,], response[x,], time=time) }
  if("VSURF_sign" %in% methods){selected[["VSURF_sign"]]<-VSURF_varselect(t(signExpr)[x,], response[x,], time=time) }
  if("mRMR_sign" %in% methods){selected[["mRMR_sign"]]<-mRMR_varselect(t(signExpr)[x,], response[x,], time=time) }
  if("sSVMifr_sign" %in% methods){selected[["sSVMifr_sign"]]<-sSVMifr_varselect(t(signExpr)[x,], response[x,], time=time) }
  
  # all methods with combined gene expression and gvsa signature expression (except DESeq2, there it is not possible to do this)
  if("coeffVar_both" %in% methods){selected[["coeffVar_both"]]<-coeffVar_varselect(combExpr[x,]) }
  if("uniCox_both" %in% methods){selected[["uniCox_both"]]<-cox_varselect(combExpr[x,],response[x,],clinical.vars = NULL)}
  if("multiCox_both" %in% methods){selected[["multiCox_both"]]<-cox_varselect(combExpr[x,],response[x,],clinical.vars = clinVar[x,])}
  if("Boruta_both" %in% methods){selected[["Boruta_both"]]<-Boruta_varselect(combExpr[x,], response[x,], time=time) }
  if("VSURF_both" %in% methods){selected[["VSURF_both"]]<-VSURF_varselect(combExpr[x,], response[x,], time=time) }
  if("mRMR_both" %in% methods){selected[["mRMR_both"]]<-mRMR_varselect(combExpr[x,], response[x,], time=time) }
  if("sSVMifr_both" %in% methods){selected[["sSVMifr_both"]]<-sSVMifr_varselect(combExpr[x,], response[x,], time=time) }
  
  
  return(selected)
}


### Define Models


### Random Forest ###

# Do 1 bootstrap for performance estimation
trainPerf.rfsrc<-function(data, mtry=5, nodesize=3){
  # message("Now in trainPerf.rfsrc")
  quietLibrary("randomForestSRC", quietly = TRUE)
  quietLibrary("survival", quietly = TRUE)
  
  indices<-base::sample(1:nrow(data),nrow(data),replace=T)
  
  mm<-rfsrc(Surv(time,status)~.,data=data[indices,],ntree=2000,mtry=mtry,nodesize=nodesize)
  #pp<-1-predict(mm,newdata=data[-indices,],outcome="test")$predicted #-> this actually uses the test-set outcome...
  pp<-1-predict(mm,newdata=data[-indices,-c(1,2)])$predicted
  
  Hmisc::rcorr.cens(1-pp, Surv(data$time, data$status)[-indices,])["C Index"]
}

# Do 10 bootstraps for a grid of parameters and return the model with best parameters.
tune.rfsrc<-function(data, number=5, tuneLength=5, showPerformance=F){
  # message("Now in tune.rfsrc")
  quietLibrary("randomForestSRC", quietly = TRUE)
  quietLibrary("survival", quietly = TRUE)
  quietLibrary("caret", quietly=T)
  
  grid <- expand.grid(mtry = caret::var_seq(p = ncol(data)-2, classification = F, len = tuneLength), 
                      nodesize = unique(floor(seq(1,6, length = tuneLength))))
  
  CI<-apply(grid,1,function(param){
    bsCI<-mean(replicate(number,{trainPerf.rfsrc(data=data,mtry=param[1],nodesize=param[2])}))
  })
  if(showPerformance==T){print(cbind(grid,CI))}
  bp<-grid[which.max(CI),]
  mm<-rfsrc(Surv(time,status)~.,data=data,ntree=2000,mtry=bp$mtry,nodesize=bp$nodesize)
  return(mm) 
}


### Supervised Principal Components
tune.superpc<-function(data, plotRatios=F){
  quietLibrary("superpc", quietly = TRUE)
  quietLibrary("survival", quietly = TRUE)
  quietLibrary("Biobase")
  
  survcols <- which(colnames(data) %in% c("time","status"))
  data.in  <- list(x=t(data[,-survcols]),y=data$time, censoring.status=data$status, featurenames=colnames(data[,-survcols]))
  train.obj<- superpc.train(data.in, type="survival")
  capture.output(cv.obj   <- superpc.cv(train.obj, data.in) )
  
  # find best settings:
  ratios<-sapply(1:3,function(i){
    x<-cv.obj$scor[i,]
    x/qchisq(0.95, i)
  })
  if(plotRatios==T){
    plot(x=cv.obj$thresholds,y=ratios[,1],type="l",ylim=c(min(ratios),max(ratios)),ylab="loglikelyhood score / threshold",xlab="threshold")
    abline(h=1)
    lines(x=cv.obj$thresholds,y=ratios[,2],col="red")
    lines(x=cv.obj$thresholds,y=ratios[,3],col="darkgreen")
    text(x=0,y=ratios[1,1],labels=1); text(x=0,y=ratios[1,2],labels=2); text(x=0,y=ratios[1,3],labels=3)
  }
  
  best.thresh<-cv.obj$thresholds[which.max(rowMax(ratios))]
  best.noComps<-which.max(rowMax(t(ratios)))
  
  return.object<-list(mm=train.obj,data.train=data.in, bestThresh=best.thresh, bestNoComp=best.noComps)
  return(return.object)
}

predictions.superpc<-function(spcm,data){
  survcols <- which(colnames(data) %in% c("time","status"))
  data.test  <- list(x=t(data[,-survcols]),y=data$time, censoring.status=data$status, featurenames=colnames(data[,-survcols]))
  
  fit.cts<- superpc.predict(spcm$mm, spcm$data.train, newdata=data.test, threshold=spcm$bestThresh, n.components=spcm$bestNoComp, prediction.type="continuous")
  pp<- fit.cts$v.pred.1df
  names(pp)<-names(fit.cts$v.pred.1df)
  return(pp)
}

### CoxNet
tune.coxnet<-function(data,tuneLength=11){
  quietLibrary("Coxnet",quietly = T)
  quietLibrary("glmnet",quietly=T)
  
  survcols <- which(colnames(data) %in% c("time","status"))
  alphas<-seq(0,1,length.out=tuneLength)
  cv.result<-t(sapply(alphas,function(a){
    cv.fit <- cv.glmnet(as.matrix(data[,-survcols]), Surv(data$time, data$status), family = "cox", maxit = 10000)
    return(c(a,min(cv.fit$cvm)))
  }))
  best.alpha<-cv.result[which.min(cv.result[,2]),1]
  cv.fit <- cv.glmnet(as.matrix(data[,-survcols]), Surv(data$time, data$status), family = "cox", maxit = 10000,alpha=best.alpha)
  return(cv.fit)
}

testModels<-function(x,vars,geneExpr.rpkm,signExpr,response, models=c("rfsrc","superpc","coxnet"), scale=T){
  quietLibrary("caret", quietly=T)

  fulldata<-t(rbind(geneExpr.rpkm,signExpr))
  if(scale){fulldata<-scale(fulldata)}
  
  res.total<-sapply(1:length(vars),function(i){
    message(names(vars)[i])
    subs<-vars[[i]]

    tr_set<-data.frame(cbind(response[x ,],fulldata[ x,subs]),stringsAsFactors = F)
    te_set<-data.frame(cbind(response[-x,],fulldata[-x,subs]),stringsAsFactors = F)
    
    res<-data.frame(matrix(nrow=length(models),ncol=3))
    rownames(res)<-models
    colnames(res)<-c("CI","Dxy","S.D.")

    for(model in models){
      tryCatch({suppressWarnings({
        message(model)
        # appearantly, caret ignores unused arguments, so it is possible to give specific arguments for the RF and NNET, and others as long as there is no overlap
        if(model=="rfsrc"){
          mm<-tune.rfsrc(tr_set)
          pp<- 1-predict(mm,newdata=te_set,outcome="test")$predicted
          res[model,1:3]<-Hmisc::rcorr.cens(pp, Surv(te_set$time, te_set$status))[1:3]
        }else if(model=="superpc"){
          spcm<-tune.superpc(tr_set)
          pp<- 1-predictions.superpc(spcm,data=te_set)
          res[model,1:3]<-Hmisc::rcorr.cens(pp, Surv(te_set$time, te_set$status))[1:3]
        }else if(model=="coxnet"){
          mm<-tune.coxnet(tr_set)
          pp<-1-predict(mm,newx=as.matrix(te_set[,-c(1,2)]),s="lambda.1se")
          res[model,1:3]<-Hmisc::rcorr.cens(pp, Surv(te_set$time, te_set$status))[1:3]
        }
      })}, error=function(e){message("ERROR :",conditionMessage(e), "\n")})
    }

    res$Model<-row.names(res)
    res$VarSelection<-names(vars)[i]
    print(res)
    return(t(res))
  })
  #require(plyr, quietly=T)
  final.result<-matrix(unlist(plyr::ldply(res.total,cbind)),ncol=5,byrow = T)
  return(final.result)

}

unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}

plotCVresults<-function(result){
  quietLibrary("tidyr", quietly=T)
  quietLibrary("dplyr", quietly=T)
  quietLibrary("gplots", quietly=T)
  my_palette <- colorRampPalette(c("black", "white"))(n = 99)
  colors <- c(seq(0,0.65,0.5),seq(0.66,1,0.01))
  my_palette <- colorRampPalette(c("black", "white"))(n = length(colors)-1)

  # CI
  CI<- as.tbl(result) %>% group_by(selection,model) %>% dplyr::summarize(CI=mean(CI,na.rm=T)) %>% spread(model,CI) %>% data.frame()
  row.names(CI)<-CI$selection
  CI<-as.matrix(CI[,-1])
  CI[!is.finite(CI)]<-0
  heatmap.2(as.matrix(CI),scale="none",dendrogram = "none",trace="none", col=my_palette, breaks=colors,cellnote=ifelse(CI==0, NA, round(CI,2)),notecol="black",key=F,margins = c(10,10), Rowv=FALSE,Colv=FALSE,main="Concordance Index")

  # Dxy
  Dxy <- as.tbl(result) %>% group_by(selection,model) %>% dplyr::summarize(Dxy=mean(Dxy,na.rm=T)) %>% spread(model,Dxy)  %>% data.frame()
  row.names(Dxy)<-Dxy$selection
  Dxy<-as.matrix(Dxy[,-1])
  Dxy[!is.finite(Dxy)]<-0
  heatmap.2(as.matrix(Dxy),scale="none",dendrogram = "none",trace="none", col=my_palette, breaks=colors,cellnote=ifelse(Dxy==0, NA, round(Dxy,2)),notecol="black",key=F,margins = c(10,10), Rowv=FALSE,Colv=FALSE,main="Dxy (Correlation)")
}

plotTopCombinations<-function(result,nr=5,returnScore="CI",oneSEM=F){
  quietLibrary("dplyr", quietly=T)
  quietLibrary("ggplot2", quietly=T)
  quietLibrary("plotrix", quietly=T)
  summ<-result %>% as.tbl() %>% group_by(selection,model) %>% dplyr::summarize(CI.av=mean(CI,na.rm=T),Dxy.av=mean(Dxy,na.rm=T),CI.sem=plotrix::std.error(CI,na.rm=T),Dxy.sem=plotrix::std.error(Dxy,na.rm=T))

  result$combination<-apply(result[,c(5,4)],1,function(x){paste0(x,collapse="_")})

  summ.CI<-summ[order(summ$CI.av,decreasing=T),]
  if(oneSEM){nr=length(summ.CI[summ.CI$CI.av>(summ.CI$CI.av[1]-summ.CI$CI.sem[1]),])}
  combi<-apply(summ.CI[1:nr,1:2],1,function(x){paste0(x,collapse="_")})
  result.sel<- result[result$combination %in% combi,]
  result.sel$combination<-factor(result.sel$combination,levels=combi)


  p<-ggplot(data=result.sel, aes(x=combination, y=CI)) +
    stat_boxplot(geom ='errorbar') +
    geom_boxplot(position=position_dodge(0.9))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0.1,1,0.1))+
    scale_fill_manual(values=c("black","lightgrey","darkgrey","grey"))+
    theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=10))+
    ylab("") + xlab("")+
    ggtitle("Concordance Index")
  plot(p)

  summ.Dxy<-summ[order(summ$Dxy.av,decreasing=T),]
  if(oneSEM){nr=length(summ.Dxy[summ.Dxy$Dxy.av>(summ.Dxy$Dxy.av[1]-summ.Dxy$Dxy.sem[1]),])}
  combi<-apply(summ.Dxy[1:nr,1:2],1,function(x){paste0(x,collapse="_")})
  result.sel<- result[result$combination %in% combi,]
  result.sel$combination<-factor(result.sel$combination,levels=combi)

  p<-ggplot(data=result.sel, aes(x=combination, y=Dxy)) +
    stat_boxplot(geom ='errorbar') +
    geom_boxplot(position=position_dodge(0.9))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0.1,1,0.1))+
    theme(axis.text.x=element_text(angle=45,vjust=1,hjust=1,size=10))+
    scale_fill_manual(values=c("black","lightgrey","darkgrey","grey"))+
    ylab("") + xlab("")+
    ggtitle("Dxy (Correlation)")
    #+geom_hline(aes(yintercept=1))
  plot(p)

  if(returnScore=="CI"){return(summ.CI[1:nr,1:2])}
  else if(returnScore=="Dxy"){return(summ.Dxy[1:nr,1:2])}
  else{message("Not returning top combinations, use returnScore='CI' or 'Dxy'.")}
}

plotCVresults2<-function(result){
  require(tidyr)
  require(dplyr)
  require(gplots)
  require(gridExtra)
  
  # AUC
  CI<- as.tbl(result) %>% group_by(selection,model) %>% dplyr::summarize(CI=mean(CI,na.rm=T),Dxy=mean(Dxy,na.rm=T),failed=n()/50) %>% data.frame()
  CI$CI[!is.finite(CI$CI)]<-0
  CI$Dxy[!is.finite(CI$Dxy)]<-0
  #CI$failed[!is.finite(CI$failed)]<-0
  
  p1<-ggplot(CI, aes(model,selection)) +
    geom_tile(aes(fill = CI), colour = "black") +
    theme(legend.position="none")+
    scale_fill_gradient(low = "black",high = "white",na.value="black")+
    geom_text(aes(label = round(CI, 2))) +
    theme(axis.text.x = element_text(size=10, angle=90),axis.text.y=element_text(size=10))+
    ggtitle("CI")
  
  p2<-ggplot(CI, aes(model,selection)) +
    geom_tile(aes(fill = Dxy), colour = "black") +
    theme(legend.position="none")+
    scale_fill_gradient(low = "black",high = "white",na.value="black")+
    geom_text(aes(label = round(Dxy, 2))) +
    theme(axis.text.x = element_text(size=10, angle=90),axis.text.y=element_text(size=10))+
    ggtitle("Dxy")
  
  p3<-ggplot(CI, aes(model,selection)) +
    geom_tile(aes(fill = failed), colour = "black") +
    theme(legend.position="none")+
    scale_fill_gradient(low = "black",high = "white",na.value="black")+
    geom_text(aes(label = round(failed, 2))) +
    
    theme(axis.text.x = element_text(size=10, angle=90),axis.text.y=element_text(size=10))+
    ggtitle("# Succesfull Attempts")
  
  grid.arrange(p1,p2,p3, ncol=3)
}


plotTopCombinations2<-function(result,nr=5,returnScore="CI",oneSD=F,diff=0,ord="mean"){
  require(gridExtra)
  
  
  require(dplyr)
  require(ggplot2)
  # Some settings for nicer plots
  theme_set(theme_bw())
  theme_update(strip.background = element_blank(),
               strip.text = element_text(size=12),
               axis.text = element_text(size=12),
               axis.text.x = element_text(angle=90,vjust=0.8,hjust=1,size=8),
               axis.ticks.length = unit(0.1, "cm"),
               #axis.text.x = element_blank(),
               legend.key = element_rect(fill = "white", color = "white"),
               #panel.grid.major = element_blank(),
               panel.grid.minor = element_blank(),
               panel.border = element_rect(colour = "white"),
               panel.background = element_blank(),
               axis.line.x = element_line(colour = "black"),
               axis.line.y = element_line(colour = "black"))
  
  if(ord=="median"){
    summ<-result %>% as.tbl() %>% group_by(selection,model) %>% dplyr::summarize(CI.av=median(CI,na.rm=T),Dxy.av=median(Dxy,na.rm=T),CI.sd=sd(CI,na.rm=T),Dxy.sd=sd(Dxy,na.rm=T))
  }else{
    summ<-result %>% as.tbl() %>% group_by(selection,model) %>% dplyr::summarize(CI.av=mean(CI,na.rm=T),Dxy.av=mean(Dxy,na.rm=T),CI.sd=sd(CI,na.rm=T),Dxy.sd=sd(Dxy,na.rm=T))
  }
  result$combination<-apply(result[,c("selection","model")],1,function(x){paste0(x,collapse="_")})
  if(returnScore=="CI"){
    summ.CI<-summ[order(summ$CI.av,decreasing=T),]
  }else if(returnScore=="Dxy"){
    summ.CI<-summ[order(summ$Dxy.av,decreasing=T),]
  }
  
  combi<-apply(summ.CI[1:nr,1:2],1,function(x){paste0(x,collapse="_")})
  result.sel<- result[result$combination %in% combi,]
  result.sel$combination<-factor(result.sel$combination,levels=combi)
  
  p1<-ggplot(data=result.sel, aes(x=combination, y=CI)) + 
    geom_boxplot(position=position_dodge(0.9))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1))+
    ylab("") + xlab("")+
    ggtitle("CI")

  p2<-ggplot(data=result.sel, aes(x=combination, y=Dxy)) + 
    geom_boxplot(position=position_dodge(0.9))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,0.1))+
    ylab("") + xlab("")+
    ggtitle("Dxy")
 
  grid.arrange(p1, p2, ncol=2)
  
  return(combi)
}



buildEnsemble<-function(geneExpr, geneExpr.raw, signExpr, response, vselmethod, model, size=20){
  # this function downsamples the dataset multiple times and makes a model each time, in the end generating an ensemble of balanced models.
  require(caret)
  # set basic caret options
  prepros<-c("center","scale")
  cont<-trainControl(method="boot",number=10, classProbs = T, returnData=F)

  # Some variable selection methods return the same every time, so no need to run them multiple times
  suppressWarnings({
    if(vselmethod %in% c("coeffVar_gene", "coeffVar_sign","coeffVar_both","uniCox_gene","uniCox_sign","uniCox_both", "multiCox_gene","multiCox_sign","multiCox_both", "DESeq2")){
      vars<-varselect(x=1:ncol(geneExpr), geneExpr=geneExpr, geneExpr.raw=geneExpr.raw, signExpr=signExpr,response=response, method=vselmethod)
      full_set<-data.frame(time=response[,1],status=response[,2],t(rbind(geneExpr.rpkm,signExpr))[,vars[[1]]])
      if(model=="rfsrc"){tryCatch({
        mm<-replicate(size,tune.rfsrc(full_set))}
        , error=function(e){message("ERROR :",conditionMessage(e), "\n")})
      }else if(model=="superpc"){tryCatch({
        mm<-replicate(size,tune.superpc(full_set))}
        , error=function(e){message("ERROR :",conditionMessage(e), "\n")})
      }else{tryCatch({
        mm<-replicate(size,tune.coxnet(full_set))}
        , error=function(e){message("ERROR :",conditionMessage(e), "\n")})
      }

    }else{

      if(model=="rfsrc"){
        mm<-replicate(size,{tryCatch({
          vars<-varselect(x=1:ncol(geneExpr), geneExpr=geneExpr, geneExpr.raw=geneExpr.raw, signExpr=signExpr,response=response, method=vselmethod)
          full_set<-data.frame(time=response[,1],status=response[,2],t(rbind(geneExpr.rpkm,signExpr))[,vars[[1]]])
          mm<-tune.rfsrc(full_set)
          return(mm)}, error=function(e){message("ERROR :",conditionMessage(e), "\n")})
        },simplify=F)
      }else if(model=="superpc"){
        mm<-replicate(size,{tryCatch({
          vars<-varselect(x=1:ncol(geneExpr), geneExpr=geneExpr, geneExpr.raw=geneExpr.raw, signExpr=signExpr,response=response, method=vselmethod)
          full_set<-data.frame(time=response[,1],status=response[,2],t(rbind(geneExpr.rpkm,signExpr))[,vars[[1]]])
          mm<-tune.superpc(full_set)
          return(mm)}, error=function(e){message("ERROR :",conditionMessage(e), "\n")})
        },simplify=F)
      }else{
        mm<-replicate(size,{tryCatch({
          vars<-varselect(x=1:ncol(geneExpr), geneExpr=geneExpr, geneExpr.raw=geneExpr.raw, signExpr=signExpr,response=response, method=vselmethod)
          full_set<-data.frame(time=response[,1],status=response[,2],t(rbind(geneExpr.rpkm,signExpr))[,vars[[1]]])
          mm<-tune.coxnet(full_set)
          return(mm)}, error=function(e){message("ERROR :",conditionMessage(e), "\n")})
        },simplify=F)
      }


    }
  })
  class(mm)<-"modelEnsemble"
  return(mm)
}
# mm.ens<-buildEnsemble(geneExpr=geneExpr.rpkm, geneExpr.raw=geenExpr.raw, signExpr=signExpr, response=response, vselmethod="mRMR_both", model="rfsrc", size=20)

predict.modelEnsemble<-function(model.ensemble,expr.set,type="class"){
  require(Biobase)
  require(randomForestSRC)
  labels=levels(predict(model.ensemble[1],expr.set)[[1]])
  preds<-sapply(labels,function(label){rowMeans(sapply(model.ensemble,function(x){predict(x, expr.set,type="prob")[,label] })) })
  preds.class<-apply(preds,1,function(i){colnames(preds)[which.max(i)] })
  if(type=="prob"){return(preds)}else{return(preds.class)}
}
# predict.SurvModelEnsemble<-function(model.ensemble,expr.set,type="class"){
#   require(Biobase)
#   require(randomForestSRC)
#   preds<-rowMeans(sapply(model.ensemble,function(x){predict(x, data.frame(expr.set))$predicted}))
#   return(preds)
# }
# predict(mm.ens,geneExpr.rpkm,type="class")


pairwise.var.test<-function(values,groups){
  #message("no p-value adjustments are implemented!")
  groups<-factor(groups)
  p.vals<-matrix(nrow=length(levels(groups))-1,ncol=length(levels(groups))-1)
  colnames(p.vals)<-paste0("col",seq(1,length(levels(groups))-1,1))
  row.names(p.vals)<-paste0("row",seq(1,length(levels(groups))-1,1))
  for(i in 1:(length(levels(groups))-1)){
    for(j in i+1:(length(levels(groups)))){
      if(j>length(levels(groups))){break}
      # message(paste0(i,"\t",j))
      x<-values[which(groups==as.character(levels(groups)[i]))]
      y<-values[which(groups==as.character(levels(groups)[j]))]
      p.vals[j-1,i]<-var.test(x,y)$p.value
      colnames(p.vals)[i]<-as.character(levels(groups)[i])
      row.names(p.vals)[j-1]<-as.character(levels(groups)[j])
    }

  }
  return(p.vals)
}

