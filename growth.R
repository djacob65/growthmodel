#------------------
# Functions used to model and display graphs for each species
#------------------

# Parameters
# * setName & Yname : set respectively the subset of data and the variable to be modelled. 
#                     See for instance metadata (https://pmb-bordeaux.fr/getdata/query/FR17AP003/metadata/?format=xml) 
#                     corresponding to the _'FR17AP003'_ dataset (Apple) managed by ODAM.
# * dataset   : an ODAM name of a dataset corresponding to a species (e.g 'FR17AP003')
# * VmaxSet   : a set of Vmax values (boundary of the slope) generating as many optimized models as there are values. 
#               A Vmax equal to zero corresponds to no bounds on the slope.
# * fneg      : a boolean indicating if the second sigmoid can have a negative sign (f parameter) 
#               i.e. when there is shrivelling of the fruit at the end of growth
# * noptim    : Number of models to calculate in order to choose the best according to the criterion
# * criterion : Criterion for choosing the final model either as the best fit (BF) or as the most probable fit (MPF)

compute_growth <- function(dataset, setName, Yname, VmaxSet, noptim=200, normalized=TRUE, fneg=FALSE, criterion="BF")
{
   cat(dataset," Getting ... ")
   dat <- get_data_from_ODAM(dataset, setName, Yname)
   cat(" OK\n\n")
    
   NV <- length(VmaxSet)
   cat("--------------------------------\n")
   cat("Model 1:\n")
   cat(dataset," ")
 # Parameter optimization
   opts <- list( model=1, noptim=noptim, vmax=0, fneg=fneg, criterion=criterion)
   fitObj1 <- fitSigmoid(dat, opts, info=TRUE)
   obj1 <- fitObj1 
 # Residuals / Parameters
   res <- obj1$data$ymodel-obj1$data$y
   if (normalized) {
       res <- (res - min(obj1$data$y))/(max(obj1$data$y)-min(obj1$data$y))
       params <- obj1$normalized_params
   } else {
       params <- obj1$params
   }
   dtrunc <- function(x,p) { trunc(x*10^p)/10^p }
 # output
   PR2 <- round(100*length(which( obj1$R2vec>dtrunc(obj1$R2,r2acc)))/noptim,3)
   M1 <- c(round(obj1$R2,4), round(obj1$adjR2,4), PR2, 100, '-', round(mean(res),6), round(sd(res),6), round(params,6) )
   names(M1) <- c('R2','Adj.R2','P(R2)%','S1final%','S2final%','Avg.Res','SD.Res', 'a','b','c','d')
   M1 <- as.matrix(t(M1))
   print_fittedParams(obj1,dataset)

   cat("--------------------------------\n")
   cat("Model 2:\n")
   M2 <- NULL
   fitObj2 <- as.list(rep(0,NV))
   for (k in 1:NV) {
      cat(dataset,"- Vmax =",VmaxSet[k],"- ")
    # Parameter optimization
      opts <- list( model=2, noptim=noptim, vmax=VmaxSet[k], fneg=fneg, criterion=criterion)
      fitObj2[[k]] <- fitSigmoid(dat, opts, info=TRUE)
      obj2 <- fitObj2[[k]]
    # Residuals / Parameters
      res <- obj2$data$ymodel-obj2$data$y
      if (normalized) {
          res <- (res - min(obj2$data$y))/(max(obj2$data$y)-min(obj2$data$y))
          params <- obj2$normalized_params
      } else { 
          params <- obj2$params
      }
    # Final Value of Models
      T <- params[1] + params[5]    # a+f
      A <- round(100*params[1]/T,2) # a/a+f
      F <- round(100*params[5]/T,2) # f/a+f
    # output
      PR2 <- round(100*length(which( obj2$R2vec>dtrunc(obj2$R2,r2acc)))/noptim,3)
      M2 <- rbind(M2, c(VmaxSet[k], round(obj2$R2,4), round(obj2$adjR2,4), PR2, A, F, round(mean(res),6), round(sd(res),6), round(params,6) ) )
      print_fittedParams(obj2,dataset)
   }
   colnames(M2) <- c('Vmax','R2','Adj.R2','P(R2)%','S1final%','S2final%','Avg.Res','SD.Res', 'a','b','c','d','f','g','h')

   list(fitObj1=fitObj1, infos1=M1, fitObj2=fitObj2, infos2=M2)
}

# Parameters
# * fitObj1 & fitObj2 : fitted models (model 1 & 2)  computed by the compute_growth routine 
# * VmaxSet   : a set of Vmax values (boundary of the slope) generating as many optimized models as there are values. 
#               A Vmax equal to zero corresponds to no bounds on the slope.
# * overlay   : a boolean indicating if graphs have to be overlaid or stacked

plot_growth <- function(fitObj1, fitObj2, VmaxSet, label='Dataset', overlay=FALSE)
{
   NV <- length(VmaxSet)
   par(mar=c(2.5,2.5,2.5,1.5))
   par(cex.main=2)
   obj1 <- fitObj1
   if (overlay) {
      par(mfrow=c(2,1))
      yintersp <- 0.25
      colors <- c("red", "cornflowerblue", "darkgreen", "blueviolet", "orange", 
                "magenta", "darkred", "coral", "mediumvioletred", "yellow4",
                "seagreen2", "lightskyblue", "darkcyan", "yellowgreen", "limegreen",
                "wheat2", "yellow4",  "violetred1", "darkorange", "cyan4")   
      legText <- c('Model 1')
      for (k in 1:NV) 
          legText <- c(legText,  gsub("Vmax=0", "No bounds", paste0('Model 2 - Vmax=',VmaxSet[k])))
      plot_fittedCurve(obj1, title=paste0(label," - Growth"))
      for (k in 1:NV) 
          add_fittedCurve(fitObj2[[k]], color=colors[k])
      legend("topleft", legend=legText, col=c("magenta", colors), 
             lty=1, bty="n", y.intersp=yintersp, pt.cex=1, cex=1.5)
      plot_RGRCurve(fitObj1, title=paste0(label," - RGR"))
      for (k in 1:NV) 
          add_RGRCurve(fitObj2[[k]], color=colors[k])
      legend("topright", legend=legText, col=c("magenta", colors), 
             lty=1, bty="n", y.intersp=yintersp, text.width = max(obj1$data$x)/6, pt.cex=1, cex=1.5)
   } else {   
      par(mfrow=c(NV,2))
      for (k in 1:NV) {
         obj2 <- fitObj2[[k]]
         vmax <- VmaxSet[k]
         title <- gsub("Vmax=0", "No bounds", paste0(label,' - Vmax=',vmax))
         plot_fittedCurve(obj1, title=title)
         add_fittedCurve(obj2)
         legend("topleft", legend=c("Model 1", "Model 2"), col=c("magenta", "cornflowerblue"), 
                lty=1, bty="n", y.intersp=0.5, pt.cex=1, cex=1.5)
         plot_RGRCurve(obj1, title=paste0(label," - RGR"))
         add_RGRCurve(obj2)
         legend("topright", legend=c("Model 1", "Model 2"), col=c("magenta", "cornflowerblue"), 
                lty=1, bty="n", y.intersp=0.5, text.width = max(obj1$data$x)/6, pt.cex=1, cex=1.5)
      }
   }
}

# Parameters
# * fitObj1 & fitObj2 : fitted models (model 1 & 2)  by the compute_growth routine 
# * VmaxSet   : a set of Vmax values (boundary of the slope) generating as many optimized models as there are values. 
# * nbval     : number of values of R-square to display in descending order

plot_R2_hist <- function(fitObj2, VmaxSet, nbval=25, label='Dataset')
{
    NV <- length(VmaxSet)
    par(mar=c(2.5,2.5,2.5,1.5))
    par(cex.main=2)
    par(mfrow=c(NV,1))
    for (k in 1:NV) {
       R2vec <-m$fitObj2[[k]]$R2vec
       V <- sort(round(R2vec,r2acc))
       U <- sort(unique(V), decreasing = TRUE)[1:nbval]
       H <- sapply(U,function(x){length(which(V==x))})
       barplot(rbind(U,H), names.arg=U,main=paste0("Distribution of R2 values - Vmax = ",VmaxSet[k]),
               xlab="R2",ylab="Count",border="blue", col="blue",density=10)
    }
}

#------------------
# Functions used to model and display graphs for all species
#------------------

compute_growth_set <- function(dataSets, setName, Yname, VmaxSet, fnegSet, model=2, noptim=200)
{
   NV <- length(dataSets)
   cat("Model",model,":\n")
   fitObj <- as.list(rep(0,NV))
   for (k in 1:NV) {
      dataset <- dataSets[k]
      cat(dataset,"- Fneg =",fnegSet[k],"- Vmax =",VmaxSet[k],"- Getting ... ")
      dat <- get_data_from_ODAM(dataset, setName, Yname)
      opts <- list( model=model, noptim=noptim, vmax=VmaxSet[k], fneg=fnegSet[k], criterion="BF")
      fitObj[[k]] <- fitSigmoid(dat, opts, info=TRUE)
      fitObj[[k]]$vmax <- VmaxSet[k]
      fitObj[[k]]$fneg <- fnegSet[k]       
   }
   fitObj
}

plot_growth_set <- function(fitObj, legend=dataNames)
{
   colors <- c("red", "cornflowerblue", "darkgreen", "blueviolet", "orange", 
             "magenta", "darkred", "coral", "mediumvioletred", "yellow4",
             "seagreen2", "lightskyblue", "darkcyan", "yellowgreen", "limegreen",
             "wheat2", "yellow4",  "violetred1", "darkorange", "cyan4")
   plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(0, 1))
   
   X <- seq(0,1,0.05)
   NV <- length(fitObj)
   for (k in 1:NV) {
      obj <- fitObj[[k]]
      Y <- obj$fsig(X,obj$normalized_params)
      Y <- Y/max(Y)
    # Increase the data size (3 fold) by data interpolating using the cubic spline method
      m <- spline(X, Y)
      lines(m$x, m$y, col=colors[k])
    }
    legend("topleft", legend=legend, col=colors[1:NV], bty="n", lty=1, lwd=2, cex=0.9)
}

plot_RGR_set <- function(fitObj, legend=dataNames)
{
   colors <- c("red", "cornflowerblue", "darkgreen", "blueviolet", "orange", 
             "magenta", "darkred", "coral", "mediumvioletred", "yellow4",
             "seagreen2", "lightskyblue", "darkcyan", "yellowgreen", "limegreen",
             "wheat2", "yellow4",  "violetred1", "darkorange", "cyan4")
   plot(1, type = "n", xlab = "", ylab = "", xlim = c(0, 1), ylim = c(-0.2, 1.1))
   abline(h=0, lty=3)
   X <- seq(0,1,0.05)
   NV <- length(fitObj)
   for (k in 1:NV) {
      obj <- fitObj[[k]]
      Y <- obj$fsig(X,obj$normalized_params)
      RGR <- obj$dfsig(X,obj$normalized_params)/Y
      RGR <- RGR/max(RGR)
    # Increase the data size (3 fold) by data interpolating using the cubic spline method
      m <- spline(X, RGR)
      lines(m$x, m$y, col=colors[k])
    }
    legend("topright", legend=legend, col=colors[1:NV], bty="n", lty=1, lwd=2, cex=0.9)
}
