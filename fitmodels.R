options(width=80)
options(warn=-1)
options(stringsAsFactors=FALSE)

if (!require("minpack.lm", quietly = TRUE))
  install.packages('minpack.lm', repos='http://cran.rstudio.com')

#-----

library(minpack.lm)

#-----
## Growth Modelisation
#-----

# n-digit precision on R-square
r2acc <- 4

# Default modelisation options
modeloptions <- list(
    model=2,        # 1 for model 1 (one sigmoïd) or 2 for model 2 (two sigmoïds)
    noptim=200,     # Number of models to calculate in order to choose the best according to the criterion
    vmax=0,         # Boundary of the slope(s); 0 means no bounds
    fneg=FALSE,     # A boolean indicating if the second sigmoid can have a negative sign (f parameter, model 2 only)
    criterion="BF"  # Choose the final model either as the best fit (BF) or as the most probable fit (MPF)
)

# Model 1 - Single Sigmoidal Model
# f(x) = d + a/(1+exp(-b*(x-c)))
sigmoid1 <- function(x,par)
{
   par[4]+par[1]/(1+exp(-par[2]*(x-par[3])))
}
deriv_sigmoid1 <- function(x,par)
{
   (par[1]*par[2]*exp(-par[2]*(x-par[3])))/(1+exp(-par[2]*(x-par[3])))^2
}

# Model 2 - Sum of two Sigmoidal Model
# f(x) = d + a/(1+exp(-b*(x-c))) + f/(1+exp(-g*(x-h)))
sigmoid2 <- function(x,par)
{
   sigmoid1(x,par[1:4]) + sigmoid1(x,c(par[5:7],0))
}
deriv_sigmoid2 <- function(x,par)
{
   deriv_sigmoid1(x,par[1:3]) + deriv_sigmoid1(x,par[5:7])
}

# Normalisation of Data before modeling
normalizeData <-function(dataInput)
{
   timeData <- dataInput$time
   timeRange <- max(timeData,na.rm = T)
   timeData <- timeData / timeRange

   intensityMin <- min(dataInput$intensity,na.rm = T)
   intensityMax <- max(dataInput$intensity,na.rm = T)
   intensityData <- dataInput$intensity - intensityMin
   intensityRange <- max(intensityData,na.rm = T)
   intensityData <- intensityData / intensityRange

   dataOutput <- data.frame(time = timeData, intensity = intensityData)
   return(list(timeIntensityData = dataOutput,
               dataScalingParameters = c(timeRange = timeRange,
                                          intensityMin = intensityMin,
                                          intensityMax = intensityMax,
                                          intensityRange = intensityRange)))
}

# Fit the sigmoidal model on dataInput (Model 1 & 2)
SigmoidalFitModel <- function(dataInput, opts=modeloptions)
{
    # The one-sigmoidal model (Model 1)
    oneSigmoid <- function(time,a,b,c) { sigmoid1(time, c(a,b,c,0)) }

    # The double-sigmoidal model (Model 2)
    sumSigmoid <- function(time,a,b,c,f,g,h) { sigmoid2(time, c(a,b,c,0,f,g,h)) }

    normalizedInput <- normalizeData(dataInput = dataInput)
    dataFrameInput <- normalizedInput$timeIntensityData

    # Modelisation options
    model <- opts$model
    noptim <- opts$noptim
    vmax <- opts$vmax
    fneg <- opts$fneg
    criterion <- opts$criterion

    # Sample size
    N <- length(dataFrameInput$intensity)

    # Number of parameters
    np <- ifelse(model==1, 3, 6)

    # Init bounds
    lowerBounds <- NULL
    upperBounds <- NULL
    if (model==1) {
       if (vmax>0) { lowerBounds <- c(a=0.3, b=0.01, c=0); upperBounds <- c(a=1, b=vmax, c=0.95); }
        formula <- intensity ~ oneSigmoid(time,a,b,c)
    } else {
       if (vmax>0) {
          lowerBounds <- c(a=0.1, b=1, c=0.1, f=ifelse(fneg, -1, 0), g=1, h=0.2)
          upperBounds <- c(a=10, b=vmax, c=1, f=1, g=vmax, h=1)
       }
       formula <- intensity ~ sumSigmoid(time,a,b,c,f,g,h)
    }

    fitCurve <- list()
    for( idx in 1:noptim ) {

       # Init Start point
       if (model==1) {
          startList <- function() {
               list(a=rnorm(1,1,0.1), b=rnorm(1,5,1), c=rnorm(1,0.5,0.1))}
       } else {
          startList <- function() {
               list(a=rnorm(1,0.5,0.25), b=rnorm(1,10,2), c=rnorm(1,0.5,0.25),
               f=ifelse(fneg, -1, 1)*rnorm(1,0.5,0.25), g=rnorm(1,10,2), h=rnorm(1,0.5,0.25)) }
       }

       # Fitting ...
       repeat {
          theFitResult <- try(minpack.lm::nlsLM(formula, dataFrameInput, 
                              start = startList(), lower = lowerBounds, upper = upperBounds, 
                              control = list(maxiter = 1000, minFactor = 1/2^20), trace=FALSE),
                              silent = TRUE)
          if (class(theFitResult)=="nls") break
       }

       # Calculate the ymodel depending on the model
       m<-as.data.frame(t(theFitResult$m$getPars()))
       o <- normalizedInput$dataScalingParameters
       if (model==1) {
          norm_p <- c( m$a, m$b, m$c, 0 )
          p <- c( m$a*o[4], m$b/o[1], m$c*o[1], o[2] )
          ymodel <- oneSigmoid(dataFrameInput$time, m$a, m$b, m$c)
       } else {
          norm_p <- c( m$a, m$b, m$c, 0, m$f, m$g, m$h )
          p <- c( m$a*o[4], m$b/o[1], m$c*o[1], o[2], m$f*o[4], m$g/o[1], m$h*o[1] )
          ymodel <- sumSigmoid(dataFrameInput$time, m$a, m$b, m$c, m$f, m$g, m$h)
       }

       # Calculate R2 and adj.R2
       SST <- sum((dataFrameInput$intensity - mean(dataFrameInput$intensity))^2)
       SSR <- sum((dataFrameInput$intensity - ymodel)^2)
       R2 <- 1 - SSR/SST
       adjR2 <- 1 - (1-R2)*(N-1)/(N-np)

       # Save fit
       fitCurve[[idx]] <- list(params=p, normalized_params=norm_p, R2=R2, adjR2=adjR2)
    }

    # Get the R2 values
    R2vec <- sapply(1:noptim, function(k){ fitCurve[[k]]$R2 })

    # Retain the best fit (BF)
    if(criterion=="BF") {
       R2.best <- 0
       idx.best <- 1
       for( idx in 1:noptim ) {
           R2.idx <- fitCurve[[idx]]$R2
           if (R2.idx>R2.best) { R2.best=R2.idx; idx.best=idx; }
       }
    }
    # Retain the most probable fit (MPF)
    if (criterion=="MPF") {
       V <- sort(round(R2vec,r2acc))
       U <- sort(unique(V), decreasing = TRUE)
       H <- sapply(U,function(x){length(which(V==x))})
       R2sel <- U[which(H==max(H))]
       R2.best <- 0
       idx.best <- 1
       for( idx in 1:noptim ) {
           R2.idx <- fitCurve[[idx]]$R2
           if (R2.idx>R2.best && round(R2.idx,4)==R2sel) { R2.best=R2.idx; idx.best=idx; }
       }
    }

    # Return the best fitting
    fitCurve[[idx.best]]$R2vec <- R2vec
    fitCurve[[idx.best]]
}

# Fit the Growth curve for dataset datasetID
fitSigmoid <- function(dat, opts=modeloptions, info=FALSE)
{
     dataInput <- data.frame(time=dat$x, intensity=dat$y)
     model <- opts$model

  if (info) cat("Fitting ... (",model,") ")
     fitObj <- SigmoidalFitModel(dataInput,opts)
  if (info) cat("R2 =",fitObj$R2," OK\n")

     if ( model==1 ) {
        fsig <- sigmoid1
        dfsig <- deriv_sigmoid1
     }
     if ( model==2 ) {
        fsig <- sigmoid2
        dfsig <- deriv_sigmoid2
     }
     par1 <- as.numeric(fitObj$params)
     norm_part1 <- as.numeric(fitObj$normalized_params)
     R2 <- fitObj$R2
     adjR2 <- fitObj$adjR2
     R2vec <- fitObj$R2vec
     yest <- fsig(dat$x,par1)
     dyest <- dfsig(dat$x,par1)
     RGR <- dyest/yest
     dat <- data.frame(x=dat$x, y=dat$y, sdev=dat$sdev,
                       ymodel=yest, dymodel=dyest, RGR=RGR)
     list(data=dat, normalized_params=norm_part1, params=par1, model=model, R2=R2, adjR2=adjR2, R2vec=R2vec, fsig=fsig, dfsig=dfsig)
}


#-----
# Plot functions
#-----

plot.with.errorbars <- function(x, y, err, ylim=NULL, xlab=NULL, ylab=NULL, ...)
{
  if (is.null(ylim))
     ylim <- c(min(y-err), max(y+err))
  plot(x, y, ylim=ylim, pch=19, col='blue', xlab=xlab, ylab=ylab, 
             cex.axis=1.5, cex.labs=1.5, cex.main=2, cex.sub=1.5, ...)
  arrows(x, y-err, x, y+err, length=0.05, angle=90, code=3, col='red')
}

plot_fittedCurve <- function(fitObj, title="fittedCurve")
{
   dat    <- fitObj$data
   params <- fitObj$params
   fsig   <- fitObj$fsig
   m <- spline(dat$x, dat$y)
   plot.with.errorbars(dat$x, dat$y, dat$sdev, main=title)
   lines(m$x, fsig(m$x,params), col='magenta')
}

plot_RGRCurve <- function(fitObj, title="RGR")
{
   dat   <- fitObj$data
   fsig  <- fitObj$fsig
   dfsig <- fitObj$dfsig
   par1  <- fitObj$params
   x <- 1:max(dat$x)
   RGR <- dfsig(x,par1)/fsig(x,par1)
   plot( x, RGR/max(RGR), ylim=c(-0.2,1), type="l", col="magenta", main=title, 
            cex.axis=1.5, cex.labs=1.5, cex.main=2, cex.sub=1.5)
}


add_fittedCurve <- function(fitObj, color='cornflowerblue')
{
   dat    <- fitObj$data
   params <- fitObj$params
   fsig   <- fitObj$fsig
   m <- spline(dat$x, dat$y)
   lines(m$x, fsig(m$x,params), col=color, lwd=2)
}

add_RGRCurve <- function(fitObj, color='cornflowerblue')
{
   dat   <- fitObj$data
   fsig  <- fitObj$fsig
   dfsig <- fitObj$dfsig
   par1  <- fitObj$params
   x <- 1:max(dat$x)
   RGR <- dfsig(x,par1)/fsig(x,par1)
   lines( x, RGR/max(RGR), col=color, lwd=2 )
}

print_fittedParams <- function(fitObj, title='Dataset')
{
   par1  <- fitObj$params
   model <- fitObj$model
   if (model==1) {
       names(par1) <- c('a','b','c','d')
       modlabel <- 'Single Sigmoid'
   }
   if (model==2) {
       modlabel <- 'Sum of two Sigmoid'
       names(par1) <- c('a','b','c','d','f','g','h')
   }
   params <- as.data.frame(t(par1))
   cat(title,' - ',modlabel,"\n")
   print(params)
   cat("\n")
}
