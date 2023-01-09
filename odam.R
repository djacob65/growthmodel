options(width=80)
options(warn=-1)
options(stringsAsFactors=FALSE)

URLWS <<- 'https://pmb-bordeaux.fr/getdata/'

if (!require("Rodam", quietly = TRUE)) {
  require(devtools)
  install_github("inrae/Rodam", dependencies = TRUE)
}


#-----

library(Rodam)

#-----
## Get Collection
#-----
get_Collection <- function()
{
  collection <- dh$getDataByName('collection')
  rownames(collection) <- 1:nrow(collection)
  collection
}

#-----
## Get Data from ODAM
#-----

get_data_from_ODAM <- function(dataset, setName, Yname)
{
     options(stringsAsFactors=FALSE)
     Yvar <- Yname
     if (Yname == 'Weight'  && dataset == 'FR17GV001') Yvar <- 'FreshMassFruit'
     if (Yname == 'Weight'  && dataset == 'FR17PP009') Yvar <- 'Weight1'

     dh <- new('odamws', URLWS, dataset)
     invisible(capture.output(
        ds1 <- dh$getSubsetByName(setName)
     ))
     M <- ds1$data[ , c('DayAfterAnthese', Yvar) ]
     M <- M[ ! is.na(M[Yvar]), ]
     M[,Yvar] <- as.numeric(M[,Yvar])

     AVG <- aggregate(. ~ DayAfterAnthese, M, mean)
     SDEV <- aggregate(. ~ DayAfterAnthese, M, sd)

     if (dataset == 'FR17GV001') AVG[,1] <- as.numeric(gsub("D","",AVG[,1]))
     dat <- data.frame(x=AVG[,1], y=AVG[,2], sdev=SDEV[,2])
     dat
}
