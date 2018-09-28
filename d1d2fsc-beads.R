library(popcycle)

offset <- function(inst, fsc, d1, d2){
      fsc <- as.numeric(fsc)
      d1<- as.numeric(d1)
      d2 <- as.numeric(d2)
      slopes <- read.csv(system.file("filter", "seaflow_filter_slopes.csv",package='popcycle'))
      slopes[slopes$ins== inst, 'notch.small.D1']
      notch.small.D1 <- slopes[slopes$ins== inst, 'notch.small.D1']
      notch.small.D2 <- slopes[slopes$ins== inst, 'notch.small.D2']
      offset.small.D1 <- round(d1 - notch.small.D1 * fsc)
      offset.small.D2 <- round(d2 - notch.small.D2 * fsc)
      notch.large.D1 <- slopes[slopes$ins== inst, 'notch.large.D1']
      notch.large.D2 <- slopes[slopes$ins== inst, 'notch.large.D2']
      offset.large.D1 <- round(d1 - notch.large.D1 * fsc)
      offset.large.D2 <- round(d2 - notch.large.D2 * fsc)

      return(c(offset.small.D1,offset.small.D2,offset.large.D1,offset.large.D2, notch.small.D1,notch.small.D2, notch.large.D1, notch.large.D2))
      }


create.filter.params <- function(inst, fsc, d1, d2, width=3000, slope.file=NULL) {
  QUANTILES <- c(2.5, 50.0, 97.5)

  # Rename to get correct dataframe headers
  beads.fsc.small <- as.numeric(fsc)
  beads.D1 <- as.numeric(d1)
  beads.D2 <- as.numeric(d2)
  width <- as.numeric(width)

  if (is.null(slope.file)) {
    slope.file <- "https://raw.githubusercontent.com/armbrustlab/seaflow-virtualcore/master/1.bead_calibration/seaflow_filter_slopes.csv"
  }
  slopes <- read.csv(slope.file)

  filter.params <- data.frame()
  headers <- c("quantile", "beads.fsc.small",
                  "beads.D1", "beads.D2", "width",
                  "notch.small.D1", "notch.small.D2",
                  "notch.large.D1", "notch.large.D2",
                  "offset.small.D1", "offset.small.D2",
                  "offset.large.D1", "offset.large.D2")
  for (quant in QUANTILES) {
    if (quant == 2.5) {
      suffix <- "_2.5"
      i <- 1
    } else if (quant == 97.5) {
      suffix <- "_97.5"
      i <- 3
    } else if (quant == 50.0) {
      suffix <- ""
      i <- 2
    }

    # Small particles
    notch.small.D1 <- beads.D1[i]/beads.fsc.small[i]
    notch.small.D2 <- beads.D2[i]/beads.fsc.small[i]
    offset.small.D1 <- offset.small.D2 <- 0

    # Large particles
    notch.large.D1 <- slopes[slopes$ins== inst, paste0('notch.large.D1', suffix)]
    notch.large.D2 <- slopes[slopes$ins== inst, paste0('notch.large.D2', suffix)]
    offset.large.D1 <- round(beads.D1[i] - notch.large.D1 * beads.fsc.small[i])
    offset.large.D2 <- round(beads.D2[i] - notch.large.D2 * beads.fsc.small[i])

      newrow <- data.frame(quant, beads.fsc.small[i],
                         beads.D1[i], beads.D2[i], width,
                         notch.small.D1, notch.small.D2,
                         notch.large.D1, notch.large.D2,
                         offset.small.D1, offset.small.D2,
                         offset.large.D1, offset.large.D2,
                         stringsAsFactors=FALSE)
    names(newrow) <- headers
    filter.params <- rbind(filter.params, newrow)
  }

  return(filter.params)
}




#path to Git repo
setwd("~/Documents/DATA/Codes/seaflow-filter")


################################################
### Copy EVT files from DAT to local machine ###
################################################
# link to EVT folder containing EVT files use with this script: dat://2b03e08a5cc4b4cc430d5a355115c6a1597ca6d6ccfc9778f1fb5b545a600deb

# Path to the raw data (DAT)
path.to.data <- "~/Documents/DATA/Codes/seaflow-filter/seaflow-filter-data/"


####################################
### CREATE concatenated EVT file ###
####################################
cruise.list <- list.files("~/Documents/DATA/Codes/seaflow-sfl/curated/", pattern='.sfl',full.names = F)

for(i in 1:length(cruise.list)){

#i <- 1
  exp <- unlist(list(strsplit(cruise.list[i],"_")))
  if(length(exp) > 2) { cruise <- paste(exp[1],exp[2],sep="_")
  } else if(length(exp) ==2) cruise <- exp[1]
    print(cruise)
  inst <-  sub(".sfl","",exp[length(exp)])
    print(inst)

evt.list <- list.files(path=paste0(path.to.data,cruise), pattern=".gz", recursive=T, full.names=T)
DF <- concatenate.evtopp(evt.list, n=100000, min.fsc = 25000, min.pe =25000, min.chl=25000, transform=F)

# # Check EVT cytograms
  # par(mfrow=c(2,3))
  # plot.cytogram(DF,'D1',"D2")
  # plot.cytogram(DF,'fsc_small',"D1")
  # plot.cytogram(DF,'fsc_small',"D2")
  # plot.cytogram(DF,'fsc_small',"pe")
  # plot.cytogram(DF,'fsc_small',"chl_small")


  write.csv(DF, paste0(cruise,"/concatenated_EVT.csv"), quote=F, row.names=F)
}

################################################################################
### GET D1, D2 and FSC coordinate of inflection point (where 1 µm beads are) ###
################################################################################
cruise.list <- list.files("~/Documents/DATA/Codes/seaflow-sfl/curated/", pattern='.sfl',full.names = F)

  #i <- 32
  exp <- unlist(list(strsplit(cruise.list[i],"_")))
  if(length(exp) > 2) { cruise <- paste(exp[1],exp[2],sep="_")
  } else if(length(exp) ==2) cruise <- exp[1]
    print(cruise)
  inst <-  sub(".sfl","",exp[length(exp)])
    print(inst)

  DF <- read.csv(paste0(cruise,"/concatenated_EVT.csv"))

  # Gates beads to find intersections of the two slopes used for OPP filtration
  ip <- inflection.point(DF)

  # check OPP filtration
  evt.list <- list.files(path=paste0(path.to.data,cruise), pattern=".gz", recursive=T, full.names=T)
  evt <- readSeaflow(evt.list[length(evt.list)/2],transform=F)
  #  evt <- readSeaflow(evt.list[34],transform=F)

  filter.params <- create.filter.params(inst, fsc=ip$fsc, d1=ip$d1, d2=ip$d2, width=5000)
    #filter.params$notch.small.D1 <- filter.params$notch.small.D2 <- 1000
  plot.filter.cytogram(evt, filter.params[2,])

    opp <- filter.notch(evt, filter.params[2,])$opp
    plot.cytogram(opp, "fsc_small", "pe")
    b <- subset(opp, pe > 50000)
    plot.cytogram(b, "fsc_small", "D1")


  # only if satisfied with filter params
  write.csv(data.frame(instrument=inst, cruise, ip), paste0(cruise,"/d1d2fsc.csv"),quote=F, row.names=F)

}


######################################################################
### CONCATENATE d1d2fsc-filterparams.csv from ALL cruises TOGETHER ###
######################################################################
csv.list <- list.files(path=".", pattern="d1d2fsc.csv", recursive=T, full.names=T)
csv.list <- csv.list[-grep("ALL-d1d2fsc", csv.list)]

DF <- NULL
for(file in csv.list){
  print(file)
   df <- read.csv(file)
   DF <- rbind(DF,df)
  }

write.csv(DF,"ALL-d1d2fsc.csv", quote=F, row.names=F)



####################################
### CREATE FILTRATION PARAMETERS ###
####################################

### DOWLOAD FILTER SLOPES
slope.file <- "https://raw.githubusercontent.com/armbrustlab/seaflow-virtualcore/master/1.bead_calibration/seaflow_filter_slopes.csv"

x <- gs_title("SeaFlow\ instrument\ log", verbose = TRUE)
list <- gs_read(x)
DF <- read.csv("ALL-d1d2fsc.csv")
DF <- merge(DF, list[,c("cruise","instrument","year")])


@@@@@@@@@@@@@@@@@@@@@@@

for( c in unique(DF$cruise)){
c <- "SCOPE_1"
df <- subset(DF, cruise  == c)
p <- create.filter.params(inst=unique(df["instrument"]), fsc=df["fsc"], d1=df["d1"], d2=df["d2"], width=5000)
params <- rbind(params, p)

}

@@@@@@@@@@@@@@@@@@@@@@@
params <- t(apply(DF,1, function(DF) create.filter.params(inst=DF["instrument"], fsc=DF["fsc"], d1=DF["d1"], d2=DF["d2"], width=5000)))
params <- do.call(rbind.data.frame, params)
params$cruise <- rep(DF$cruise, each=3)

write.csv(params,"ALL-filterparams.csv", quote=F, row.names=F)






##############################
### PLOT BEADS COORDINATES ###
##############################
library(scales)
cols <- colorRampPalette(c("blue4","royalblue4","deepskyblue3", "seagreen3", "yellow", "orangered2","darkred"))
DF <- DF[order(DF$year),]

df <- as.data.frame(t(apply(DF,1, function(DF) offset(inst=DF["instrument"], fsc=DF["fsc"], d1=DF["d1"], d2=DF["d2"]))))
colnames(df) <- c("offset.small.D1","offset.small.D2","offset.large.D1","offset.large.D2", "notch.small.D1","notch.small.D2", "notch.large.D1", "notch.large.D2")
df$cruise <- as.vector(DF$cruise)

png("ALL-d1d2fsc.png",width=12, height=12, unit='in', res=400)

plot(c(0,2^16), c(0,2^16), pch=NA, xlab="fsc_small", ylab="D1 & D2",las=1, main='Bead coordinates')
  for(i in 1:nrow(df)){
      abline(b=mean(df$notch.small.D1[i],df$notch.small.D2[i]), a=mean(df$offset.small.D1[i],df$notch.small.D2[i]),col='grey',lty=2)
      abline(b=mean(df$notch.large.D1[i],df$notch.large.D2[i]), a=mean(df$offset.large.D1[i],df$notch.large.D2[i]),col='grey',lty=2)
      points(DF[i,"fsc"],DF[i,"d1"],bg=alpha(cols(nrow(df))[i],0.5),pch=21,cex=2)
    }
    legend('topleft',legend=DF$cruise,pch=21, pt.bg=alpha(cols(nrow(df)),0.5), bty='n', ncol=2)

dev.off()

id <- match(c("SCOPE_16","SCOPE_6","MBARI_1","DeepDOM"),DF$cruise)
points(DF[id,"fsc"],DF[id,"d1"],pch=1,cex=3,col=2,lwd=2)
