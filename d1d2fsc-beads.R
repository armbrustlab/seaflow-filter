library(popcycle)

#' Plot helpful cytograms for estimating the D1, D2 and FSC coordinates of the inflection point (corresponds to location of 1µm beads).
#'
#' @param dataframe containing EVT data.
#' @return D1, D2 and fsc values of presumed 1 µm beads
#' @export
inflection.point <- function(DF){
  QUANTILES <- c(2.5, 50.0, 97.5)

  def.par <- par(no.readonly = TRUE) # save default, for resetting...
  par(mfrow=c(1,3),pty='s')

  plot_cyt(DF, "fsc_small", "pe")

    poly.beads <- splancs::getpoly(quiet=TRUE)
    b <- subset(DF,splancs::inout(DF[,c("fsc_small", "pe")],poly=poly.beads, bound=TRUE, quiet=TRUE))

  plot_cyt(b, "fsc_small", "D1")
      polyd1 <- splancs::getpoly(quiet=TRUE)
      opp.d1 <- subset(b,splancs::inout(b[,c("fsc_small", "D1")],poly=polyd1, bound=TRUE, quiet=TRUE))

  plot_cyt(b, "fsc_small", "D2")
      polyd2 <- splancs::getpoly(quiet=TRUE)
      opp.d2 <- subset(b,splancs::inout(b[,c("fsc_small", "D2")],poly=polyd2, bound=TRUE, quiet=TRUE))
  par(def.par)

      FSC <- round(summary(c(opp.d1$fsc_small, opp.d2$fsc_small)))
      D1 <- round(summary(opp.d1$D1))
      D2 <- round(summary(opp.d2$D2))

      inflection <- data.frame()
      for (quant in QUANTILES) {
        if (quant == 2.5) { i <- 2; j <- 5
        } else if (quant == 50.0) { i <- j <- 3
        } else if (quant == 97.5) { i <- 5; j <- 2
          }
        fsc <- as.vector(FSC[i])
        d1 <- as.vector(D1[j])
        d2 <- as.vector(D2[j])
        newrow <- data.frame(quantile=quant, fsc, d1, d2, stringsAsFactors=FALSE)
        inflection <- rbind(inflection, newrow)
        }

  return(inflection)
}


create.filter.params <- function(inst, fsc, d1, d2, min.d1, min.d2, width=3000, slope.file=NULL) {
  QUANTILES <- c(2.5, 50.0, 97.5)

  # Rename to get correct dataframe headers
  beads.fsc.small <- as.numeric(fsc)
  beads.D1 <- as.numeric(d1)
  beads.D2 <- as.numeric(d2)
  min.D1 <- as.numeric(min.d1)
  min.D2 <- as.numeric(min.d2)

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
    offset.small.D1 <- min.D1
    offset.small.D2 <- min.D2
    notch.small.D1 <- round((beads.D1[i]-min.D1)/beads.fsc.small[i],3)
    notch.small.D2 <- round((beads.D2[i]-min.D2)/beads.fsc.small[i],3)

    # Large particles
    notch.large.D1 <- round(slopes[slopes$ins== inst, paste0('notch.large.D1', suffix)],3)
    notch.large.D2 <- round(slopes[slopes$ins== inst, paste0('notch.large.D2', suffix)],3)
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

#i <- 6
  exp <- unlist(list(strsplit(cruise.list[i],"_")))
  if(length(exp) > 2) { cruise <- paste(exp[1],exp[2],sep="_")
  } else if(length(exp) ==2) cruise <- exp[1]
    print(cruise)
  inst <-  sub(".sfl","",exp[length(exp)])
    print(inst)

evt.list <- list.files(path=paste0(path.to.data,cruise), pattern=".gz", recursive=T, full.names=T)
DF <- concatenate.evtopp(evt.list, n=100000, min.fsc = 0, min.pe =25000, min.chl=25000, transform=F)

#Check EVT cytograms
  plot_cytogram(DF,'D1',"D2", bins=300, transform=F)
  plot_cytogram(DF,'fsc_small',"D1", transform=F)
  plot_cytogram(DF,'fsc_small',"D2", transform=F)
  plot_cytogram(DF,'fsc_small',"pe", transform=F)
  plot_cytogram(DF,'fsc_small',"chl_small", transform=F)


  write.csv(DF, paste0(cruise,"/concatenated_EVT.csv"), quote=F, row.names=F)
}

################################################################################
### GET D1, D2 and FSC coordinate of inflection point (where 1 µm beads are) ###
################################################################################
cruise.list <- list.files("~/Documents/DATA/Codes/seaflow-sfl/curated/", pattern='.sfl',full.names = F)

  i <- 6
  exp <- unlist(list(strsplit(cruise.list[i],"_")))
  if(length(exp) > 2) { cruise <- paste(exp[1],exp[2],sep="_")
  } else if(length(exp) ==2) cruise <- exp[1]
    print(cruise)
  inst <-  sub(".sfl","",exp[length(exp)])
    print(inst)

  evt.list <- list.files(path=paste0(path.to.data,cruise), pattern="00", recursive=T, full.names=T)

  DF <- read.csv(paste0(cruise,"/concatenated_EVT.csv"))

  # Gates beads to find intersections of the two slopes used for OPP filtration
  ip <- inflection.point(DF)

  filter.params <- create.filter.params(inst, fsc=ip$fsc, d1=ip$d1, d2=ip$d2, min.d1 =min(DF$D1), min.d2 = min(DF$D2), width=5000)
    #filter.params$notch.small.D1 <- filter.params$notch.small.D2 <- 1000

  # check OPP filtration
  evt <- readSeaflow(evt.list[length(evt.list)/2],transform=F)
  #  evt <- readSeaflow(evt.list[1],transform=F)
  plot_filter_cytogram(evt, filter.params[2,])

    opp <- filter.notch(evt, filter.params[1,])$opp
    plot_cyt(opp, "fsc_small", "pe")
    b <- subset(opp, pe > 40000)
    plot_cyt(b, "fsc_small", "D2")
    plot_cyt(b, "fsc_small", "D2")


  # only if satisfied with filter params
  write.csv(data.frame(instrument=inst, cruise, filter.params), paste0(cruise,"/filterparams.csv"),quote=F, row.names=F)



### Combine results
csv.list <- list.files(path=".", pattern="filterparams.csv", recursive=T, full.names=T)
csv.list <- csv.list[-grep("ALL-filterparams.csv", csv.list)]
DF <- do.call(rbind, lapply(csv.list, function(x) read_csv(x)))

write.csv(DF,"ALL-filterparams.csv", quote=F, row.names=F)







##############################
### PLOT BEADS COORDINATES ###
##############################
library(scales)
library(plotrix)
library(viridis)

DF <- read_csv("ALL-filterparams.csv")

# Get official cruise ID
seaflow.meta <- gs_read(gs_title("SeaFlow\ instrument\ log", verbose = FALSE))
id <- match(DF$cruise,seaflow.meta$cruise)
DF$cruise <- seaflow.meta$"Cruise ID"[id]
DF$time <- as.POSIXlt(paste(1,seaflow.meta$"month"[id],seaflow.meta$"year"[id]), format="%d %B %Y")
DF <- DF[order(DF$time),]


# split by quantile
DF1 <- subset(DF, quantile == 50.0)
DF2 <- subset(DF, quantile == 2.5)
DF3 <- subset(DF, quantile == 97.5)



slope <- read.csv("https://raw.githubusercontent.com/armbrustlab/seaflow-virtualcore/master/1.bead_calibration/seaflow_filter_slopes.csv")

png("ALL-filterparams.png",width=12, height=12, unit='in', res=400)

plot(DF1$"beads.fsc.small", DF1$"beads.D1", pch=21,cex=2.5, bg=alpha(viridis(nrow(DF1)),0.5), col=1, xlab="fsc_small", ylab="D1 & D2",las=1, main='Bead coordinates', xlim=c(0,2^16), ylim=c(0,2^16))
  segments(x0=DF2$"beads.fsc.small", y0=DF2$"beads.D1",x1=DF1$"beads.fsc.small", y1=DF1$"beads.D1",  col=alpha(viridis(nrow(DF2)),0.5),pch=21,cex=2, lwd=4)
  segments(x0=DF3$"beads.fsc.small", y0=DF3$"beads.D1",x1=DF1$"beads.fsc.small", y1=DF1$"beads.D1",  col=alpha(viridis(nrow(DF3)),0.5),pch=21,cex=2, lwd=4)

  legend('topleft',legend=DF1$cruise,pch=21, pt.bg=alpha(viridis(nrow(DF1)),0.5), bty='n', ncol=2, pt.cex=1.5)
  abline(b=mean(c(slope$notch.small.D1, slope$notch.small.D2)), a=0, lty=2, col='grey',lwd=2)
  abline(b=mean(c(slope$notch.large.D1, slope$notch.large.D2)), a=-44500, lty=2, col='grey',lwd=2)

    draw.circle(44500,29000,2000, lwd=2, border='red3', col=alpha('grey',0.5))
    abline(h=29000, lwd=2, col='red3')
    abline(v=44500, lwd=2, col='red3')

dev.off()
