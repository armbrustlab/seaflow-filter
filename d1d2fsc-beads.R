library(popcycle)
library(arrow)

# popcycle now contains functions for getting filter parameters.

# Use subsampled EVT 


###################
### Directories ###
###################

## Francois
setwd("~/Documents/Codes/seaflow-filter") #path to Git repo
path.to.data <- "~/Desktop/"
cruise.list <- list.files("~/Documents/Codes/seaflow-sfl/curated/", pattern='.sfl',full.names = F)

## Annette
setwd("~/SeaFlow/Forks/seaflow-filter") #path to Git repo
path.to.data <- "/Users/annettehynes/Library/CloudStorage/GoogleDrive-ahynes@uw.edu/Shared drives/SeaFlow-VCT/snakemake/"
cruise.list <- list.files("~/SeaFlow/Forks/seaflow-sfl/curated/", pattern='.sfl',full.names = F)

print(cruise.list)


####################################
### CREATE concatenated EVT file ###
####################################
i <- 4; print(cruise.list[i])

# metadata
exp <- unlist(list(strsplit(cruise.list[i],"_")))
if(length(exp) > 2) { cruise <- paste(exp[1],exp[2],sep="_")
} else if(length(exp) == 2) cruise <- exp[1]
print(cruise)
inst <-  sub(".sfl","",exp[length(exp)])
print(inst)

# Concatenation
#evt.list <- list.files(path=paste0(path.to.data,cruise), pattern=".gz", recursive=T, full.names=T)
#DF <- concatenate.evt(evt.list,evt.dir=paste0(path.to.data,cruise), n=100000, min.fsc = 0, min.pe =25000, min.chl=25000, transform=F)

evt_sample <- paste0(path.to.data, cruise, "/sample/", cruise, ".beads-sample.parquet")
evt <- arrow::read_parquet(evt_sample)
DF <- dplyr::slice_sample(evt, n = 100000)

#Check EVT cytograms
popcycle::plot_cyt(DF,'fsc_small',"pe", transform=F)

# save
system(paste('mkdir', cruise))
write.csv(DF, paste0(cruise,"/concatenated_EVT.csv"), quote=F, row.names=F)


################################################################################
### GET D1, D2 and FSC coordinate of inflection point (where 1 µm beads are) ###
################################################################################


# Gates beads to find intersections of the two slopes used for OPP filtration
ip <- popcycle::inflection_point(DF)


# create filter parameters
filter.params <- popcycle::create_filter_params(inst, fsc=ip$fsc, d1=ip$d1, d2=ip$d2, 
                                      min_d1 = -5000, min_d2 = -5000, width=15000) 

# check OPP filtration
full_sample <- paste0(path.to.data, cruise, "/sample/", cruise, ".unfiltered-sample.parquet")
evt <- arrow::read_parquet(full_sample)

popcycle::plot_filter_cytogram(evt, filter.params)

# only if satisfied with filter params
write.csv(data.frame(instrument=inst, cruise, filter.params), paste0(cruise,"/filterparams.csv"),quote=F, row.names=F)

### Combine results ###
csv.list <- list.files(path=".", pattern="filterparams.csv", recursive=T, full.names=T)
csv.list <- csv.list[-grep("ALL-filterparams.csv", csv.list)]
DF <- do.call(rbind, lapply(csv.list, function(x) read.csv(x)))

write.csv(DF,"ALL-filterparams.csv", quote=F, row.names=F)

##############################
### PLOT BEADS COORDINATES ###
##############################
library(scales)
library(plotrix)
library(viridis)
library(googlesheets4)
library(tidyverse)

DF <- read_csv("ALL-filterparams.csv")

# Get official cruise ID
seaflow.meta <- googlesheets4::read_sheet("https://docs.google.com/spreadsheets/d/1Tsi7OWIZWfCQJqLDpId2aG_i-8Cp-p63PYjjvDkOtH4")


id <- match(DF$cruise,seaflow.meta$cruise)
DF$cruise <- seaflow.meta$"Cruise ID"[id]
DF$time <- as.POSIXlt(paste("01",seaflow.meta$"Month"[id],seaflow.meta$"Year"[id]), format="%d %B %Y", tz="GMT")
DF <- DF[order(DF$time),]


# split by quantile
DF1 <- subset(DF, quantile == 50.0)
DF2 <- subset(DF, quantile == 2.5)
DF3 <- subset(DF, quantile == 97.5)



slope <- read_csv("https://raw.githubusercontent.com/armbrustlab/seaflow-virtualcore/master/1.bead_calibration/seaflow_filter_slopes.csv")

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

















#################
### RT cruise ###
#################
inst <- 751; cruise <- "SR1917"
path <- "~/Downloads"
evt.list <- list.files(path=path, pattern=".gz", recursive=T, full.names=T)
evt <- readSeaflow(evt.list[1],transform=T)

plot_cytogram(evt, "fsc_small", "D1")
plot_cytogram(evt, "D1", "D2")
plot_cytogram(evt, "fsc_small", "chl_small")

ip <- inflection.point(evt)
filter.params <- create.filter.params(inst, fsc=ip$fsc, d1=ip$d1, d2=ip$d2, min.d1 =4000, min.d2 = 3000, width=5000)


plot_filter_cytogram(evt, filter.params)

  par(mfrow=c(2,2))
    opp <- filter.notch(evt, filter.params)
    plot_cyt(opp, "fsc_small", "chl_small")
    plot_cyt(opp, "fsc_small", "pe"); abline(v=filter.params[,'beads.fsc.small'], lty=2,col=2)
    b <- subset(opp, pe > 40000)
    plot_cyt(b, "fsc_small", "D1"); abline(h=filter.params[,'beads.D1'], lty=2,col=2)
    plot_cyt(b, "fsc_small", "D2"); abline(h=filter.params[,'beads.D2'], lty=2,col=2)


write.csv(data.frame(instrument=inst, cruise, filter.params), paste0("~/Desktop/filterparams.csv"),quote=F, row.names=F)

