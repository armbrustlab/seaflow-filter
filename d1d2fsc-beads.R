library(popcycle)
library(aws.s3)
library(googlesheets)



setwd("~/Documents/DATA/Codes/seaflow-filter")


################################################
### Copy EVT files from DAT to local machine ###
################################################
# link to EVT folder containing EVT files use with this script:

dat://2b03e08a5cc4b4cc430d5a355115c6a1597ca6d6ccfc9778f1fb5b545a600deb


####################################
### CREATE concatenated EVT file ###
####################################

x <- gs_title("SeaFlow\ instrument\ log", verbose = TRUE)
list <- gs_read(x)
cruise.list <- list$cruise


for(cruise in cruise.list){
  print(cruise)
  evt.list <- list.files(path=paste0("evt/",cruise), pattern=".gz", recursive=T, full.names=T)
  DF <- concatenate.evtopp(evt.list, n=100000, min.fsc = 2, min.pe =5, min.chl=0, transform=T)
  write.csv(DF, paste0(cruise,"/concatenated_EVT.csv"), quote=F, row.names=F)

}


# Check EVT cytograms
  DF <- read.csv(paste0(cruise,"/concatenated_EVT.csv"))
  par(mfrow=c(2,3))
  plot.cytogram(DF,'D1',"D2")
  plot.cytogram(DF,'fsc_small',"D1")
  plot.cytogram(DF,'fsc_small',"D2")
  plot.cytogram(DF,'fsc_small',"pe")
  plot.cytogram(DF,'fsc_small',"chl_small")



#########################################################
### GET D1, D2 and FSC coordinate of inflection point ###
#########################################################
cruise <- cruise.list[49]
  print(cruise)
  DF <- read.csv(paste0(cruise,"/concatenated_EVT.csv"))

  # Gates beads to find intersections of the two slopes used for OPP filtration
  ip <- inflection.point(DF)

  write.csv(data.frame(cruise, ip), paste0(cruise,"/d1d2fsc-filterparams.csv"),quote=F, row.names=F)





######################################################################
### CONCATENATE d1d2fsc-filterparams.csv from ALL cruises TOGETHER ###
######################################################################
csv.list <- list.files(path=".", pattern="d1d2fsc-filterparams.csv", recursive=T, full.names=T)

DF <- NULL
for(file in csv.list){
  print(file)
   df <- read.csv(file)
   DF <- rbind(DF,df)
  }

write.csv(DF,"ALL-filterparams.csv", quote=F, row.names=F)

png("ALL-filterparams.png",width=20, height=30, unit='in', res=100)

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
par(pty='m', oma=c(0,5,0,0),cex=1.4)
barplot(DF$fsc, names.arg=DF$cruise, horiz=T, las=1, main="fsc")
barplot(DF$d1, names.arg=DF$cruise, horiz=T, las=1, main="D1 & D2")
barplot(DF$d2, horiz=T, las=1, col='lightgreen',density=50, add=T)
#plot(DF$cruise, 0.5*c(DF$d2+DF$d1)/DF$fsc, ylab="D/FSC")
barplot(0.5*c(DF$d2+DF$d1)/DF$fsc, names.arg=DF$cruise, horiz=T, las=1, col='lightblue',density=50, main="D/FSC")

dev.off()
