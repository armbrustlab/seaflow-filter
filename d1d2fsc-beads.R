library(popcycle)
library(aws.s3)
library(googlesheets)


#path to Git repo
setwd("~/Documents/DATA/Codes/seaflow-filter")


################################################
### Copy EVT files from DAT to local machine ###
################################################
# link to EVT folder containing EVT files use with this script:

dat://2b03e08a5cc4b4cc430d5a355115c6a1597ca6d6ccfc9778f1fb5b545a600deb

#Path to the raw data (DAT)
path.to.data <- "~/Documents/DATA/Codes/seaflow-filter/seaflow-filter-data/"


####################################
### CREATE concatenated EVT file ###
####################################
x <- gs_title("SeaFlow\ instrument\ log", verbose = TRUE)
list <- gs_read(x)
cruise.list <- list$cruise


cruise <- cruise.list[46]
print(cruise)
evt.list <- list.files(path=paste0(path.to.data,cruise), pattern=".gz", recursive=T, full.names=T)
DF <- concatenate.evtopp(evt.list, n=100000, min.fsc = 2, min.pe =25000, min.chl=0, transform=F)
write.csv(DF, paste0(cruise,"/concatenated_EVT.csv"), quote=F, row.names=F)




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

QUANTILES <- c(2.5, 50.0, 97.5) # NEED TO FIX THIS BUG...

cruise <- cruise.list[52]
  print(cruise)
  DF <- read.csv(paste0(cruise,"/concatenated_EVT.csv"))

  # Gates beads to find intersections of the two slopes used for OPP filtration
  ip <- inflection.point(DF)

  # check OPP filtration
  inst <- unlist(list(list[which(list$cruise == cruise),'instrument']))
  filter.params <- create.filter.params(inst, fsc=ip$fsc, d1=ip$d1, d2=ip$d2)[2,]
  evt <- readSeaflow(evt.list[10],transform=F)
  plot.filter.cytogram(evt, filter.params)


  # only if satisfied with params
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

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
par(pty='m', oma=c(0,5,0,0),cex=1.4)
barplot(DF$fsc, names.arg=DF$cruise, horiz=T, las=1, main="fsc")
barplot(DF$d2, names.arg=DF$cruise, horiz=T, las=1, main="D1 & D2")
barplot(DF$d1, horiz=T, las=1, col='lightgreen',density=50, add=T)
barplot(0.5*c(DF$d2+DF$d1)/DF$fsc, names.arg=DF$cruise, horiz=T, las=1, col='lightblue',density=50, main="D/FSC")







########################################
### CHECK OFFSET for small particles ###
########################################
cols <- colorRampPalette(c("blue4","royalblue4","deepskyblue3", "seagreen3", "yellow", "orangered2","darkred"))

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

x <- gs_title("SeaFlow\ instrument\ log", verbose = TRUE)
list <- gs_read(x)
DF <- read.csv("ALL-filterparams.csv")
DF <- merge(DF, list[,c("cruise","instrument","year")])
DF <- DF[order(DF$year),]

df <- as.data.frame(t(apply(DF,1, function(DF) offset(inst=DF["instrument"], fsc=DF["fsc"], d1=DF["d1"], d2=DF["d2"]))))
colnames(df) <- c("offset.small.D1","offset.small.D2","offset.large.D1","offset.large.D2", "notch.small.D1","notch.small.D2", "notch.large.D1", "notch.large.D2")
df$cruise <- as.vector(DF$cruise)

png("ALL-filterparams.png",width=12, height=12, unit='in', res=400)

plot(c(0,2^16), c(0,2^16), pch=NA, xlab="fsc_small", ylab="D1 & D2",las=1, main='Bead coordinates')
  for(i in 1:nrow(df)){
      abline(b=mean(df$notch.small.D1[i],df$notch.small.D2[i]), a=mean(df$offset.small.D1[i],df$notch.small.D2[i]),col='grey',lty=2)
      abline(b=mean(df$notch.large.D1[i],df$notch.large.D2[i]), a=mean(df$offset.large.D1[i],df$notch.large.D2[i]),col='grey',lty=2)
      points(DF[i,"fsc"],DF[i,"d1"],bg=alpha(cols(nrow(df))[i],0.5),pch=21,cex=2)
    }
    legend('topleft',legend=DF$cruise,pch=21, pt.bg=alpha(cols(nrow(df)),0.5), bty='n', ncol=2)

dev.off()
