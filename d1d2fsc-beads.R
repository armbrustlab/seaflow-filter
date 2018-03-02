library(popcycle)
library(aws.s3)
library(googlesheets)



setwd("~/Documents/DATA/Codes/")

################################################
### Copy EVT files from AWS to local machine ###
################################################

system("aws s3 ls s3://armbrustlab.seaflow/")


#list of SEAFLOW cruise from official LOG
 x <- gs_title("SeaFlow\ instrument\ log", verbose = TRUE)
list <- gs_read(x)
cruise.list <- list$cruise


### COPY FROM AWS TO LOCAL
for(cruise in cruise.list){
  print(cruise)

  # Get the list of ALL EVT files from a specifc cruise
  files <- get_bucket(bucket = 'armbrustlab.seaflow', prefix=paste0(cruise,"/"), max=Inf) # if it doesn't work, try opening S3 on a browser

  if(length(files) == 0){
    print(paste0("no EVT files found for cruise ", cruise, ", check that cruise name in official LOG matches the one in AWS"))
    next
    }
    ##########################
    # select every 100 files #
    ##########################
    id <- seq(30, length(files), by=100)

    # COPY evt files from AWS to local
  for(i in id) system(paste0("aws s3 cp s3://armbrustlab.seaflow/",files[i]$Contents$Key,"seaflow-filter",files[i]$Contents$Key))

}



####################################
### CREATE concatenated EVT file ###
####################################

x <- gs_title("SeaFlow\ instrument\ log", verbose = TRUE)
list <- gs_read(x)
cruise.list <- list$cruise


for(cruise in cruise.list){
  print(cruise)
  evt.list <- list.files(path=paste0("~/Desktop/EVT_allcruises/",cruise), pattern=".gz", recursive=T, full.names=T)


DF <- concatenate.evtopp(evt.list, n=100000, min.fsc = 2, min.pe =5, min.chl=0, transform=T)


write.csv(DF, paste0("seaflow-filter",cruise,"/concatenated_EVT.csv"), quote=F, row.names=F)

}


# Check EVT cytograms
  DF <- read.csv(paste0("seaflow-filter",cruise,"/concatenated_EVT.csv"))
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
  DF <- read.csv(paste0("seaflow-filter",cruise,"/concatenated_EVT.csv"))
  ip <- inflection.point(DF)

  write.csv(data.frame(cruise, ip), paste0("seaflow-filter",cruise,"/d1d2fsc-filterparams.csv"),quote=F, row.names=F)





######################################################################
### CONCATENATE d1d2fsc-filterparams.csv from ALL cruises TOGETHER ###
######################################################################
csv.list <- list.files(path="seaflow-filter", pattern="d1d2fsc-filterparams.csv", recursive=T, full.names=T)

DF <- NULL
for(file in csv.list){
  print(file)
   df <- read.csv(file)
   DF <- rbind(DF,df)
  }

write.csv(DF,"seaflow-filter/ALL-filterparams.csv", quote=F, row.names=F)

png("seaflow-filter/ALL-filterparams.png",width=20, height=30, unit='in', res=100)

layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
par(pty='m', oma=c(0,5,0,0),cex=1.4)
barplot(DF$fsc, names.arg=DF$cruise, horiz=T, las=1, main="fsc")
barplot(DF$d1, names.arg=DF$cruise, horiz=T, las=1, main="D1 & D2")
barplot(DF$d2, horiz=T, las=1, col='lightgreen',density=50, add=T)
#plot(DF$cruise, 0.5*c(DF$d2+DF$d1)/DF$fsc, ylab="D/FSC")
barplot(0.5*c(DF$d2+DF$d1)/DF$fsc, names.arg=DF$cruise, horiz=T, las=1, col='lightblue',density=50, main="D/FSC")

dev.off()
