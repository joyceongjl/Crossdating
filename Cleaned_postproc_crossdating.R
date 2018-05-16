#Scripts from lines 7 to 56 deal with post-processing after getting the individual csv files as an output from Adam Rountrey's plugin using Image J.
#If using a different plugin or software to generate increment widths, skip to line 58, noting the format of the files

#Load dplR program and set working directory to the folder that contains the increment width measurements (in csv format).
#Note that this is specific to the output from Adam Rountrey's plugin.
library(dplR)
setwd("D:/Pzonatus Chronology TIF V3/inc measurements")
#to get masterdata file with all increment widths
# inc width files are in D:\Pzonatus Chronology TIF V3\inc measurements

list.files("D:/Pzonatus Chronology TIF V3/inc measurements") # checking of list.files function
# how to get number of objects in folder = number of columns
fnames<-list.files("D:/Pzonatus Chronology TIF V3/inc measurements")
fnames1<-gsub("T", "0", fnames)#replace T with 0, just to help make things tidier in later lines of code
ncols<-length(fnames1)

#set up loop to load all data
masterdata<-matrix(data=NA,nrow=61,ncol=ncols)#nrow=no years
#column names to file names
colnames(masterdata)<-c( sub("\\.csv","",fnames1))
head(masterdata) # data check
rownames(masterdata) <- seq(1950,2010,1)
head(masterdata) #data check

#set up masterlist with all csv files but only column 2 and 5, without first row. For ALL DATA
list<-lapply(fnames, function(x) read.csv(x, header=TRUE)[-1,])
list<-lapply(list, function(z) z[c(2,5)])

#set up loop to load all csv files into masterdata
#column header in list is "IncNum.v1.3c"
head(masterdata)
for (i in 1:length(fnames)) {
  masterdata[,i]<-list[[i]][,2][match(rownames(masterdata),list[[i]][,1])]
}
write.csv(masterdata,"D:/Pzonatus Chronology TIF V3/data files/masterdata_rawinc_20180224.csv") #this is just to save a copy
length(masterdata[,1]) # data check

#save masterdata as ALL, with date. Then run masterdata for adults
alldat<-masterdata

#to remove first 7 years of increment widths, in my case it was to get rid of the juvenile years as those are highly variable.
for (k in 1:dim(alldat)[2]){
  alldat[which(!is.na(alldat[,k]))[c(1:7)],k] <- NA
}
write.csv(alldat,"D:/Pzonatus Chronology TIF V3/data files/masterdata_rawinc_1st7excld_20180224.csv") #saving data

#sometimes the format of the files seem weird, so I usually read in the file again to make sure everything is ok
alldat<-read.csv(file="D:/Pzonatus Chronology TIF V3/data files/masterdata_rawinc_1st7excld_20180224.csv")
str(alldat)
rownames(alldat)<-alldat[,1]
alldat<-alldat[,-1]#to get rid of rownames
head(alldat)
rownames(alldat)#contains years
colnames(alldat)#transects
dat1<-as.data.frame(alldat[11:59,])#because data frame is needed to change data into rwl file and also to get rid of years that only have NAs
str(dat1)#data check

#if starting from this line, ensure that dat1 is a data frame with rownames=years, colnames=individual fish/transects names, eg. PZ00301 for first transect of fish#3
write.rwl(dat1, fname=file.choose())#save as filename.rwl somewhere, make sure '.rwl' is included when saving the file
dat2<-read.rwl(file.choose())#read in the rwl file
rwl.report(dat2)#shows number of series, measurements, avg length, etc.
plot(dat2, plot.type="spag")#quick plot of rwl object

#detrending series to do crossdating on
#note that for adult fish, I have been using a modified negative exponential curve (ModNegExp) to detrend because that fits the normal growth curve of a fish
#there are other detrending methods such as a spline, choose the right one for your species
#we have kept the spline rigidity parameter at 22 years, following Black et al. 2005, as this works for long-lived species <30 years old. This should be tested for each species.
dat2negexp<-detrend(dat2, method="ModNegExp", nyrs=22)#does this make sense if the lowest #yrs i hav is 10?
str(dat2negexp)#data check
colMeans(dat2negexp, na.rm=TRUE)# data check, to ensure the column means are close to 1
rwl.report(dat2negexp)#shows number of series, measurements, avg length, etc.
dat2stats<-rwl.stats(dat2)#shows number of years in each series, mean, median, stdev, etc.
max(dat2stats$year)#maximum age of oldest fish in dataset

#####crossdating using detrended series with segment length of 14, we have historically used 15 year segments lengths for crossdating in COFECHA, but this package requires even lengths
rwl.14.negexp<-corr.rwl.seg(dat2negexp, seg.length=14, n=NULL, bin.floor=10,
                            prewhiten=FALSE, pcrit=0.05, biweight=TRUE)#this comes up with a plot of the good segments and those that might have errors, as well as 6 attributes to help in crossdating
str(rwl.14.negexp)#shows the 6 attributes, see help in package for details
rwl.14.negexp$flags
#checking each series that has low correlations to show where potential error might be
corr.series.seg(rwl=dat2negexp, series="PZ0040", seg.length=14, n=NULL, bin.floor=10,
                prewhiten=FALSE, pcrit=0.05, biweight=TRUE)
#once flagged series and rough years are known, go back to images and check the raw increments. 
#edit original increment width data if there is indeed an error. Then redo whole process of detrending and crossdating.


#additional statistics on final detrended series
#note that for long-lived fishes, a mean interseries correlation of ~0.2 is decent, mean rbar should be positive and eps>0.8 is considered good (although eps>0.5 is sometimes alright)
rwi.stats(dat2negexp, ids=NULL, period="max", method="spearman")#shows total rbar (rbar.tot), eps values (eps), etc.
dat2intser<-interseries.cor(dat2negexp, prewhiten=FALSE, method="spearman")#read about prewhiten and methods to choose the right one
mean(dat2intser[,1])#mean interseries correlation
dat2negexp.stats<-rwi.stats.running(dat2negexp, ids=NULL, n=NULL, prewhiten=FALSE, running.window=TRUE, window.length=15, window.overlap=14, min.corr.overlap=14)# length of window and overlap depends on fish longevity and number of years in the data
str(dat2negexp.stats)#data check
write.csv(dat2negexp.stats, "D:/Pzonatus Chronology TIF V3/data files/dat2negexp_stats_15yr_20180306.csv")
plot(dat2negexp.stats$start.year,dat2negexp.stats$eps, ylim = c(0,1), ylab="EPS", xlab="15yr mid-point")#plot of eps values over running windows

#calculation of master chronology using final series and years that have good rbar and eps values
exprwi<-dat2negexp[10:48,]#selection of data that are considered 'good quality'
head(exprwi)#data check
rwi.stats(exprwi, ids=NULL, period="max", method="spearman")#total rbar and eps values for the 'good quality' series
expintser<-interseries.cor(exprwi, prewhiten=FALSE, method="spearman")
mean(expintser[,1])#mean interseries correlation
expcrn<-chron(exprwi, prefix="exp")#builds a mean value chronology
exp.stats<-rwi.stats.running(exprwi, ids=NULL, n=NULL, prewhiten=FALSE, running.window=TRUE, window.length=15, window.overlap=14, min.corr.overlap=14)
plot(exp.stats$mid.year,exp.stats$eps, ylim = c(0,1), ylab="EPS", xlab="15yr mid-point")
