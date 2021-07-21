###################################################################################################
# Multifrequency Indicator formulas
# calclulate multifrequency indicators from Trenkel and Berger (2013, Ecological Indicators)
#
# jech
#
# setwd("/home/mjech/NOAA_Gdrive/sonarpros/R_programs/Classification")
# source("EK60_MFI.R")

# install package gtools--various R programming Toosl
library(gtools)
library(plot.matrix)

# set the margins a bit larger
opar = par()
par(mar=c(5.1, 4.1, 4.1, 4.1))

### start with a clean slate
rm(list=ls(all=TRUE))

#import/export data
#datapath = "x:/Butterfish"
#datapath = paste0('/home/mjech/NOAA_Gdrive/NOAA_SAIP_Butterfish-Mackerel/',
 #                 'Victoria Work/Exported Data/Archive/HB0905')

datapath = paste0('/Volumes/GeringSSD/HB0905_CSV')
#outpath = "x:/Butterfish"

### read in the data and create the Sv array
dfiles = c('HB0905_20090916_140733-20090916_140829_Sv_18kHz.csv',
           'HB0905_20090916_140733-20090916_140829_Sv_38kHz.csv',
           'HB0905_20090916_140733-20090916_140829_Sv_70kHz.csv',
           'HB0905_20090916_140733-20090916_140829_Sv_120kHz.csv',
           'HB0905_20090916_140733-20090916_140829_Sv_200kHz.csv')
# Frequencies in kHz
# for now these match the data files and in order. in the future, this should
# be done dynamically
f=c(18,38,70,120,200)
#f=c(38,125,200, 455)

# read in the first data set, create the Sv array, then fill with the first
# data set
Svdata = read.csv(file=paste(datapath, dfiles[1], sep='/'), header=TRUE, 
                  sep=',')
nd = nrow(Svdata) # number of depth bins
np = ncol(Svdata)-1 # number of pings
nf = length(dfiles) # number of frequencies
Sv_f = array(0, dim=(c(nd, np, nf)))
for (p in 1:np) {
  Sv_f[,p,1] = Svdata[,p+1]
}

# plot the "echogram"
plot(Sv_f[,,1], breaks=c(-100, -10), col=topo.colors, main=paste(f[1], ' kHz'))

# do the remaining frequencies
for (fq in 2:nf) {
  Svdata = read.csv(file=paste(datapath, dfiles[fq], sep='/'), header=TRUE,
                    sep=',')
  for (p in 1:np) {
    Sv_f[,p,fq] = Svdata[,p+1]
  }
  plot(Sv_f[,,fq], breaks=c(-100, -10), col=topo.colors, main=paste(f[fq], ' kHz'))
}

### Calculate the MFI parameters
# Inverse frequencies
e_f=1/f
for (i in 1:length(f)) {
  print(sprintf('f(kHz)=%d, 1/f=%1.4f', f[i], e_f[i]))
}

# Delta paramater
delta=40

# Distance function
#find the unique combinations of frequency pairs
#combinations(n=5, r=2, v=f, set=F, repeats.allowed=F)
#index of unique combinations of frequency pairs
f_idx = combinations(n=length(f), r=2, set=F, repeats.allowed=F)
d_f = cbind(f_idx, 0)

# calculate the "frequency distance"
for (i in 1:length(f_idx[,1])) {
  d_f[i,3] = 1-exp(-abs(f[f_idx[i,1]]-f[f_idx[i,2]])/delta)
  print(sprintf('f1,f2: %d, %d; d=%1.4f', f[f_idx[i,1]], f[f_idx[i,2]], d_f[i,3]))
}

# Convert Sv to sv (linear)
sv_f = 10^(Sv_f/10)

# there are a number of different ways to scale the data
D_f = sv_f
# scale each frequency to its maximum
#scale_max = 0
# scale to a global maximum
scale_max = 1

if (scale_max == 0) {
  for (fq in 1:nf) {
    D_f[,,fq] = sv_f[,,fq]/max(sv_f[,,fq], na.rm=TRUE)
  }
}
if (scale_max == 1) {
  D_f = sv_f/max(sv_f, na.rm=TRUE)
}


# calculate the MFI values
MFI = array(0, dim=(c(nd, np)))
num = array(0, dim=(c(nd, np)))
den = array(0, dim=(c(nd, np)))
for (i in 1:(nf-1)) {
  for (l in (i+1):nf) {
    #print(sprintf('f1,f2: %d, %d', f[i], f[l]))
    num = num+
      d_f[which(d_f[,1] == i & d_f[,2] == l),3]*D_f[,,i]*D_f[,,l]*e_f[i]*e_f[l]
    den = den+D_f[,,i]*D_f[,,l]*e_f[i]*e_f[l]
  }
}
MFI = ((num/den)-0.4)/0.6

# plot the MFI
plot(MFI, breaks=c(0, 1), col=topo.colors, main='MFI')
plot(MFI, breaks=c(0, 0.4, 0.6, 0.7, 0.8, 1.0), 
     col=c('blue', 'green', 'yellow', 'orange', 'red'), main='MFI')

# end of main
