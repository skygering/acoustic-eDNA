###################################################################################################
# Multifrequency Indicator formulas
# calclulate multifrequency indicators from Trenkel and Berger (2013, Ecological Indicators)
#
# jech
#
# setwd("/home/jjech/NOAA_Gdrive/sonarpros/R_programs/Classification")
# source("MFI.R")

# install package gtools--various R programming Toosl
library(gtools)

### start with a clean slate
rm(list=ls(all=TRUE))

#import/export data
datapath = "x:/Butterfish"
outpath = "x:/Butterfish"

# Frequencies in kHz
#f=c(18,38,70,120,200)
f=c(38,125,200, 455)

# Inverse frequencies
e_f=1/f

for (i in 1:length(f)) {
  print(sprintf('f(kHz)=%d, 1/f=%1.4f', f[i], e_f[i]))
}

# Volume backscatter for each frequency
Sv_f=c(-76,-75,-75,-71,-65)

# Convert Sv to sv (linear)
sv_f=10^(Sv_f/10)

# Scale to maximum sv
D_f=sv_f/max(sv_f)

# Delta paramater
delta=40

# Distance function
#find the unique combinations of frequency pairs
#combinations(n=5, r=2, v=f, set=F, repeats.allowed=F)
#index of unique combinations of frequency pairs
f_idx=combinations(n=length(f), r=2, set=F, repeats.allowed=F)
d_f=seq(length(f_idx[,1]))*0

# calculate the "frequency distance"
for (i in 1:length(f_idx[,1])) {
  #cat("f1,f2:",i,",",f[f_idx[i,1]],f[f_idx[i,2]],"\n")
  d_f[i]=1-exp(-abs(f[f_idx[i,1]]-f[f_idx[i,2]])/delta)
  print(sprintf('f1,f2: %d, %d; d=%1.4f', f[f_idx[i,1]], f[f_idx[i,2]], d_f[i]))
}
#d_f=1-exp(-abs(f_i-f_l)/delta)


#Multiplying variables to sum (pairing)
#row1=D_f[1]*D_f[2]*e_f[1]*e_f[2]*d_f[1]
#row2=D_f[1]*D_f[3]*e_f[1]*e_f[3]*d_f[2]
#row3=D_f[1]*D_f[4]*e_f[1]*e_f[4]*d_f[3]
#row4=D_f[1]*D_f[5]*e_f[1]*e_f[5]*d_f[4]
#row5=D_f[2]*D_f[3]*e_f[2]*e_f[3]*d_f[5]
#row6=D_f[2]*D_f[4]*e_f[2]*e_f[4]*d_f[6]
#row7=D_f[2]*D_f[5]*e_f[2]*e_f[5]*d_f[7]
#row8=D_f[3]*D_f[4]*e_f[3]*e_f[4]*d_f[8]
#row9=D_f[3]*D_f[5]*e_f[3]*e_f[5]*d_f[9]
#row10=D_f[4]*D_f[5]*e_f[4]*e_f[5]*d_f[10]


#MF Indicator
#I=(((sum(row1, row2, row3, row4, row5, row6, row7, row8, row9,row10))/(sum(row11, row12, row13, row14, row15, row16, row17,row18, row19, row20))-0.4)/0.6)





