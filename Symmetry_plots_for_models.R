#Model 1 - single pulse before EU_AS split
#f1 = 0.025, f2/f3/f4 = 0
model_1 <- read.table("/Users/fernandovillanea/Documents/Neanderthal_admix/Neanderthal_admix/outfile_map_chr1_1f.txt",sep="\t",header=TRUE)
EU1 = model_1$frequency_EU
AS1 = model_1$frequency_AS
mean(EU1)/170
mean(AS1)/394

#2D heatgraph matrix
EU_AS1 = matrix(0,nrow=170+1,ncol=394+1) # +1 to include fixed alleles
for (i in 1:nrow(model_1)) {
  EU1_freq = model_1[i,2]
  AS1_freq = model_1[i,3]
  EU_AS1[(EU1_freq+1), (AS1_freq+1)] = EU_AS1[(EU1_freq+1),(AS1_freq+1)]+1 #+1 is because R is one offset
}

#Project Down
projectDown = function(d,m) {
  n = length(d)-1
  l = 0:n
  res = numeric(m+1)
  for (i in 0:m) {
    res[i+1] = sum(d*exp(lchoose(l,i)+lchoose(n-l,m-i)-lchoose(n,m)))
  }
  res
}

EU_AS_d1 <-matrix(0,nrow=101,ncol=394)
for (i in 1:394){
  EU_AS_d1[,i] = projectDown(EU_AS1[,i],100)
}

EU_AS_PD1 <-matrix(0,ncol=101,nrow=101)
for (i in 1:101){
  EU_AS_PD1[i,] = projectDown(EU_AS_d1[i,],100)
}
EU_AS_PD1[1,1] = 0

#Mask counts less than 1
EU_AS_PD1_masked = EU_AS_PD1
EU_AS_PD1_masked[EU_AS_PD1_masked<.5] = 0
log_EU_AS_PD1_masked = log(EU_AS_PD1_masked)
image(log_EU_AS_PD1_masked, xlim = c(0, 1),
      ylim = c(0, 1),xlab = "EUR", ylab = "ASIA", main = "Freq. of archaic human introgression (1 pulse model)")

#Symmetry statistic
sym_stat1 = c()
for (i in 1:101){
  sym_stat1 = c(sym_stat1, sum((EU_AS_PD1[i,] - EU_AS_PD1[,i]))/sum((EU_AS_PD1[i,] + EU_AS_PD1[,i]+1)))
}

plot(sym_stat1, main = "1 pulse model", ylim = c(-0.2, 0.2))
abline(h=0,lty=2)
legend("topright", legend=c("EU +", "AS -"))

#Model 2 - two pulses into AS
#f1 = 0.02, f2 = 0.002, f3/f4 = 0
model_2 <- read.table("/Users/fernandovillanea/Documents/Neanderthal_admix/Neanderthal_admix/outfile_all_2f.txt",sep="\t",header=TRUE)
EU2 = model_2$frequency_EU
AS2 = model_2$frequency_AS
mean(EU2)/170
mean(AS2)/394

#2D heatgraph matrix
EU_AS2 = matrix(0,nrow=170+1,ncol=394+1) # +1 to include fixed alleles
for (i in 1:nrow(model_2)) {
  EU2_freq = model_2[i,2]
  AS2_freq = model_2[i,3]
  EU_AS2[(EU2_freq+1), (AS2_freq+1)] = EU_AS2[(EU2_freq+1),(AS2_freq+1)]+1 #+1 is because R is one offset
}

#Project Down
projectDown = function(d,m) {
  n = length(d)-1
  l = 0:n
  res = numeric(m+1)
  for (i in 0:m) {
    res[i+1] = sum(d*exp(lchoose(l,i)+lchoose(n-l,m-i)-lchoose(n,m)))
  }
  res
}

EU_AS_d2 <-matrix(0,nrow=101,ncol=394)
for (i in 1:394){
  EU_AS_d2[,i] = projectDown(EU_AS2[,i],100)
}

EU_AS_PD2 <-matrix(0,ncol=101,nrow=101)
for (i in 1:101){
  EU_AS_PD2[i,] = projectDown(EU_AS_d2[i,],100)
}
EU_AS_PD2[1,1] = 0

#Mask counts less than 1
EU_AS_PD2_masked = EU_AS_PD2
EU_AS_PD2_masked[EU_AS_PD2_masked<.5] = 0
log_EU_AS_PD2_masked = log(EU_AS_PD2_masked)
image(log_EU_AS_PD2_masked, xlim = c(0, 1),
      ylim = c(0, 1),xlab = "EUR", ylab = "ASIA", main = "Freq. of archaic human introgression (2 pulse model)")

#Symmetry statistic
sym_stat2 = c()
for (i in 1:101){
  sym_stat2 = c(sym_stat2, sum((EU_AS_PD2[i,] - EU_AS_PD2[,i]))/sum((EU_AS_PD2[i,] + EU_AS_PD2[,i]+1)))
}
plot(sym_stat2, main = "2 pulse model", ylim = c(-0.2, 0.2))
abline(h=0,lty=2)
legend("topright", legend=c("EU +", "AS -"))

#Model 3 - One pulse into EU_AS and dilution into EU 
#f1 = 0.025, f2/f3 = 0, f4 = 0.10
model_3 <- read.table("/Users/fernandovillanea/Documents/Neanderthal_admix/Neanderthal_admix/outfile_map_wholegen_dil.txt",sep="\t",header=TRUE)
EU3 = model_3$frequency_EU
AS3 = model_3$frequency_AS
mean(EU3)/170
mean(AS3)/394

#2D heatgraph matrix
EU_AS3 = matrix(0,nrow=170+1,ncol=394+1) # +1 to include fixed alleles
for (i in 1:nrow(model_3)) {
  EU3_freq = model_3[i,2]
  AS3_freq = model_3[i,3]
  EU_AS3[(EU3_freq+1), (AS3_freq+1)] = EU_AS3[(EU3_freq+1),(AS3_freq+1)]+1 #+1 is because R is one offset
}

#Project Down
projectDown = function(d,m) {
  n = length(d)-1
  l = 0:n
  res = numeric(m+1)
  for (i in 0:m) {
    res[i+1] = sum(d*exp(lchoose(l,i)+lchoose(n-l,m-i)-lchoose(n,m)))
  }
  res
}

EU_AS_d3 <-matrix(0,nrow=101,ncol=394)
for (i in 1:394){
  EU_AS_d3[,i] = projectDown(EU_AS3[,i],100)
}

EU_AS_PD3 <-matrix(0,ncol=101,nrow=101)
for (i in 1:101){
  EU_AS_PD3[i,] = projectDown(EU_AS_d3[i,],100)
}
EU_AS_PD3[1,1] = 0

#Mask counts less than 1
EU_AS_PD3_masked = EU_AS_PD3
EU_AS_PD3_masked[EU_AS_PD3_masked<.5] = 0
log_EU_AS_PD3_masked = log(EU_AS_PD3_masked)
image(log_EU_AS_PD3_masked, xlim = c(0, 1),
      ylim = c(0, 1),xlab = "EUR", ylab = "ASIA", main = "Freq. of archaic human introgression (dilution model)")

#Symmetry statistic
sym_stat3 = c()
for (i in 1:101){
  sym_stat3 = c(sym_stat3, sum((EU_AS_PD3[i,] - EU_AS_PD3[,i]))/sum((EU_AS_PD3[i,] + EU_AS_PD3[,i]+1)))
}
plot(sym_stat3, main = "dilution model", ylim = c(-0.2, 0.2))
abline(h=0,lty=2)
legend("topright", legend=c("EU +", "AS -"))

#ONE PLOT TO RULE THEM ALL
par(mfcol=c(2,2))

#empirical
plot(sym_stat, main = "Empirical data (PD_100)", ylim = c(-0.2, 0.2))
abline(h=0,lty=2)
legend("topright", legend=c("EU +", "AS -"))
#1 pulse
plot(sym_stat1, main = "1 pulse model", ylim = c(-0.2, 0.2))
abline(h=0,lty=2)
legend("topright", legend=c("EU +", "AS -"))
#2 pulse
plot(sym_stat2, main = "2 pulse model", ylim = c(-0.2, 0.2))
abline(h=0,lty=2)
legend("topright", legend=c("EU +", "AS -"))
#dilution
plot(sym_stat3, main = "dilution model", ylim = c(-0.2, 0.2))
abline(h=0,lty=2)
legend("topright", legend=c("EU +", "AS -"))