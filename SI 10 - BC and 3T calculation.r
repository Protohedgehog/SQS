#Read in data
# Obtain data file as occurrence matrix from Fossilworks.
# Perform a standard search, select the generate diversity curves option, and download the 'presences.txt' file

occMatrix<-read.delim("presences.txt")

#Remove genus name and total occurrences columns
occMatrix<-occMatrix[,3:ncol(occMatrix)]

genera<-nrow(occMatrix)
bins<-ncol(occMatrix)
rangeFreq<-rep(0,bins)

# Replace all NAs with 0s
for (t in 1:nrow(occMatrix)){
  for (b in 1:bins){
    if(is.na(occMatrix[t,b])){
      occMatrix[t,b]=0
    }
  }
}

bins<-ncol(occMatrix)
twoTimers<-rep(0,10000)
threeTimers<-rep(0,10000)
partTimers<-rep(0,10000)
gapFillers<-rep(0,10000)
firstSeen<-rep(0,10000)
lastSeen<-rep(0,10000)

for (t in 1:nrow(occMatrix)){
  for (b in 1:bins){
    if (occMatrix[t,b]>0){
      if (firstSeen[t]==0){
        firstSeen[t]=b
      }
      lastSeen[t]=b
    }
    if (b>1){
      if ((occMatrix[t,b-1]>0) & (b<bins)){
        if (occMatrix[t,b]>0){
          twoTimers[b]<-twoTimers[b]+1
        }
        if ((occMatrix[t,b]>0) & (occMatrix[t,b+1]>0)){
          threeTimers[b]<-threeTimers[b]+1
        }
        else if((occMatrix[t,b]==0) & (occMatrix[t,b+1]>0)){
          partTimers[b]<-partTimers[b]+1
        }
      }
      if(b>1 & (occMatrix[t,b-1])>0 & b+1<bins){  
        if((occMatrix[t,b+1]==0) & (occMatrix[t,b+2]>0)){
          gapFillers[b]<-gapFillers[b]+1
        }
      }
    }
  }
}

#plot(1:bins,log(gapFillers[1:13]),type="o",ylim=c(0,5),xlab="bins",ylab="log count")
#points(1:bins,log(partTimers[1:13]),col="green",type="o")
#points(1:bins,log(twoTimers[1:13]),col="red",type="o")
#points(1:bins,log(threeTimers[1:13]),col="blue",type="o")

#Calculate Foote rates (BC)
Nbt<-rep(0,14) # Crosses both boundaries
Nb<-rep(0,14) # Crosses bottom boundary
Nt<-rep(0,14) # Crosses top boundary 

for (t in 1:nrow(occMatrix)){
  if(firstSeen[t]+1<lastSeen[t]){
    for (b in firstSeen[t]+1:lastSeen[t]-1){
      Nbt[b]<-Nbt[b]+1
    }
  }
  if (firstSeen[t]<lastSeen[t]){
    Nb[lastSeen[t]]<-Nb[lastSeen[t]]+1
    Nt[firstSeen[t]]<-Nt[firstSeen[t]]+1
  }
}

Nb
Nt
Nbt<-Nbt[1:14]

(Fext<-log((Nb+Nbt)/Nbt)) #Extinction rate
(Forig<-log((Nt+Nbt)/Nbt)) #Origination rate


#Calculate gap-filler estimates
twoTimers<-twoTimers[1:13]
partTimers<-partTimers[1:13]
threeTimers<-threeTimers[1:13]
gapFillers<-gapFillers[1:13]
GFest<-log((twoTimers+partTimers)/(threeTimers+partTimers+gapFillers))
GFest
median(GFest,na.rm=TRUE)

#Calculate three timer (3T) extinction rate
output=NULL
for (i in 1:bins){
  ThreeExt<-log(twoTimers[i]/threeTimers[i])+log(threeTimers[i+1]/(threeTimers[i+1]+partTimers[i+1]))
  output<-cbind(ThreeExt,output)  
}
output
median(output,na.rm=TRUE)


# 3T origination rates
output=NULL
for (i in 1:bins){
  ThreeExt<-log(twoTimers[i+1]/threeTimers[i])+log(threeTimers[i+1]/(threeTimers[i+1]+partTimers[i+1]))
  output<-cbind(ThreeExt,output)  
}
output

par(mfrow=c(2,2))
plot(Fext,type="o",xlab="bins")
plot(Forig, type="o",xlab="bins")
plot(gapFillers,type="o",xlab="bins")
plot(1:length(output),output,type="o",ylab="3T")
