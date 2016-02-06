# Shareholder Quorum Subsampling, version 3.3
# Author: John Alroy
# Modified by Jon Tennant
# SQS function from http://bio.mq.edu.au/~jalroy/SQS.html
# Note that the file is no longer available from source

#Read in data
data<-read.csv("filename",header=TRUE)
                     
#Get number of unique genus occurrences (i.e., raw taxonomic diversity)
data.diversity<-length(unique(Hettangian[,1])) # Change column number to whichever has the taxon name list

#Convert occurrences into an array with counts of occurrences by taxon name
data.occs<-table(Hettangian[,1]) #change column number to whichever has your taxon name list

# Read in John's sqs function

sqs<-function(ab,q,trials,method,ignore.singletons,dominant)  {
  
  params <- array(data=NA,dim=0,dimnames=c("raw richness"))
  if (missing(trials))  {
    trials <- 100
  }
  if (missing(method))	{
    method <- ""
  } else if (method != "" && method != "rarefaction" && method != "CR")	{
    return(print('If the method is rarefaction enter method="rarefaction" or "CR"',quote=F))
  }
  if ((q <= 0 || q >= 1) && method != "rarefaction" && method != "CR")	{
    
    return(print("If the method is SQS the quota must be greater than zero and less than one",quote=F))
  } else if (q < 1 && (method == "rarefaction" || method == "CR"))	{
    
    return(print("If the method is rarefaction the quota must be an integer",quote=F))
  }
  ignore <- 0
  if (! missing(ignore.singletons) && (ignore.singletons == T || ignore.singletons == "yes" || ignore.singletons == "y"))	{
    ignore <- 1
  }
  if (missing(dominant))	{
    dominant <- 0
  } else if (dominant != "" && dominant != "exclude" && dominant != "no")	{
    return(print('To exclude the dominant taxon, enter dominant="exclude" or "no"',quote=F))
  }
  
  # compute basic statistics
  specimens <- sum(ab)
  singletons <- 0
  doubletons <- 0
  highest <- 0
  for (i in 1:length(ab))	{
    if (ab[i] == 1)	{
      singletons <- singletons + 1
    } else if (ab[i] == 2)	{
      doubletons <- doubletons + 1
    }
    if (ab[i] > highest)	{
      highest <- ab[i]
      mostfrequent <- i
    }
  }
  u <- 1 - singletons / specimens
  # optionally, exclude the dominant taxon
  if (dominant == "exclude" || dominant == "no")	{
    u <- 1 - singletons / (specimens - highest)
  }
  
  if (u == 0)	{
    return(print("Coverage is zero because all taxa are singletons",quote=F))
  }
  
  # compute raw taxon frequencies (temporarily)
  freq <- ab / specimens
  
  # standard recursive equation for Fisher's alpha
  alpha <- 10
  oldalpha <- 0
  while (abs(alpha - oldalpha) > 0.0000001)	{
    oldalpha <- alpha
    alpha <- length(ab) / log(1 + specimens/alpha)
  }
  
  params["raw richness"] <- length(ab)
  params["Good's u"] <- u
  params["subsampled richness"] <- NA
  params["subsampled u"] <- NA
  params["Chao 1"] <- length(ab) + singletons**2/(2* doubletons)
  params["subsampled Chao 1"] <- NA
  # governing parameter of the geometric series distribution
  params["k"] <- abs(lm(log(sort(freq)) ~ c(1:length(freq)))$coefficients[2])
  params["Fisher's alpha"] <- alpha
  params["Shannon's H"] <- -1 * sum(freq * log(freq))
  params["Hurlbert's PIE"] <- (1 - sum(freq**2)) * length(ab) / (length(ab) - 1)
  params["dominance"] <- highest / specimens
  params["specimens"] <- specimens
  params["singletons"] <- singletons
  params["doubletons"] <- doubletons
  params["specimens drawn"] <- 0
  
  if (dominant != "exclude" && dominant != "no")	{
    highest <- 0
    mostfrequent <- 0
  }
  
  # return if the rarefaction quota is equal to or higher than the
  #  specimen count
  if (method == "rarefaction" && q >= specimens - highest)	{
    return(params)
  }
  
  # compute taxon frequencies (tweak added in version 3.1)
  freq <- ab - (singletons + doubletons / 2) / length(ab)
  freq <- freq / (specimens - highest)
  
  # return if the quorum target is higher than estimated coverage
  if ((q > sum(freq) && method != "rarefaction" && method != "CR") || (q >= sum(ab)))	{
    return(params)
  }
  
  # create an array in which each cell corresponds to one specimen
  ids <- array()
  n <- 0
  for (i in 1:length(ab))	{
    for (j in 1:ab[i])	{
      n <- n + 1
      ids[n] <- i
    }
  }
  
  # subsampling trial loop
  # s will be the subsampled taxon count
  s <- array(rep(0,trials))
  drawn <- array(rep(0,trials))
  mostfrequentdrawn <- array(rep(0,trials))
  subsingle <- array(rep(0,trials))
  subdouble <- array(rep(0,trials))
  subchao <- array(rep(0,trials))
  
  for (trial in 1:trials)	{
    pool <- ids
    left <- length(pool)
    seen <- array(data=rep(0,length(ab)))
    subfreq <- array(rep(0,length(ab)))
    if (method != "rarefaction" && method != "CR")	{
      udrawn <- 0
      # equation new to version 3.0
      # the exponent corrects for downwards bias
      while (udrawn < q)	{
        # draw a specimen
        x <- floor(runif(1,min=1,max=left+1))
        # add to frequency and taxon sums if species has
        #  not been drawn previously
        subfreq[pool[x]] <- subfreq[pool[x]] + 1
        if (seen[pool[x]] == 0)	{
          if (pool[x] != mostfrequent && (ignore == 0 || ab[pool[x]] > 1))	{
            udrawn <- udrawn + freq[pool[x]]
          }
          seen[pool[x]] <- 1
          # randomly throw back some draws that put the sum over q
          #  (an even better algorithm added in version 3.0)
          if (runif(1) <= freq[pool[x]] || udrawn < q)	{
            s[trial] <- s[trial] + 1
          } else	{
            subfreq[pool[x]] <- subfreq[pool[x]] - 1
          }
        }
        # decrease pool of specimens not yet drawn
        pool[x] <- pool[left]
        left <- left - 1
      }
    } else	{
      i <- 0
      draws <- 0
      while (i < q)	{
        draws <- draws + 1
        x <- floor(runif(1,min=1,max=length(ids)-draws+2))
        subfreq[pool[x]] <- subfreq[pool[x]] + 1
        if (pool[x] != mostfrequent)	{
          i <- i + 1
        }
        if (seen[pool[x]] == 0)	{
          seen[pool[x]] <- 1
          s[trial] <- s[trial] + 1
        }
        pool[x] <- pool[length(ids)-draws+1]
      }
    }
    for (i in 1:length(ab))	{
      if (subfreq[i] == 1 && i != mostfrequent)	{
        subsingle[trial] <- subsingle[trial] + 1
      } else if (subfreq[i] == 2 && i != mostfrequent)	{
        subdouble[trial] <- subdouble[trial] + 1
      }
    }
    if (subsingle[trial] > 0 && subdouble[trial] > 0)	{
      subchao[trial] <- s[trial] + subsingle[trial]**2/(2*subdouble[trial])
    } else	{
      subchao[trial] <- s[trial]
    }
    drawn[trial] <- sum(subfreq)
    if (mostfrequent != 0)	{
      mostfrequentdrawn[trial] <- subfreq[mostfrequent]
    }
  }
  
  # compute vectors of non-zero counts
  options(warn=-1)
  s2 <- sort(sqrt(s-1))^2+1
  d2 <- sort(sqrt(drawn-1))^2+1
  m2 <- sort(sqrt(mostfrequentdrawn-1))^2+1
  ss2 <- sort(sqrt(subsingle-1))^2+1
  options(warn=0)
  
  # compute geometric means
  params["subsampled richness"] <- exp(mean(log(s2))) * length(s2)/length(s)
  params["specimens drawn"] <- exp(mean(log(d2))) * length(d2)/length(drawn)
  meanmost <- 0
  if (sum(mostfrequentdrawn) > 0)	{
    meanmost <- exp(mean(log(m2))) * length(m2)/length(mostfrequentdrawn)
  }
  meansubsingle <- exp(mean(log(ss2))) * length(ss2)/length(subsingle)
  
  params["subsampled u"] <- 1 - meansubsingle / (params["specimens drawn"] - meanmost)
  params["subsampled Chao 1"] <- exp(mean(log(subchao)))
  return(params)
  
}

# Use on occurrence data
# (Occurrences,quorum,trials)
# Calculate Good's u for each time interval with a quorum of 0.1
sqs(data.occs,0.1,1)["Good's u"]

#Perform analysis using quorum selected from above (Good's u minus a little bit), and a number of trials
sqs(data.occs,0.4,1000)["subsampled richness"] 

