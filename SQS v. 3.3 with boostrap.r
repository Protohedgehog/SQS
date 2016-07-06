# sqs version 3.2 by John Alroy
# performs shareholder quorum subsampling on an array of specimen counts
# can be used to perform classical rarefaction instead of SQS
# written 29 July 2010; version 2.0 completed 14 February 2011; versions 3.0,
#  3.1, 3.2, and 3.3 written 3 June, 12 July, 19-23 August, and 2 December 2011
# change in version 3.3: fixed bug in computation of maximum estimated coverage
# changes in version 3.2: specimens drawn and subsampled U are based on
#  geometric means of counts; ignore.singletons option added
# change in version 3.1: an even better frequency adjustment involving
#  singletons and doubletons
# changes in version 3.0: an even better subsampling algorithm involving a new
#  adjustment to u combined with a new throwback criterion
# changes in version 2.0: improved subsampling algorithm; including the dominant 
#  taxon is now the default; improved reporting of errors and basic statistics
# warning: do not use this program with taxonomic occurrence data drawn from
#  multiple published references because it is not designed to count
#  single-reference taxa or adjust for long taxonomic lists
# warning: version 1.0 yields estimates that are downwards-biased when q < 0.6
#  and abundance distributions are highly uneven

# Read in data
data<-read.csv("file.csv",header=TRUE)

#Convert occurrence list to occurrence table

data.occs<-table(data[,"genus"])


sqs<-function(ab,q,trials,method,ignore.singletons,dominant)  {
  
  params <- array(data=NA,dim=0,dimnames=c("raw richness"))
  if (missing(trials))  {
    trials <- 100
  }
  if (missing(method))  {
    method <- ""
  } else if (method != "" && method != "rarefaction" && method != "CR")  {
    return(print('If the method is rarefaction enter method="rarefaction" or "CR"',quote=F))
  }
  if ((q <= 0 || q >= 1) && method != "rarefaction" && method != "CR")  {
    
    return(print("If the method is SQS the quota must be greater than zero and less than one",quote=F))
  } else if (q < 1 && (method == "rarefaction" || method == "CR"))  {
    
    return(print("If the method is rarefaction the quota must be an integer",quote=F))
  }
  ignore <- 0
  if (! missing(ignore.singletons) && (ignore.singletons == T || ignore.singletons == "yes" || ignore.singletons == "y"))  {
    ignore <- 1
  }
  if (missing(dominant))  {
    dominant <- 0
  } else if (dominant != "" && dominant != "exclude" && dominant != "no")  {
    return(print('To exclude the dominant taxon, enter dominant="exclude" or "no"',quote=F))
  }
  
  # compute basic statistics
  specimens <- sum(ab)
  singletons <- 0
  doubletons <- 0
  highest <- 0
  for (i in 1:length(ab))  {
    if (ab[i] == 1)  {
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
  params["subsampled richness"] <- u
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

# Check real value with a quorum of 0.4
sqs(data.occs,0.4)['subsampled richness']

# Bootstrap
s<-rep(NA,1000)
for(i in 1:length(s)){
  bootstrapped<-data.occs[sample(1:length(data.occs),replace=TRUE)]
  s[i]<-sqs(bootstrapped,0.4)['subsampled richness']
}

# To remove any NAs (bootstrapping can produce a lot of pseudo-singletons)
s<-na.omit(s)
s<-as.numeric(s)

# Returns the median and 95% confidence intervals
quantile(s,probs=c(0.5,0.025,0.975))
hist(s)

