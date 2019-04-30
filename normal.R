install.packages("ggplot2")
install.packages("plyr")

library(ggplot2)
library(parallel)
library(plyr)
options(mc.cores=12)


#Generate data
norm_data <- rnorm(50,0,2)
xbar <- mean(norm_data)
sigma <- var(norm_data)^(1/2)


#Rejection algorithm
accepted <- c()
batches <- 10^1
bsize=10^7
start <- Sys.time()
for (i in 1:batches){
  mu <- rnorm(bsize,0,1)
  sam <- mclapply(mu,sampling)
  sam <- plyr::compact(sam[])
  #sam <- sam[!is.na(sam)]
  sam <- as.numeric(sam)  
  accepted <- c(accepted,sam)
  print(i)
}
end <- Sys.time()
rej <- accepted

sampling <- function(mu){
  samp <- rnorm(50,mu,2)
  sbar <- mean(samp)
  if ((abs(sbar-xbar))<0.00005){
    return(mu)
  }
}

hist(sam,breaks=150,freq=FALSE)

###MCMC
#Normal transition kernel, sd=0.2
runs <- 10^7
mcmc_accept <- rep(NA,runs)
mcmc_accept[1] <- 0
mu <- 0
samples <- c()

start_mcmc <- Sys.time()
for(i in 1:runs){
  sample <- rnorm(1,mu,1)
  sbar <- sampling(sample)
  if (!is.null(sbar)){
    alpha <- dnorm(sbar,mean=0,sd=1)/dnorm(mu,mean=0,sd=1)
    alpha <- min(1,alpha)
    compare <- runif(1)
    if (alpha > compare){
      mu <- sample
    }
  }
  mcmc_accept[i] <- mu
  print(i)
}
end_mcmc <- Sys.time()

#Set burn in period
mcmc_burn <- mcmc_accept[5000:length(mcmc_accept)]

#Plot 
plot(mcmc_accept,type="l")

#Create histogram of unique samples
mcmc_unique <- unique(mcmc_burn)
hist(mcmc_unique,breaks=100,freq=FALSE)

#Create histogram of all samples
plot(samples)
hist(samples,freq=FALSE)




##SMC
#Uniform pertubation kernel
n <- 5*10^4
epsilons <- c(1,0.5,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001,0.00005)
t <- length(epsilons)
sd <- c(2.5,2,1.7,1.5,1,1.2,1,0.7,0.5,0.1)
params <- array(999,dim=c(t+1,n,2))

initial <- rnorm(n,0,1)
params[1,,1] <- initial
initial_weights <- rep(1/n,n)
params[1,,2] <- initial_weights


smc_start <- Sys.time()
for (i in 1:t){
  samp <- sample(params[i,,1],size=n,replace=TRUE,prob=params[i,,2])   #Generate n samples
  sample_array <- mclapply(samp,smc_sampling,i,mc.cores=12)           #Pass samples through sampling function to make use of multicore capabilities
  sample_array <- simplify2array(sample_array)                        #Convert to array
  total_weight <- sum(sample_array[2,])                               #Normalise weights (this step may not be needed, I'm not entirely sure how R handles sampling)
  sample_array[2,] <- sample_array[2,]/total_weight
  params[i+1,,1]<- sample_array[1,]                                   #Store values in array
  params[i+1,,2]<- sample_array[2,]
  print(i)                                                            #Progress indicator
}
smc_end <- Sys.time()


#Sampling function
smc_sampling <- function(mu_samp,i){
  sbar <- abs(xbar) + epsilons[i] + 400
  while ((abs(sbar-xbar))>epsilons[i]){
    mu_samp <- sample(params[i,,1],size=1,replace=TRUE,prob=params[i,,2])
    new_mu <- runif(1,mu_samp-sd[i],mu_samp+sd[i])
    samp <- rnorm(50,new_mu,2)
    sbar <- mean(samp)
  }
  weight <- dnorm(new_mu,mean=0,sd=1)
  return( c(new_mu,weight))
}




##Actual posterior
mu_0 <- 0
sigma_0 <- 1
n <- 50
sigma_1 <- ((n/sigma^2)+(1/sigma_0^2))^(-1/2)
mu_1 <- (sigma_1)^2*((n/sigma^2)*xbar)

#Plot actual posterior as overlay
lines(x_axis,dnorm(x_axis,mean=mu_1,sd=sigma_1),lwd=2) 

#par(mfrow = c(2, 2))



#Rejection histogram
hist(accepted,breaks=100,freq=FALSE,ylim=c(0,3),xlab="Accepted values of μ",main="Histogram of accepted values of μ")


#SMC histograms
hist(params[2,,1],breaks=100,freq=FALSE,ylim=c(0,2),xlab="Accepted values of μ",main="Histogram of accepted values of μ for epsilon = 1")
hist(params[3,,1],breaks=100,freq=FALSE,ylim=c(0,2),xlab="Accepted values of μ",main="Histogram of accepted values of μ for epsilon = 0.5")
hist(params[4,,1],breaks=100,freq=FALSE,ylim=c(0,2),xlab="Accepted values of μ",main="Histogram of accepted values of μ for epsilon = 0.1")
hist(params[10,,1],breaks=100,freq=FALSE,ylim=c(0,2),xlab="Accepted values of μ",main="Histogram of accepted values of μ for epsilon = 0.00005")


