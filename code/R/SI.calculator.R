# (c) 2016 John Huber, Alex Perkins
# jhuber3@nd.edu, taperkins@nd.edu
# CRAPL v0.1 license - http://matt.might.net/articles/crapl/
# See README.md in git repository for more info.
# This code calculates the GT/SI under specific contexts

# set working directory 
setwd('../../malaria_serial_interval/code')

# load appropriate packages 

library(pracma)

# load appropriate data

gametocytemia.U = as.matrix(read.csv('../../malaria_serial_interval/data/Gametocytemia.csv',header=F))
IDP.A = hist(unlist(sapply(1:nrow(gametocytemia.U),function(ii) #Infection to Detection Period (Asymptomatic)
  which(gametocytemia.U[ii,]>20))),breaks=seq(.5,801.5,1),plot=F)$density[1:365]
IDP.S = as.numeric(unlist(read.csv('../../malaria_serial_interval/data/bite_to_clinic.csv'))) #Infection to Detection Symptomatic
GT.UT.data = as.matrix(read.csv('../../malaria_serial_interval/data/secondary_probabilities_general.csv')) #Untreated GT
GT.T.data = as.matrix(read.csv('../../malaria_serial_interval/data/new_drug_secondary_probabilities_general.csv')) #Treated GT
gamma.shape = c(27.37, 5.38)
gamma.scale = c(1.73, 15.65)

# normalize data 

normalize = function(x){
  x = x / sum(x)
}
GT.UT.data = normalize(GT.UT.data)
GT.T.data = normalize(GT.T.data)

# approximate functions for the treated and untreated GT

GT.UT.fun = approxfun(x = 1:863, y = GT.UT.data)
GT.T.fun = approxfun(x = 1:863, y = GT.T.data)

# calculate the characteristic equation for a specific case history.
# argument case.history: a vector of 1's and 2's where 1 signifies "treated" or symptomatic" and 2 signifies "untreated" or "asymptomatic"
# argument SI.status: a boolean that indicates whether the serial interval will be calculated. 

characteristic.equation = function(case.history, SI.status){
  number.cases = length(case.history)
  number.generations = number.cases - 1

  shape = vector(mode = "double", length = number.generations)
  scale = vector(mode = "double", length = number.generations)
  
  for(ii in 1:number.generations){
    case.index = case.history[ii]
    shape[ii] = gamma.shape[as.numeric(case.index)]
    scale[ii] = gamma.scale[as.numeric(case.index)]
  }
  
  char_eq = eval(parse(text = paste('(fft(dgamma(1:863,shape = ', shape, ', scale = ', scale, ')))', sep = " ", collapse  = '*')))
  GT.mean = Re(ifft(char_eq))  
    
  if(SI.status == TRUE){
    if(case.history[1] == 1 && case.history[number.cases] == 1){
      serial.interval = rep(0, length(IDP.S) + length(GT.mean) + length(IDP.S))
      timing = seq(-length(IDP.S)+1, length(GT.mean)+length(IDP.S),1)
      for(ii in 1:length(IDP.S)){
        for(jj in 1:length(GT.mean)){
          for(kk in 1:length(IDP.S)){
            serial.interval[length(IDP.S) - ii + jj + kk] =
              serial.interval[length(IDP.S) - ii + jj + kk] + 
              IDP.S[ii] * GT.mean[jj] * IDP.S[kk]
          }
        }
      }
    }
    if(case.history[1] == 1 && case.history[number.cases] == 2){
      serialinterval = rep(0, length(IDP.S) + length(GT.mean) + length(IDP.A))
      timing = seq(-length(IDP.S)+1, length(GT.mean) + length(IDP.A), 1)
      for(ii in 1:length(IDP.S)){
        for(jj in 1:length(GT.mean)){
          for(kk in 1:length(IDP.A)){
            serial.interval[length(IDP.S) - ii + jj + kk] = 
              serial.interval[length(IDP.S) - ii + jj + kk] +
              IDP.S[ii] * GT.mean[jj] * IDP.A[kk]
          }
        }
      }
    }
    if(case.history[1] == 2 && case.history[number.cases] == 1){
      serial.interval = rep(0, length(IDP.A) + length(GT.mean) + length(IDP.S))
      timing = seq(-length(IDP.S)+1, length(GT.mean) + length(IDP.A), 1)
      for(ii in 1:length(IDP.A)){
        for(jj in 1:length(GT.mean)){
          for(kk in 1:length(IDP.S)){
            serial.interval[length(IDP.A) - ii + jj + kk] =
              serial.interval[length(IDP.A) - ii + jj + kk] + 
              IDP.A[ii] * GT.mean[jj] * IDP.S[kk]
          }
        }
      }
    }
    if(case.history[1] == 2 && case.history[number.cases] == 2){
      serial.interval = rep(0, length(IDP.A) + length(GT.mean) + length(IDP.A))
      timing = seq(-length(IDP.A)+1, length(GT.mean)+length(IDP.A), 1)
      for(ii in 1:length(IDP.A)){
        for(jj in 1:length(GT.mean)){
          for(kk in 1:length(IDP.A)){
            serial.interval[length(IDP.A) - ii + jj + kk] =
              serial.interval[length(IDP.A) - ii + jj + kk] +
              IDP.A[ii] * GT.mean[jj] * IDP.A[kk]
          }
        }
      }
    }
    SI.dataframe = data.frame(timing, serial.interval)
    colnames(SI.dataframe) = c("t", "p")
    return(SI.dataframe)
    #return(serial.interval)
  }
  else
    timing = seq(1,length(GT.mean), 1)
    GT.dataframe = data.frame(timing, GT.mean)
    colnames(GT.dataframe) = c("t", "p")
    return(GT.dataframe)
    #return(GT.mean)
}

# function to compute output the GI or SI based on user input. 
# argument casehistory: vector of 1's and 2's where 1 signifies "treated" or "symptomatic" and 2 signifies "untreated' or "asymptomatic"
# argument output.type: "Gamma" (Gamma distribution parameters), "NB" (Negative Binomial Parameters), "PDF" (Probability Density Function)
# argument SI.boolean: boolean which determines whether a serial interval will be calculated. passed to characteristic.equation function

SI.calculator  = function(casehistory, output.type, SI.boolean){
  # call characteristic.equation function. Returns output based on user input
  characteristic.eq = characteristic.equation(case.history = casehistory, SI.status = SI.boolean)
  
  # normalize vector
  generated.output = characteristic.eq
  
  # evaluates by method of least squares. computes gamma parameters shape and scale 
  
  if(output.type == "Gamma"){
    leastsquares.gamma = function(par){
      x = min(characteristic.eq$t):max(characteristic.eq$t)
      values = dgamma(x, shape = par[1], scale = par[2])
      sum((values - characteristic.eq$p)^2)
    }
    generated.output = optim(par = c(5,20), leastsquares.gamma)$par
  }
  
  #evaluates by method of least squares. computes negative binomial parameters n and probability
  if(output.type == "NB"){
    leastsquares.nbinom = function(par){
      x = min(characteristic.eq$Timing):max(characteristic.eq$Timing)
      par = exp(par) + 10^-10
      values = dnbinom(x, size = par[1], prob = par[2])
      sum((values - characteristic.eq$p)^2)
    }
    optimized = optim(par=c(log(10), log(0.1)), leastsquares.nbinom)
    generated.output = exp(optimized$par)
  }
  
  # returns generated.output. If output.type == "PDF", a vector of densities will be returned
  # otherwise, generated.output is a set of fitted parameters
  generated.output
}


