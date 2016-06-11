# define function to normalize a vector
normalize = function(x){x / sum(x)}


# load or calculate components of the transmission cycle
LEP = c(rep(0,5),1,rep(0,5))

HMTP.T = normalize(as.matrix(read.csv('../../malaria_serial_interval/data/infectivities_treated.csv')))
HMTP.U = normalize(scan('../../malaria_serial_interval/data/infectivities_untreated.csv',sep=','))

EIP.sd = 2.472235
EIP.mean = normalize(dnorm(1:ceiling(qnorm(.99,10.475,EIP.sd)),10.475,EIP.sd))
EIP.kankiya = normalize(dnorm(1:ceiling(qnorm(.99,10.3,EIP.sd)),10.3,EIP.sd))
EIP.kaduna = normalize(dnorm(1:ceiling(qnorm(.99,11.6,EIP.sd)),11.6,EIP.sd))
EIP.namawala = normalize(dnorm(1:ceiling(qnorm(.99,11.1,EIP.sd)),11.1,EIP.sd))
EIP.butelgut = normalize(dnorm(1:ceiling(qnorm(.99,8.9,EIP.sd)),8.9,EIP.sd))

MHTP.mean = normalize(dexp(1:ceiling(qexp(.99,.1175)),.1175))
MHTP.kankiya = normalize(dexp(1:ceiling(qexp(.99,.06)),.06))
MHTP.kaduna = normalize(dexp(1:ceiling(qexp(.99,.1)),.1))
MHTP.namawala = normalize(dexp(1:ceiling(qexp(.99,.17)),.17))
MHTP.butelgut = normalize(dexp(1:ceiling(qexp(.99,.14)),.14))


# calculate generation interval distributions
GI.T.mean = rep(0, which(LEP==1) + length(HMTP.T) + length(EIP.mean) + length(MHTP.mean))
for(ii in 1:length(HMTP.T)){
  for(jj in 1:length(EIP.mean)){
    for(kk in 1:length(MHTP.mean)){
      GI.T.mean[which(LEP==1) + ii + jj + kk] =
        GI.T.mean[which(LEP==1) + ii + jj + kk] +
        HMTP.T[ii] * EIP.mean[jj] * MHTP.mean[kk]
    }
  }
}

GI.U.mean = rep(0, which(LEP==1) + length(HMTP.U) + length(EIP.mean) + length(MHTP.mean))
for(ii in 1:length(HMTP.U)){
  for(jj in 1:length(EIP.mean)){
    for(kk in 1:length(MHTP.mean)){
      GI.U.mean[which(LEP==1) + ii + jj + kk] =
        GI.U.mean[which(LEP==1) + ii + jj + kk] +
        HMTP.U[ii] * EIP.mean[jj] * MHTP.mean[kk]
    }
  }
}

GI.T.kankiya = rep(0, which(LEP==1) + length(HMTP.T) + length(EIP.kankiya) + length(MHTP.kankiya))
for(ii in 1:length(HMTP.T)){
  for(jj in 1:length(EIP.kankiya)){
    for(kk in 1:length(MHTP.kankiya)){
      GI.T.kankiya[which(LEP==1) + ii + jj + kk] =
        GI.T.kankiya[which(LEP==1) + ii + jj + kk] +
        HMTP.T[ii] * EIP.kankiya[jj] * MHTP.kankiya[kk]
    }
  }
}

GI.U.kankiya = rep(0, which(LEP==1) + length(HMTP.U) + length(EIP.kankiya) + length(MHTP.kankiya))
for(ii in 1:length(HMTP.U)){
  for(jj in 1:length(EIP.kankiya)){
    for(kk in 1:length(MHTP.kankiya)){
      GI.U.kankiya[which(LEP==1) + ii + jj + kk] =
        GI.U.kankiya[which(LEP==1) + ii + jj + kk] +
        HMTP.U[ii] * EIP.kankiya[jj] * MHTP.kankiya[kk]
    }
  }
}

GI.T.kaduna = rep(0, which(LEP==1) + length(HMTP.T) + length(EIP.kaduna) + length(MHTP.kaduna))
for(ii in 1:length(HMTP.T)){
  for(jj in 1:length(EIP.kaduna)){
    for(kk in 1:length(MHTP.kaduna)){
      GI.T.kaduna[which(LEP==1) + ii + jj + kk] =
        GI.T.kaduna[which(LEP==1) + ii + jj + kk] +
        HMTP.T[ii] * EIP.kaduna[jj] * MHTP.kaduna[kk]
    }
  }
}

GI.U.kaduna = rep(0, which(LEP==1) + length(HMTP.U) + length(EIP.kaduna) + length(MHTP.kaduna))
for(ii in 1:length(HMTP.U)){
  for(jj in 1:length(EIP.kaduna)){
    for(kk in 1:length(MHTP.kaduna)){
      GI.U.kaduna[which(LEP==1) + ii + jj + kk] =
        GI.U.kaduna[which(LEP==1) + ii + jj + kk] +
        HMTP.U[ii] * EIP.kaduna[jj] * MHTP.kaduna[kk]
    }
  }
}

GI.T.namawala = rep(0, which(LEP==1) + length(HMTP.T) + length(EIP.namawala) + length(MHTP.namawala))
for(ii in 1:length(HMTP.T)){
  for(jj in 1:length(EIP.namawala)){
    for(kk in 1:length(MHTP.namawala)){
      GI.T.namawala[which(LEP==1) + ii + jj + kk] =
        GI.T.namawala[which(LEP==1) + ii + jj + kk] +
        HMTP.T[ii] * EIP.namawala[jj] * MHTP.namawala[kk]
    }
  }
}

GI.U.namawala = rep(0, which(LEP==1) + length(HMTP.U) + length(EIP.namawala) + length(MHTP.namawala))
for(ii in 1:length(HMTP.U)){
  for(jj in 1:length(EIP.namawala)){
    for(kk in 1:length(MHTP.namawala)){
      GI.U.namawala[which(LEP==1) + ii + jj + kk] =
        GI.U.namawala[which(LEP==1) + ii + jj + kk] +
        HMTP.U[ii] * EIP.namawala[jj] * MHTP.namawala[kk]
    }
  }
}

GI.T.butelgut = rep(0, which(LEP==1) + length(HMTP.T) + length(EIP.butelgut) + length(MHTP.butelgut))
for(ii in 1:length(HMTP.T)){
  for(jj in 1:length(EIP.butelgut)){
    for(kk in 1:length(MHTP.butelgut)){
      GI.T.butelgut[which(LEP==1) + ii + jj + kk] =
        GI.T.butelgut[which(LEP==1) + ii + jj + kk] +
        HMTP.T[ii] * EIP.butelgut[jj] * MHTP.butelgut[kk]
    }
  }
}

GI.U.butelgut = rep(0, which(LEP==1) + length(HMTP.U) + length(EIP.butelgut) + length(MHTP.butelgut))
for(ii in 1:length(HMTP.U)){
  for(jj in 1:length(EIP.butelgut)){
    for(kk in 1:length(MHTP.butelgut)){
      GI.U.butelgut[which(LEP==1) + ii + jj + kk] =
        GI.U.butelgut[which(LEP==1) + ii + jj + kk] +
        HMTP.U[ii] * EIP.butelgut[jj] * MHTP.butelgut[kk]
    }
  }
}


# save all variables to file
save(list=ls(),file='../../malaria_serial_interval/output/Fig1.RData')
