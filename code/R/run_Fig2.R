load('../../malaria_serial_interval/output/Fig1.RData')

gametocytemia.U = as.matrix(read.csv('../../malaria_serial_interval/data/total_asexual_parasitemia.csv',header=F))
IDP.A = hist(unlist(sapply(1:nrow(gametocytemia.U),function(ii)
which(gametocytemia.U[ii,]>50))),breaks=seq(.5,801.5,1),plot=F)$density[1:365]
IDP.S = as.numeric(unlist(read.csv('../../malaria_serial_interval/data/bite_to_clinic.csv')))


SI.TT = rep(0, length(IDP.S) + length(GI.T.mean) + length(IDP.S))
for(ii in 1:length(IDP.S)){
  for(jj in 1:length(GI.T.mean)){
    for(kk in 1:length(IDP.S)){
      SI.TT[length(IDP.S) - ii + jj + kk] =
      SI.TT[length(IDP.S) - ii + jj + kk] +
      IDP.S[ii] * GI.T.mean[jj] * IDP.S[kk]
    }
  }
}

SI.UT = rep(0, length(IDP.A) + length(GI.T.mean) + length(IDP.S))
for(ii in 1:length(IDP.A)){
  for(jj in 1:length(GI.U.mean)){
    for(kk in 1:length(IDP.S)){
      SI.UT[length(IDP.A) - ii + jj + kk] =
        SI.UT[length(IDP.A) - ii + jj + kk] +
        IDP.A[ii] * GI.U.mean[jj] * IDP.S[kk]
    }
  }
}

SI.TU = rep(0, length(IDP.S) + length(GI.T.mean) + length(IDP.A))
for(ii in 1:length(IDP.S)){
  for(jj in 1:length(GI.T.mean)){
    for(kk in 1:length(IDP.A)){
      SI.TU[length(IDP.S) - ii + jj + kk] =
        SI.TU[length(IDP.S) - ii + jj + kk] +
        IDP.S[ii] * GI.T.mean[jj] * IDP.A[kk]
    }
  }
}

SI.UU = rep(0, length(IDP.A) + length(GI.T.mean) + length(IDP.A))
for(ii in 1:length(IDP.A)){
  for(jj in 1:length(GI.U.mean)){
    for(kk in 1:length(IDP.A)){
      SI.UU[length(IDP.A) - ii + jj + kk] =
        SI.UU[length(IDP.A) - ii + jj + kk] +
        IDP.A[ii] * GI.U.mean[jj] * IDP.A[kk]
    }
  }
}

# save all variables to file
save(list=ls(),file='../../malaria_serial_interval/output/Fig2.RData')
