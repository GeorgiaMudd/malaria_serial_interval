load('../../malaria_serial_interval/output/Fig1.RData')

gametocytemia.U = as.matrix(read.csv('../../malaria_serial_interval/data/total_asexual_parasitemia.csv',header=F))
IDP.A.100 = hist(unlist(sapply(1:nrow(gametocytemia.U),function(ii)
  which(gametocytemia.U[ii,]>100))),breaks=seq(.5,801.5,1),plot=F)$density[1:365]
IDP.A.50 = hist(unlist(sapply(1:nrow(gametocytemia.U),function(ii)
  which(gametocytemia.U[ii,]>50))),breaks=seq(.5,801.5,1),plot=F)$density[1:365]
IDP.A.20 = hist(unlist(sapply(1:nrow(gametocytemia.U),function(ii)
  which(gametocytemia.U[ii,]>20))),breaks=seq(.5,801.5,1),plot=F)$density[1:365]
IDP.A.10 = hist(unlist(sapply(1:nrow(gametocytemia.U),function(ii)
  which(gametocytemia.U[ii,]>10))),breaks=seq(.5,801.5,1),plot=F)$density[1:365]
IDP.A.1 = hist(unlist(sapply(1:nrow(gametocytemia.U),function(ii)
  which(gametocytemia.U[ii,]>1))),breaks=seq(.5,801.5,1),plot=F)$density[1:365]
IDP.A.0.1 = hist(unlist(sapply(1:nrow(gametocytemia.U),function(ii)
  which(gametocytemia.U[ii,]>0.1))),breaks=seq(.5,801.5,1),plot=F)$density[1:365]


SI.TU.100 = rep(0, length(IDP.S) + length(GI.T.mean) + length(IDP.A.100))
for(ii in 1:length(IDP.S)){
  for(jj in 1:length(GI.T.mean)){
    for(kk in 1:length(IDP.A.100)){
      SI.TU.100[length(IDP.S) - ii + jj + kk] =
        SI.TU.100[length(IDP.S) - ii + jj + kk] +
        IDP.S[ii] * GI.T.mean[jj] * IDP.A.100[kk]
    }
  }
}

SI.TU.50 = rep(0, length(IDP.S) + length(GI.T.mean) + length(IDP.A.50))
for(ii in 1:length(IDP.S)){
  for(jj in 1:length(GI.T.mean)){
    for(kk in 1:length(IDP.A.50)){
      SI.TU.50[length(IDP.S) - ii + jj + kk] =
        SI.TU.50[length(IDP.S) - ii + jj + kk] +
        IDP.S[ii] * GI.T.mean[jj] * IDP.A.50[kk]
    }
  }
}

SI.TU.20 = rep(0, length(IDP.S) + length(GI.T.mean) + length(IDP.A.20))
for(ii in 1:length(IDP.S)){
  for(jj in 1:length(GI.T.mean)){
    for(kk in 1:length(IDP.A.20)){
      SI.TU.20[length(IDP.S) - ii + jj + kk] =
        SI.TU.20[length(IDP.S) - ii + jj + kk] +
        IDP.S[ii] * GI.T.mean[jj] * IDP.A.20[kk]
    }
  }
}

SI.TU.10 = rep(0, length(IDP.S) + length(GI.T.mean) + length(IDP.A.10))
for(ii in 1:length(IDP.S)){
  for(jj in 1:length(GI.T.mean)){
    for(kk in 1:length(IDP.A.10)){
      SI.TU.10[length(IDP.S) - ii + jj + kk] =
        SI.TU.10[length(IDP.S) - ii + jj + kk] +
        IDP.S[ii] * GI.T.mean[jj] * IDP.A.10[kk]
    }
  }
}

SI.TU.1 = rep(0, length(IDP.S) + length(GI.T.mean) + length(IDP.A.1))
for(ii in 1:length(IDP.S)){
  for(jj in 1:length(GI.T.mean)){
    for(kk in 1:length(IDP.A.1)){
      SI.TU.1[length(IDP.S) - ii + jj + kk] =
        SI.TU.1[length(IDP.S) - ii + jj + kk] +
        IDP.S[ii] * GI.T.mean[jj] * IDP.A.1[kk]
    }
  }
}

SI.TU.0.1 = rep(0, length(IDP.S) + length(GI.T.mean) + length(IDP.A.0.1))
for(ii in 1:length(IDP.S)){
  for(jj in 1:length(GI.T.mean)){
    for(kk in 1:length(IDP.A.0.1)){
      SI.TU.0.1[length(IDP.S) - ii + jj + kk] =
        SI.TU.0.1[length(IDP.S) - ii + jj + kk] +
        IDP.S[ii] * GI.T.mean[jj] * IDP.A.0.1[kk]
    }
  }
}

# save all variables to file
save(list=ls(),file='../../malaria_serial_interval/output/Fig3.RData')