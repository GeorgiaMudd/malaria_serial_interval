# load data
IND.UT = read.csv('../../malaria_serial_interval/data/total_untreated.csv')
IND.T = read.csv('../../malaria_serial_interval/data/total_within_treated.csv')

# generate matrices
cummulativesum.T = matrix(nrow = 1000, ncol = 365)
cummulativesum.UT = matrix(nrow = 1000, ncol = 365)
quantiles.T = matrix(data = rep(0, 5*365), nrow = 5, ncol = 365)
quantiles.UT = matrix(data = rep(0, 5*365), nrow = 5, ncol = 365)
quantiles.range = c(0.025,0.25,0.50,0.75,0.975)
mean.T = matrix(nrow = 365, ncol = 1)
mean.UT = matrix(nrow = 365, ncol = 1)

# function to normalize vector 
normalize = function(x){x / sum(x)}

#generate mean curves
mean.T = normalize(rowMeans(IND.T)[1:365])
mean.UT = normalize(rowMeans(IND.UT)[1:365])

# remove NaN's 
mask = apply(IND.T, 2, is.nan)
IND.T[mask] = 0

# fill matrices 
for(i in 1:1000){
  cummulativesum.T[i,] = cumsum(normalize(IND.T[1:365,i]))
  cummulativesum.UT[i,] = cumsum(normalize(IND.UT[1:365,i]))
}

for(i in 1:5){
  for(j in 1:365){
    quantiles.T[i,j] = quantile((cummulativesum.T[1:1000,j]), quantiles.range[i], names = FALSE, na.rm=T)
    quantiles.UT[i,j] = quantile((cummulativesum.UT[1:1000,j]), quantiles.range[i], names = FALSE, na.rm=T)
  }
}