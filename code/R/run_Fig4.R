# load data
t.1.14 = as.matrix(read.csv('../../malaria_serial_interval/data/treated_seasonal_1_14.csv'))
t.9.14 = as.matrix(read.csv('../../malaria_serial_interval/data/treated_seasonal_9_14.csv'))
t.1.120 = as.matrix(read.csv('../../malaria_serial_interval/data/treated_seasonal_1_120.csv'))
t.9.120 = as.matrix(read.csv('../../malaria_serial_interval/data/treated_seasonal_9_120.csv'))
u.1.14 = as.matrix(read.csv('../../malaria_serial_interval/data/untreated_seasonal_1_14.csv'))
u.9.14 = as.matrix(read.csv('../../malaria_serial_interval/data/untreated_seasonal_9_14.csv'))
u.1.120 = as.matrix(read.csv('../../malaria_serial_interval/data/untreated_seasonal_1_120.csv'))
u.9.120 = as.matrix(read.csv('../../malaria_serial_interval/data/untreated_seasonal_9_120.csv'))

# make a 365 x 365 matrix to store in
new.t.1.14 = matrix(0,nrow=730,ncol=365)
new.t.9.14 = matrix(0,nrow=730,ncol=365)
new.t.1.120 = matrix(0,nrow=730,ncol=365)
new.t.9.120 = matrix(0,nrow=730,ncol=365)
new.u.1.14 = matrix(0,nrow=730,ncol=365)
new.u.9.14 = matrix(0,nrow=730,ncol=365)
new.u.1.120 = matrix(0,nrow=730,ncol=365)
new.u.9.120 = matrix(0,nrow=730,ncol=365)

# populate the first day
new.t.1.14[1:365,1] = t.1.14[1:365,1]
new.t.9.14[1:365,1] = t.9.14[1:365,1]
new.t.1.120[1:365,1] = t.1.120[1:365,1]
new.t.9.120[1:365,1] = t.9.120[1:365,1]
new.u.1.14[1:365,1] = u.1.14[1:365,1]
new.u.9.14[1:365,1] = u.9.14[1:365,1]
new.u.1.120[1:365,1] = u.1.120[1:365,1]
new.u.9.120[1:365,1] = u.9.120[1:365,1]

# populate all other days to actually start on their day
for(ii in 2:365){
  new.t.1.14[ii:(ii+364),ii] = t.1.14[1:365,ii]
  new.t.9.14[ii:(ii+364),ii] = t.9.14[1:365,ii]
  new.t.1.120[ii:(ii+364),ii] = t.1.120[1:365,ii]
  new.t.9.120[ii:(ii+364),ii] = t.9.120[1:365,ii]
  new.u.1.14[ii:(ii+364),ii] = u.1.14[1:365,ii]
  new.u.9.14[ii:(ii+364),ii] = u.9.14[1:365,ii]
  new.u.1.120[ii:(ii+364),ii] = u.1.120[1:365,ii]
  new.u.9.120[ii:(ii+364),ii] = u.9.120[1:365,ii]
}

# normalize
for(ii in 1:365){
  new.t.1.14[,ii] = new.t.1.14[,ii] / sum(new.t.1.14[,ii])
  new.t.9.14[,ii] = new.t.9.14[,ii] / sum(new.t.9.14[,ii])
  new.t.1.120[,ii] = new.t.1.120[,ii] / sum(new.t.1.120[,ii])
  new.t.9.120[,ii] = new.t.9.120[,ii] / sum(new.t.9.120[,ii])
  new.u.1.14[,ii] = new.u.1.14[,ii] / sum(new.u.1.14[,ii])
  new.u.9.14[,ii] = new.u.9.14[,ii] / sum(new.u.9.14[,ii])
  new.u.1.120[,ii] = new.u.1.120[,ii] / sum(new.u.1.120[,ii])
  new.u.9.120[,ii] = new.u.9.120[,ii] / sum(new.u.9.120[,ii])
}

mod = function(x, y){
  x%%y
}