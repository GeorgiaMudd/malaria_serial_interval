load('../../malaria_serial_interval/output/Fig2.RData')
load('../../malaria_serial_interval/output/Fig3.RData')

tiff('../../malaria_serial_interval/output/Fig3.tiff',width=6.5,height=6.5,units='in',res=150,compression='lzw')

layout(matrix(1:2,2,1,byrow=T))

par(oma=c(1,2,0,0),mar=c(3.5,2.5,2,1.5))

yrange = c(0, max(IDP.A.100, IDP.A.50, IDP.A.20, IDP.A.10, IDP.A.0.1))
plot(1:365, IDP.A.100, axes = F, type = "l", col = "red", bty = "n", xlab = "", ylab = "", ylim = yrange)
lines(1:365, IDP.A.50, col = "green")
lines(1:365, IDP.A.20, col = "blue")
lines(1:365, IDP.A.10, col = "orange")
lines(1:365, IDP.A.0.1, col = "purple")
axis(1, pos = 0, at=seq(0,360,60), cex.axis = 0.7)
axis(2, pos = 0, cex.axis = 0.7)
mtext('A',3,at=0)
mtext('Pr(case detected on this date)',2,line=2.6,cex=.7)
mtext('Days since asymptomatic case infected',1,line=2.25,cex=.7,at=180)
legend("topright", title = "Threshold", legend = c("100", "50","20","10", "0.1"), col = 
         c("red", "green", "blue", "orange", "purple"), lty = 1, bty = "n", cex=0.7)

x = seq(-length(IDP.S)+1,length(GI.T.mean)+length(IDP.A),1)
y = SI.TU.100
yrange2 = c(0, max(SI.TU.100,SI.TU.50, SI.TU.20, SI.TU.10, SI.TU.0.1))
plot(x[44:409],y[44:409],xlim=c(0,365),type='l',axes = F,bty='n',ylim=yrange2, col = "red",
     xlab = "", ylab = "")
lines(x[44:409], SI.TU.50[44:409], col = "green")
lines(x[44:409], SI.TU.20[44:409], col = "blue")
lines(x[44:409], SI.TU.10[44:409], col = "orange")
lines(x[44:409], SI.TU.0.1[44:409], col = "purple")

axis(1, pos = 0, at = seq(-120, 360, 60), cex.axis = 0.7)
axis(2, pos = 0, cex.axis = 0.7)

mtext('B',3,at=0)
mtext('Pr(secondary case detected on this date)',2,line=2.6,cex=.7)
mtext('Days since primary case detected',1,line=2.25,cex=.7,at=180)

dev.off()

