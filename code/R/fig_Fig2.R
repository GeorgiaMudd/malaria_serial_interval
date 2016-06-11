load('../../malaria_serial_interval/output/Fig2.RData')

tiff('../../malaria_serial_interval/output/Fig2.tiff',width=6.5,height=6.5,units='in',res=150,compression='lzw')

layout(matrix(1:6,3,2,byrow=T))

par(oma=c(1,2,0,0),mar=c(3.5,2.5,2,1.5))

barplot(IDP.S,xaxs='i',yaxs='i',space=0,border=NA,col=rgb(0,0,0,.4))
legend('topright',' Treated \n Sympto.',bty='n')
axis(1,at=seq(0,50,10))
mtext('A',3,at=0)
mtext('Pr(case detected on this date)',2,line=2.6,cex=.7)

barplot(IDP.A,xaxs='i',yaxs='i',space=0,border=NA,col=rgb(0,0,0,.4))
legend('topright','Untreated \n Asympto.',bty='n')
axis(1,at=seq(0,360,60))
mtext('B',3,at=0)
mtext('Days since case infected',1,line=2.25,cex=.7,at=-40)

par(mar=c(2.5,2.5,2,1.5))
x = seq(-length(IDP.S)+1,length(GI.T.mean)+length(IDP.S),1)
y = SI.TT
plot(x,y,xlim=c(-120,365),type='l',xaxs='i',yaxs='i',bty='n',ylim=c(0,1.05*max(y)))
polygon(c(x,x[1]),c(y,y[1]),col=rgb(0,0,1,.4))
lines(x,y,col=4,lwd=3)
legend(x=180,y=.039,expression('1' * degree * ' = Treated'),fill=rgb(0,0,1,.4),border=rgb(0,0,1,.8),bty='n')
legend(x=158,y=.034,expression('2' * degree * ' = Sympto.'),lwd=2,col=4,bty='n')
abline(v=0,lty=2)
mtext('Pr(secondary case detected on this date)',2,line=2.6,at=-.007,cex=.7)
mtext('C',3,at=-120)

x = seq(-length(IDP.S)+1,length(GI.T.mean)+length(IDP.A),1)
y = SI.TU
plot(x,y,xlim=c(-120,365),type='l',xaxs='i',yaxs='i',bty='n',ylim=c(0,1.05*max(y)))
polygon(c(x,x[1]),c(y,y[1]),col=rgb(0,0,1,.4))
lines(x,y,col=2,lwd=3)
legend(x=180,y=.009,expression('1' * degree * ' = Treated'),fill=rgb(0,0,1,.4),border=rgb(0,0,1,.8),bty='n')
legend(x=158,y=.0079,expression('2' * degree * ' = Asympto.'),lwd=2,col=2,bty='n')
abline(v=0,lty=2)
mtext('D',3,at=-120)

x = seq(-length(IDP.A)+1,length(GI.T.mean)+length(IDP.S),1)
y = SI.UT
plot(x,y,xlim=c(-120,365),type='l',xaxs='i',yaxs='i',bty='n',ylim=c(0,1.05*max(y)))
polygon(c(x,x[1]),c(y,y[1]),col=rgb(1,0,0,.4))
lines(x,y,col=4,lwd=3)
legend(x=180,y=.007,expression('1' * degree * ' = Untreated'),fill=rgb(1,0,0,.4),border=rgb(1,0,0,.8),bty='n')
legend(x=158,y=.0061,expression('2' * degree * ' = Sympto.'),lwd=2,col=4,bty='n')
abline(v=0,lty=2)
mtext('E',3,at=-120)

x = seq(-length(IDP.A)+1,length(GI.T.mean)+length(IDP.A),1)
y = SI.UU
plot(x,y,xlim=c(-120,365),type='l',xaxs='i',yaxs='i',bty='n',ylim=c(0,1.05*max(y)))
polygon(c(x,x[1]),c(y,y[1]),col=rgb(1,0,0,.4))
lines(x,y,col=2,lwd=3)
legend(x=180,y=.005,expression('1' * degree * ' = Untreated'),fill=rgb(1,0,0,.4),border=rgb(1,0,0,.8),bty='n')
legend(x=158,y=.0044,expression('2' * degree * ' = Asympto.'),lwd=2,col=2,bty='n')
abline(v=0,lty=2)
mtext('F',3,at=-120)
mtext('Days since primary case detected',1,line=2.25,cex=.7,at=-175)

dev.off()
