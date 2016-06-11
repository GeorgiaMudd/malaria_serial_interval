tiff(filename = '../../malaria_serial_interval/output/Fig4.tiff', width = 6.5, height = 9, units = 'in', res = 400)
#pdf(file = '../../malaria_serial_interval/output/Fig3.pdf', width = 5, height = 6.5)

layout(matrix(1:8,4,2))
par(oma=c(0,10,3,0),mar=c(1,0,0,0))

amp = 1; sigma = 14
s = list(x=1:730, y=1:365, z=new.t.1.14)
psp = persp(s,theta=0,phi=30,col=rgb(1,1,1,0),border=NA,ylab='',xlab='Time',zlab='Density',axes=F)
text(-.36,.27,expression('Day (1' * degree * ')'),cex=1,srt=45)
text(-.4,-.14,'Density',cex=1,srt=104)
mtext(expression('Day (2' * degree * ')'),1,cex=.7,line=-.1)
moz = amp * dnorm(mod(1:730,365),180,sigma) / dnorm(180,180,sigma) + 1
tr = trans3d(x=c(1:730,730,1),y=rep(365,732),z=c(moz,0,0)/max(moz)*max(s$z),psp)
polygon(tr,col=rgb(0,0,0,.175),border=NA)
for(ii in rev(seq(1,360,by=30))){
  tr = trans3d(x=ii:730,y=rep(ii,730-ii+1),z=s$z[ii:730,ii],psp)
  lines(tr,col=rainbow(360,alpha=.75)[ii],lwd=2)
}
mtext('Treated',3,line=1)
mtext('Low amplitude',2,las=1,at=0,adj=0,line=9)
mtext('Very seasonal',2,las=1,at=-.1,adj=0,line=9)

amp = 9; sigma = 14
s = list(x=1:730, y=1:365, z=new.t.9.14)
psp = persp(s,theta=0,phi=30,col=rgb(1,1,1,0),border=NA,ylab='',xlab='Time',zlab='Density',axes=F)
text(-.36,.27,expression('Day (1' * degree * ')'),cex=1,srt=45)
text(-.4,-.14,'Density',cex=1,srt=104)
mtext(expression('Day (2' * degree * ')'),1,cex=.7,line=-.1)
moz = amp * dnorm(mod(1:730,365),180,sigma) / dnorm(180,180,sigma) + 1
tr = trans3d(x=c(1:730,730,1),y=rep(365,732),z=c(moz,0,0)/max(moz)*max(s$z),psp)
polygon(tr,col=rgb(0,0,0,.175),border=NA)
for(ii in rev(seq(1,360,by=30))){
  tr = trans3d(x=ii:730,y=rep(ii,730-ii+1),z=s$z[ii:730,ii],psp)
  lines(tr,col=rainbow(360,alpha=.75)[ii],lwd=2)
}
mtext('High amplitude',2,las=1,at=0,adj=0,line=9)
mtext('Very seasonal',2,las=1,at=-.1,adj=0,line=9)

amp = 1; sigma = 120
s = list(x=1:730, y=1:365, z=new.t.1.120)
psp = persp(s,theta=0,phi=30,col=rgb(1,1,1,0),border=NA,ylab='',xlab='Time',zlab='Density',axes=F)
text(-.36,.27,expression('Day (1' * degree * ')'),cex=1,srt=45)
text(-.4,-.14,'Density',cex=1,srt=104)
mtext(expression('Day (2' * degree * ')'),1,cex=.7,line=-.1)
moz = amp * dnorm(mod(1:730,365),180,sigma) / dnorm(180,180,sigma) + 1
tr = trans3d(x=c(1:730,730,1),y=rep(365,732),z=c(moz,0,0)/max(moz)*max(s$z),psp)
polygon(tr,col=rgb(0,0,0,.175),border=NA)
for(ii in rev(seq(1,360,by=30))){
  tr = trans3d(x=ii:730,y=rep(ii,730-ii+1),z=s$z[ii:730,ii],psp)
  lines(tr,col=rainbow(360,alpha=.75)[ii],lwd=2)
}
mtext('Low amplitude',2,las=1,at=0,adj=0,line=9)
mtext('Less seasonal',2,las=1,at=-.1,adj=0,line=9)

amp = 9; sigma = 120
s = list(x=1:730, y=1:365, z=new.t.9.120)
psp = persp(s,theta=0,phi=30,col=rgb(1,1,1,0),border=NA,ylab='',xlab='Time',zlab='Density',axes=F)
text(-.36,.27,expression('Day (1' * degree * ')'),cex=1,srt=45)
text(-.4,-.14,'Density',cex=1,srt=104)
mtext(expression('Day (2' * degree * ')'),1,cex=.7,line=-.1)
moz = amp * dnorm(mod(1:730,365),180,sigma) / dnorm(180,180,sigma) + 1
tr = trans3d(x=c(1:730,730,1),y=rep(365,732),z=c(moz,0,0)/max(moz)*max(s$z),psp)
polygon(tr,col=rgb(0,0,0,.175),border=NA)
for(ii in rev(seq(1,360,by=30))){
  tr = trans3d(x=ii:730,y=rep(ii,730-ii+1),z=s$z[ii:730,ii],psp)
  lines(tr,col=rainbow(360,alpha=.75)[ii],lwd=2)
}
mtext('High amplitude',2,las=1,at=0,adj=0,line=9)
mtext('Less seasonal',2,las=1,at=-.1,adj=0,line=9)

amp = 1; sigma = 14
s = list(x=1:730, y=1:365, z=new.u.1.14)
psp = persp(s,theta=0,phi=30,col=rgb(1,1,1,0),border=NA,ylab='',xlab='Time',zlab='Density',axes=F)
text(-.36,.27,expression('Day (1' * degree * ')'),cex=1,srt=45)
text(-.4,-.14,'Density',cex=1,srt=104)
mtext(expression('Day (2' * degree * ')'),1,cex=.7,line=-.1)
moz = amp * dnorm(mod(1:730,365),180,sigma) / dnorm(180,180,sigma) + 1
tr = trans3d(x=c(1:730,730,1),y=rep(365,732),z=c(moz,0,0)/max(moz)*max(s$z),psp)
polygon(tr,col=rgb(0,0,0,.175),border=NA)
for(ii in rev(seq(1,360,by=30))){
  tr = trans3d(x=ii:730,y=rep(ii,730-ii+1),z=s$z[ii:730,ii],psp)
  lines(tr,col=rainbow(360,alpha=.75)[ii],lwd=2)
}
mtext('Untreated',3,line=1)

amp = 9; sigma = 14
s = list(x=1:730, y=1:365, z=new.u.9.14)
psp = persp(s,theta=0,phi=30,col=rgb(1,1,1,0),border=NA,ylab='',xlab='Time',zlab='Density',axes=F)
text(-.36,.27,expression('Day (1' * degree * ')'),cex=1,srt=45)
text(-.4,-.14,'Density',cex=1,srt=104)
mtext(expression('Day (2' * degree * ')'),1,cex=.7,line=-.1)
moz = amp * dnorm(mod(1:730,365),180,sigma) / dnorm(180,180,sigma) + 1
tr = trans3d(x=c(1:730,730,1),y=rep(365,732),z=c(moz,0,0)/max(moz)*max(s$z),psp)
polygon(tr,col=rgb(0,0,0,.175),border=NA)
for(ii in rev(seq(1,360,by=30))){
  tr = trans3d(x=ii:730,y=rep(ii,730-ii+1),z=s$z[ii:730,ii],psp)
  lines(tr,col=rainbow(360,alpha=.75)[ii],lwd=2)
}

amp = 1; sigma = 120
s = list(x=1:730, y=1:365, z=new.u.1.120)
psp = persp(s,theta=0,phi=30,col=rgb(1,1,1,0),border=NA,ylab='',xlab='Time',zlab='Density',axes=F)
text(-.36,.27,expression('Day (1' * degree * ')'),cex=1,srt=45)
text(-.4,-.14,'Density',cex=1,srt=104)
mtext(expression('Day (2' * degree * ')'),1,cex=.7,line=-.1)
moz = amp * dnorm(mod(1:730,365),180,sigma) / dnorm(180,180,sigma) + 1
tr = trans3d(x=c(1:730,730,1),y=rep(365,732),z=c(moz,0,0)/max(moz)*max(s$z),psp)
polygon(tr,col=rgb(0,0,0,.175),border=NA)
for(ii in rev(seq(1,360,by=30))){
  tr = trans3d(x=ii:730,y=rep(ii,730-ii+1),z=s$z[ii:730,ii],psp)
  lines(tr,col=rainbow(360,alpha=.75)[ii],lwd=2)
}

amp = 9; sigma = 120
s = list(x=1:730, y=1:365, z=new.u.9.120)
psp = persp(s,theta=0,phi=30,col=rgb(1,1,1,0),border=NA,ylab='',xlab='Time',zlab='Density',axes=F)
text(-.36,.27,expression('Day (1' * degree * ')'),cex=1,srt=45)
text(-.4,-.14,'Density',cex=1,srt=104)
mtext(expression('Day (2' * degree * ')'),1,cex=.7,line=-.1)
moz = amp * dnorm(mod(1:730,365),180,sigma) / dnorm(180,180,sigma) + 1
tr = trans3d(x=c(1:730,730,1),y=rep(365,732),z=c(moz,0,0)/max(moz)*max(s$z),psp)
polygon(tr,col=rgb(0,0,0,.175),border=NA)
for(ii in rev(seq(1,360,by=30))){
  tr = trans3d(x=ii:730,y=rep(ii,730-ii+1),z=s$z[ii:730,ii],psp)
  lines(tr,col=rainbow(360,alpha=.75)[ii],lwd=2)
}
dev.off()
