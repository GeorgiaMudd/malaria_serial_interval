tiff('../../malaria_serial_interval/output/Fig5.tiff', width = 6.5, height = 6.5, units = 'in', res = 150, compression = 'lzw')

par(mfrow = c(2,2), mar = c(2.5,4,2,1), oma = c(1,0,0,0))

plot(1:365, normalize(IND.UT[1:365,1]), type = "l", col = "blue", xlim = c(0,365), ylim = c(0,0.07), xlab = "", ylab = "Density", las = 1, bty = 'l', lwd = 0.5, yaxs = "i", xaxs = "i")
#polygon(IND.UT[1:365,1], col = rgb(0,0,1,0.5), border = NA)
for(i in 2:15){
  lines(normalize(IND.UT[1:365,i]), col = "blue", lwd = 0.5)
  #polygon(IND.UT[1:365,i], col = rgb(0,0,1,0.5), border = NA)
}
lines(mean.UT, col = "black", lwd = 1.5)
mtext('A',3,at=0,line=.5)
legend("bottomright", c("Mean", "Individual"), col = c("black", "blue"), lwd = 2, bty = "n", cex = 0.8)
mtext('Untreated',side=3,cex=0.8,line=0)

plot(1:365, normalize(IND.T[1:365,1]), type = "l", col = "red", xlim = c(0,365), ylim = c(0,0.07), xlab = "", ylab = "", las = 1, bty = 'l', lwd = 0.5, yaxs = "i", xaxs = "i")
#polygon(IND.T[1:365,1], col = rgb(1,0,0,0.5), border = NA)
for(i in 2:15){
  lines(normalize(IND.T[1:365,i]), col = "red", lwd = 0.5)
  #polygon(IND.T[1:365,i], col = rgb(1,0,0,0.5), border = NA)
}
lines(mean.T, col = "black", lwd = 1.5)
legend("bottomright", c("Mean", "Individual"), col = c("black", "red"), lwd = 2, bty = "n", cex = 0.8)
mtext('B',3,at=0,line=.5)
mtext('Treated',side=3,cex=.8,line=0)

x = c(1:365)
plot(1:365, quantiles.UT[1,], type = "l", col = rgb(0,0,1,0.5), xlim = c(0,365), ylim = c(0,1.0), xlab  = "", ylab = "Cumulative density", las = 1, bty = "l", yaxs = "i", xaxs = "i")
lines(quantiles.UT[2,], col = rgb(1,0,0,0.5))
lines(quantiles.UT[4,], col = rgb(1,0,0,0.5))
lines(quantiles.UT[5,], col = rgb(0,0,1,0.5))
polygon(c(x, rev(x)), c(quantiles.UT[5,],rev(quantiles.UT[1,])), col = rgb(0,0,1,0.5), border = NA)
polygon(c(x, rev(x)), c(quantiles.UT[4,],rev(quantiles.UT[2,])), col = rgb(1,0,0,0.5), border = NA)
lines(quantiles.UT[3,])
legend("bottomright", c("0.025 - 0.975 quantile", "0.25-0.75 quantile", "Median"), col = c(rgb(0,0,1,0.5),rgb(1,0,0,0.5), "black"), lty = 1, bty = "n", cex = 0.8, lwd=2)
mtext('C',3,at=0,line=.5)

plot(1:365, quantiles.T[1,], type = "l", col = rgb(0,0,1,0.5), xlim = c(0,365), ylim = c(0,1.0), xlab = "", ylab = "", las = 1, bty = "l", yaxs = "i", xaxs = "i")
lines(quantiles.T[2,], col = rgb(1,0,0,0.5))
lines(quantiles.T[4,], col = rgb(1,0,0,0.5))
lines(quantiles.T[5,], col = rgb(1,0,0,0.5))
polygon(c(x, rev(x)), c(quantiles.T[5,], rev(quantiles.T[1,])), col = rgb(0,0,1,0.5), border = NA)
polygon(c(x, rev(x)), c(quantiles.T[4,], rev(quantiles.T[2,])), col = rgb(1,0,0,0.5), border = NA)
lines(quantiles.T[3,])
legend("bottomright", c("0.025 - 0.975 quantile", "0.25-0.75 quantile", "Median"), col = c(rgb(0,0,1,0.5),rgb(1,0,0,0.5), "black"), lty = 1, bty = "n", cex = 0.8, lwd=2)
mtext('D',3,at=0,line=.5)
mtext('Days since primary infection',side=1,line=2.25,cex=.8, at = -70)

dev.off()