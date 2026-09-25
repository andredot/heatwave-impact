### Ambient temperatures ###
ET <- read.table(file=paste0(indir,"/temperature.txt"), header = T, sep =";", dec=".")
ET <- ET[with(ET, order(year, week)),]
ET$t = ifelse((ET$week %% 2) == 0, min(ET$temp, ET$ptmin, na.rm = TRUE), max(ET$temp, ET$ptmax, na.rm = TRUE))
ET$temp = ifelse(is.na(ET$temp), ET$ptemp, ET$temp)
ET$Uexecc <- with(ET, pmax(ptmax,temp))
ET$Lexecc <- with(ET, pmin(ptmin,temp))
ET$wk <- 1:nrow(ET)
ET$yw <- sprintf("%04iw%02d", ET$year, ET$week)

png(paste0(outdir,"/temperature_v4R.png"), width=1200, height=800, res=130)
par(mar=c(9,4,3,2))
with(ET, {
  plot(x=wk, y=t, type="n", xlab=NA, ylab="Degree Celsius", 
      main="Ambient temperatures", bty="l", xaxt="n")
  axis(1, at=wk[as.integer(substr(yw,6,7)) %in% c(1,27)], 
      labels=yw[as.integer(substr(yw,6,7)) %in% c(1,27)], las=2, cex.axis=0.9)
  polygon(x=c(wk,rev(wk)), y=c(Uexecc, rev(ptmax)),
    border=do.call(rgb, c(as.list(col2rgb("red")/255), alpha=0.5)),
    col=do.call(rgb, c(as.list(col2rgb("red")/255), alpha=0.5)))
  polygon(x=c(wk,rev(wk)), y=c(Lexecc, rev(ptmin)),
    border=do.call(rgb, c(as.list(col2rgb("blue")/255), alpha=0.5)),
    col=do.call(rgb, c(as.list(col2rgb("blue")/255), alpha=0.5)))
  points(x=wk, y=ptemp, type="l", col="black")
  points(x=wk, y=ptmin, type="l", col="blue")
  points(x=wk, y=ptmax, type="l", col="red")
  points(x=wk, y=temp, pch=19, cex=0.7)
  mtext("Year-week", side=1, line=4.5)
  legend("bottom", c("expected temperature", "max expected temperature", "min expected temperature"),
    lty="solid", col=c("black", "red", "blue"), lwd=2,
    horiz=TRUE, inset=c(0,-0.39), xpd=NA, bty="n", cex=1)
  legend("bottom", c("observed temperature", "extreme heat", "extreme cold"),
    col=c("black", do.call(rgb, c(as.list(col2rgb("red")/255), alpha=0.5)), do.call(rgb, c(as.list(col2rgb("blue")/255), alpha=0.5))), pch=c(19,15,15),
    horiz=TRUE, inset=c(0,-0.46), xpd=NA, bty="n", cex=1, pt.cex=c(1,2,2))  
})
dev.off()
rm(ET)


### Influenza Activity ###

results4 <- read.table(
    file = paste0(outdir, "/", country, "_output_v4", ifelse(IArest, "_IArestricted", ""), ".txt"), 
    header = TRUE, sep =";", dec=".") 
results4 <- results4[with(results4, order(agegrp, year, week)),]
results4$wk <- unlist(lapply(table(results4$agegrp), seq))
results4$agegrp.labels <- ordered(results4$agegrp, 
    labels = c("0-4 years", "5-14 years", "15-64 years", "Aged 65", "Total"))
results4$yw <- sprintf("%04iw%02d", results4$year, results4$week)


png(paste0(outdir,"/IA_agegroups_v4R.png"), width=1200, height=1600, res=130)
par(mar=c(3,4,2,2), oma=c(5,2,2,0), mfrow=c(5,1))
for (a in 0:4)
  with(results4[(results4$agegrp==a),], {
    plot(x=wk, y=IA, xaxt="n", ylab=NA, xlab=NA, type="n", bty="l", 
        main=unique(results4$agegrp.labels)[a+1])
    if (a==4) {
      axis(1, at=wk[as.integer(substr(yw,6,7)) %in% c(1,27)], 
          labels=yw[as.integer(substr(yw,6,7)) %in% c(1,27)], las=2, cex.axis=0.9)
    }
    grid(nx=NA, ny=NULL)
    abline(v=wk[as.integer(substr(yw,6,7)) %in% c(1,27)], col="lightgray", lty="dotted")
    points(x=wk, y=IA, col="red", lwd=2, type="l")
  })
mtext("Influenza Activity per week", side=2, outer=TRUE)
mtext("Influenza Activity by age group", side=3, outer=TRUE, adj=0, cex=1.2)
mtext("year-week", side=1, line=5)
dev.off()


for (a in 0:4) {
  png(paste0(outdir,"/IA_agegroup_",a,"_v4R.png"), width=1200, height=800, res=130)
  par(mar=c(6,4,4,2))
    with(results4[(results4$agegrp==a),], {
      plot(x=wk, y=IA, xaxt="n", ylab="Influenza Activity per week", xlab=NA, type="n", bty="l")
      axis(1, at=wk[as.integer(substr(yw,6,7)) %in% c(1,27)], 
          labels=yw[as.integer(substr(yw,6,7)) %in% c(1,27)], las=2, cex.axis=0.9)
      grid(nx=NA, ny=NULL)
      abline(v=wk[as.integer(substr(yw,6,7)) %in% c(1,27)], col="lightgray", lty="dotted")
      mtext(paste("Influenza Activity, age group:",unique(results4$agegrp.labels)[a+1]), 
          side=3, line=1, adj=0)
      points(x=wk, y=IA, col="red", lwd=2, type="l")
      mtext("year-week", side=1, line=4.5)
    })
  dev.off()
}
rm(a)