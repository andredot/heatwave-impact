### Graphs - calendar time ###

if (IArest) {
  results4 <- read.table(file=paste0(outdir,"/",country,"_output_v4_IArestricted.txt"), header = TRUE, sep =";", dec=".")
  results4$EIA <- results4$EB + pmax(0,results4$EdIA)
  results4$EIAET <- results4$EB + pmax(0,results4$EdIA) + results4$EdET
  note <- "Note: Negative effects of IA are not shown"
} else {
  results4 <- read.table(file=paste0(outdir,"/",country,"_output_v4.txt"), header = TRUE, sep =";", dec=".")
  results4$EIA <- results4$EB + results4$EdIA
  results4$EIAET <- results4$EB + results4$EdIA + results4$EdET
  note <- ""
}
results4 <- results4[with(results4, order(agegrp, year, week)),]
results4$wk <- unlist(lapply(table(results4$agegrp), seq))
results4$agegrp.labels <- ordered(results4$agegrp, 
    labels = c("0-4 years","5-14 years","15-64 years","Aged 65","Total"))
results4$yw <- sprintf("%04dw%02d", results4$year, results4$week)


plotDeaths <- function(a, mortality=FALSE, annot=TRUE) {
  if (mortality) {
    results4[,c("deaths","EB","EB_95L","EB_95U","EIAET","EIA")] <- 
      100000 * results4[,c("deaths","EB","EB_95L","EB_95U","EIAET","EIA")] / results4[,"N"]
  }
  with(results4[(results4$agegrp==a),], {
    plot(wk, EB, type="n", xaxt="n", bty="l", xlab=NA, ylab=NA, 
        ylim=range(pretty(range(c(EIAET, EIA, EB_95U, EB_95L)))))
    if (annot) {
      axis(1, at=wk[as.integer(substr(yw,6,7)) %in% c(1,27)], 
          labels=yw[as.integer(substr(yw,6,7)) %in% c(1,27)], las=2, cex.axis=0.9)
    }
    grid(nx=NA, ny=NULL)
    abline(v=wk[as.integer(substr(yw,6,7)) %in% c(1,27)], col="lightgray", lty="dotted")
    polygon(x=c(wk,rev(wk)), y=c(deaths,rev(EB)), 
      col=do.call(rgb, c(as.list(col2rgb("gray")/255), alpha=0.7)),
      border=do.call(rgb, c(as.list(col2rgb("gray")/255), alpha=0.7)))
    points(wk, EIAET, col="darkgreen", lwd=2, type="l")
    points(wk, EIA, col="red", lwd=2, type="l")
    points(wk, EB, lwd=2, type="l")
    points(wk, EB_95L, lty="dashed", lwd=2, type="l")
    points(wk, EB_95U, lty="dashed", lwd=2, type="l")
  })
  if (annot) {
    mtext(sprintf("FluMOMO v4 - week %s, %s", end_week, end_year), side=3, line=1, adj=1, cex=0.8)
    mtext("year-week", side=1, line=4.5)
    if (mortality) {
      mtext(paste("Mortality by age group:",unique(results4$agegrp.labels)[a+1]), side=3, line=1, adj=0, cex=1.2)
      mtext("Deaths / 100.000 per week", side=2, line=2.5)
    } else {
      mtext(paste("Number of deaths, age group:",unique(results4$agegrp.labels)[a+1]), side=3, line=1, adj=0, cex=1.2)
      mtext("Number of deaths per week", side=2, line=2.5)
    }
    legend("bottom", c("Effect of Influenza Activity (IA)",
        "Effect of excess temperature (on top of IA)", "Baseline", "95% confidence interval"), 
        lwd=2, col=c("red","darkgreen","black","black"), lty=c(rep("solid",3),"dashed"), 
        bty="n", ncol=2, xpd=NA, inset=c(0,-0.5), seg.len=4)
  }
}


plotDeathsMultiple <- function(mortality=FALSE) {
  par(mar=c(3,2,3,2), mfrow=c(5,1), oma=c(8,3.5,3,0))
  for(a in 0:4) {
    plotDeaths(a, mortality=mortality, annot=FALSE)
    mtext(unique(results4$agegrp.labels)[a+1], side=3, line=1)
    if (a==4) {
      with(results4[(results4$agegrp==a),], 
        axis(1, at=wk[as.integer(substr(yw,6,7)) %in% c(1,27)], 
          labels=yw[as.integer(substr(yw,6,7)) %in% c(1,27)], las=2))
    }
  }
  mtext(sprintf("FluMOMO v4 - week %s, %s", end_week, end_year), side=3, line=0, adj=0.9, cex=0.8, outer=TRUE)
  mtext("year-week", side=1, line=5)
  if (mortality) {
    mtext("Mortality by age group", side=3, line=0, adj=0, cex=1.2, outer=TRUE)
    mtext("Deaths / 100.000 per week", side=2.5, line=1, outer=TRUE)
  } else {
    mtext("Number of deaths by age group", side=3, line=0, adj=0, cex=1.2, outer=TRUE)
    mtext("Number of deaths per week", side=2.5, line=1, outer=TRUE)
  }
  legend("bottom", c("Effect of Influenza Activity (IA)",
      "Effect of excess temperature (on top of IA)", "Baseline", "95% confidence interval"), 
      lwd=2, col=c("red","darkgreen","black","black"), lty=c(rep("solid",3),"dashed"), 
      bty="n", ncol=2, xpd=NA, inset=c(0,-0.9), seg.len=4)
}


png(paste0(outdir,"/deaths_agegroups_v4", ifelse(IArest, "_IArestricted", ""), ".png"), width=1200, height=1600, res=130)
plotDeathsMultiple(mortality=FALSE)
dev.off()
for (a in 0:4) {
  png(paste0(outdir, "/deaths_agegroup_", a, "_v4", ifelse(IArest, "_IArestricted", ""), ".png"), width=1200, height=800, res=130)
  par(mar=c(10,4,3,2))
  plotDeaths(a, mortality=FALSE)
  dev.off()
}

if (population) {
  png(paste0(outdir,"/mr_agegroups_v4", ifelse(IArest, "_IArestricted", ""), ".png"), width=1200, height=1600, res=130)
  plotDeathsMultiple(mortality=TRUE)
  dev.off()  
  for (a in 0:4) {
    png(paste0(outdir,"/mr_agegroup", a, "_v4", ifelse(IArest, "_IArestricted", ""), ".png"), width=1200, height=800, res=130)
    par(mar=c(10,4,3,2))
    plotDeaths(a, mortality=TRUE)
    dev.off()
  }
}
rm(plotDeaths, plotDeathsMultiple, results4, a, note)
