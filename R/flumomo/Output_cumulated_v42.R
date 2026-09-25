### Graphs - cumulated ###

results4 <- read.table(
      paste0(outdir, "/", country, "_output_v4", ifelse(IArest, "_IArestricted", ""), ".txt"), 
      header = TRUE, sep =";", dec=".")
note <- ifelse(IArest, "Note: Negative effects of IA are not shown", "")
results4$agegrp.labels <- factor(results4$agegrp, 
    labels = c("0-4 years", "5-14 years", "15-64 years", "Aged 65", "Total"))


pal <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF")

plotCumulSingle <- function(a, mortality=FALSE, lwd=2, ylab=NULL) {
  if (mortality) {
    yvar <- "cm"
    ylb <- "Deaths / 100,000"
    title <- paste0("Cumulated ",s," mortality attributable to ",t,", age group: ",unique(results4$agegrp.labels)[a+1])
  } else {
    yvar <- "cd"
    ylb <- "Number of deaths"
    title <- paste0("Cumulated ",s," number of deaths attributable to ",t,", age group: ",unique(results4$agegrp.labels)[a+1])
  }
  if (!is.null(ylab)) ylb <- ylab
  xlb <- "Week number"
  if (is.na(ylb)) xlb <- NA
  
  with(out[(out$agegrp==a),], {
    plot(x=wk, y=get(yvar), type="n", xaxt="n", bty="l", 
      xlab=xlb, ylab=ylb,
      col=pal[1:length(levels(season))][as.integer(season)]
    )
    if (!is.na(ylb)) {
      mtext(title, side=3, line=2, adj=0)
      mtext(sprintf("FluMOMO v4 - week %s, %s", end_week, end_year), side=3, line=2, adj=1, cex=0.8)
    }
    axis(1, at=seq(1, max(wk), 4), labels = unique(floor(week))[seq(1, max(wk), 4)])
    grid(nx=NA, ny=NULL)
    abline(v=seq(1, max(wk), 4), col="lightgray", lty="dotted")
    mapply(function(s, col, l) {
      points(x=wk[season==s], y=get(yvar)[season==s], col=col, type="l", lty=l, lwd=lwd)
    }, s=levels(season), col=rev(pal[1:length(levels(season))]), 
    l=c(rep("dotted", length(levels(season))-1), "solid"))
  })
}


plotCumulMultiple <- function(mortality=FALSE, lwd=2, ylab=NULL) {
  par(mfrow=c(5,1), oma=c(6,3,2,0))
  for (a in 0:4) {
    par(mar=c(3,2,3,2))
    plotCumulSingle(a, mortality, ylab=NA)
    legend("topleft", levels(results4$agegrp.labels)[a+1], bty="n", xpd=NA, cex=1.3, inset=c(0,-0.2))
  }
  mtext(paste0("Cumulated ",s, ifelse(mortality, " mortality", " number of deaths"), " attributable to ",t,", age group: ",unique(results4$agegrp.labels)[a+1]), side=3, line=0, adj=0, outer=TRUE)
  mtext(sprintf("FluMOMO v4 - week %s, %s", end_week, end_year), side=3, line=0, adj=0.95, cex=0.8, outer=TRUE)
  mtext(ifelse(mortality, "Deaths / 100,000", "Number of deaths"), side=2, line=1, outer=TRUE)
  mtext("Week number", side=1, line=3, cex=0.8)
  legend("bottom", legend=rev(levels(out$season)),
    lty = c("solid", rep("dotted", length(levels(out$season))-1)),
    col = pal[1:length(levels(out$season))], lwd=2, cex=1.2,
    bty="n", horiz=TRUE, xpd=NA, inset=c(0,-0.59))
}


for (s in c("summer","winter","year")) {
  for (t in c("IA")) {
    out <- results4[!is.na(results4[,s]), 
        c("agegrp","agegrp.labels","year","week",paste0("cEd",t,"_",s),"N")]
    colnames(out) <- c("agegrp","agegrp.labels","year","week","cd","N")
    out$cm <- 100000 * out$cd / out$N
    out[out$week==53,"week"] <- 52.5
    if (s == "winter") {
      out$wk <- with(out, 1 + (week-40)*(week>=40) + (week+12)*(week<=20))
      out$season <- factor(paste0(as.character(with(out, year*(week>=40) + (year-1)*(week<=20))),"/",substr(as.character(with(out, year*(week>=40) + (year-1)*(week<=20))+1),3,4)))
    }
    if (s == "summer") {
      out$wk <- out$week-20
      out$season <- factor(out$year)
    }
    if (s == "year") {
      out$wk <- out$week
      out$season <- factor(out$year)
    }
    out <- out[with(out, order(agegrp, wk)),]
    png(paste0(outdir,"/cumulated_",t,"_deaths_agegroups_",s,"_v4", ifelse(IArest, "_IArestricted", ""), ".png"), width=1200, height=1600, res=130)
    plotCumulMultiple()
    dev.off()
    if (population) {
      png(paste0(outdir,"/cumulated_",t,"_mr_agegroups_",s,"_v4", ifelse(IArest, "_IArestricted", ""), ".png"), width=1200, height=1600, res=130)
      plotCumulMultiple(mortality=TRUE)
      dev.off()
    }
    for (a in 0:4) {
      png(paste0(outdir,"/cumulated_",t,"_deaths_agegroup_",a,"_",s,"_v4", ifelse(IArest, "_IArestricted", ""), ".png"), width=1200, height=800, res=130)
      par(oma=c(3,0,0,0))
      plotCumulSingle(a)
      legend("bottom", legend=rev(levels(out$season)),
        lty = c("solid", rep("dotted", length(levels(out$season))-1)),
        col = pal[1:length(levels(out$season))], lwd=2,
        bty="n", horiz=TRUE, xpd=NA, inset=c(0,-0.35))
      dev.off()
      if (population) {
        png(paste0(outdir,"/cumulated_",t,"_mr_agegroup_",a,"_",s,"_v4", ifelse(IArest, "_IArestricted", ""), ".png"), width=1200, height=800, res=130)
        par(oma=c(3,0,0,0))
        plotCumulSingle(a, mortality=TRUE)
        legend("bottom", legend=rev(levels(out$season)),
          lty = c("solid", rep("dotted", length(levels(out$season))-1)),
          col = pal[1:length(levels(out$season))], lwd=2,
          bty="n", horiz=TRUE, xpd=NA, inset=c(0,-0.35))
        dev.off()
      }
    }
  }
}
rm(plotCumulSingle, plotCumulMultiple, results4, a, s, t, out, note, pal)
