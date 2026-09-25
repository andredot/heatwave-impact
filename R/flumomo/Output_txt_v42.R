### ;-separated txt files ###

if (IArest) { 
  results4 <- read.table(file=paste0(outdir,"/",country,"_output_v4_IArestricted.txt"), 
      header = TRUE, sep =";", dec=".")
} else { 
  results4 <- read.table(file=paste0(outdir,"/",country,"_output_v4.txt"), 
      header = TRUE, sep =";", dec=".")
}

csvout <- function(x, s) {
  y <- do.call("rbind", lapply(sort(unique(x[,s])), function(i) {
    cbind(unique(x[,"agegrp"]), min(x[(x[,s]==i),"yw"]), max(x[(x[,s]==i),"yw"]), 
        round(tail(
          x[(x[,s]==i), 
              c(paste0("cEdIA_",s), paste0("cEdIA_",s,"_95L"), paste0("cEdIA_",s,"_95U"), 
              paste0("cEdET_",s), paste0("cEdET_",s,"_95L"), paste0("cEdET_",s,"_95U"), 
              paste0("cexcess_",s), paste0("cuexcess_",s), "N")],
          n=1)))
    }))
  colnames(y) <- c("AgeGroup", "Start", "End", "cIA", "cIA_95L", "cIA_95U", 
                    "cET", "cET_95L", "cET_95U", "cexcess", "cuexcess", "N")
  return(y)
}

results4$yw <- sprintf("%04iw%02i",results4$year, results4$week)

for (s in c("summer","winter","year")) {
  out <- do.call("rbind", lapply(0:4, function(a) {
    csvout(results4[(!is.na(results4[,s])) & (results4$agegrp==a),], s)
  }))
  out$AgeGroup <- factor(out$AgeGroup, levels = c(0,1,2,3,4), 
        labels = c("0-4 years","5-14 years","15-64 years","Aged 65","Total"))
  write.table(out[,c("AgeGroup", "Start", "End", "cIA", "cIA_95L", "cIA_95U", 
              "cET", "cET_95L", "cET_95U", "cexcess", "cuexcess")], 
      file = paste0(outdir, "/", s, "_deaths_v4_IA", ifelse(IArest, "restricted", ""), ".txt" ), 
      row.names = FALSE, quote = FALSE, sep =";", dec=".", na="")
  if (population) {
    out[,c("cIA", "cIA_95L", "cIA_95U", "cET", "cET_95L", "cET_95U", "cexcess", "cuexcess")] <-
      round(100000 * out[,c("cIA", "cIA_95L", "cIA_95U", "cET", "cET_95L", "cET_95U", "cexcess", "cuexcess")] / out[,"N"], 2)
    colnames(out) <- c("AgeGroup", "Start", "End", "cmrIA", "cmrIA_95L", "cmrIA_95U", "cmrET", "cmrET_95L", "cmrET_95U", "cexcess", "cuexcess", "N")
    write.table(out[,c("AgeGroup", "Start", "End", "cmrIA", "cmrIA_95L", "cmrIA_95U", 
              "cmrET", "cmrET_95L", "cmrET_95U", "cexcess", "cuexcess")],
      file=paste0(outdir,"/",s,"_mr_v4_IA", ifelse(IArest, "restricted", ""), ".txt"), 
      row.names = FALSE, quote = FALSE, sep =";", dec=".")
  }
}

rm(out, s, csvout, results4)
