load("./data/LHD_CD2.rda")
load("./data/LHD_MD2.rda")
load("./data/UD_CD2.rda")
load("./data/UD_MD2.rda")
DesignQuery <- function(n,s,q,crit="CD2", ShowCrit = TRUE)
{
  if (q==n) {
    if (crit=="CD2") DataX = LHD_CD2
    if (crit=="MD2") DataX = LHD_MD2
  } else {
    if (crit=="CD2") DataX = UD_CD2
    if (crit=="MD2") DataX = UD_MD2
  }
  idx = which(DataX$n == n & DataX$s == s & DataX$q == q)
  if (length(idx) == 0) {
    warning("No design found.")
    return(NULL)
  } else{
    D = DataX$Design[idx]
    D = matrix(as.numeric(strsplit(D, ",")[[1]]), n, s+1)[,-1]
    if(ShowCrit) cat("CD2 =", DesignEval(D, "CD2"),
                     "MD2 =", DesignEval(D, "MD2"),
                     "Maximin =", DesignEval(D, "maximin"), "\n")
    return (as.data.frame(D))
  }
}

panel.bar <- function(x, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  qx = table(x)
  par(usr = c(0, length(qx), 0, max(qx)*1.5))
  barplot(qx, width=1, space = 0, col=5, axes=FALSE, add = TRUE)
}

panel.scatter <- function(x, y, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  tmp = as.data.frame(table(x,y))
  points(tmp$x, tmp$y, cex=tmp$Freq/max(tmp$Freq), pch=19, col=4, xpd=TRUE)
  grid()
}

panel.heatmap <- function(x, y, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  tmp = table(x,y)
  adj = 1/2/(dim(tmp)-1)
  par(usr = c(-adj[1], 1+adj[1], -adj[2], 1+adj[2]))
  colmap = adjustcolor(terrain.colors(length(unique(tmp))), alpha.f=0.8)
  if (length(colmap)==1) colmap = adjustcolor(terrain.colors(3), alpha.f=0.8)[2]
  image(tmp, col=rev(colmap), add=TRUE)
}

DesignPairPlot <- function(D, Diag=FALSE) {
  if (Diag==TRUE) { pairs(D, lower.panel = panel.scatter, upper.panel = panel.heatmap, diag.panel = panel.bar, cex.labels=1.2) }
  else {
    pairs(D, lower.panel = panel.scatter, upper.panel = panel.heatmap, cex.labels=1.2)
  }
}

UniTracePlot <- function(output, skip=0) {
  x = seq(skip+1, length(output$criterion_list))
  y = output$criterion_list[x]
  plot(x, y, type="l",
       xlab="Iteration", ylab="Criterion")
  abline(v=x[which.min(y)], lty=2, col=2)
  abline(h=min(y), lty=2, col=2)
  title(main=paste("Best Value = ", round(min(y),6), " in ", round(output$time_consumed,2), " sec", sep=""), cex.main=1)
}

GenAUD_MS <- function(X0, n, crit="MD2", maxiter=30, nshoot = 5, vis=FALSE)
{
  crit_list = c()
  shoot_idx = c()
  bestcrit = 1e10
  for (i in 1:nshoot)
  {
    list0 = GenAUD(X0=X0, n=n, crit=crit, maxiter = maxiter)
    crit_list = c(crit_list, list0$criterion_list)
    shoot_idx = c(shoot_idx, length(list0$criterion_list))
    tmp = DesignEval(list0$final_design, crit = crit)
    if (tmp < bestcrit) {
      bestcrit = tmp
      bestdesign = list0$final_design
    }
  }
  if(vis) {
    par(mar=rep(2,4))
    plot(crit_list, type = "l", xlab = "", ylab = "")
    abline(v = which.min(crit_list), col=2, lty = 2)
    # abline(h = min(crit_list), col=2, lty = 2)
    if (nshoot>1) abline(v = cumsum(shoot_idx)[1:(nshoot-1)], col=4, lty=2)
    }
  return(as.data.frame(bestdesign))
}

GenUD_MS <- function(n, s, q, crit="MD2", maxiter=30, nshoot = 5, vis=FALSE)
{
  crit_list = c()
  shoot_idx = c()
  bestcrit = 1e10
  for (i in 1:nshoot)
  {
    list0 = GenUD(n=n,s=s,q=q, crit=crit, maxiter = maxiter)
    crit_list = c(crit_list, list0$criterion_list)
    shoot_idx = c(shoot_idx, length(list0$criterion_list))
    tmp = DesignEval(list0$final_design, crit = crit)
    if (tmp < bestcrit) {
      bestcrit = tmp
      bestdesign = list0$final_design
    }
  }
  if(vis) {
    par(mar=rep(2,4))
    plot(crit_list, type = "l", xlab = "", ylab = "")
    abline(v = which.min(crit_list), col=2, lty = 2)
    # abline(h = min(crit_list), col=2, lty = 2)
    if (nshoot>1)  abline(v = cumsum(shoot_idx)[1:(nshoot-1)], col=4, lty=2)
  }
  return(as.data.frame(bestdesign))
}
