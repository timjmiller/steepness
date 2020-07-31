library(Hmisc)
library(RColorBrewer)
library(sfsmisc)
temp = function(mod,y) cbind(mod$rep$log_FMSY, mod$rep$log_SSB_MSY, mod$rep$log_MSY, mod$rep$log_SPR0, mod$rep$log_SPR_MSY, mod$rep$log_SR_R0, mod$rep$logit_SR_h,
  mod$rep$log_SR_a, mod$rep$log_SR_b)
temp = temp(bh0)
temp[,c(1:6,8:9)] = exp(temp[,c(1:6,8:9)])
temp[,7] = 0.2 + 0.8/(1 + exp(-temp[,7]))
rownames(temp) = c(bh0$years)
colnames(temp) = c(
  "$F_{\\text{MSY}}$",
  "$S_{\\text{MSY}}$ (mt)",
  "MSY (mt)",
  "$\\phi_0$",
  "$\\phi_{F_{\\text{MSY}}}$",
  "$R_0 ~ (10^3)$",
  "$h$",
  "$\\alpha$",
  "$\\beta$")
temp = as.data.frame(temp)
temp[[9]] = pretty10exp(temp[[9]], digits=3, lab.type = "latex", lab.sep = "times")
temp = temp[,c(4,6,7,1:3,5,8:9)]
temp[,7]/temp[,1]
x = format.df(temp[,1:7], big.mark = " ", cdec = c(2,0,2,2,0,0,2), numeric.dollar = FALSE)
x = latex(x, file = 'BH_SR_MSY.tex', rowlabel = '', table.env = FALSE, rowlabel.just = "l", col.just = rep("r", 7))

temp = function(mod,y) cbind(mod$rep$log_FMSY, mod$rep$log_SSB_MSY, mod$rep$log_MSY, mod$rep$log_SPR0, mod$rep$log_SPR_MSY, mod$rep$log_SR_R0, mod$rep$logit_SR_h,
  mod$rep$log_SR_a, mod$rep$log_SR_b)
temp = temp(rick0)
temp[,c(1:6,8:9)] = exp(temp[,c(1:6,8:9)])
temp[,7] = 0.2 + exp(temp[,7])
rownames(temp) = c(rick0$years)
colnames(temp) = c(
  "$F_{\\text{MSY}}$",
  "$S_{\\text{MSY}}$ (mt)",
  "MSY (mt)",
  "$\\phi_0$",
  "$\\phi_{F_{\\text{MSY}}}$",
  "$R_0 ~ (10^3)$",
  "$h$",
  "$\\alpha$",
  "$\\beta$")
temp = as.data.frame(temp)
temp[[9]] = pretty10exp(temp[[9]], digits=3, lab.type = "latex", lab.sep = "times")
temp = temp[,c(4,6,7,1:3,5,8:9)]
temp[,7]/temp[,1]
x = format.df(temp[,1:7], big.mark = " ", cdec = c(2,0,2,2,0,0,2), numeric.dollar = FALSE)
x = latex(x, file = 'Ricker_SR_MSY.tex', rowlabel = '', table.env = FALSE, rowlabel.just = "l", col.just = rep("r", 7))
  
palette.fn <- function(n) colorRampPalette(c("dodgerblue","green","red"), space = "Lab")(n)
  
cairo_pdf("snemaytf_biology.pdf", family = "Times", height = 10, width = 10)
par(mfcol = c(2,2), mar = c(1,5,1,1), oma = c(4,0,1,0))
plot(1:6, bh0$env$data$mature[1,], type = 'n', lwd = 2, ylab = "Proportion mature", xlab = "", axes = FALSE, cex.lab = 1.3)
axis(1, labels = FALSE, lwd= 2, cex.axis = 1.3)
axis(2, lwd = 2, cex.axis = 1.3)
box(lwd = 2)
grid(col = gray(0.7))
lines(1:6, bh0$env$data$mature[1,], lwd = 2)

plot(1:6, bh0$rep$MAA[1,], type = 'n', lwd = 2, ylab = bquote(italic(M)), xlab = "", axes = FALSE, cex.lab = 1.3)
axis(1, at = 1:6, labels = bh0$ages.lab, lwd= 2, cex.axis = 1.3)
axis(2, lwd = 2, cex.axis = 1.3)
box(lwd = 2)
grid(col = gray(0.7))
lines(1:6, bh0$rep$MAA[1,], lwd = 2)
mtext(side = 1, "Age (years)", line = 2.5, cex = 1.3)


matplot(bh0$years, bh0$env$data$waa[bh0$env$data$waa_pointer_ssb,,], type = 'n', lty = 1, col = palette.fn(6), xlab = "", ylab = "Mass (kg)", axes = FALSE, cex.lab = 1.3)
axis(1, labels = FALSE, lwd= 2, cex.axis = 1.3)
axis(2, lwd = 2, cex.axis = 1.3)
box(lwd = 2)
grid(col = gray(0.7))
matplot(bh0$years, bh0$env$data$waa[bh0$env$data$waa_pointer_ssb,,], type = 'l', lty = 1, col = palette.fn(6), add = TRUE)
legend("topright", col = palette.fn(6), title = "Age (years)", legend = bh0$ages.lab,  lty = 1, box.lwd = 2) 

plot(bh0$years, exp(bh0$rep$log_SPR0), type = 'n', lwd = 2, ylab = bquote(italic("\u3D5")[0]), xlab = "", axes = FALSE, cex.lab = 1.3)
axis(1, lwd= 2, cex.axis = 1.3)
axis(2, lwd = 2, cex.axis = 1.3)
box(lwd = 2)
grid(col = gray(0.7))
lines(bh0$years, exp(bh0$rep$log_SPR0), lwd = 2)
mtext(side = 1, "Year", line = 2.5, cex = 1.3)
dev.off()

cairo_pdf("snemaytf_selectivity.pdf", family = "Times", height = 10, width = 10)
x = bh0$rep$FAA_tot/apply(bh0$rep$FAA_tot,1,max)
par(cex.axis = 1.3, cex.lab = 1.3, lwd = 2)
plot(1:6, x[1,], type = "l", axes = FALSE, ylab = "Selectivity", xlab = "Age (years)", ylim = c(0,1))
axis(1, at = 1:6, labels = bh0$ages.lab, lwd= 2)
axis(2, lwd = 2)
box()
grid(col = gray(0.7))
dev.off()

SR.fn = function(bh=bh0, rick=rick0, bhy=1, ricky=1)
{
  #Beverton-Holt
  a.bh = exp(bh$rep$log_SR_a[bhy])
  b.bh = exp(bh$rep$log_SR_b[bhy])
  S.bh = seq(0,max(bh$rep$SSB),10)
  phi0 = exp(bh$rep$log_SPR0)
  R.bh = a.bh * S.bh/(1 + b.bh * S.bh)

  #Ricker
  a.r = exp(rick$rep$log_SR_a[ricky])
  b.r = exp(rick$rep$log_SR_b[ricky])
  S.r = seq(0,max(rick$rep$SSB),10)
  R.r = a.r * S.r * exp(-b.r*S.r)

  ymax = max(R.bh, R.r, bh$rep$NAA[-1,1], rick$rep$NAA[-1,1])
  
  R0 = (a.bh-1/phi0)/b.bh
  S0 = R0 * phi0
  h = a.bh*phi0/(4 + a.bh*phi0)
  n = length(h)
  par(mfrow = c(1,2), mar = c(2,1,1,1), oma = c(3,4,1,0))
  plot(S.bh,R.bh, type = 'l', xlab = "", ylab = "", ylim = c(0, ymax), cex.lab= 1.5, lwd = 2, axes = FALSE)
  grid(col = gray(0.7))
  x = bquote(paste(italic("\u3B1") == .(round(a.bh,2)) ~ "and" ~ italic("\u3B2")))
  y = pretty10exp(b.bh, digits=3)
  z = bquote(paste(.(x[[2]]) == .(y[[1]])))
  legend("topright", legend = z, cex = 1.5, bg = "white", box.col = "transparent")
  axis(1, lwd= 2, cex.axis = 1.3)
  axis(2, lwd = 2, cex.axis = 1.3)
  box(lwd = 2)
  palette.fn <- colorRampPalette(c("dodgerblue","green","red"), space = "Lab")
  mypalette = function(n) palette.fn(n)
  
  cols = mypalette(length(phi0))
  points(S0, R0, pch = 19, col = cols)
  points(bh$rep$SSB[-length(phi0)], bh$rep$NAA[-1,1], pch = 17, col = cols[-length(phi0)])
  ind = seq(1,length(phi0),10)
  legend("bottomright", col = cols[ind], legend = bh$years[ind], pch = 19, box.lwd = 2) 
  mtext("Beverton-Holt", side = 3, outer = FALSE, line = 0, cex = 1.5)
  
  
  R0 = log(a.r*phi0)/(b.r*phi0)
  S0 = R0 * phi0
  h = 0.2 * (a.r * phi0)^(0.8)
  n = length(h)
  plot(S.r,R.r, type = 'l', xlab = "", ylab = "", ylim = c(0, ymax), cex.lab= 1.5, lwd = 2, axes = FALSE)
  grid(col = gray(0.7))
  x = bquote(paste(italic("\u3B1") == .(round(a.r,2)) ~ "and" ~ italic("\u3B2")))
  y = pretty10exp(b.r, digits=3)
  z = bquote(paste(.(x[[2]]) == .(y[[1]])))
  legend("topright", legend = z, cex = 1.5, bg = "white", box.col = "transparent")
  axis(1, lwd= 2, cex.axis = 1.3)
  axis(2, labels = FALSE, lwd = 2, cex.axis = 1.3)
  box(lwd = 2)
  
  cols = mypalette(length(phi0))
  points(S0, R0, pch = 19, col = cols)
  points(rick$rep$SSB[-length(phi0)], rick$rep$NAA[-1,1], pch = 17, col = cols[-length(phi0)])
  ind = seq(1,length(phi0),10)
  legend("bottomright", col = cols[ind], legend = rick$years[ind], pch = 19, box.lwd = 2) 
  mtext("Ricker", side = 3, outer = FALSE, line = 0, cex = 1.5)
  
  mtext(side = 2, bquote(paste(italic(R), " (", 10^3, ")")), line = 2, cex = 1.5, outer = TRUE)
  mtext("Spawning Biomass (mt)", side = 1, outer = TRUE, line = 1.5, cex = 1.5)
}
cairo_pdf("annual_R0_S0.pdf", family = "Times", height = 8, width = 12)
SR.fn()
dev.off()

get_SPR = function(F, M, sel, mat, waassb, fracyrssb, at.age = FALSE)
{
  n_ages = length(sel)
  SPR = numeric()
  n = 1
  F = F * sel
  Z = F + M
  for(a in 1:(n_ages-1))
  {
    SPR[a] = n[a] * mat[a] * waassb[a] * exp(-fracyrssb * Z[a])
    n[a+1] = n[a] * exp(-Z[a])
  }
  n[n_ages] = n[n_ages]/(1-exp(-Z[n_ages]))
  SPR[n_ages] = n[n_ages] * mat[n_ages] * waassb[n_ages] * exp(-fracyrssb * Z[n_ages])
  if(at.age) return(SPR)
  else return(sum(SPR))
}

#-------Y/R -----------------------------
get_YPR = function(F, M, sel, waacatch, at.age = FALSE)
{
  n_ages = length(sel)
  YPR = numeric()
  n = 1
  F = F * sel
  Z = F + M
  for(a in 1:(n_ages-1))
  {
    YPR[a] = n[a] * F[a] * waacatch[a] * (1.0 - exp(-Z[a]))/Z[a]
    n[a+1] = n[a] * exp(-Z[a])
  }
  n[n_ages] = n[n_ages]/(1 - exp(-Z[n_ages]))
  YPR[n_ages] = n[n_ages] * F[n_ages] * waacatch[n_ages] * (1.0 - exp(-Z[n_ages]))/Z[n_ages]
  if(at.age) return(YPR)
  else return(sum(YPR))
}

f.of.phi0.fn = function(phi0 = seq(0.5,2.5,0.01))
{
  phi0y = exp(bh0$rep$log_SPR0)
  a = exp(bh0$rep$log_SR_a[1])
  b = exp(bh0$rep$log_SR_b[1])
  h = a*phi0/(4 + a*phi0)
  hy = a*phi0y/(4 + a*phi0y)
  R0 = (a-1/phi0)/b
  R0y = (a-1/phi0y)/b

  a.r = exp(rick0$rep$log_SR_a[1])
  b.r = exp(rick0$rep$log_SR_b[1])
  R0.r = log(a.r*phi0)/(b.r*phi0)
  R0y.r = log(a.r*phi0y)/(b.r*phi0y)
  h.r = 0.2 * (a.r * phi0)^(0.8)
  hy.r = 0.2 * (a.r * phi0y)^(0.8)
  range_h = range(c(h,hy,h.r, hy.r))
  range_R0 = range(c(R0,R0y,R0.r, R0y.r))
  cols <- colorRampPalette(c("dodgerblue","green","red"), space = "Lab")(length(phi0y))
  cols = col2rgb(cols)
  cols = rgb(cols[1,],cols[2,],cols[3,], alpha = 150, maxColorValue = 255)
  
  par(mfcol = c(2,2), mar = c(1,1,1,1), oma = c(4,4,1,0))
  plot(phi0, h, type = 'n', xlab = "", ylab = "", cex.lab = 1.5, axes = FALSE, ylim = range_h)
  grid(col = gray(0.7))
  x = bquote(paste(italic("\u3B1") == .(round(a,2)) ~ "and" ~ italic("\u3B2")))
  y = pretty10exp(b, digits=3)
  z = bquote(paste(.(x[[2]]) == .(y[[1]])))
  legend("topleft", legend = z, cex = 1.5, bg = "white", box.col = "transparent")
  box(lwd = 2)
  axis(1, labels =FALSE, lwd = 2, cex.axis = 1.5)
  axis(2, lwd = 2, cex.axis = 1.5)
  lines(phi0,h, lwd  = 2)
  points(phi0y, hy, pch = 19, col = cols)
  mtext(side = 2, bquote(italic(h)), cex = 1.5, line = 2.5)
  #text(par()$usr[1]*1.05, par()$usr[4]*0.98, z, cex = 1.5, adj = c(0,1), bg = "white")
  mtext("Beverton-Holt", side = 3, outer = FALSE, line = 0, cex = 1.5)
  
  plot(phi0, R0, type = 'n', xlab = "", ylab = "", cex.lab = 1.5, axes = FALSE, ylim = range_R0)
  grid(col = gray(0.7))
  box(lwd = 2)
  axis(1, lwd = 2, cex.axis = 1.5)
  axis(2, lwd = 2, cex.axis = 1.5)
  lines(phi0,R0, lwd  = 2)
  points(phi0y, R0y, pch = 19, col = cols)
  ind = seq(1,length(phi0y),10)
  legend("bottomright", col = cols[ind], legend = bh0$years[ind], pch = 19, box.lwd = 2, cex = 1.5, bg = "white") 
  mtext(side = 2, bquote(italic(R)[0] ~ (10^3)), cex = 1.5, line = 2.5)

  plot(phi0, h.r, type = 'n', xlab = "", ylab = "", cex.lab = 1.5, axes = FALSE, ylim = range_h)
  grid(col = gray(0.7))
  x = bquote(paste(italic("\u3B1") == .(round(a.r,2)) ~ "and" ~ italic("\u3B2")))
  y = pretty10exp(b.r, digits=3)
  z = bquote(paste(.(x[[2]]) == .(y[[1]])))
  legend("topleft", legend = z, cex = 1.5, bg = "white", box.col = "transparent")
  box(lwd = 2)
  axis(1, labels =FALSE, lwd = 2, cex.axis = 1.5)
  axis(2, labels =FALSE, lwd = 2, cex.axis = 1.5)
  lines(phi0,h.r, lwd  = 2)
  points(phi0y, hy.r, pch = 19, col = cols)
  mtext("Ricker", side = 3, outer = FALSE, line = 0, cex = 1.5)
  
  plot(phi0, R0.r, type = 'n', xlab = "", ylab = "", cex.lab = 1.5, axes = FALSE, ylim = range_R0)
  grid(col = gray(0.7))
  box(lwd = 2)
  axis(1, lwd = 2, cex.axis = 1.5)
  axis(2, labels =FALSE, lwd = 2, cex.axis = 1.5)
  lines(phi0,R0.r, lwd  = 2)
  points(phi0y, R0y.r, pch = 19, col = cols)
  ind = seq(1,length(phi0y),10)
  mtext(side = 1, bquote(italic("\u3D5")[0]), cex= 1.5, outer = TRUE, line = 2)
  
}
cairo_pdf("steepness_R0_phi0.pdf", family = "Times", height = 8, width = 8)
f.of.phi0.fn()
dev.off()

S0_R0_curves.fn = function(years, bhmod = bh0, rickmod = rick0)
{
  mod = bhmod
  phi0 = exp(mod$rep$log_SPR0)
  n = length(phi0)
  if(missing(years)) years = 1:n
  palette.fn <- colorRampPalette(c("dodgerblue","green","red"), space = "Lab")
  mypalette = function(n) palette.fn(n)
  cols = mypalette(length(phi0))
  
  h = 0.2 + 0.8/(1+ exp(-mod$rep$logit_SR_h[n]))
  R0a = exp(mod$rep$log_SR_R0[n])
  S0a = R0a * phi0
  aa = 4*h/((1-h)*phi0)
  ba = (5*h -1)/(R0a*phi0*(1-h))
  R0b = exp(mod$rep$log_SR_R0)
  S0b = R0b[n] * phi0[n]
  ab = 4*h/((1-h)*phi0)
  bb = (5*h -1)/(S0b*(1-h))
  S = seq(0,max(mod$rep$SSB, S0a, S0b),10)
  Ra = sapply(years, function(x) aa[x]* S/(1 + ba[x] * S))
  Rb = sapply(years, function(x) ab[x]* S/(1 + bb * S))
  maxR = max(Ra,Rb)
  maxS = max(S)
  print(range(R0a))
  print(range(R0b))
  print(range(S0a))
  print(range(S0b))
  
  mod = rickmod
  h.r = 0.2 + exp(mod$rep$logit_SR_h[n])
  R0a.r = exp(mod$rep$log_SR_R0[n])
  S0a.r = R0a.r * phi0
  aa.r = (5*h.r)^(5/4)/phi0
  ba.r =  (5/4)*log(5*h.r)/(R0a.r*phi0)
  R0b.r = exp(mod$rep$log_SR_R0)
  S0b.r = R0b.r[n] * phi0[n]
  ab.r = (5*h.r)^(5/4)/phi0
  bb.r = (5/4)*log(5*h.r)/S0b.r
  S.r = seq(0,max(mod$rep$SSB, S0a.r, S0b.r),10)
  Ra.r = sapply(years, function(x) aa.r[x]* S.r * exp(-ba.r[x] * S.r))
  Rb.r = sapply(years, function(x) ab.r[x]* S.r * exp(-bb.r * S.r))
  maxR = max(maxR,Ra.r,Rb.r)
  maxS = max(maxS,S.r)
  print(range(R0a.r))
  print(range(R0b.r))
  print(range(S0a.r))
  print(range(S0b.r))
  
  par(mfcol = c(2,2), mar = c(1,1,1,1), oma = c(4,4,1,1))
  plot(S,Ra[,1], xlab = "", ylab = "", ylim = c(0,maxR), xlim = c(0,maxS), col = cols[1], type = 'n', axes = FALSE)
  grid(col = gray(0.7))
  x = bquote(paste(italic(h) == .(round(h,2)) ~ "and" ~ italic(R)[0]))
  y = pretty10exp(R0a*1e3, digits=3)
  z = bquote(paste(.(x[[2]]) == .(y[[1]])))
  legend("topleft", legend = z, cex = 1, bg = "white", box.col = "transparent", xjust = 0, yjust = 0.5)
  axis(1, labels = FALSE, lwd = 2)
  axis(2, lwd = 2)
  box(lwd = 2)
  mtext(bquote(paste(italic(R), " (", 10^3, ")")), side = 2, outer = TRUE, cex = 1.5, line = 2)
  sapply(1:length(years), function(x) lines(S, Ra[,x], col = cols[years[x]]))
  points(S0a, rep(R0a,n), pch = 19, col = cols)
  mtext("Beverton-Holt", side = 3, line = 0, cex = 1.5, outer = FALSE)

  plot(S,Rb[,1], xlab = "", ylab = "", ylim = c(0,maxR), xlim = c(0,maxS), col = cols[1], type = 'n', axes = FALSE)
  grid(col = gray(0.7))
  x = bquote(paste(italic(h) == .(round(h,2)) ~ "and" ~ italic(S)[0]))
  y = pretty10exp(S0b, digits=3)
  z = bquote(paste(.(x[[2]]) == .(y[[1]]) ~ "(mt)"))
  legend("topleft", legend = z, cex = 1, bg = "white", box.col = "transparent", xjust = 0, yjust = 0.5)
  box(lwd = 2)
  axis(1, lwd = 2)
  axis(2, lwd = 2)
  sapply(1:length(years), function(x) lines(S, Rb[,x], col = cols[years[x]]))
  points(rep(S0b,n), S0b/phi0, pch = 19, col = cols)
  ind = seq(1,length(phi0),10)
  legend("bottomright", col = cols[ind], legend = mod$years[ind], pch = 19, box.lwd = 2) 

  plot(S.r,Ra.r[,1], xlab = "", ylab = "", ylim = c(0,maxR), xlim = c(0,maxS), col = cols[1], type = 'n', axes = FALSE)
  grid(col = gray(0.7))
  x = bquote(paste(italic(h) == .(round(h.r,2)) ~ "and" ~ italic(R)[0]))
  y = pretty10exp(R0a.r*1e3, digits=3)
  z = bquote(paste(.(x[[2]]) == .(y[[1]])))
  legend("topright", legend = z, cex = 1, bg = "white", box.col = "transparent", xjust = 1, yjust = 0.5)
  axis(1, labels = FALSE, lwd = 2)
  axis(2, labels = FALSE, lwd = 2)
  box(lwd = 2)
  mtext(bquote(paste(italic(R), " (", 10^3, ")")), side = 2, outer = TRUE, cex = 1.5, line = 2)
  sapply(1:length(years), function(x) lines(S.r, Ra.r[,x], col = cols[years[x]]))
  points(S0a.r, rep(R0a.r,n), pch = 19, col = cols)
  mtext("Ricker", side = 3, line = 0, cex = 1.5, outer = FALSE)

  plot(S.r,Rb.r[,1], xlab = "", ylab = "", ylim = c(0,maxR), xlim = c(0,maxS), col = cols[1], type = 'n', axes = FALSE)
  grid(col = gray(0.7))
  x = bquote(paste(italic(h) == .(round(h.r,2)) ~ "and" ~ italic(S)[0]))
  y = pretty10exp(S0b.r, digits=3)
  z = bquote(paste(.(x[[2]]) == .(y[[1]]) ~ "(mt)"))
  legend("topright", legend = z, cex = 1, bg = "white", box.col = "transparent", xjust = 1, yjust = 0.5)
  box(lwd = 2)
  axis(1, lwd = 2)
  axis(2, labels = FALSE, lwd = 2)
  sapply(1:length(years), function(x) lines(S.r, Rb.r[,x], col = cols[years[x]]))
  points(rep(S0b.r,n), S0b.r/phi0, pch = 19, col = cols)
  mtext("Spawning Biomass (mt)", side = 1, outer = TRUE, cex = 1.5, line = 2)
  ind = seq(1,length(phi0),10)
}
cairo_pdf("steepness_R0_S0_curves.pdf", family = "Times", height = 8, width = 8)
S0_R0_curves.fn()
dev.off()

SR.fn = function(phi0y = 1973, hy = 1999, mod = bh0, mod_2 = bh4, maxS = 30000, maxR = 30000, fn = "bh")
{
  hy = which(mod$years == hy) #1999
  phi0y = which(mod$years == phi0y) #1973
  phi0 = exp(mod$rep$log_SPR0)
  R0 = exp(mod$rep$log_SR_R0[phi0y]) #correct for 1973.
  RMSY = exp(mod$rep$log_R_MSY[phi0y]) #correct for 1973
  phiMSY = exp(mod$rep$log_SPR_MSY[phi0y]) #correct for 1973
  SSBMSY = exp(mod$rep$log_SSB_MSY[phi0y]) #correct for 1973
  a = exp(mod$rep$log_SR_a[1]) #all the same and estimated.
  b = exp(mod$rep$log_SR_b[1]) #all the same and estimated.
  if(fn == "bh") h = 0.2 + 0.8/(1+ exp(-mod$rep$logit_SR_h[hy])) #the assumed steepness
  if(fn == "ricker") h = 0.2 + exp(mod$rep$logit_SR_h[hy])
  print(h)
  if(fn == "bh") h = 0.2 + 0.8/(1+ exp(-mod_2$rep$logit_SR_h[phi0y])) #all the same, unestimated, and equal to the assumed steepness
  if(fn == "ricker") h = 0.2 + exp(mod_2$rep$logit_SR_h[phi0y])
  print(h)
  R0. = exp(mod_2$rep$log_SR_R0[phi0y]) #all the same and estimated.
  RMSY. = exp(mod_2$rep$log_R_MSY[phi0y]) #all the same and estimated.
  phiMSY. = exp(mod_2$rep$log_SPR_MSY[phi0y]) #all the same and estimated.
  SSBMSY. = exp(mod_2$rep$log_SSB_MSY[phi0y])
  a. = 4*h/((1-h)*phi0[phi0y])
  print(a.)
  print(exp(mod_2$rep$log_SR_a[phi0y]))
  a. = exp(mod_2$rep$log_SR_a[phi0y])
  print(a.)
  b. = (5*h -1)/(R0.*phi0[phi0y]*(1-h))
  print(b.)
  b. = exp(mod_2$rep$log_SR_b[phi0y])
  print(b.)
  S = seq(0,maxS,10)
  if(fn == "bh") 
  {
    R = a * S/(1 + b * S)
    R. = a. * S/(1 + b. * S)
  }
  if(fn == "ricker") 
  {
    R = a * S * exp(-b * S)
    R. = a. * S * exp(-b. * S)
  }
  plot(S,R, type = 'n', axes = FALSE, ylab = "", xlab = "", xlim = c(0,maxS), ylim = c(0,maxR))#c(0, max(R,R.)))
  grid(col = gray(0.7))
  lines(S,R, lwd = 2)
  points(R0*phi0[phi0y], R0, pch = 19)
  points(RMSY*phiMSY, RMSY, pch = 17)
  lines(S,R., col = "red", lwd = 2)
  points(R0.*phi0[phi0y], R0., col = "red",, pch = 19)
  points(RMSY.*phiMSY., RMSY., col = "red",, pch = 17)
  abline(0, 1/phi0[phi0y], col = gray(0.7), lwd = 2)
  abline(0, 1/phiMSY, lwd = 2, lty = 2)
  abline(0, 1/phiMSY., col = "red",, lwd = 2, lty = 2)
  print(c(R0,R0*phi0[phi0y]))
  print(c(R0.,R0.*phi0[phi0y]))
  print(c(RMSY,RMSY*phiMSY))
  print(c(RMSY.,RMSY.*phiMSY.))
  text(par()$xaxp[2]*0.9,par()$yaxp[1]*1.1, bquote(paste(italic(h) == .(round(h,2)))), cex = 1.5)
}

cairo_pdf("steepness_assumptions_plots.pdf", family = "Times", height = 8, width = 12)
par(mfrow = c(1,2), oma = c(2,5,2,0), mar = c(3,1,1,1))
SR.fn(maxR = 40000, maxS = 50000)
  axis(1, lwd = 2, cex.lab = 1.5)
  axis(2, lwd = 2, cex.lab = 1.5)
  box(lwd = 2)
mtext(side = 3, "Beverton-Holt", cex = 1.5, line = 0)
SR.fn(mod = rick0, mod_2 = rick4, fn = "ricker", maxR = 40000, maxS = 50000)
  axis(1, lwd = 2, cex.lab = 1.5)
  axis(2, labels = FALSE, lwd = 2, cex.lab = 1.5)
  box(lwd = 2)
mtext(side = 3, "Ricker", cex = 1.5, line = 0)
mtext(side = 2, bquote(paste(italic(R), " (", 10^3, ")")), line = 2, cex = 1.5, outer = TRUE)
mtext(side = 1, "Spawning Biomass (mt)", line = 0, cex = 1.5, outer = TRUE)
dev.off()

temp = function(mod,y) c(mod$rep$log_FMSY[y], mod$rep$log_SSB_MSY[y], mod$rep$log_MSY[y], mod$rep$log_SPR0[y], mod$rep$log_SR_R0[y], mod$rep$logit_SR_h[y],
  mod$rep$log_SR_a[y], mod$rep$log_SR_b[y])
temp = cbind(temp(bh0,1), temp(bh4,1), temp(rick0,1), temp(rick4,1))
temp[c(1:5,7:8),] = exp(temp[c(1:5,7:8),])
temp[6,1:2] = 0.2 + 0.8/(1 + exp(-temp[6,1:2]))
temp[6,3:4] = 0.2 + exp(temp[6,3:4])
temp = t(temp)
rownames(temp) = c(
  "Beverton-Holt ($\\alpha$,$\\beta$)", 
  "Beverton-Holt ($h$ from 1999 with $\\phi_0$ from 1973)",
  "Ricker ($\\alpha$,$\\beta$)", 
  "Ricker ($h$ from 1999 with $\\phi_0$ from 1973)")
colnames(temp) = c(
  "$F_{\\text{MSY}}$",
  "$S_{\\text{MSY}}$ (mt)",
  "MSY (mt)",
  "$\\phi_0$",
  "$R_0 ~ (10^3)$",
  "$h$",
  "$\\alpha$",
  "$\\beta$")
temp = as.data.frame(temp)
temp[[8]] = pretty10exp(temp[[8]], digits=3, lab.type = "latex", lab.sep = "times")
x = format.df(temp, big.mark = " ", cdec = c(2,0,0,2,0,2,2), numeric.dollar = FALSE)
x = latex(x, file = 'SR_BRP_table.tex', rowlabel = '', table.env = FALSE, rowlabel.just = "l", col.just = rep("r", 8))


dynB0.fn = function(alpha, beta, N1, mat, M, wt, fracyrssb, recdevs, bh=FALSE)
{
  ny = length(recdevs)
  na = length(N1)
  SSB = sum(N1*mat[1,]*wt[1,]*exp(-M[1,]*fracyrssb[1]))
  NAA = matrix(NA, ny, na)
  NAA[1,] = N1
  for(i in 1:(ny-1)) 
  {
    if(bh) NAA[i+1,1] = exp(recdevs[i+1]) * alpha * SSB[i]/(1 + beta*SSB[i])
    else NAA[i+1,1] = exp(recdevs[i+1]) * alpha * SSB[i] * exp(-beta * SSB[i])
    NAA[i+1,-1] = NAA[i,1:(na-1)] * exp(-M[i, 1:(na-1)])
    NAA[i+1,na] = NAA[i+1,na] + NAA[i,na]*exp(-M[i,na])
    SSB[i+1] = sum(NAA[i+1,]*mat[i+1,]*wt[i+1,]*exp(-M[i+1,]*fracyrssb[i+1]))
  }
  return(list(NAA=NAA, SSB = SSB))
}

make.dynb0.res = function(mod = bh2)
{
  h_y = mod$env$data$steepness_year
  bh = mod$env$data$recruit_model == 3
  res = list()
  waa = mod$env$data$waa
  for(i in 1:3) waa[i,,] = rep(waa[i,h_y,], each = length(mod$years))
  #first just see what happens when there are no rec devs, and constant waa
  res$nodev_constwaa = dynB0.fn(
    alpha = exp(mod$rep$log_SR_a[h_y]),
    beta = exp(mod$rep$log_SR_b[h_y]),
    N1 = mod$rep$NAA[1,], 
    mat = mod$env$data$mature,
    M = mod$rep$MAA,
    wt = waa[mod$env$data$waa_pointer_ssb,,],
    fracyrssb = mod$env$data$fracyr_SSB,
    recdevs = rep(0, length(mod$years)),
    bh = bh)
#now with realized recdevs
  res$constwaa = dynB0.fn(
    alpha = exp(mod$rep$log_SR_a[h_y]),
    beta = exp(mod$rep$log_SR_b[h_y]),
    N1 = mod$rep$NAA[1,], 
    mat = mod$env$data$mature,
    M = mod$rep$MAA,
    wt = waa[mod$env$data$waa_pointer_ssb,,],
    fracyrssb = mod$env$data$fracyr_SSB,
    recdevs = log(mod$rep$NAA[,1]) - log(mod$rep$pred_NAA[,1]) - 0.5*exp(2*mod$parList$log_NAA_sigma[1]),
    bh = bh)
  #now no rec devs, but variable waa
  res$nodev = dynB0.fn(
    alpha = exp(mod$rep$log_SR_a[h_y]),
    beta = exp(mod$rep$log_SR_b[h_y]),
    N1 = mod$rep$NAA[1,], 
    mat = mod$env$data$mature,
    M = mod$rep$MAA,
    wt = mod$env$data$waa[mod$env$data$waa_pointer_ssb,,],
    fracyrssb = mod$env$data$fracyr_SSB,
    recdevs = rep(0, length(mod$years)),
    bh = bh)
  #now rec devs and variable waa
  res$all = dynB0.fn(
    alpha = exp(mod$rep$log_SR_a[h_y]),
    beta = exp(mod$rep$log_SR_b[h_y]),
    N1 = mod$rep$NAA[1,], 
    mat = mod$env$data$mature,
    M = mod$rep$MAA,
    wt = mod$env$data$waa[mod$env$data$waa_pointer_ssb,,],
    fracyrssb = mod$env$data$fracyr_SSB,
    recdevs = log(mod$rep$NAA[,1]) - log(mod$rep$pred_NAA[,1]) - 0.5*exp(2*mod$parList$log_NAA_sigma[1]),
    bh = bh)
  return(res)
}

dynb0.bh2 = make.dynb0.res(bh2)
dynb0.bh1 = make.dynb0.res(bh1)
dynb0.rick2 = make.dynb0.res(rick2)
dynb0.rick1 = make.dynb0.res(rick1)

#what if we used 1973?
cairo_pdf("dynamic_B0_deterministic.pdf", family = "Times", height = 10, width = 10)
par(mfrow = c(2,2), mar = c(1,5,1,1), oma = c(4,0,1,3))
max.xy = c(max(dynb0.bh1[[1]]$SSB, dynb0.bh2[[1]]$SSB), max(dynb0.bh2[[1]]$NAA[,1], dynb0.bh1[[1]]$NAA[,1]))/1000
plot(bh2$years, dynb0.bh2[[1]]$SSB, axes = FALSE, ann = FALSE, type = "n", ylim = c(0, max.xy[1]))
axis(1, labels = FALSE, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7))
abline(h = exp(bh2$rep$log_SR_R0 + bh2$rep$log_SPR0)[bh2$env$data$steepness_year]/1000, lwd = 2)
abline(h = exp(bh1$rep$log_SR_R0 + bh1$rep$log_SPR0)[bh1$env$data$steepness_year]/1000, lwd = 2, col = "gray")
points(bh2$years, dynb0.bh2[[1]]$SSB/1000, col = "black", pch = 19)
points(bh2$years, dynb0.bh1[[1]]$SSB/1000, col = "gray", pch = 19)
mtext("Spawning Biomass (kmt)", side = 2, line = 2.5, cex = 1.5)

plot(bh2$years, dynb0.bh2[[1]]$NAA[,1]/1000, axes = FALSE, ann = FALSE, type = "n", ylim = c(0, max.xy[2]))
axis(1, labels = FALSE, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7))
abline(h = exp(bh2$rep$log_SR_R0)[bh2$env$data$steepness_year]/1000, lwd = 2)
abline(h = exp(bh1$rep$log_SR_R0)[bh1$env$data$steepness_year]/1000, lwd = 2, col = "gray")
points(bh2$years, dynb0.bh2[[1]]$NAA[,1]/1000, col = "black", pch = 19)
points(bh2$years, dynb0.bh1[[1]]$NAA[,1]/1000, col = "gray", pch = 19)
mtext(bquote("Recruits " ~ (10^3)), side = 2, line = 2.5, cex = 1.5)
mtext("Beverton-Holt", side = 4, line = 1, cex = 1.5)

max.xy = c(max(dynb0.rick1[[1]]$SSB, dynb0.rick2[[1]]$SSB), max(dynb0.rick2[[1]]$NAA[,1], dynb0.rick1[[1]]$NAA[,1]))/1000
plot(rick2$years, dynb0.rick2[[1]]$SSB, axes = FALSE, ann = FALSE, type = "n", ylim = c(0, max.xy[1]))
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7))
abline(h = exp(rick2$rep$log_SR_R0 + rick2$rep$log_SPR0)[rick2$env$data$steepness_year]/1000, lwd = 2)
abline(h = exp(rick1$rep$log_SR_R0 + rick1$rep$log_SPR0)[rick1$env$data$steepness_year]/1000, lwd = 2, col = "gray")
points(rick2$years, dynb0.rick2[[1]]$SSB/1000, col = "black", pch = 19)
points(rick2$years, dynb0.rick1[[1]]$SSB/1000, col = "gray", pch = 19)
mtext("Spawning Biomass (kmt)", side = 2, line = 2.5, cex = 1.5)

plot(rick2$years, dynb0.rick2[[1]]$NAA[,1]/1000, axes = FALSE, ann = FALSE, type = "n", ylim = c(0, max.xy[2]))
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7))
abline(h = exp(rick2$rep$log_SR_R0)[rick2$env$data$steepness_year]/1000, lwd = 2)
abline(h = exp(rick1$rep$log_SR_R0)[rick1$env$data$steepness_year]/1000, lwd = 2, col = "gray")
points(rick2$years, dynb0.rick2[[1]]$NAA[,1]/1000, col = "black", pch = 19)
points(rick2$years, dynb0.rick1[[1]]$NAA[,1]/1000, col = "gray", pch = 19)
mtext(bquote("Recruits " ~ (10^3)), side = 2, line = 2.5, cex = 1.5)
mtext("Ricker", side = 4, line = 1, cex = 1.5)
mtext("Year", side = 1, outer = TRUE, line = 2, cex = 1.5)
dev.off()

plot(bh1$years, bh0$rep$SSB/dynb0.bh1[[4]]$SSB)
plot(bh1$years, bh0$rep$SSB/exp(bh0$rep$log_SR_R0 + bh0$rep$log_SPR0))
plot(bh1$years, (bh0$rep$SSB/dynb0.bh1[[4]]$SSB)/(bh0$rep$SSB/exp(bh0$rep$log_SR_R0 + bh0$rep$log_SPR0)))
plot(bh1$years, dynb0.bh1[[4]]$SSB/exp(bh0$rep$log_SR_R0 + bh0$rep$log_SPR0))
plot(bh1$years, exp(bh0$rep$log_SR_R0 + bh0$rep$log_SPR0))
plot(bh1$years, dynb0.rick1[[4]]$SSB/exp(rick0$rep$log_SR_R0 + rick0$rep$log_SPR0))
plot(bh1$years, exp(rick0$rep$log_SR_R0 + rick0$rep$log_SPR0))

cairo_pdf("../tex/dynamic_B0_constwaa.pdf", family = "Times", height = 10, width = 10)
par(mfrow = c(2,2), mar = c(1,5,1,1), oma = c(4,0,1,3))
max.xy = c(max(dynb0.bh1[[2]]$SSB, dynb0.bh2[[2]]$SSB), max(dynb0.bh2[[2]]$NAA[,1], dynb0.bh1[[2]]$NAA[,1]))/1000
plot(bh2$years, dynb0.bh2[[2]]$SSB, axes = FALSE, ann = FALSE, type = "n", ylim = c(0, max.xy[1]))
axis(1, labels = FALSE, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7))
abline(h = exp(bh2$rep$log_SR_R0 + bh2$rep$log_SPR0)[bh2$env$data$steepness_year]/1000, lwd = 2)
abline(h = exp(bh1$rep$log_SR_R0 + bh1$rep$log_SPR0)[bh1$env$data$steepness_year]/1000, lwd = 2, col = "gray")
abline(h = median(dynb0.bh2[[2]]$SSB[-(1:6)])/1000, lwd = 2, lty = 2)
abline(h = median(dynb0.bh1[[2]]$SSB[-(1:6)])/1000, lwd = 2, lty = 2, col = "gray")
points(bh2$years, dynb0.bh2[[2]]$SSB/1000, col = "black", pch = 19)
points(bh2$years, dynb0.bh1[[2]]$SSB/1000, col = "gray", pch = 19)
mtext("Spawning Biomass (kmt)", side = 2, line = 2.5, cex = 1.5)

plot(bh2$years, dynb0.bh2[[2]]$NAA[,1]/1000, axes = FALSE, ann = FALSE, type = "n", ylim = c(0, max.xy[2]))
axis(1, labels = FALSE, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7))
abline(h = exp(bh2$rep$log_SR_R0)[bh2$env$data$steepness_year]/1000, lwd = 2)
abline(h = exp(bh1$rep$log_SR_R0)[bh1$env$data$steepness_year]/1000, lwd = 2, col = "gray")
abline(h = median(dynb0.bh2[[2]]$NAA[-(1:6),1])/1000, lwd = 2, lty = 2)
abline(h = median(dynb0.bh1[[2]]$NAA[-(1:6),1])/1000, lwd = 2, lty = 2, col = "gray")
points(bh2$years, dynb0.bh2[[2]]$NAA[,1]/1000, col = "black", pch = 19)
points(bh2$years, dynb0.bh1[[2]]$NAA[,1]/1000, col = "gray", pch = 19)
mtext(bquote("Recruits " ~ (10^3)), side = 2, line = 2.5, cex = 1.5)
mtext("Beverton-Holt", side = 4, line = 1, cex = 1.5)

max.xy = c(max(dynb0.rick1[[2]]$SSB, dynb0.rick2[[2]]$SSB), max(dynb0.rick2[[2]]$NAA[,1], dynb0.rick1[[2]]$NAA[,1]))/1000
plot(rick2$years, dynb0.rick2[[2]]$SSB, axes = FALSE, ann = FALSE, type = "n", ylim = c(0, max.xy[1]))
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7))
abline(h = exp(rick2$rep$log_SR_R0 + rick2$rep$log_SPR0)[rick2$env$data$steepness_year]/1000, lwd = 2)
abline(h = exp(rick1$rep$log_SR_R0 + rick1$rep$log_SPR0)[rick1$env$data$steepness_year]/1000, lwd = 2, col = "gray")
abline(h = median(dynb0.rick2[[2]]$SSB[-(1:6)])/1000, lwd = 2, lty = 2)
abline(h = median(dynb0.rick1[[2]]$SSB[-(1:6)])/1000, lwd = 2, lty = 2, col = "gray")
points(rick2$years, dynb0.rick2[[2]]$SSB/1000, col = "black", pch = 19)
points(rick2$years, dynb0.rick1[[2]]$SSB/1000, col = "gray", pch = 19)
mtext("Spawning Biomass (kmt)", side = 2, line = 2.5, cex = 1.5)

plot(rick2$years, dynb0.rick2[[2]]$NAA[,1]/1000, axes = FALSE, ann = FALSE, type = "n", ylim = c(0, max.xy[2]))
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7))
abline(h = exp(rick2$rep$log_SR_R0)[rick2$env$data$steepness_year]/1000, lwd = 2)
abline(h = exp(rick1$rep$log_SR_R0)[rick1$env$data$steepness_year]/1000, lwd = 2, col = "gray")
abline(h = median(dynb0.rick2[[2]]$NAA[-(1:6),1])/1000, lwd = 2, lty = 2)
abline(h = median(dynb0.rick1[[2]]$NAA[-(1:6),1])/1000, lwd = 2, lty = 2, col = "gray")
points(rick2$years, dynb0.rick2[[2]]$NAA[,1]/1000, col = "black", pch = 19)
points(rick2$years, dynb0.rick1[[2]]$NAA[,1]/1000, col = "gray", pch = 19)
mtext(bquote("Recruits " ~ (10^3)), side = 2, line = 2.5, cex = 1.5)
mtext("Ricker", side = 4, line = 1, cex = 1.5)
dev.off()

cairo_pdf("../tex/dynamic_B0_nodev.pdf", family = "Times", height = 10, width = 10)
par(mfrow = c(2,2), mar = c(1,5,1,1), oma = c(4,0,1,3))
max.xy = c(max(dynb0.bh1[[3]]$SSB, dynb0.bh2[[3]]$SSB), max(dynb0.bh2[[3]]$NAA[,1], dynb0.bh1[[3]]$NAA[,1]))/1000
plot(bh2$years, dynb0.bh2[[3]]$SSB, axes = FALSE, ann = FALSE, type = "n", ylim = c(0, max.xy[1]))
axis(1, labels = FALSE, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7))
abline(h = median(dynb0.bh2[[3]]$SSB[-(1:6)])/1000, lwd = 2, lty = 2)
abline(h = median(exp(bh0$rep$log_SR_R0 + bh0$rep$log_SPR0)[-(1:6)])/1000, lwd = 2, lty = 2, col = "gray")
points(bh2$years, dynb0.bh2[[3]]$SSB/1000, col = "black", pch = 19)
mtext("Spawning Biomass (kmt)", side = 2, line = 2.5, cex = 1.5)

plot(bh2$years, dynb0.bh2[[3]]$NAA[,1]/1000, axes = FALSE, ann = FALSE, type = "n", ylim = c(0, max.xy[2]))
axis(1, labels = FALSE, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7))
abline(h = median(dynb0.bh2[[3]]$NAA[-(1:6),1])/1000, lwd = 2, lty = 2)
abline(h = median(exp(bh0$rep$log_SR_R0[-(1:6)]))/1000, lwd = 2, lty = 2, col = "gray")
points(bh2$years, dynb0.bh2[[3]]$NAA[,1]/1000, col = "black", pch = 19)
mtext(bquote("Recruits " ~ (10^3)), side = 2, line = 2.5, cex = 1.5)
mtext("Beverton-Holt", side = 4, line = 1, cex = 1.5)

max.xy = c(max(dynb0.rick1[[3]]$SSB, dynb0.rick2[[3]]$SSB), max(dynb0.rick2[[3]]$NAA[,1], dynb0.rick1[[3]]$NAA[,1]))/1000
plot(rick2$years, dynb0.rick2[[3]]$SSB, axes = FALSE, ann = FALSE, type = "n", ylim = c(0, max.xy[1]))
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7))
abline(h = median(dynb0.rick2[[3]]$SSB[-(1:6)])/1000, lwd = 2, lty = 2)
abline(h = median(exp(rick0$rep$log_SR_R0 + rick0$rep$log_SPR0)[-(1:6)])/1000, lwd = 2, lty = 2, col = "gray")
points(rick2$years, dynb0.rick2[[3]]$SSB/1000, col = "black", pch = 19)
mtext("Spawning Biomass (kmt)", side = 2, line = 2.5, cex = 1.5)

plot(rick2$years, dynb0.rick2[[3]]$NAA[,1]/1000, axes = FALSE, ann = FALSE, type = "n", ylim = c(0, max.xy[2]))
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7))
abline(h = median(dynb0.rick2[[3]]$NAA[-(1:6),1])/1000, lwd = 2, lty = 2)
abline(h = median(exp(rick0$rep$log_SR_R0[-(1:6)]))/1000, lwd = 2, lty = 2, col = "gray")
points(rick2$years, dynb0.rick2[[3]]$NAA[,1]/1000, col = "black", pch = 19)
mtext(bquote("Recruits " ~ (10^3)), side = 2, line = 2.5, cex = 1.5)
mtext("Ricker", side = 4, line = 1, cex = 1.5)
mtext("Year", side = 1, outer = TRUE, line = 2, cex = 1.5)
dev.off()

cairo_pdf("../tex/dynamic_B0_all.pdf", family = "Times", height = 10, width = 10)
par(mfrow = c(2,2), mar = c(1,5,1,1), oma = c(4,0,1,3))
max.xy = c(max(dynb0.bh1[[4]]$SSB, dynb0.bh2[[4]]$SSB), max(dynb0.bh2[[4]]$NAA[,1], dynb0.bh1[[4]]$NAA[,1]))/1000
plot(bh2$years, dynb0.bh2[[4]]$SSB, axes = FALSE, ann = FALSE, type = "n", ylim = c(0, max.xy[1]))
axis(1, labels = FALSE, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7))
abline(h = median(dynb0.bh2[[4]]$SSB[-(1:6)])/1000, lwd = 2, lty = 2)
abline(h = median(exp(bh0$rep$log_SR_R0 + bh0$rep$log_SPR0)[-(1:6)])/1000, lwd = 2, lty = 2, col = "gray")
points(bh2$years, dynb0.bh2[[4]]$SSB/1000, col = "black", pch = 19)
mtext("Spawning Biomass (kmt)", side = 2, line = 2.5, cex = 1.5)

plot(bh2$years, dynb0.bh2[[4]]$NAA[,1]/1000, axes = FALSE, ann = FALSE, type = "n", ylim = c(0, max.xy[2]))
axis(1, labels = FALSE, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7))
abline(h = median(dynb0.bh2[[4]]$NAA[-(1:6),1])/1000, lwd = 2, lty = 2)
abline(h = median(exp(bh0$rep$log_SR_R0[-(1:6)]))/1000, lwd = 2, lty = 2, col = "gray")
points(bh2$years, dynb0.bh2[[4]]$NAA[,1]/1000, col = "black", pch = 19)
mtext(bquote("Recruits " ~ (10^3)), side = 2, line = 2.5, cex = 1.5)
mtext("Beverton-Holt", side = 4, line = 1, cex = 1.5)

max.xy = c(max(dynb0.rick1[[4]]$SSB, dynb0.rick2[[4]]$SSB), max(dynb0.rick2[[4]]$NAA[,1], dynb0.rick1[[4]]$NAA[,1]))/1000
plot(rick2$years, dynb0.rick2[[4]]$SSB, axes = FALSE, ann = FALSE, type = "n", ylim = c(0, max.xy[1]))
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7))
abline(h = median(dynb0.rick2[[4]]$SSB[-(1:6)])/1000, lwd = 2, lty = 2)
abline(h = median(exp(rick0$rep$log_SR_R0 + rick0$rep$log_SPR0)[-(1:6)])/1000, lwd = 2, lty = 2, col = "gray")
points(rick2$years, dynb0.rick2[[4]]$SSB/1000, col = "black", pch = 19)
mtext("Spawning Biomass (kmt)", side = 2, line = 2.5, cex = 1.5)

plot(rick2$years, dynb0.rick2[[4]]$NAA[,1]/1000, axes = FALSE, ann = FALSE, type = "n", ylim = c(0, max.xy[2]))
axis(1, lwd = 2, cex.axis = 1.5)
axis(2, lwd = 2, cex.axis = 1.5)
box(lwd = 2)
grid(col = gray(0.7))
abline(h = median(dynb0.rick2[[4]]$NAA[-(1:6),1])/1000, lwd = 2, lty = 2)
abline(h = median(exp(rick0$rep$log_SR_R0[-(1:6)]))/1000, lwd = 2, lty = 2, col = "gray")
points(rick2$years, dynb0.rick2[[4]]$NAA[,1]/1000, col = "black", pch = 19)
mtext(bquote("Recruits " ~ (10^3)), side = 2, line = 2.5, cex = 1.5)
mtext("Ricker", side = 4, line = 1, cex = 1.5)
mtext("Year", side = 1, outer = TRUE, line = 2, cex = 1.5)
dev.off()
