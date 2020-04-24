library(TMB)
#changed WHAM, to allow alpha and beta corresponding to a particular year of steepness for estimating other things like R0, MSY
compile("wham_alt.cpp")
dyn.load(dynlib("wham_alt"))
#Southern New England - Mid-Atlantic yellowtail flounder 1973-2016
load("wham_v1.RData")
fit.tmb.fn = function(model, n.newton=3, do.sdrep = TRUE)
{
  require(TMB)
  model$opt <- nlminb(model$par, model$fn,model$gr, control = list(iter.max = 1000, eval.max = 1000))
  if(n.newton) for(i in 1:n.newton) { # Take a few extra newton steps 
    g <- as.numeric(model$gr(model$opt$par))
    h <- optimHess(model$opt$par, model$fn, model$gr)
    model$opt$par <- model$opt$par - solve(h, g)
    model$opt$objective <- model$fn(model$opt$par)
  }
  model$date = Sys.time()
  model$dir = getwd()
  model$rep <- model$report()
  model$TMB_version = packageVersion("TMB")
  model$final_gradient = model$gr()
  if(do.sdrep) 
  {
    model$sdrep <- try(sdreport(model))
    model$is_sdrep = !is.character(model$sdrep)
    if(model$is_sdrep) model$na_sdrep = any(is.na(summary(model$sdrep)[,2]))
    else model$na_sdrep = NA
  }
  if(do.sdrep & model$is_sdrep) model$parList = as.list(model$sdrep, "Est") 
  else model$parList = model$env$parList()
  return(model)
}

fit.wham.fn = function(input, n.newton = 3, do.sdrep = TRUE, do.retro = FALSE, n.peels = 7)
{
  mod <- MakeADFun(input$data,input$par,DLL="wham_alt", random = input$random, map = input$map)
  mod = fit.tmb.fn(mod, n.newton = n.newton, do.sdrep = do.sdrep)
  if(do.retro) mod$peels = retro.fn(mod, ran = unique(names(mod$env$par[mod$env$random])), n.peels= n.peels)
  mod$years = input$years
  mod$ages.lab = input$ages.lab
  mod$model_name = input$model_name
  return(mod)
}
fix.all.pars = function(input)
{
  fe = names(input$par)[!(names(input$par) %in% input$random)]
  map = lapply(fe, function(x) factor(rep(NA, length(input$par[[x]]))))
  names(map) = fe
  return(map)
}

#B-H: not using steepness parameterization
temp = input
temp$data$steepness_year = which(input$years == 1973)
temp$data$R0_is_S0 = 0
temp$data$FMSY_startval = 0.2
bh0 = fit.wham.fn(temp,do.retro = FALSE)

#Ricker: not using steepness parameterization
temp = input
temp$data$steepness_year = which(input$years == 1973)
temp$data$R0_is_S0 = 0
temp$data$FMSY_startval = 0.8
temp$data$recruit_model = 4
temp$par$mean_rec_pars[2] = -5
rick0 = fit.wham.fn(temp,do.retro = FALSE)

#B-H: Estimate R0 and steepness using phi0 from 1973 (low steepness)
temp = input
temp$data$use_steepness = 1 #Beverton-Holt S-R
temp$data$steepness_year = which(input$years == 1973)
temp$data$R0_is_S0 = 0
temp$data$FMSY_startval = 0.2
bh1 = fit.wham.fn(temp,do.retro = FALSE)

#Ricker: Estimate R0 and steepness using phi0 from 1973 (low steepness)
temp = input
temp$data$use_steepness = 1 
temp$data$steepness_year = which(input$years == 1973)
temp$data$R0_is_S0 = 0
temp$data$FMSY_startval = 0.8
temp$data$recruit_model = 4
temp$par$mean_rec_pars[2] = 10
rick1 = fit.wham.fn(temp,do.retro = FALSE)

#B-H: Estimate R0 and steepness using phi0 from 1999 (high steepness)
temp = input
temp$data$use_steepness = 1
temp$data$steepness_year = which(input$years == 1999)
temp$data$R0_is_S0 = 0
temp$data$FMSY_startval = 0.2
bh2 = fit.wham.fn(temp,do.retro = FALSE)

#Ricker: Estimate R0 and steepness using phi0 from 1999 (high steepness)
temp = input
temp$data$use_steepness = 1 
temp$data$steepness_year = which(input$years == 1999)
temp$data$R0_is_S0 = 0
temp$data$FMSY_startval = 0.8
temp$data$recruit_model = 4
temp$par$mean_rec_pars[2] = 10
rick2 = fit.wham.fn(temp,do.retro = FALSE)

#B-H: Use phi0 from 1999 and steepness estimated from 1973
temp = input
temp$data$use_steepness = 1
temp$data$steepness_year = which(input$years == 1999)
temp$data$R0_is_S0 = 0
temp$data$FMSY_startval = 0.2
temp$par$mean_rec_pars = bh1$parList$mean_rec_pars
temp$map$mean_rec_pars = factor(c(NA,1)) #fix steepness
bh3_free = fit.wham.fn(temp,do.retro = FALSE) #estimate all parameters EXCEPT 
temp$par = bh1$parList #steepness in 1973
temp$map = fix.all.pars(temp)
temp$map$log_NAA = factor(rep(NA,length(temp$par$log_NAA)))
temp$map$mean_rec_pars = factor(c(NA,1)) #fix steepness
#temp$map$log_NAA_sigma = factor(c(1,NA))
temp$random = NULL
bh3 = fit.wham.fn(temp,do.retro = FALSE)
#temp = fit.wham.fn(temp,do.retro = FALSE)
0.2 + 0.8/(1 + exp(-bh3$rep$logit_SR_h))
0.2 + 0.8/(1 + exp(-bh0$rep$logit_SR_h[temp$data$steepness_year]))
exp(bh3$rep$log_SR_R0)
exp(bh0$rep$log_SR_R0[temp$data$steepness_year])
exp(bh3$rep$log_FMSY[temp$data$steepness_year])
exp(bh0$rep$log_FMSY[temp$data$steepness_year])
exp(bh3$rep$log_SSB_MSY[temp$data$steepness_year])
exp(bh0$rep$log_SSB_MSY[temp$data$steepness_year])


#Ricker: Use phi0 from 1999 and steepness estimated from 1973
temp = input
temp$data$use_steepness = 1 
temp$data$steepness_year = which(input$years == 1999)
temp$data$R0_is_S0 = 0
temp$data$FMSY_startval = 0.8
temp$data$recruit_model = 4
temp$par$mean_rec_pars = rick1$parList$mean_rec_pars
temp$map$mean_rec_pars = factor(c(NA,1)) #fix steepness
rick3_free = fit.wham.fn(temp,do.retro = FALSE) #estimate all parameters EXCEPT steepness
temp$par = rick1$parList #steepness in 1973
temp$map = fix.all.pars(temp)
temp$map$log_NAA = factor(rep(NA,length(temp$par$log_NAA)))
temp$map$mean_rec_pars = factor(c(NA,1)) #fix steepness
#temp$map$log_NAA_sigma = factor(c(1,NA))
temp$random = NULL
#temp$par$mean_rec_pars[2] = 10
rick3 = fit.wham.fn(temp,do.retro = FALSE)
0.2 + exp(rick3$rep$logit_SR_h)
0.2 + exp(rick0$rep$logit_SR_h[temp$data$steepness_year])
exp(rick3$rep$log_SR_R0)
exp(rick0$rep$log_SR_R0[temp$data$steepness_year])
exp(rick3$rep$log_FMSY[temp$data$steepness_year])
exp(rick0$rep$log_FMSY[temp$data$steepness_year])
exp(rick3$rep$log_SSB_MSY[temp$data$steepness_year])
exp(rick0$rep$log_SSB_MSY[temp$data$steepness_year])


#B-H: Use phi0 from 1973 and steepness estimated from 1999
temp = input
temp$data$use_steepness = 1
temp$data$steepness_year = which(input$years == 1973)
temp$data$R0_is_S0 = 0
temp$data$FMSY_startval = 0.5
temp$par = bh2$parList #steepness in 1999
temp$map = fix.all.pars(temp)
temp$map$log_NAA = factor(rep(NA,length(temp$par$log_NAA)))
temp$map$mean_rec_pars = factor(c(NA,1)) #fix steepness
#temp$map$log_NAA_sigma = factor(c(1,NA))
temp$random = NULL
bh4 = fit.wham.fn(temp,do.retro = FALSE)
0.2 + 0.8/(1 + exp(-bh4$rep$logit_SR_h))
0.2 + 0.8/(1 + exp(-bh0$rep$logit_SR_h[temp$data$steepness_year]))
exp(bh4$rep$log_SR_R0)
exp(bh0$rep$log_SR_R0[temp$data$steepness_year])
exp(bh4$rep$log_FMSY[temp$data$steepness_year])
exp(bh0$rep$log_FMSY[temp$data$steepness_year])
exp(bh4$rep$log_SSB_MSY[temp$data$steepness_year])
exp(bh0$rep$log_SSB_MSY[temp$data$steepness_year])

#Ricker: Use phi0 from 1973 and steepness estimated from 1999
temp = input
temp$data$use_steepness = 1 
temp$data$steepness_year = which(input$years == 1973)
temp$data$R0_is_S0 = 0
temp$data$FMSY_startval = 0.8
temp$data$recruit_model = 4
temp$par = rick2$parList #steepness in 1999
temp$map = fix.all.pars(temp)
temp$map$log_NAA = factor(rep(NA,length(temp$par$log_NAA)))
temp$map$mean_rec_pars = factor(c(NA,1)) #fix steepness
temp$random = NULL
rick4 = fit.wham.fn(temp,do.retro = FALSE)
0.2 + exp(rick4$rep$logit_SR_h)
0.2 + exp(rick0$rep$logit_SR_h[temp$data$steepness_year])
exp(rick4$rep$log_SR_R0)
exp(rick0$rep$log_SR_R0[temp$data$steepness_year])
exp(rick4$rep$log_FMSY[temp$data$steepness_year])
exp(rick0$rep$log_FMSY[temp$data$steepness_year])
exp(rick4$rep$log_SSB_MSY[temp$data$steepness_year])
exp(rick0$rep$log_SSB_MSY[temp$data$steepness_year])
exp(rick4$rep$log_MSY[temp$data$steepness_year])
exp(rick0$rep$log_MSY[temp$data$steepness_year])


save(bh0,bh1,bh2,bh3,bh4, file = "bh_models.RData")
save(rick0,rick1,rick2,rick3,rick4, file = "rick_models.RData")

plot(bh0$years, exp(bh0$rep$log_FMSY), type = 'l', ylim = c(0, exp(max(bh0$rep$log_FMSY, bh3$rep$log_FMSY, bh4$rep$log_FMSY))))
lines(bh0$years, exp(bh3$rep$log_FMSY), lty = 2, col = "red")
lines(bh0$years, exp(bh4$rep$log_FMSY), lty = 2, col = "blue")
abline(v = 1973, col = "blue")
abline(v = 1999, col = "red")


#dynamic B0?
load("bh_models.RData")
temp = input
temp$data$use_steepness = 1
#R0 is estimated for phi0 in 1999 (highest value)
temp$data$steepness_year = which(input$years == 1999)
temp$data$R0_is_S0 = 0
temp$data$FMSY_startval = 0.2
#bh2 = fit.wham.fn(temp,do.retro = FALSE)
0.2 + 0.8/(1 + exp(-bh2$rep$logit_SR_h))
0.2 + 0.8/(1 + exp(-bh0$rep$logit_SR_h[temp$data$steepness_year]))
exp(bh2$rep$log_SR_R0)
exp(bh0$rep$log_SR_R0[temp$data$steepness_year])

#first just see what happens when there are no rec devs
set.seed(034582)
#temp$data$simulate_state = 0
temp$par = bh2$parList
temp$par$log_F1[] = -100
temp$par$F_devs[] = 0
temp$par$log_NAA_sigma[] = -100 #var(devs) = 0

x = MakeADFun(temp$data, temp$par, map = temp$map, random = temp$random, DLL = "wham_alt")
y = x$simulate(complete = TRUE)
plot(bh2$years, y$SSB)
plot(bh2$years, y$R)
plot(y$SSB[-44],y$R[-1]) #all on the curve
#if we don't make post-recruit productivity constant, things just never stabilize

for(i in 1:3) temp$data$waa[i,,] = rep(temp$data$waa[i,temp$data$steepness_year,], each = length(bh2$years))
set.seed(034582)
x = MakeADFun(temp$data, temp$par, map = temp$map, random = temp$random, DLL = "wham_alt")
y = x$simulate(complete = TRUE)
plot(bh2$years, y$SSB)
plot(bh2$years, y$NAA[,1])
plot(y$SSB[-44],y$NAA[-1,1], col = palette.fn(44), pch = 19) #all on the curve

#now add variability in R
set.seed(034582)
temp$par = bh2$parList
rdevs = log(bh2$rep$NAA[,1]) - log(bh2$rep$pred_NAA[,1])
temp$par$log_F1[] = -100
temp$par$F_devs[] = 0
temp$par$log_NAA_sigma[2] = -100 #var(devs) = 0
temp$data$bias_correct_pe = 1
x = MakeADFun(temp$data, temp$par, map = temp$map, random = temp$random, DLL = "wham_alt")
y = x$simulate(complete = TRUE)
plot(bh2$years, y$SSB)
plot(bh2$years, y$NAA[,1])
plot(y$SSB[-44],y$NAA[-1,1], col = palette.fn(44), pch = 19) #all on the curve

dynB0.fn = function(alpha, beta, N1, mat, M, wt, fracyrssb, recdevs)
{
  ny = length(recdevs)
  na = length(N1)
  SSB = sum(N1*mat[1,]*wt[1,]*exp(-M[1,]*fracyrssb[1]))
  NAA = matrix(NA, ny, na)
  NAA[1,] = N1
  for(i in 1:(ny-1)) 
  {
    NAA[i+1,1] = exp(recdevs[i+1]) * alpha * SSB[i]/(1 + beta*SSB[i])
    NAA[i+1,-1] = NAA[i,1:(na-1)] * exp(-M[i, 1:(na-1)])
    NAA[i+1,na] = NAA[i+1,na] + NAA[i,na]*exp(-M[i,na])
    SSB[i+1] = sum(NAA[i+1,]*mat[i+1,]*wt[i+1,]*exp(-M[i+1,]*fracyrssb[i+1]))
  }
  return(list(NAA=NAA, SSB = SSB))
}
z = dynB0.fn(alpha = exp(bh2$rep$log_SR_a[bh2$env$data$steepness_year]),
  beta = exp(bh2$rep$log_SR_b[bh2$env$data$steepness_year]),
  N1 = bh2$rep$NAA[1,], 
  mat = bh2$env$data$mature,
  M = bh2$rep$MAA,
  wt = bh2$env$data$waa[bh2$env$data$waa_pointer_ssb,,],
  fracyrssb = bh2$env$data$fracyr_SSB,
  #recdevs = rep(0, length(bh2$years)))$SSB
  recdevs = log(bh2$rep$NAA[,1]) - log(bh2$rep$pred_NAA[,1]))
plot(z$SSB)  

#  
temp = input
temp$data$use_steepness = 1
#R0 is estimated for phi0 in 1999 (highest value)
temp$data$steepness_year = which(input$years == 1999)
temp$data$R0_is_S0 = 0
temp$data$FMSY_startval = 0.2
#temp$data$simulate_state = 0
temp$par = bh2$parList
temp$par$log_F1[] = -100
temp$par$F_devs[] = 0
temp$par$log_NAA_sigma[] = -100 #var(devs) = 0
for(i in 1:3) temp$data$waa[i,,] = rep(temp$data$waa[i,temp$data$steepness_year,], each = length(bh2$years))
set.seed(034582)
  
  
  
