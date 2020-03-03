if(!require("here")) {
  install.packages("here")
  library("here")
}

if(!require("AICcmodavg")) {
  install.packages("AICcmodavg")
  library("AICcmodavg")
}

if(!require("readr")) {
  install.packages("readr")
  library("readr")
}

if(!require("plyr")) {
  install.packages("plyr")
  library("plyr")
}


datadir <- here("data")
figsdir <- here("figures")  

model.sel.table <- function(model_objects,lifestage){
  
  if(lifestage == "spawner_smolt"){
    mod.sel <- data.frame(model = c(rep("Ricker",4),rep("Bev Holt",4),rep("DI",4)),
                          cov = rep(c("winter flow + pink esc","winter flow","pink esc","NA"),3),
                          k = sapply(model_objects,function(x) length(coef(x))),
                          AICc = sapply(model_objects,function(x) AICc(x)))
  }  
  
  if(lifestage == "smolt_nearshore"|lifestage == "nearshore_adult"){
    mod.sel <- data.frame(model = c(rep("Ricker",5),rep("Bev Holt",5),rep("DI",5)),
                          cov = rep(c("spring flow + pink esc + NPGO","spring flow","pink esc","NPGO","NA"),3),
                          k = sapply(model_objects,function(x) length(coef(x))),
                          AICc = sapply(model_objects,function(x) AICc(x)))
  }  
  
  
  
  # calculate delta-AICc
  mod.sel$delta.AICc <- mod.sel$AICc - min(mod.sel$AICc)
  
  # calculate Akaike weights
  wt <- exp(-0.5*mod.sel$delta.AICc)
  mod.sel$Ak.wt <- wt/sum(wt)
  
  # sort results
  AICc.tbl <- mod.sel[order(mod.sel$AICc),]
  
  return(AICc.tbl)
  
}

##data
dd_data <- read_csv(file.path(datadir,"skagit_bay_cpue.csv"))
dd_data <- data.frame(dd_data)


model <- lm(nearshore_chum_fry_growth~chum_smolt_abundance,data = dd_data)
summary(model)

logit_growth <- 1/(1 + exp(-dd_data$nearshore_chum_fry_growth))
model <- glm(logit_growth~dd_data$chum_smolt_abundance, family = "binomial")
summary(model)

model <- lm(logit_growth~dd_data$chum_smolt_abundance*dd_dat)
summary(model)
plot(nearshore_chum_fry_growth~chum_smolt_abundance,data = dd_data)
lines(exp(pred)~seq(0,3e7,1e5))

pred <- predict(model,newdata = list(chum_smolt_abundance = seq(0,3e7,1e5)))

##
esc_full_ts <- read_csv(file.path(datadir,"skagit_chum_esc.csv"))
catch_full_ts <- read_csv(file.path(datadir,"skagit_chum_catch.csv"))


dat_sar <- read.csv(file.path(datadir,"ps_hatchery_chum_return_rates.csv"))

covars <- data.frame(scale(dd_data[,31:34]))

##spawner_smolt
 
y <- log(dd_data$chum_smolt_abundance)
x <- dd_data$escapement
pink_esc <- covars$northSound_Pink_esc  
flow_winter <- covars$winter_flow
data <- data.frame(y,x,pink_esc,flow_winter)

##Ricker
rk_spawner_smolt <- nls(y~ln_Rkr_a + log(x) - beta*x + gamma1*flow_winter + gamma2*pink_esc,data = data, start = list(ln_Rkr_a = 1.6, beta = 1e-05,gamma1 = 0,gamma2 = 0))
rk_spawner_smolt_1 <- nls(y~ln_Rkr_a + log(x) - beta*x + gamma1*flow_winter,data = data, start = list(ln_Rkr_a = 1.6, beta = 1e-05,gamma1 = 0))
rk_spawner_smolt_2 <- nls(y~ln_Rkr_a + log(x) - beta*x + gamma2*pink_esc,data = data, start = list(ln_Rkr_a = 1.6, beta = 1e-05,gamma2 = 0))
rk_spawner_smolt_3 <- nls(y~ln_Rkr_a + log(x) - beta*x,data = data, start = list(ln_Rkr_a = 1.6, beta = 1e-05))


##BH
bh_spawner_smolt <- nls(y~ln_BH_a + log(x) - log(1 + beta*x) + gamma1*flow_winter + gamma2*pink_esc,data = data, start = list(ln_BH_a = 1.6, beta = 1e-05,gamma1 = 0,gamma2 = 0))
bh_spawner_smolt_1 <- nls(y~ln_BH_a + log(x) - log(1 + beta*x) + gamma1*flow_winter,data = data, start = list(ln_BH_a = 1.6, beta = 1e-05,gamma1 = 0))
bh_spawner_smolt_2 <- nls(y~ln_BH_a + log(x) - log(1 + beta*x) + gamma2*pink_esc,data = data, start = list(ln_BH_a = 1.6, beta = 1e-05,gamma2 = 0))
bh_spawner_smolt_3 <- nls(y~ln_BH_a + log(x) - log(1 + beta*x),data = data, start = list(ln_BH_a = 1.6, beta = 1e-05))

##linear
di_spawner_smolt <- lm(y~log(x) + flow_winter + pink_esc,data = data)
di_spawner_smolt_1 <- lm(y~log(x) + flow_winter,data = data)
di_spawner_smolt_2 <- lm(y~log(x) + pink_esc,data = data)
di_spawner_smolt_3 <- lm(y~log(x),data = data)

## model selection
mod_fits_spawner_smolt <- vector("list", 12)
mod_fits_spawner_smolt[[1]] <- rk_spawner_smolt;mod_fits_spawner_smolt[[2]] <- rk_spawner_smolt_1;mod_fits_spawner_smolt[[3]] <- rk_spawner_smolt_2;mod_fits_spawner_smolt[[4]] <- rk_spawner_smolt_3
mod_fits_spawner_smolt[[5]] <- bh_spawner_smolt;mod_fits_spawner_smolt[[6]] <- bh_spawner_smolt_1;mod_fits_spawner_smolt[[7]] <- bh_spawner_smolt_2;mod_fits_spawner_smolt[[8]] <- bh_spawner_smolt_3
mod_fits_spawner_smolt[[9]] <- di_spawner_smolt;mod_fits_spawner_smolt[[10]] <- di_spawner_smolt_1;mod_fits_spawner_smolt[[11]] <- di_spawner_smolt_2;mod_fits_spawner_smolt[[12]] <- di_spawner_smolt_3


model.sel.table_spawner_smolt <- model.sel.table(model_objects = mod_fits_spawner_smolt,lifestage = "spawner_smolt")
plot(residuals(rk_spawner_smolt),type = "b")




## smolt_nearshore

y <- log(dd_data$CH.0.)
x <- (dd_data$chum_smolt_abundance)
pink_esc <- covars$northSound_Pink_esc  
flow_spring <- covars$spring_flow
NPGO <- covars$NPGO


data <- data.frame(y,x,pink_esc,flow_spring,NPGO)

plot(dd_data$CH.0.~dd_data$NPGO,pch = 16,col = "black")
model <- lm(y~NPGO)
summary(model)

##Ricker
rk_smolt_nearshore <- nls(y~ln_Rkr_a + log(x) - beta*x + gamma1*flow_spring + gamma2*pink_esc + gamma3*NPGO,data = data, start = list(ln_Rkr_a = 1.6, beta = 1e-05,gamma1 = 0,gamma2 = 0,gamma3 = 0))
rk_smolt_nearshore_1 <- nls(y~ln_Rkr_a + log(x) - beta*x + gamma1*flow_spring,data = data, start = list(ln_Rkr_a = 1.6, beta = 1e-05,gamma1 = 0))
rk_smolt_nearshore_2 <- nls(y~ln_Rkr_a + log(x) - beta*x + gamma2*pink_esc,data = data, start = list(ln_Rkr_a = 1.6, beta = 1e-05,gamma2 = 0))
rk_smolt_nearshore_3 <- nls(y~ln_Rkr_a + log(x) - beta*x + gamma3*NPGO,data = data, start = list(ln_Rkr_a = 1.6, beta = 1e-05,gamma3 = 0))
rk_smolt_nearshore_4 <- nls(y~ln_Rkr_a + log(x) - beta*x,data = data, start = list(ln_Rkr_a = 1.6, beta = 1e-05))



##BH
bh_smolt_nearshore <- nls(y~ln_BH_a + log(x) - log(1 + beta*x) + gamma1*flow_spring + gamma2*pink_esc + gamma3*NPGO,data = data, start = list(ln_BH_a = 0.5, beta = 1e-08,gamma1 = 0,gamma2 = 0,gamma3 = 0.01))
bh_smolt_nearshore_1 <- nls(y~ln_BH_a + log(x) - log(1 + beta*x) + gamma1*flow_spring,data = data, start = list(ln_BH_a = 0.5, beta = 1e-08,gamma1 = 0))
bh_smolt_nearshore_2 <- nls(y~ln_BH_a + log(x) - log(1 + beta*x) + gamma2*pink_esc,data = data, start = list(ln_BH_a = 0.5, beta = 1e-08,gamma2 = 0))
bh_smolt_nearshore_3 <- nls(y~ln_BH_a + log(x) - log(1 + beta*x) + gamma3*NPGO,data = data, start = list(ln_BH_a = 0.5, beta = 1e-08,gamma3 = 0))
bh_smolt_nearshore_4 <- nls(y~ln_BH_a + log(x) - log(1 + beta*x),data = data, start = list(ln_BH_a = 0.5, beta = 1e-08))


##linear
di_smolt_nearshore <- lm(y~log(x) + flow_spring + pink_esc + NPGO,data = data)
di_smolt_nearshore_1 <- lm(y~log(x) + flow_spring,data = data)
di_smolt_nearshore_2 <- lm(y~log(x) + pink_esc,data = data)
di_smolt_nearshore_3 <- lm(y~log(x) + NPGO,data = data)
di_smolt_nearshore_4 <- lm(y~log(x),data = data)


mod_fits_smolt_nearshore <- vector("list", 15)
mod_fits_smolt_nearshore[[1]] <- rk_smolt_nearshore;mod_fits_smolt_nearshore[[2]] <- rk_smolt_nearshore_1;mod_fits_smolt_nearshore[[3]] <- rk_smolt_nearshore_2;mod_fits_smolt_nearshore[[4]] <- rk_smolt_nearshore_3;mod_fits_smolt_nearshore[[5]] <- rk_smolt_nearshore_4
mod_fits_smolt_nearshore[[6]] <- bh_smolt_nearshore;mod_fits_smolt_nearshore[[7]] <- bh_smolt_nearshore_1;mod_fits_smolt_nearshore[[8]] <- bh_smolt_nearshore_2;mod_fits_smolt_nearshore[[9]] <- bh_smolt_nearshore_3;mod_fits_smolt_nearshore[[10]] <- bh_smolt_nearshore_4
mod_fits_smolt_nearshore[[11]] <- di_smolt_nearshore;mod_fits_smolt_nearshore[[12]] <- di_smolt_nearshore_1;mod_fits_smolt_nearshore[[13]] <- di_smolt_nearshore_2;mod_fits_smolt_nearshore[[14]] <- di_smolt_nearshore_3;mod_fits_smolt_nearshore[[15]] <- di_smolt_nearshore_4

model.sel.table_smolt_nearshore <- model.sel.table(model_objects = mod_fits_smolt_nearshore,lifestage = "smolt_nearshore")


plot(residuals(di_smolt_nearshore_3),type = "b")


## nearshore_adult
y <- log(dd_data$adult_chum_recruits)
x <- dd_data$CH.0.
pink_esc <- covars$northSound_Pink_esc  
flow_spring <- covars$spring_flow
NPGO <- covars$NPGO

data <- data.frame(y,x,flow_spring,NPGO,pink_esc)

##Ricker
rk_nearshore_adult <- nls(y~ln_Rkr_a + log(x) - beta*x + gamma1*flow_spring + gamma2*pink_esc + gamma3*NPGO,data = data, start = list(ln_Rkr_a = 1.6, beta = 1e-05,gamma1 = 0,gamma2 = 0,gamma3 = 0))
rk_nearshore_adult_1 <- nls(y~ln_Rkr_a + log(x) - beta*x + gamma1*flow_spring,data = data, start = list(ln_Rkr_a = 1.6, beta = 1e-05,gamma1 = 0))
rk_nearshore_adult_2 <- nls(y~ln_Rkr_a + log(x) - beta*x + gamma2*pink_esc,data = data, start = list(ln_Rkr_a = 1.6, beta = 1e-05,gamma2 = 0))
rk_nearshore_adult_3 <- nls(y~ln_Rkr_a + log(x) - beta*x + gamma3*NPGO,data = data, start = list(ln_Rkr_a = 1.6, beta = 1e-05,gamma3 = 0))
rk_nearshore_adult_4 <- nls(y~ln_Rkr_a + log(x) - beta*x,data = data, start = list(ln_Rkr_a = 1.6, beta = 1e-05))


##BH
bh_nearshore_adult <- nls(y~ln_BH_a + log(x) - log(1 + beta*x) + gamma1*flow_spring + gamma2*pink_esc + gamma3*NPGO,data = data, start = list(ln_BH_a = 0.5, beta = 1e-08,gamma1 = 0,gamma2 = 0,gamma3 = 0.01))
bh_nearshore_adult_1 <- nls(y~ln_BH_a + log(x) - log(1 + beta*x) + gamma1*flow_spring,data = data, start = list(ln_BH_a = 0.5, beta = 1e-08,gamma1 = 0))
bh_nearshore_adult_2 <- nls(y~ln_BH_a + log(x) - log(1 + beta*x) + gamma2*pink_esc,data = data, start = list(ln_BH_a = 0.5, beta = 1e-08,gamma2 = 0))
bh_nearshore_adult_3 <- nls(y~ln_BH_a + log(x) - log(1 + beta*x) + gamma3*NPGO,data = data, start = list(ln_BH_a = 0.5, beta = 1e-08,gamma3 = 0))
bh_nearshore_adult_4 <- nls(y~ln_BH_a + log(x) - log(1 + beta*x),data = data, start = list(ln_BH_a = 0.5, beta = 1e-08))


##linear
di_nearshore_adult <- lm(y~log(x) + flow_spring  + pink_esc + NPGO,data = data)
di_nearshore_adult_1 <- lm(y~log(x) + flow_spring,data = data)
di_nearshore_adult_2 <- lm(y~log(x) + pink_esc,data = data)
di_nearshore_adult_3 <- lm(y~log(x) + NPGO,data = data)
di_nearshore_adult_4 <- lm(y~log(x),data = data)

mod_fits_nearshore_adult <- vector("list", 15)
mod_fits_nearshore_adult[[1]] <- rk_nearshore_adult;mod_fits_nearshore_adult[[2]] <- rk_nearshore_adult_1;mod_fits_nearshore_adult[[3]] <- rk_nearshore_adult_2;mod_fits_nearshore_adult[[4]] <- rk_nearshore_adult_3;mod_fits_nearshore_adult[[5]] <- rk_nearshore_adult_4
mod_fits_nearshore_adult[[6]] <- bh_nearshore_adult;mod_fits_nearshore_adult[[7]] <- bh_nearshore_adult_1;mod_fits_nearshore_adult[[8]] <- bh_nearshore_adult_2;mod_fits_nearshore_adult[[9]] <- bh_nearshore_adult_3;mod_fits_nearshore_adult[[10]] <- bh_nearshore_adult_4
mod_fits_nearshore_adult[[11]] <- di_nearshore_adult;mod_fits_nearshore_adult[[12]] <- di_nearshore_adult_1;mod_fits_nearshore_adult[[13]] <- di_nearshore_adult_2;mod_fits_nearshore_adult[[14]] <- di_nearshore_adult_3;mod_fits_nearshore_adult[[15]] <- di_nearshore_adult_4


model.sel.table_nearshore_adult <- model.sel.table(model_objects = mod_fits_nearshore_adult,lifestage = "nearshore_adult")


data.frame(residuals(di_nearshore_adult_2))

## output candidate models to table
outFile="Density Dependent Model Results.csv"
write.table(model.sel.table_spawner_smolt, outFile, quote=FALSE, sep=",", append=TRUE, col.names=NA)
write.table(model.sel.table_smolt_nearshore, outFile, quote=FALSE, sep=",", append=TRUE, col.names=NA)
write.table(model.sel.table_nearshore_adult, outFile, quote=FALSE, sep=",", append=TRUE, col.names=NA)

##
png(file.path(figsdir, "life_stage_transitions.png"),
    height = 8, width = 8, units = "in", res = 500)
par(mfrow=c(3,3), mai=c(0.6,0.6,0.3,0.1), omi=c(0.2,0.5,0,0))

xoffSet <- 0.9
yoffSet <- 0.03
x <- seq(min(dd_data$escapement),max(dd_data$escapement),1000)
pred <- exp(predict(rk_spawner_smolt,newdata = list(x = x,flow_winter = rep(flow_winter[1],length(x)),pink_esc = rep(pink_esc[1],length(x))))) 

pred_spawner_smolt <- matrix(NA,length(x),length(flow_winter))
for(i in 1:dim(pred_spawner_smolt)[2]){
  
  pred_spawner_smolt[,i] <- exp(predict(rk_spawner_smolt,newdata = list(x = x,flow_winter = rep(flow_winter[i],length(x)),pink_esc = rep(pink_esc[i],length(x))))) 
  
}

plot(dd_data$chum_smolt_abundance~dd_data$escapement,pch = 16, xlab = "Spawners", 
     ylab = "Smolt abundance", cex.lab = 1.2, bty = "L", ylim = c(0,3.0e7),type = "n")
for(i in 1:dim(pred_spawner_smolt)[2]){
  
  lines(pred_spawner_smolt[,i]~x,col = "darkgray")
  
}

points(dd_data$chum_smolt_abundance~dd_data$escapement,pch = 16)
text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
     y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
     paste0("(",letters[1],")"),
     cex = 1.2)

plot(dd_data$chum_smolt_abundance~dd_data$year,pch = 16, xlab = "Year", 
     ylab = "Smolt abundance", cex.lab = 1.2, bty = "L", ylim = c(0,3.0e7))
lines(exp(predict(rk_spawner_smolt))~dd_data$year,lwd = 1.2,col = "blue3")
text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
     y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
     paste0("(",letters[2],")"),
     cex = 1.2)

plot(residuals(rk_spawner_smolt)~dd_data$year,type = "b",pch = 16, xlab = "Year", 
     ylab = "Residuals", cex.lab = 1.2, bty = "L")
text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
     y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
     paste0("(",letters[3],")"),
     cex = 1.2)
##
x <- seq(min(dd_data$chum_smolt_abundance),max(dd_data$chum_smolt_abundance),50000)

pred_smolt_nearshore <- matrix(NA,length(x),length(NPGO))
for(i in 1:dim(pred_spawner_smolt)[2]){
  
  pred_smolt_nearshore[,i] <- exp(predict(di_smolt_nearshore_3,newdata = list(x = x,NPGO = rep(NPGO[i],length(x))))) 
  
}

plot(dd_data$CH.0.~dd_data$chum_smolt_abundance,pch = 16, xlab = "Smolt abundance", 
     ylab = "Nearshore fry cpue", cex.lab = 1.2, bty = "L",type = "n")
for(i in 1:dim(pred_smolt_nearshore)[2]){
  
  lines(pred_smolt_nearshore[,i]~x,col = "darkgray")
  
}

points(dd_data$CH.0.~dd_data$chum_smolt_abundance,pch = 16)
text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
     y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
     paste0("(",letters[4],")"),
     cex = 1.2)

plot(dd_data$CH.0.~dd_data$year,pch = 16, xlab = "Year", 
     ylab = "Nearshore fry cpue", cex.lab = 1.2, bty = "L")
lines(exp(predict(di_smolt_nearshore_3))~dd_data$year,col = "blue3")
text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
     y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
     paste0("(",letters[5],")"),
     cex = 1.2)


plot(residuals(di_smolt_nearshore_3)~dd_data$year,type = "b",pch = 16, xlab = "Year", 
     ylab = "Residuals", cex.lab = 1.2, bty = "L")
text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
     y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
     paste0("(",letters[6],")"),
     cex = 1.2)
##
x <- seq(min(dd_data$CH.0.),max(dd_data$CH.0.),20)

pred_nearshore_adult <- matrix(NA,length(x),length(NPGO))
for(i in 1:dim(pred_nearshore_adult)[2]){
  
  pred_nearshore_adult[,i] <- exp(predict(di_nearshore_adult_2,newdata = list(x = x,pink_esc = rep(pink_esc[i],length(x))))) 
  
}

plot(dd_data$adult_chum_recruits~dd_data$CH.0.,pch = 16, xlab = "Nearshore fry cpue", 
     ylab = "Adult recruits", cex.lab = 1.2, bty = "L",type = "n")

for(i in 1:dim(pred_smolt_nearshore)[2]){
  
  lines(pred_nearshore_adult[,i]~x,col = "darkgrey")
  
}

points(dd_data$adult_chum_recruits~dd_data$CH.0.,pch = 16)
text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
     y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
     paste0("(",letters[7],")"),
     cex = 1.2)

plot(dd_data$adult_chum_recruits[1:20]~dd_data$year[1:20],pch = 16, xlab = "Year", 
     ylab = "Adult recruits", cex.lab = 1.2, bty = "L")
lines(exp(predict(di_nearshore_adult_2))~dd_data$year[1:20],col = "blue3")
text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
     y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
     paste0("(",letters[8],")"),
     cex = 1.2)

plot(residuals(di_nearshore_adult_2)~dd_data$year[1:20],type = "b", pch = 16,xlab = "Year", 
     ylab = "Residuals", cex.lab = 1.2, bty = "L")
text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
     y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
     paste0("(",letters[9],")"),
     cex = 1.2)

dev.off()

##smolt-adult
y <- log(dd_data$adult_chum_recruits)
x <- dd_data$chum_smolt_abundance
pink_esc <- covars$northSound_Pink_esc  
flow_spring <- covars$spring_flow
NPGO <- covars$NPGO
data <- data.frame(y,x)

di_smolt_adult <- lm(y~log(x) + pink_esc,data = data)
summary(di_smolt_adult)

plot(residuals(di_smolt_adult),type = "b")

data <- data.frame(y,x)

  png(file.path(figsdir, "life_stage_transitions.png"),
      height = 8, width = 7, units = "in", res = 500)
  par(mfrow=c(3,2), mai=c(0.6,0.6,0.3,0.1), omi=c(0.2,0.5,0,0))
  
  xoffSet <- 0.9
  yoffSet <- 0.03
  #smolts per spawner
  plot(dd_data$chum_smolt_abundance~dd_data$escapement,pch = 16, xlab = "spawners", 
       ylab = "smolt abundance", cex.lab = 1.2, bty = "L", ylim = c(0,3.0e7))
  
  legend("topleft",bty = "n",legend = c("Beverton Holt", "Ricker", "Linear"),lty = c(1,1,1),col = c("blue","red","black"))
  lines(exp(predict(bh_spawner_smolt,newdata = list(x = seq(1,400000,100))))~seq(1,400000,100),col = "blue")
  lines(exp(predict(rk_spawner_smolt,newdata = list(x = seq(1,400000,100))))~seq(1,400000,100),col = "red")
  lines(exp(predict(di_spawner_smolt,newdata = list(x = seq(1,400000,100))))~seq(1,400000,100),col = "black")
  text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
       y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
       paste0("(",letters[1],")"),
       cex = 1.2)
  plot(predict(bh_spawner_smolt,newdata = list(x = seq(1,400000,100)))~seq(1,400000,100),type =)
  x <- predict(bh_spawner_smolt,newdata = list(x = seq(1,400000,100)))~seq(1,400000,100)
  plot(residuals(di_spawner_smolt)~dd_data$year,type = "b", pch = 16,xlab = "Outmigration year",
       ylab = "Model residuals",cex.lab = 1.2,bty = "L")
  
  text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
       y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
       paste0("(",letters[4],")"),
       cex = 1.2)
  
  #nearshore fry per smolt
  plot(CH.0.~chum_smolt_abundance,data = dd_data, pch = 16, xlab = "Chum smolt abundance", 
       ylab = "Nearshore chum fry cpue", cex.lab = 1.2, bty = "L")
  lines(exp(predict(bh_smolt_nearshore,newdata = list(x = seq(1,4e7,1e6))))~seq(1,4e7,1e6),col = "blue")
  lines(exp(predict(rk_smolt_nearshore,newdata = list(x = seq(1,4e7,1e6))))~seq(1,4e7,1e6), col = "red")
  lines(exp(predict(di_smolt_nearshore,newdata = list(x = seq(1,4e7,1e6))))~seq(1,4e7,1e6), col = "black")
  text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
       y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
       paste0("(",letters[2],")"),
       cex = 1.2)
  plot(residuals(di_smolt_nearshore)~dd_data$year,type = "b", pch = 16,xlab = "Outmigration year",
       ylab = "Model residuals",cex.lab = 1.2,bty = "L")
  text(x = par()$usr[1] + diff(par()$usr[1:2]) * xoffSet,
       y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
       paste0("(",letters[5],")"),
       cex = 1.2)
  
  
  plot(dd_data$adult_chum_recruits~dd_data$CH.0.,pch = 16, xlab = "Nearshore chum fry cpue", 
       ylab = "Adult recruitment", cex.lab = 1.2, bty = "L")
  lines(exp(predict(bh_nearshore_adult,newdata = list(x = seq(1,2000,100))))~seq(1,2000,100), col = "blue")
  lines(exp(predict(rk_nearshore_adult,newdata = list(x = seq(1,2000,100))))~seq(1,2000,100), col = "red")
  lines(exp(predict(di_nearshore_adult,newdata = list(x = seq(1,2000,100))))~seq(1,2000,100), col = "black")
  text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
       y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
       paste0("(",letters[3],")"),
       cex = 1.2)
  plot(residuals(di_nearshore_adult)~dd_data$year[1:(dim(dd_data)[1]-2)], type = "b", pch = 16,xlab = "Outmigration year",
       ylab = "Model residuals",cex.lab = 1.2,bty = "L")
  text(x = par()$usr[1] + diff(par()$usr[1:2]) * xoffSet,
       y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
       paste0("(",letters[6],")"),
       cex = 1.2)
  
  dev.off()
  
  ##
  png(file.path(figsdir, "life_stage_abundance.png"),
      height = 8, width =6, units = "in", res = 500)
  par(mfrow=c(3,1), mai=c(0.6,0.6,0.1,0.1), omi=c(0.2,0.5,0,0))
  
  plot(esc_full_ts$escapement~esc_full_ts$year,type = "b",bty = "L",pch = 16,xlab = "", 
       ylab = "Chum spawner abundance",cex = 1.2, cex.lab = 1.2)
  text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
       y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
       paste0("(",letters[1],")"),
       cex = 1.2)
  
  plot(dd_data$chum_smolt_abundance~(dd_data$year-1),xlim = c(min(esc_full_ts$year),max(esc_full_ts$year)),type = "b",bty = "L",pch = 16,xlab = "", 
       ylab = "Chum smolt abundance", cex = 1.2, cex.lab = 1.2)
  text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
       y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
       paste0("(",letters[2],")"),
       cex = 1.2)
  
  plot(dd_data$CH.0.~(dd_data$year-1),xlim = c(min(esc_full_ts$year),max(esc_full_ts$year)),type = "b",bty = "L",pch = 16,xlab = "Brood year", 
       ylab = "Nearshore chum fry cpue", cex = 1.2, cex.lab = 1.2)
  text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
       y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
       paste0("(",letters[3],")"),
       cex = 1.2)
  
  dev.off()
  
  png(file.path(figsdir, "terminalRunSize.png"),
      height = 4, width =6, units = "in", res = 500)
  par(mai = c(0.8,0.8,0.1,0.1), omi = c(0.5,0.2,0.1,0.2))
  trs <- esc_full_ts$escapement + catch_full_ts$catch
  plot(trs~esc_full_ts$year,type = "l",bty = "L",pch = 16,col = "black",lwd = 2,
       xlab = "Year",ylab = "Chum terminal run size",cex =1.2,cex.lab = 1.2)
  points(trs~esc_full_ts$year,pch = 16)
  dev.off()

  
  ##growth model
  y <- dd_data$nearshore_chum_fry_growth
  x <- dd_data$chum_smolt_abundance
  
  model_growth <- lm((nearshore_chum_fry_growth)~chum_smolt_abundance,data = dd_data)
  summary(model_growth)
  
  pred_linear <- 
  
  # new fit
  fitnls <- nls(y ~SSasymp(x, Asym, r0, lrc)) 
  summary(fitnls)
  
  predict <- seq(min(dd_data$chum_smolt_abundance),max(dd_data$chum_smolt_abundance),100000)
  pred <- predict(fitnls,newdata = list(x = predict))
  
 
  

  ## data
  fl.data <- "chum_lengths.csv" 
  
  
  fl <- read_csv(file.path(datadir,fl.data))
  
  
  months <- c(3,4,5,6,7,8)
  years <- seq(1997,2018,1)
  

  
  fl <- fl[which(fl$Month %in% months & fl$Year %in% years),]
  fl_average <- aggregate(length.mm~Month,FUN = "mean",data = fl)
  
  data_chum_growth <- matrix(NA,nrow = 1, ncol = length(months),dimnames = list("",months[order(months)]))
  a <- NULL
  r <- NULL
  
  for(i in 1:length(years)){
    #i<-1
    data.temp <- fl[fl$Year %in% years[i],]
    
    fl_data <- data.frame(Month = data.temp$Month,length.mm = data.temp$length.mm)
    
    #matplot(cum_prop[,12:19],type = "b",pch = 16)
    
    
    growth_model <- nls(length.mm~a*(1+r)^Month,start = list(a = 2,r = 2),data = fl_data)
    pred <- predict(growth_model,newdata = list(Month = months))
    temp_growth <- data.frame(rbind(months,pred))
    colnames(temp_growth) <- temp_growth[1,]
    temp_growth <- temp_growth[-1,]
    data_chum_growth <- rbind.fill.matrix(data_chum_growth,temp_growth)
    a[i] <- coef(growth_model)[1]
    r[i] <- coef(growth_model)[2]
    
  }
  
  
  
  data_chum_growth <- data_chum_growth[-1,]
  

  
  png(file.path(figsdir, "chumGrowth.png"),
      height = 7, width =7, units = "in", res = 500)
  
  layout(matrix(c(1,2,3,0),2,2),c(3,3),c(3,3))
  xoffSet <- 0.05
  yoffSet <- 0.03
  par(mai = c(0.8,0.8,0.1,0.1), omi = c(0,0,0,0))
  
 
  matplot(t(data_chum_growth),type = "l",col = "grey",bty = "L",lty = 1,xaxt = "n",
          xlab = "Month",ylab = "Chum fork length (mm)",cex.lab = 1.2)
  axis(1,at = c(1,2,3,4,5,6),labels = c(3,4,5,6,7,8))
  text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
       y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
       paste0("(",letters[1],")"),
       cex = 1.2)
  
  plot(r~years,type = "l",col = "black",bty = "L",lty = 1,
       xlab = "Year",ylab = "growth rate",cex.lab = 1.2, ylim = c(0.08,0.22))
  text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
       y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
       paste0("(",letters[2],")"),
       cex = 1.2)
  
  plot((nearshore_chum_fry_growth)~(chum_smolt_abundance),pch = 16,bty = "L",data= dd_data,
       xlab = "chum smolt abundance",ylab = "growth rate",cex.lab = 1.2)
  lines(pred~predict,type = "l")
  text(x = par()$usr[1] + diff(par()$usr[1:2])*xoffSet,
       y = par()$usr[4] - diff(par()$usr[3:4]) * yoffSet,
       paste0("(",letters[3],")"),
       cex = 1.2)
  
  dev.off()
  
  
  ##timing model
  timing_model <- lm(d50_Chum~CH.0.*spring_flow,data = dd_data)
  summary(timing_model)
  