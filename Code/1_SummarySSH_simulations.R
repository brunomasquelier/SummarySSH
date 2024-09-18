rm(list=ls(all=TRUE))

  # setwd("C:\\SIMULATIONS\\RAMMPS") 

  # load packages
  library(demography, warn.conflicts = F)
  library(RColorBrewer, warn.conflicts = F)
  library(foreign, warn.conflicts = F)
  library(survival, warn.conflicts = F)
  library(pspline, warn.conflicts = F)
  library(gdata)
  library(demogsurv)
  library(xtable)
  library(stargazer)
  library(demogR)

  # Set of colors
  palette_G = rep(brewer.pal(9,"Greens")[3:9], 10)
	palette_B = rep(brewer.pal(9,"Blues")[3:9], 10)
	palette_R = rep(brewer.pal(9,"Reds")[3:9], 10)
	palette_Gr = rep(brewer.pal(9,"Greys")[3:9], 10)
	palette_Pp = rep(brewer.pal(9,"Purples")[3:9], 10)
	palette_O = rep(brewer.pal(9,"Oranges")[3:9], 10)
	
	# Set of symbols 
	pch_type = rep(c(15:20,21:25), 3)
	
	# Useful functions
	myseq <- function(x){unlist(lapply(x,function(i){1:i}))}
	zna<-function(x){return(ifelse(x==0,NA,x))} 
	"%w/o%" <- function(x,y) x[!x %in% y] #--  x without y
	
	# Standard age groups
	age_ref = c(0,1:4, seq(5,120, by = 1))*12
	age_ref5 = c(0,1, seq(5,120, by = 5))*12
	
	# set seed
  set.seed(3)

  # Other parameters  
  max_age = 100
  
  # Standards, pre-formated set of rates, etc.
  death <- read.delim2("death_100", as.is = TRUE )
  logitsBrass <- read.delim2("LogitsBrass.txt", as.is = TRUE)
  Boothstandard <- read.xls("Booth standard.xls")
  load("taux_m.rda")
  
  coeff_TZ2001 = data.frame(Age_group = c('20-24', '25-29', '30-34', '35-39', '40-44', '45-49'),
                            n = seq(25, 50, 5),
                            a = c(-0.0003, -0.1546, -0.1645,-0.1388, -0.1140,-0.1018),
                            b = c(1.0011,1.1560,1.1660,1.1406,1.1168, 1.1066))
  location_TZ2001 = data.frame(Age_group = c('20-24', '25-29', '30-34', '35-39', '40-44', '45-49'),
                               n = seq(25, 50, 5),
                               c = c(3.23,5.46,7.52,9.38,11.00, 12.32),
                               d = c(1.12,1.95,2.78,3.62,4.45,5.28))
  
# Re-create the simulations for each value of growth rate?
redosimgrowth = TRUE
if(redosimgrowth){
for(growthrate in c(1,3)){
  r = growthrate/100; print(paste("**** SIM", growthrate))
  if(!dir.exists(paste0(getwd(), "\\TZA0", growthrate))){
  dir.create(paste0(getwd(), "\\TZA0", growthrate))}
  if(!dir.exists(paste0(getwd(), "\\Results_TZA0", growthrate))){
    dir.create(paste0(getwd(), "\\Results_TZA0", growthrate))}
  dirrates = paste0("TZA0", growthrate)
  
  # Mortality parameters (Brass logits, General Standard)
  alpha_mort = c(0.2, -0.2, -0.6, -1.0)
  beta_mort =  c(0.7, 1.1)
  n_mort = length(alpha_mort)*length(beta_mort)
  mort_param <- expand.grid(alpha_mort=alpha_mort, beta_mort=beta_mort)
  mort_param$n_mort = 1:n_mort
  
  # Fertility parameters (Brass relational model, Booth Standard)
  alpha_fert =  c(-0.5, -0.2,0.1,0.4)
  if(growthrate == 1){beta_fert =  c(1.15, 1.4, 1.8)}
  if(growthrate == 3){beta_fert =  c(1.0, 1.15, 1.4)}
  n_fert = length(alpha_fert)*length(beta_fert)
  
  fert_param <- expand.grid(alpha_fert=alpha_fert, beta_fert=beta_fert)
  fert_param$n_fert = 1:n_fert
  
  all_param <- expand.grid(n_fert=1:n_fert, n_mort=1:n_mort)
  all_param = merge(all_param, fert_param, by = "n_fert")
  all_param = merge(all_param, mort_param, by = "n_mort")
  n_sim = n_mort * n_fert
  
  # ___________________________________________________________________________________________________________
  # 1. Life tables
  F_life_tables <- array(0, dim = c(101,1 + 15,n_fert, n_mort)); dimnames(F_life_tables) = list(c(0:max_age),c("npx", "nqx", "nax", "nx", "tx", "lx", "dx", "Lx", "Tx", "ex", "expaLx", "1ca", "age", "age5", "fx", "nb_naiss"), as.character(1:n_fert), as.character(1:n_mort)) 
  M_life_tables <- array(0, dim = c(101,1 + 15,n_fert, n_mort)); dimnames(M_life_tables) = list(c(0:max_age),c("npx", "nqx", "nax", "nx", "tx", "lx", "dx", "Lx", "Tx", "ex", "expaLx", "1ca", "age", "age5", "fx", "nb_naiss"), as.character(1:n_fert), as.character(1:n_mort)) 
  
  # Exact ages
  F_life_tables[,13,,] = c(0:(dim(F_life_tables)[1] - 1)) 
  M_life_tables[,13,,] = c(0:(dim(M_life_tables)[1] - 1)) 
  
  # Survivors
  F_life_tables[1,6,,] = 10000; M_life_tables[1,6,,] = 10000 
  for(i in 2:(nrow(F_life_tables))) {
    for(j in 1:dim(F_life_tables)[4]) {
      F_life_tables[i,6,,j] = (1/(1+exp(2*mort_param[j,1]+2*mort_param[j,2]*as.numeric(logitsBrass[,]))))[i-1]*10000
      M_life_tables[i,6,,j] = (1/(1+exp(2*mort_param[j,1]+2*mort_param[j,2]*as.numeric(logitsBrass[,]))))[i-1]*10000
    }}
  
  # Deaths
  for(i in 1:(nrow(F_life_tables)-1)) {
    F_life_tables[i,7,,] = F_life_tables[i,6,,] - F_life_tables[i+1,6,,]
    M_life_tables[i,7,,] = M_life_tables[i,6,,] - M_life_tables[i+1,6,,]
  }
  F_life_tables[nrow(F_life_tables),7,,] = F_life_tables[nrow(F_life_tables),6,,]
  M_life_tables[nrow(M_life_tables),7,,] = M_life_tables[nrow(M_life_tables),6,,]
  
  # Survival probabilities
  F_life_tables[,1,,] = (F_life_tables[,6,,] - F_life_tables[,7,,]) / F_life_tables[,6,,]
  M_life_tables[,1,,] = (M_life_tables[,6,,] - M_life_tables[,7,,]) / M_life_tables[,6,,]
  
  # Smoothing at old ages
  for(i in 1:dim(F_life_tables)[3]) {
    for(j in 1:dim(F_life_tables)[4]) {
      F_life_tables[60:101,1,i,j] = spline(c(60:90, 101), c(F_life_tables[60:90,1,i,j], 0), xout = 60: 101)$y
      M_life_tables[60:101,1,i,j] = spline(c(60:90, 101), c(M_life_tables[60:90,1,i,j], 0), xout = 60: 101)$y
    }}
  
  # nqx
  F_life_tables[,2,,] = 1-F_life_tables[,1,,]
  M_life_tables[,2,,] = 1-M_life_tables[,1,,]
  
  # revise to account for smoothing
  M_life_tables[,6,i,j] = M_life_tables[1,6,i,j] * cumprod(c(1,1-M_life_tables[,2,i,j]))[1:101]
  F_life_tables[,6,i,j] = F_life_tables[1,6,i,j] * cumprod(c(1,1-F_life_tables[,2,i,j]))[1:101]
  
  for(i in 1:(nrow(F_life_tables)-1)) {
    F_life_tables[i,7,,] = F_life_tables[i,6,,] - F_life_tables[i+1,6,,]
    M_life_tables[i,7,,] = M_life_tables[i,6,,] - M_life_tables[i+1,6,,]
  }
  
  F_life_tables[nrow(F_life_tables),7,,] = F_life_tables[nrow(F_life_tables),6,,]
  M_life_tables[nrow(M_life_tables),7,,] = M_life_tables[nrow(M_life_tables),6,,]
  
  # Age groups
  F_life_tables[,14,,] = c(0, rep(1,4), rep(seq(5, 95, by = 5), each = 5), 100)
  M_life_tables[,14,,] = c(0, rep(1,4), rep(seq(5, 95, by = 5), each = 5), 100)
  
  # nAx
  for(i in 2:(nrow(F_life_tables)-1)) {
    F_life_tables[i,3,,] = ((-1/24)*F_life_tables[i-1,7,,] + 0.5*F_life_tables[i,7,,] + (1/24)*F_life_tables[i+1,7,,])/F_life_tables[i,7,,]
    M_life_tables[i,3,,] = ((-1/24)*M_life_tables[i-1,7,,] + 0.5*M_life_tables[i,7,,] + (1/24)*M_life_tables[i+1,7,,])/M_life_tables[i,7,,]
  }
  
  # length of interval
  F_life_tables[,4,,] = 1; M_life_tables[,4,,] = 1
  
  # ntx
  F_life_tables[,5,,] = F_life_tables[,2,,] / (F_life_tables[,4,,]-F_life_tables[,2,,]*(F_life_tables[,4,,]-F_life_tables[,3,,]))
  M_life_tables[,5,,] = M_life_tables[,2,,] / (M_life_tables[,4,,]-M_life_tables[,2,,]*(M_life_tables[,4,,]-M_life_tables[,3,,]))
  
  F_life_tables[,5,,][F_life_tables[,5,,] > 1]  = 1; F_life_tables[101,5,,] = 1 
  M_life_tables[,5,,][M_life_tables[,5,,] > 1]  = 1; M_life_tables[101,5,,] = 1
  
  
  #nLx
  for(i in 1:(nrow(F_life_tables)-1)) {
    F_life_tables[i,8,,] = F_life_tables[i,4,,]*F_life_tables[i+1,6,,] + F_life_tables[i,3,,]*F_life_tables[i,7,,]
    M_life_tables[i,8,,] = M_life_tables[i,4,,]*M_life_tables[i+1,6,,] + M_life_tables[i,3,,]*M_life_tables[i,7,,]
  }
  
  F_life_tables[nrow(F_life_tables),8,,] = F_life_tables[nrow(F_life_tables),3,,]*F_life_tables[nrow(F_life_tables),7,,]
  M_life_tables[nrow(M_life_tables),8,,] = M_life_tables[nrow(M_life_tables),3,,]*M_life_tables[nrow(M_life_tables),7,,]
  
  # Tx
  for(i in 1:dim(F_life_tables)[3]) {
    for(j in 1:dim(F_life_tables)[4]) {
      F_life_tables[1,9,i,j] = sum(F_life_tables[,8,i,j]) 
      M_life_tables[1,9,i,j] = sum(M_life_tables[,8,i,j]) 
    }}
  
  for(i in 1:dim(F_life_tables)[3]) {
    for(j in 1:dim(F_life_tables)[4]) {
      F_life_tables[2,9,i,j] = sum(F_life_tables[,8,i,j]) - F_life_tables[1,8,i,j]
      M_life_tables[2,9,i,j] = sum(M_life_tables[,8,i,j]) - M_life_tables[1,8,i,j]
    }}
  
  for(i in 3:(nrow(F_life_tables))) {
    for(j in 1:dim(F_life_tables)[4]) {
      F_life_tables[i,9,,j] = apply(F_life_tables[,8,,j], c(2), sum) - apply(F_life_tables[1:(i-1),8,,j], c(2), sum)			
      M_life_tables[i,9,,j] = apply(M_life_tables[,8,,j], c(2), sum) - apply(M_life_tables[1:(i-1),8,,j], c(2), sum)
    }}
  
  # Ex
  F_life_tables[,10,,] = F_life_tables[,9,,]/F_life_tables[,6,,]
  M_life_tables[,10,,] = M_life_tables[,9,,]/M_life_tables[,6,,]
  
  #____________________________________________________________________________________________________________________________________________________________________________
  # 2. ASFR 
  # Estimation of F/TF = exp(-exp(-Y)) with Y = alpha + beta*Ys
  for(i in 1:dim(F_life_tables)[3]) {
    F_life_tables[12:50,15,i,]  =  exp(-exp(-(fert_param[i,1] + fert_param[i,2] * Boothstandard[,2])))
  }
  
  F_life_tables[12:49,15,,]  =  F_life_tables[13:50,15,,1] - F_life_tables[12:49,15,,1]
  F_life_tables[50,15,,] = 1 - F_life_tables[50,15,,]
  
  # Mean age at childbearing
  for(i in 1:n_fert) {		
    fert_param$MAC[i] = weighted.mean(F_life_tables[,13,i,1]+0.5, F_life_tables[,15,i,1])
  }
  
  #____________________________________________________________________________________________________________________________________________________________________________
  # 3. Stable-equivalent age structure
  # Crude birth rate
  both = rep(NA, length(1:dim(F_life_tables)[4]))
  for(i in 1:dim(F_life_tables)[4]) {
    both[i] = 1/sum(exp(-r * (F_life_tables[,13,1,i] + 0.5)) * (F_life_tables[,8,1,i] + F_life_tables[,8,1,i])/(F_life_tables[1,6,1,i] + F_life_tables[1,6,1,i]))
  }
  for(k in 1:dim(F_life_tables)[4]) {
    F_life_tables[,12,,k] = (both[k]/1) * exp(-r * (F_life_tables[,13,,k] + 0.5)) * F_life_tables[,8,,k]/F_life_tables[1,6,,k]
    M_life_tables[,12,,k] = (both[k]/1) * exp(-r * (M_life_tables[,13,,k] + 0.5)) * M_life_tables[,8,,k]/M_life_tables[1,6,,k]
  }
  for(i in 1:dim(F_life_tables)[4]) {
    F_life_tables[,16,,i] =  F_life_tables[,15,,i] *F_life_tables[,12,,i]
  }
  for(i in 1:n_sim) {		
    all_param$M[i] = weighted.mean(F_life_tables[,13,all_param$n_fert[i],all_param$n_mort[i]]+0.5, F_life_tables[,16,all_param$n_fert[i],all_param$n_mort[i]])
  }
  round(range( all_param$M), 1); round(mean( all_param$M), 1)				
  
  # TFR is defined by mortality rates and age structure of fertility
  for(k in 1:dim(F_life_tables)[4]) {
    stopifnot(F_life_tables[,"1ca",,k] == (both[k]) * exp(-r * (F_life_tables[,"age",,k] + 0.5)) * F_life_tables[,"Lx",,k]/F_life_tables[1,"lx",,k])
  }
  
  for(k in 1:dim(F_life_tables)[4]) {
    for(j in 1:dim(F_life_tables)[3]) {
      up = both[k]/sum(F_life_tables[,"fx",j,k]*F_life_tables[,"1ca",j,k]/2)
      F_life_tables[,"fx",j,k] =  F_life_tables[,"fx",j,k]*up
    }}
  
  check =matrix(NA,  dim(F_life_tables)[3], dim(F_life_tables)[4])
  for(k in 1:dim(F_life_tables)[4]) {
    for(j in 1:dim(F_life_tables)[3]) {
      check[j,k] =  sum(exp(-r * (F_life_tables[12:50,"age",j,k] + 0.5)) * (F_life_tables[12:50,"Lx",j,k]/F_life_tables[1,"lx",j,k])  *  F_life_tables[12:50,"fx",j,k]/2)
    }}
  stopifnot(unique(round(check, 9) == 1))
  
  # TFR and life expectancy
  for(i in 1:n_sim) {		
    all_param$e0[i] = (F_life_tables[1,'ex',all_param$n_fert[i],all_param$n_mort[i]])
    all_param$TFR[i] = sum(F_life_tables[,'fx',all_param$n_fert[i],all_param$n_mort[i]])
  }
  
  #_____________________________________________________________________________________________
  # 4. Transfert to SOCSIM
  # Mortality
  taux_TM = cbind(
    c(c(0, 1), 2:100),
    c(c(1,0), rep(0, 99)))
  birth_f = cbind(c(9:51, 100), rep(0, 44))
   
  RATES_M_M = as.list(1:nrow(all_param)); RATES_M_F = as.list(1:nrow(all_param)); MORT_SEG = as.list(1:nrow(all_param))
  for(i in 1:nrow(all_param)) {
    RATES_M_M[[i]] = cbind(taux_TM[2:nrow(taux_TM),], format(rbind(matrix(data = c(#all_param$nmrFq[i], 1-((1-all_param$pnmrFq[i])^(1/11)), 
      round(1-((1-M_life_tables[1:99,"nqx",all_param$n_fert[i],all_param$n_mort[i]])^(1/12)),12), 1-(0.01^(1/12))), nrow = (max_age), ncol = 1))[,1], scientific = F))
    RATES_M_F[[i]] = cbind(taux_TM[2:nrow(taux_TM),], format(rbind(matrix(data = c(#all_param$nmrFq[i], 1-((1-all_param$pnmrFq[i])^(1/11)), 
      round(1-((1-F_life_tables[1:99,"nqx",all_param$n_fert[i],all_param$n_mort[i]])^(1/12)),12), 1-(0.01^(1/12))), nrow = (max_age), ncol = 1))[,1], scientific = F))
    
    MORT_SEG[[i]]= 	
      MORT_SEG[[i]] =rbind(
        cbind("death ", 1, " M single"),
        RATES_M_M[[i]],
        cbind("death ", 1, " F single"),
        RATES_M_F[[i]])
  }
  
  # Fertility
  FERT_SEG = as.list(1:nrow(all_param))		
  for(i in 1:nrow(all_param)) {
    j=1
    FERT_SEG[[i]] = rbind(c("birth   ", j, paste("F single", 0, sep = ' ')), cbind(birth_f[,1:2], format(c(round(F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F)), 
                          c("birth   ", j, paste("F married", 0, sep = ' ')),  cbind(birth_f[,1:2], format(c(round(F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F)), 
                          c("birth   ", j, paste("F divorced", 0, sep = ' ')), cbind(birth_f[,1:2], format(c(round(F_life_tables[9:51,"fx",all_param$n_fert[i],all_param$n_mort[i]]/12,12), 0), sc = F)))
  }
  
  # Stable population: 
  starting_pop =  round(40000/exp(r*150), 0)/2
  
  pop = as.list(1:nrow(all_param))
  for(i in 1:nrow(all_param)) {
    pop[[i]] = cbind(append(round(starting_pop * M_life_tables[,12,all_param$n_fert[i],all_param$n_mort[i]], 0), round(starting_pop * F_life_tables[,12,all_param$n_fert[i],all_param$n_mort[i]], 0)), c(rep(0,101), rep(1,101)), c(rep(c(1), 202)), rep(3,202), rep(0:100,2), rep(99,202), rep(0,202), rep(0,202), rep(0,202), rep(0,202), rep(0,202), rep(0,202), rep(1,202), rep(0,202), rep(0,202))
    pop[[i]][,6] = (100 - pop[[i]][,5])*12
    pop[[i]] = pop[[i]][,-c(5)]
    pop[[i]] = as.data.frame(cbind(c(1:sum(pop[[i]][,1])),pop[[i]][rep(1:dim(pop[[i]])[1], times = pop[[i]][,1]),c(-1)][order(pop[[i]][rep(1:dim(pop[[i]])[1], times = pop[[i]][,1]),c(-1)][,4]),]))
    colnames(pop[[i]]) = c("person_id", "sex", "group", "next", "month_d_o_b", "mother", "father", "e_s_mom", "e_s_dad", "lborn", "last_marr", "mstatus", "month_d_o_d", "fmult")
  } 
  # Simulations will end in january 2020 (which is 1441 in CMC code)
  # according to socim, it will start in the month 1201 and end at the end of 3000
  # so to convert to cmc, deduce 1560 to each month in socsim
  # CMC = (Year-1900 ) * 12 + Month
  
  # Universal marriage
  NUPT_SEG =      rbind(
    cbind("marriage ", 1, " M single"),
    cbind(taux_m[c(1,10,100),1:2], c(0, 0, 0.99)),
    cbind("marriage ", 1, " F single"),
    cbind(taux_m[c(1,10,100),1:2], c(0, 0, 0.99)),
    cbind("marriage ", 1, " M widowed"),
    cbind(taux_m[c(1,10,100),1:2], c(0, 0, 0.99)),
    cbind("marriage ", 1, " F widowed"),
    cbind(taux_m[c(1,10,100),1:2], c(0, 0, 0.99)))
  
  Rates = as.list(1:nrow(all_param))
  
  for(i in 1:n_sim) {		
    # Print starting pop
    cat(paste("write.table(as.data.frame(pop[[", i, "]]), file = \"", dirrates, "/", "Pop_", dirrates, "_", paste("sim", i, sep = "_"), ".opop\", quote = F, row.names=F, col.names=F); cat(NULL, file = \"", dirrates, "/", "Pop_", dirrates, "_", paste("sim", i, sep = "_"), ".opox\"); cat(NULL, file = \"", dirrates, "/", "Pop_", dirrates, "_", paste("sim", i, sep = "_"), ".omar\")",  sep = ""), file = paste(dirrates, "/", "Generate_Pop_", dirrates, "_", paste("sim", i, sep = "_"),".R", sep = ""))				
    source(paste(dirrates, "/", "Generate_Pop_", dirrates, "_", paste("sim", i, sep = "_"),".R", sep = ""))
    unlink(paste(dirrates, "/", "Generate_Pop_", dirrates, "_", paste("sim", i, sep = "_"),".R", sep = ""))
    # Create sup files
    SEGMENT= c(paste("\nduration", 150*12, sep = " "),
               "\nbint 12", 
               paste("\nsex_ratio", 0.5, sep = ' '),
               "\nhetfert 0", 
               "\nalpha  0",
               "\nbeta 1",
               paste("\ninclude ",  dirrates, "/", "Rates_", dirrates, "_", paste("sim", i, sep = "_"), sep = ""),
               "\n*birthtarget  1 958", 
               "\nrun")
    cat("segments", 1,
        paste("\ninput_file ",  dirrates,  "/", "Pop_", dirrates, "_", paste("sim", i, sep = "_"), sep = ""),
        paste("\noutput_file ", "Results_", dirrates, "/", "Pop_", dirrates, "_", paste("sim", i, sep = "_"), "_out", sep = ""),
        SEGMENT,
        file = paste(dirrates, "/", dirrates, "_", paste("sim", i, sep = "_"), ".sup", sep = ""))
    # RATES : segment 1 
    Rates[[i]] = rbind(c("*" , "SIM", i), 
                       MORT_SEG[[i]], 
                       unlist(FERT_SEG[[i]]),
                       unlist(NUPT_SEG))
    
    
    cat(paste("write.table(as.data.frame(Rates[[", i, "]]), file = \"", dirrates, "/", "Rates_", dirrates, "_", paste("sim", i, sep = "_"), "\", quote = F, row.names=F, col.names=F)",  sep = ""), 
        file = paste("TZA01", "/", "GenerateRate_", dirrates, "_", paste("sim", i, sep = "_"), ".R", sep = ""))		
    source(paste(dirrates, "/", "GenerateRate_", dirrates, "_", paste("sim", i, sep = "_"), ".R", sep = ""))
    unlink(paste(dirrates, "/", "GenerateRate_", dirrates, "_", paste("sim", i, sep = "_"), ".R", sep = ""))
  }
  
  if(FALSE){      
    for (isim in 1:n_sim) {
      # Launch simulations
      system(paste('fromdos ', dirrates, '/', dirrates, '_', paste("sim", isim, sep = "_"), '.sup ', sep = ''), wait = T) 
      system(paste('fromdos ', dirrates, '/', 'Rates_', dirrates, '_', paste("sim", isim, sep = "_"), sep = ''), wait = T) 
      system(paste('socsim ', dirrates, '/', dirrates, '_', paste("sim", isim, sep = "_"), '.sup ', isim, sep = ''), wait = T) 
      # pop
      final_pop <- data.frame(read.table(file = paste("Results_", dirrates, '/Pop_', dirrates, '_sim_', isim, "_out.opop", sep = ''))); colnames(final_pop) = c("pid", "sex", "group", "next", "month_dob", "mother", "father", "nesibm", "e_s_dad", "lborn", "last_marr", "mstatus", "month_dod", "fmult")
      save(final_pop, file = paste("Results_", dirrates, '/Pop_', dirrates, '_sim_', isim, "_out.rda", sep = ''), ascii = F)
      # marriage
      marr_file <- data.frame(read.table(file = paste("Results_", dirrates, '/Pop_', dirrates, '_sim_', isim, "_out.omar", sep = '')))
      colnames(marr_file) = c("union_id", "wife_id", "husband_id", "start", "end", "reason", "wife_prior", "husband_prior")
      save(marr_file, file = paste("Results_", dirrates, '/Marr_', dirrates, '_sim_', isim, "_out.rda", sep = ''), ascii = F)
      # 
      system(paste("rm ", "Results_", dirrates,  "/*.opop", sep = ''))
      system(paste("rm ", "Results_", dirrates, "/*.opox", sep = ''))
      system(paste("rm ", "Results_", dirrates, "/*.pyr", sep = ''))
      system(paste("rm ", "Results_", dirrates, "/*.otx", sep = ''))
    }
  }
  
  # Re-estimate  fertility and mortality rates from simulations
  #____________________________________________________________________________________________________________________________________________________________
  all_param$growth= NA; all_param$livingpop = NA; all_param$totalpop = NA
  for(isim in 1:n_sim) {
      load(paste("C:\\SIMULATIONS\\Results_", dirrates, "\\", "Pop_", dirrates, "_sim_", isim, "_out.rda", sep = ""))
      file = final_pop;
      names(file) = sub("e_s_mom", "nesibm", names(file)); names(file) = sub("d_o_b", "dob", names(file)); names(file) = sub("d_o_d", "dod", names(file)); names(file) = sub("person_id", "pid", names(file))
      names(file) = sub("last_marr", "union_id", names(file))
      stopifnot(nrow(file[file$mother == 0 & file$nesibm !=0,]) == 0)
      cat("\nSIM:", isim)
      cat("\nLiving pop:", nrow(file[file$month_dod == 0,]))
      cat("\nGrowth rate:", 100*(log(nrow(file[file$month_dod == 0,])/nrow(file[file$mother == 0,])))/150)
      all_param$livingpop[isim] = nrow(file[file$month_dod == 0,])
      all_param$growth[isim] = 100*(log(nrow(file[file$month_dod == 0,])/nrow(file[file$mother == 0,])))/150
      all_param$totalpop[isim] = nrow(file)
    }
    plot(density(all_param$growth/100), main = "Density of r", xlim = c(1.90/100, 2.20/100), ylab = "Density")					
    rug(all_param$growth/100, col = 'grey50')
  
  for(isim in 1:n_sim) {
      cat("\n", "** Start -", paste("sim", isim, sep = "_"), "**", "\n")
      
     # 1. Load file, compute d_o_d and d_o_b
      load(file = paste("Results_", dirrates, '/Pop_', dirrates,'_sim_', isim, "_out.rda", sep = ''))
      file = final_pop
      names(file) = sub("e_s_mom", "nesibm", names(file)); names(file) = sub("d_o_b", "dob", names(file)); names(file) = sub("d_o_d", "dod", names(file)); names(file) = sub("person_id", "pid", names(file))
      names(file) = sub("last_marr", "union_id", names(file))
      
      stopifnot(nrow(file[file$mother == 0 & file$nesibm !=0,]) == 0)
      cat("Living pop:", nrow(file[file$month_dod == 0,]))
      cat("Growth rate:", 100*(log(nrow(file[file$month_dod == 0,])/nrow(file[file$mother == 0,])))/150)
      
      file$month_dod[file$month_dod == 0] = NA
      
      start_sim = max(file$month_dob[file$mother == 0]+1, na.rm = T); end_sim = max(c(file$month_dod, file$month_dob),na.rm = TRUE) + 1; length_sim = end_sim - start_sim
      time.cut = seq(start_sim, end_sim + 12, by = 12); time.cut.labels = as.character(time.cut[1:(length(time.cut)-1)])
      
      w = file[,c(grep("pid|dob|mother|father|dod|sex|union_id|nesibm", colnames(file)))]
      
      # 1. DOB and DOD
      w$journ = trunc(runif(nrow(w), min =0 , max = 30)); w$journ[w$month_dob ==  start_sim] = trunc(runif(nrow(w[w$month_dob ==  start_sim,]), min = 1 , max = 30))
      w$dob = w$month_dob + w$journ/31
      w$month_dod[w$month_dod == 0] = NA
      w$jourd = trunc(runif(nrow(w), min =0 , max = 30)); w$jourd[!is.na(w$month_dod) & w$month_dod ==  start_sim] = trunc(runif(nrow(w[!is.na(w$month_dod) & w$month_dod ==  start_sim,]), min = 1 , max = 30))
      w$jourd[which(w$month_dod == w$month_dob)] = trunc(runif(nrow(w[which(w$month_dod == w$month_dob),]), min =w$journ[which(w$month_dod == w$month_dob)]+1, max = 30))
      w$dod = w$month_dod + w$jourd/31
      w$dob[!is.na(w$dod) & w$dob == w$dod] = w$dob[!is.na(w$dod) & w$dob == w$dod] - 1/31
      
      w = w[order(w$pid),]; stopifnot(which(rownames(w) != w$pid) == 0)
      
      # Slight change if age of mothers at birth is a multiple of 12 months 
      w$dob_m = c(rep(NA,length(w$mother[w$mother ==0])), w$dob[w$mother[(length(w$mother[w$mother == 0]) +1):nrow(w)]])
      w$dob[which(abs(floor((w$dob - w$dob_m)/12) - (w$dob - w$dob_m)/12) == 0)] = w$dob[which(abs(floor((w$dob - w$dob_m)/12) - (w$dob - w$dob_m)/12) == 0)] + 1/31 
      
      # Check that children are not born after their mother's death and get mother's dob
      w$dod_m = c(rep(NA,length(w$mother[w$mother ==0])), w$dod[w$mother[(length(w$mother[w$mother == 0]) +1):nrow(w)]])
      w$dob[!is.na(w$dod_m) & w$dob >= w$dod_m] = w$dob[!is.na(w$dod_m) & w$dob >= w$dod_m] - 1
      w$dob_m = c(rep(NA,length(w$mother[w$mother ==0])), w$dob[w$mother[length(w$mother[w$mother == 0]):nrow(w)]])
      
      # Date of birth of next eldest sibling
      w$dob_nesibm = w$dob[zna(w$nesibm)]
      while( nrow(w[!is.na(w$dob_nesibm) & w$dob <= w$dob_nesibm,]) != 0) {
        w$dob[!is.na(w$dob_nesibm) & w$dob <= w$dob_nesibm] = w$dob_nesibm[!is.na(w$dob_nesibm) & w$dob <= w$dob_nesibm] + 1/31
        w$dob_nesibm = w$dob[zna(w$nesibm)]
      }
      w$dod[!is.na(w$dod) & w$dob >= w$dod] = w$dob[!is.na(w$dod) & w$dob >= w$dod] + 1/31
      w$dob[w$mother > 0 & w$dob < start_sim] = start_sim + 0.01 
      w$dob[w$dob > end_sim] = end_sim - 0.01
      
      # 2.  Estimation of fertility rates 
      bio.fert = w[,-c(grep("dod_m|month|jour|bornfromHIV|partum|dob_m", colnames(w)))]
      
      bio.fert$rang[bio.fert$mother > 0] = ave(bio.fert$dob[bio.fert$mother > 0], bio.fert$mother[bio.fert$mother > 0], FUN = rank)
      
      # Locate births into mothers' reproductive periods and recover information about siblings
      gen_hist <- bio.fert[!is.na(bio.fert$rang), grep("^dob$|^dod$|^mother|rang", colnames(bio.fert))]; names(gen_hist) = sub("mother", "pid", sub("rang", "parite", sub("dob", "child_dob", sub("dod", "child_dod", names(gen_hist)))))
      gen_hist = cbind(gen_hist, 1); names(gen_hist) = sub("^1$", "event", names(gen_hist))
      gen_hist$parite = as.numeric(gen_hist$parite); stopifnot(nrow(file[file$mother > 0,]) == nrow(gen_hist))
      
      bio.fert = bio.fert[w$sex == 1,-c(grep("mother|rang", names(bio.fert)))]
      
      gen_hist$DF = ave(gen_hist$parite, gen_hist$pid, FUN = max)
      
      last_born = gen_hist[,-c(grep("event|DF|dod", names(gen_hist)))]; last_born$parite = last_born$parite + 1; names(last_born)[grep("child_dob", names(last_born))] = "last_born"
      gen_hist = merge(gen_hist, last_born, by=c("pid","parite"), all.x = TRUE, all.y = TRUE); gen_hist$parite = gen_hist$parite - 1
      
      gen_hist <- merge(gen_hist, bio.fert[,c(grep("pid|dob|dod|cmc.inf", names(bio.fert)))], by="pid", all.x = TRUE, all.y = TRUE)
      gen_hist$parite[is.na(gen_hist$parite)] = 0; gen_hist$DF[is.na(gen_hist$DF)] = gen_hist$parite[is.na(gen_hist$DF)]; gen_hist$event[is.na(gen_hist$event)] = 0
      
      gen_hist$start_sim = start_sim; gen_hist$end_sim = end_sim
      
      gen_hist = gen_hist[order(gen_hist$DF, decreasing = T),]
      
      # Follow-up time
      gen_hist$start.time = gen_hist$dob
      gen_hist$start.time[gen_hist$dob < gen_hist$start_sim] = gen_hist$start_sim[gen_hist$dob < gen_hist$start_sim]
      gen_hist$start.time[!is.na(gen_hist$parite) & gen_hist$parite >= 1] = gen_hist$last_born[!is.na(gen_hist$parite) & gen_hist$parite >= 1]
      
      gen_hist$end.time = gen_hist$dod 
      gen_hist$end.time[!is.na(gen_hist$parite)] = gen_hist$child_dob[!is.na(gen_hist$parite)] 
      gen_hist$end.time[!is.na(gen_hist$parite) & gen_hist$parite == gen_hist$DF] = gen_hist$dod[!is.na(gen_hist$parite) & gen_hist$parite == gen_hist$DF] 
      gen_hist$end.time[is.na(gen_hist$end.time)] = end_sim
      
      gen_hist$fu.time = gen_hist$end.time - gen_hist$start.time
      
      # Survival object	
      gen_hist_Surv = with(gen_hist, Surv(time = gen_hist$fu.time, event = gen_hist$event, type = "right"))
      
      annee <- tcut(gen_hist$start.time, time.cut, labels=time.cut.labels)
      
      gen_hist$startage = (gen_hist$start.time - gen_hist$dob)/12
      age <- tcut(gen_hist$startage*12, age_ref, labels=as.character(age_ref[1:(length(age_ref)-1)]))
      
      gen_hist$parite[is.na(gen_hist$parite)] = 0
      
      gen_hist_pp <- pyears(gen_hist_Surv ~ annee + age + parite, gen_hist, scale = 12, data.frame=TRUE)
      
      # Check that offtable  = 0
      stopifnot(round(sum(gen_hist$fu.time)/12 - sum(gen_hist_pp$data$pyears), 4) == 0)
      
      fert_pp <- gen_hist_pp$data; stopifnot(nrow(fert_pp[fert_pp$event > 1 & fert_pp$pyears == 0,]) == 0)
      # Check that we have no births before age 10
      stopifnot(min(fert_pp$event[as.numeric(as.character(fert_pp$age))/12 < 10 | as.numeric(as.character(fert_pp$age))/12 >= 50]) == 0)
      
      fert_pp = fert_pp[as.numeric(as.character(fert_pp$age))/12 >= 10 & as.numeric(as.character(fert_pp$age))/12 < 50,]
      
      fert_pp$age <- as.factor(as.numeric(as.character(fert_pp$age))/12)
      fert_pp$annee<- (as.numeric(as.character(fert_pp$annee))-start_sim)/12 + 1870
      fert_pp$annee = as.factor(fert_pp$annee)

      fert_pp$sim.id = isim
      save(fert_pp, file = paste(getwd(),  "\\Results_", dirrates, "\\fert_pp_",  paste("sim", isim, sep = "_"),   sep = ""), ascii = FALSE)
      
      # The number of births must be equal to the number of individuals generated by the simulation 
      stopifnot(abs(sum(fert_pp$event) - nrow(w[w$mother != 0,])) < 4) 
      
      # 3. Estimation of mortality rates 
      bio.mort = w[,-c(grep("month|jour|dod_m|dob_m|bornfromHIV|intrapartum", colnames(w)))]
      # Follow-up time for the deceased and the survivors
      bio.mort$fu.time = (bio.mort$dod - bio.mort$dob)
      bio.mort$fu.time[is.na(bio.mort$fu.time)] = (end_sim - bio.mort$dob[is.na(bio.mort$fu.time)])
      # Folluw-up time for deceased persons born before start_sim
      bio.mort$fu.time[!is.na(bio.mort$dod) &  (bio.mort$dob < start_sim)] = (bio.mort$dod[!is.na(bio.mort$dod) &  (bio.mort$dob < start_sim)] - start_sim)
      stopifnot(nrow(bio.mort[is.na(bio.mort$fu.time),]) == 0)
      bio.mort$event = as.numeric(!is.na(bio.mort$dod))
      # Survival object	
      bio.mort_Surv = with(bio.mort, Surv(time = bio.mort$fu.time, event = bio.mort$event, type = "right"))
      bio.mort$startage = NA
      bio.mort$startage[which(bio.mort$dob >= start_sim)] = 0 
      bio.mort$startage[which(bio.mort$dob < start_sim)] = (start_sim - bio.mort$dob[which(bio.mort$dob < start_sim)])/12
      stopifnot(nrow(bio.mort[is.na(bio.mort$startage),]) == 0)
      
      age <- tcut(bio.mort$startage*12, age_ref, labels=as.character(age_ref[1:(length(age_ref)-1)]))
      
      time.cut = seq(start_sim, end_sim + 12, by = 12)
      time.cut.labels = as.character(time.cut[1:(length(time.cut)-1)])
      
      annee <- tcut(bio.mort$dob, time.cut, labels=time.cut.labels)
      annee[bio.mort$dob < start_sim] <- tcut(start_sim, time.cut, labels=time.cut.labels)
      
      bio.mort_pp <- pyears(bio.mort_Surv ~ annee + age + sex , bio.mort, scale = 12, data.frame=TRUE)
      stopifnot(round(sum(bio.mort$fu.time)/12 - sum(bio.mort_pp$data$pyears), 4) == 0)
      mort_pp <- bio.mort_pp$data; stopifnot(nrow(mort_pp[mort_pp$event > 1 & mort_pp$pyears == 0,]) == 0)
      mort_pp$age <- as.numeric(as.character(mort_pp$age))/12
      mort_pp$age[mort_pp$age > 80] = 80
      mort_pp$agegrp = trunc(mort_pp$age/5)*5
      mort_pp$annee<- (as.numeric(as.character(mort_pp$annee))-start_sim)/12 + 1870
      
      # Check : Total female exposure = fertility exposure
      stopifnot(round(sum(gen_hist$fu.time)/12  - sum(mort_pp$pyears[mort_pp$sex == 1]), 5) == 0)
      
      mort_pp$sim.id = isim
      save(mort_pp, file = paste(getwd(),  "\\Results_", dirrates, "\\mort_pp_",  paste("sim", isim, sep = "_"),   sep = ""), ascii = FALSE)
      
      all = w; stopifnot(nrow(file) == nrow(all))
      save(all, file = paste(getwd(),  "\\Results_", dirrates, "\\all_",  paste("sim", isim, sep = "_"),   sep = ""), ascii = FALSE)
      
      warnings()
      cat("\n", "** End -", "sim ", isim, "**", "\n")
    }
  
  # Background mortality 
  sexes = c("Both", "female", "male")
  par(mai=c(0.5,1,0.5,0.5))
  plot(F_life_tables[,"tx",all_param$n_fert[1],all_param$n_mort[1]], type = 'l', lwd = 1, log = 'y', ylim = c(0.0001,1), xlim = c(0, 75), ylab = "tx", xlab = "Age")
  for (i in 1:12) {
    load(file = paste(getwd(),  "\\Results_", dirrates, "\\mort_pp_", paste("sim", isim, sep = "_"),   sep =""))          	   
    mort_pp$sex = as.numeric(as.character(mort_pp$sex))	
    points(tapply(mort_pp$event[mort_pp$annee >= 1900 & mort_pp$sex < 2], list(mort_pp$age[mort_pp$annee >= 1900 & mort_pp$sex < 2]), sum)/tapply(mort_pp$pyears[mort_pp$annee >= 1900 & mort_pp$sex < 2], list(mort_pp$age[mort_pp$annee >= 1900 & mort_pp$sex < 2]), sum), type = 'l', col = palette_G[7], ylab = 'ntx')
  }
  points(F_life_tables[,"tx",all_param$n_fert[1],all_param$n_mort[1]], type = 'l', lwd = 1)
  
  for (i in 13:24) {
    load(file = paste(getwd(),  "\\Results_", dirrates, "\\mort_pp_", paste("sim", isim, sep = "_"),   sep =""))          	   
    mort_pp$sex = as.numeric(as.character(mort_pp$sex))	
    points(tapply(mort_pp$event[mort_pp$annee >= 1900 & mort_pp$sex < 2], list(mort_pp$age[mort_pp$annee >= 1900 & mort_pp$sex < 2]), sum)/tapply(mort_pp$pyears[mort_pp$annee >= 1900 & mort_pp$sex < 2], list(mort_pp$age[mort_pp$annee >= 1900 & mort_pp$sex < 2]), sum), type = 'l', col = palette_G[5], ylab = 'ntx')
  }
  points(F_life_tables[,"tx",all_param$n_fert[13],all_param$n_mort[13]], type = 'l', lwd = 1)
  
  for (i in 25:36) {
    load(file = paste(getwd(),  "\\Results_", dirrates, "\\mort_pp_", paste("sim", isim, sep = "_"),   sep =""))          	   
    mort_pp$sex = as.numeric(as.character(mort_pp$sex))	
    points(tapply(mort_pp$event[mort_pp$annee >= 1900 & mort_pp$sex < 2], list(mort_pp$age[mort_pp$annee >= 1900 & mort_pp$sex < 2]), sum)/tapply(mort_pp$pyears[mort_pp$annee >= 1900 & mort_pp$sex < 2], list(mort_pp$age[mort_pp$annee >= 1900 & mort_pp$sex < 2]), sum), type = 'l',  col = palette_G[3], ylab = 'ntx')
  }
  points(F_life_tables[,"tx",all_param$n_fert[25],all_param$n_mort[25]], type = 'l', lwd = 1)
  
  for (i in 37:48) {
    load(file = paste(getwd(),  "\\Results_", dirrates, "\\mort_pp_", paste("sim", isim, sep = "_"),   sep =""))          	   
    mort_pp$sex = as.numeric(as.character(mort_pp$sex))	
    points(tapply(mort_pp$event[mort_pp$annee >= 1900 & mort_pp$sex < 2], list(mort_pp$age[mort_pp$annee >= 1900 & mort_pp$sex < 2]), sum)/tapply(mort_pp$pyears[mort_pp$annee >= 1900 & mort_pp$sex < 2], list(mort_pp$age[mort_pp$annee >= 1900 & mort_pp$sex < 2]), sum), type = 'l',  col = palette_G[1], ylab = 'ntx')
  }
  points(F_life_tables[,"tx",all_param$n_fert[37],all_param$n_mort[37]], type = 'l', lwd = 1)
  
  legend("bottomright", expression(paste(alpha ==0.2), paste(alpha == -0.2), paste(alpha == -0.6), paste(alpha == -1)), bty= 'n', lty = c(1, 1, 1,1), col = palette_G[-1 +2*(4:1)], cex = 1, lwd = c(2,2,2,2))
  
  # Format mort_pp to be used with the function life table from the demography package
  for(isim in 1:n_sim) {
    
    load(file = paste(getwd(),  "\\Results_", dirrates, "\\mort_pp_", paste("sim", isim, sep = "_"),   sep =""))          	   
    
    # Re-estimation of rates from lifetable (after 1950)
    exposure_M = tapply(mort_pp$pyears[mort_pp$annee >= 1950 & mort_pp$annee < 2050], list(mort_pp$age[mort_pp$annee >= 1950 & mort_pp$annee < 2050], mort_pp$annee[mort_pp$annee >= 1950 & mort_pp$annee < 2050], mort_pp$sex[mort_pp$annee >= 1950 & mort_pp$annee < 2050]), sum)[,,1]
    exposure_F = tapply(mort_pp$pyears[mort_pp$annee >= 1950 & mort_pp$annee < 2050], list(mort_pp$age[mort_pp$annee >= 1950 & mort_pp$annee < 2050], mort_pp$annee[mort_pp$annee >= 1950 & mort_pp$annee < 2050], mort_pp$sex[mort_pp$annee >= 1950 & mort_pp$annee < 2050]), sum)[,,2]
    
    exp_serie_F = cbind(rep(1950:2019, each = length(unique(mort_pp$age))), rep(sort(unique(mort_pp$age)), length(1950:2019)), as.vector(exposure_F), as.vector(exposure_F))
    colnames(exp_serie_F) = c("Year", "Age", "Female", "Male")
    exp_serie_F[,3] = ifelse(exp_serie_F[,3]  < 0.0001 | is.na(exp_serie_F[,3]), 0.0001, exp_serie_F[,3]); exp_serie_F[,4] = ifelse(exp_serie_F[,4]  < 0.0001 | is.na(exp_serie_F[,4]), 0.0001, exp_serie_F[,4])
    
    exp_serie_M = cbind(rep(1950:2019, each = length(unique(mort_pp$age))), rep(sort(unique(mort_pp$age)), length(1950:2019)), as.vector(exposure_M), as.vector(exposure_M))
    colnames(exp_serie_M) = c("Year", "Age", "Female", "Male")
    exp_serie_M[,3] = ifelse(exp_serie_M[,3]  < 0.0001 | is.na(exp_serie_M[,3]), 0.0001, exp_serie_M[,3]); exp_serie_M[,4] = ifelse(exp_serie_M[,4]  < 0.0001 | is.na(exp_serie_M[,4]), 0.0001, exp_serie_M[,4])
    
    taux_M = (tapply(mort_pp$event[mort_pp$annee >= 1950 & mort_pp$annee < 2050], list(mort_pp$age[mort_pp$annee >= 1950 & mort_pp$annee < 2050], mort_pp$annee[mort_pp$annee >= 1950 & mort_pp$annee < 2050], mort_pp$sex[mort_pp$annee >= 1950 & mort_pp$annee < 2050]), sum)/tapply(mort_pp$pyears[mort_pp$annee >= 1950 & mort_pp$annee < 2050], list(mort_pp$age[mort_pp$annee >= 1950 & mort_pp$annee < 2050], mort_pp$annee[mort_pp$annee >= 1950 & mort_pp$annee < 2050], mort_pp$sex[mort_pp$annee >= 1950 & mort_pp$annee < 2050]), sum))[,,1]
    taux_F = (tapply(mort_pp$event[mort_pp$annee >= 1950 & mort_pp$annee < 2050], list(mort_pp$age[mort_pp$annee >= 1950 & mort_pp$annee < 2050], mort_pp$annee[mort_pp$annee >= 1950 & mort_pp$annee < 2050], mort_pp$sex[mort_pp$annee >= 1950 & mort_pp$annee < 2050]), sum)/tapply(mort_pp$pyears[mort_pp$annee >= 1950 & mort_pp$annee < 2050], list(mort_pp$age[mort_pp$annee >= 1950 & mort_pp$annee < 2050], mort_pp$annee[mort_pp$annee >= 1950 & mort_pp$annee < 2050], mort_pp$sex[mort_pp$annee >= 1950 & mort_pp$annee < 2050]), sum))[,,2]
    
    taux_serie_F = cbind(rep(1950:2019, each = length(unique(mort_pp$age))), rep(sort(unique(mort_pp$age)), length(1950:2019)), as.vector(taux_F), as.vector(taux_F))
    colnames(taux_serie_F) = c("Year", "Age", "Female", "Male")
    taux_serie_F[,3] = ifelse(taux_serie_F[,3]  < 0.0001 | is.na(taux_serie_F[,3]), 0.0001, taux_serie_F[,3]); taux_serie_F[,4] = ifelse(taux_serie_F[,4]  < 0.0001 | is.na(taux_serie_F[,4]), 0.0001, taux_serie_F[,4])
    
    taux_serie_M = cbind(rep(1950:2019, each = length(unique(mort_pp$age))), rep(sort(unique(mort_pp$age)), length(1950:2019)), as.vector(taux_M), as.vector(taux_M))
    colnames(taux_serie_M) = c("Year", "Age", "Female", "Male")
    taux_serie_M[,3] = ifelse(taux_serie_M[,3]  < 0.0001 | is.na(taux_serie_M[,3]), 0.0001, taux_serie_M[,3]); taux_serie_M[,4] = ifelse(taux_serie_M[,4]  < 0.0001 | is.na(taux_serie_M[,4]), 0.0001, taux_serie_M[,4])
    
    write.table(taux_serie_F, file = paste(getwd(),  "\\Results_", dirrates, "\\taux_F_",    paste("sim", isim, sep = "_"), ".txt", sep = ""), quote = FALSE, row.names=FALSE, col.names=TRUE)
    write.table(exp_serie_F,  file = paste(getwd(),  "\\Results_", dirrates, "\\exposure_F_", paste("sim", isim, sep = "_"), ".txt", sep = ""), quote = FALSE, row.names=FALSE, col.names=TRUE)
    write.table(taux_serie_M, file = paste(getwd(),  "\\Results_", dirrates, "\\taux_M_",    paste("sim", isim, sep = "_"), ".txt", sep = ""), quote = FALSE, row.names=FALSE, col.names=TRUE)
    write.table(exp_serie_M,  file = paste(getwd(),  "\\Results_", dirrates, "\\exposure_M_", paste("sim", isim, sep = "_"), ".txt", sep = ""), quote = FALSE, row.names=FALSE, col.names=TRUE)
    
  }
  
    F_fromsim = read.demogdata(paste(getwd(),  "\\Results_", dirrates, "\\taux_F_",  paste("sim", isim, sep = "_"), ".txt", sep = ""), paste("C:\\SIMULATIONS\\Results_TZA01\\", "exposure_F_",  paste("sim", isim, sep = "_"), ".txt", sep = ""), type="mortality", skip = 0, popskip = 0, label= paste("sim", isim, sep = "_"))
    names(F_fromsim$rate) = c("female", "male")	 
    M_fromsim = read.demogdata(paste(getwd(),  "\\Results_", dirrates, "\\taux_M_",  paste("sim", isim, sep = "_"), ".txt", sep = ""), paste("C:\\SIMULATIONS\\Results_TZA01\\", "exposure_M_",  paste("sim", isim, sep = "_"), ".txt", sep = ""), type="mortality", skip = 0, popskip = 0, label= paste("sim", isim, sep = "_"))
    names(M_fromsim$rate) = c("female", "male")	 
    
    x_lim = c(1950, 2020)
    # Life expectancy
    plot(1950:2019, lifetable(F_fromsim, series = "female", type = c("period"))$ex[1,], col = palette_G[6], type = 'l', lwd = 2, lty = 1, ylim = c(35, 75), xlim = x_lim, main= i, xlab = "", ylab = expression(paste("e"[0])))
    points(1950:2019, lifetable(M_fromsim, series = "male", type = c("period"))$ex[1,], col = palette_B[6], type = 'l', lwd = 2, lty = 1, xlim = c(1980, 2010), ylim = c(35, 75), xlab = "", ylab = expression(paste("e"[0])))
    abline(h = F_life_tables[1,"ex",all_param$n_fert[isim],all_param$n_mort[isim]], col = palette_O[4], lwd = 2) 
    # 45q15
    plot(1950:2019, 1-lifetable(F_fromsim, series = "female", type = c("period"))$lx[grep("^60$", rownames(   lifetable(F_fromsim, series = "female", type = c("period"))$lx)),]/lifetable(F_fromsim, series = "female", type = c("period"))$lx[grep("^15$", rownames(   lifetable(F_fromsim, series = "female", type = c("period"))$lx)),], col = palette_G[6], type = 'l', lwd = 2, lty = 1, ylim = c(0, 0.8), xlim = x_lim,  xlab = "", ylab = expression(paste(" "[45], "q"["15"])))
    points(1950:2019, 1-lifetable(M_fromsim, series = "male", type = c("period"))$lx[grep("^60$", rownames(   lifetable(M_fromsim, series = "male", type = c("period"))$lx)),]/lifetable(M_fromsim, series = "male", type = c("period"))$lx[grep("^15$", rownames(   lifetable(M_fromsim, series = "male", type = c("period"))$lx)),], col = palette_B[6], type = 'l', lwd = 2, lty = 1, ylim = c(0, 0.8), xlim = x_lim,  xlab = "", ylab = expression(paste(" "[45], "q"["15"])))
    abline(h =  1- F_life_tables[61,"lx",all_param$n_fert[isim],all_param$n_mort[isim]]/F_life_tables[16,"lx",all_param$n_fert[isim],all_param$n_mort[isim]], col = palette_O[4], lwd = 2) 
    # 5q0
    plot(1950:2019, 1-lifetable(F_fromsim, series = "female", type = c("period"))$lx[grep("^5$", rownames(   lifetable(F_fromsim, series = "female", type = c("period"))$lx)),]/lifetable(F_fromsim, series = "female", type = c("period"))$lx[grep("^0$", rownames(   lifetable(F_fromsim, series = "female", type = c("period"))$lx)),], col = palette_G[6], type = 'l', lwd = 2, lty = 1, ylim = c(0, 0.5), xlim = x_lim, main= "", xlab = "", ylab = expression(paste(" "[5], "q"["0"])))
    points(1950:2019, 1-lifetable(M_fromsim, series = "male", type = c("period"))$lx[grep("^5$", rownames(   lifetable(M_fromsim, series = "male", type = c("period"))$lx)),]/lifetable(M_fromsim, series = "male", type = c("period"))$lx[grep("^0$", rownames(   lifetable(M_fromsim, series = "male", type = c("period"))$lx)),], col = palette_B[6], type = 'l', lwd = 2, lty = 1, ylim = c(0, 0.3), xlim = x_lim,  xlab = "", ylab = expression(paste(" "[5], "q"["0"])))
    abline(h =  1- F_life_tables[6,"lx",all_param$n_fert[isim],all_param$n_mort[isim]]/F_life_tables[1,"lx",all_param$n_fert[isim],all_param$n_mort[isim]], col = palette_O[4], lwd = 2) 
    # nqx 
    plot(apply(lifetable(M_fromsim, series = "male", type = c("period"))$qx, 1, mean), type = 'l', log = 'y')
    points(F_life_tables[,"nqx",all_param$n_fert[isim],all_param$n_mort[isim]], type = 'l', lwd = 2, col = 'grey50')
    # lx 
    plot(0:80, apply(lifetable(M_fromsim, series = "male", type = c("period"))$lx[,which(lifetable(M_fromsim, series = "male", type = c("period"))$ex[1,] < 80)], 1, median), type = 'l', xlab = 'Age', ylab = 'lx')
    points(0:100, F_life_tables[,"lx",all_param$n_fert[isim],all_param$n_mort[isim]]/10000, type = 'l', lwd = 2, col = 'grey50')
    
  # Check on ASFRs
  for (isim in 1:n_sim) {
    load(file = paste(getwd(),  "\\Results_", dirrates,  "\\fert_pp_", paste("sim", isim, sep = "_"),   sep =""))          	   
    fert_pp$annee = as.numeric(as.character(fert_pp$annee))
    fert_pp$age = as.numeric(as.character(fert_pp$age))
    plot(sort(unique(fert_pp$age)), tapply(fert_pp$event[fert_pp$annee >= 1950 ], list(fert_pp$age[fert_pp$annee >= 1950 ]), sum)/tapply(fert_pp$pyears[fert_pp$annee >= 1950 ], list(fert_pp$age[fert_pp$annee >= 1950 ]), sum), type = 'l', col = 3, ylab = 'ntx', xlim = c(10, 50), lwd=2, main = i)
    points(0:100, F_life_tables[,"fx",all_param$n_fert[isim],all_param$n_mort[isim]], type = 'l', lwd = 2, col = 'grey50')
  }
  # Mean age at birth of children
  for(isim in 1:n_sim) {		
    all_param$M[isim] = weighted.mean(F_life_tables[,13,all_param$n_fert[isim],all_param$n_mort[isim]]+0.5, F_life_tables[,"fx",all_param$n_fert[isim],all_param$n_mort[isim]]*F_life_tables[,"1ca",all_param$n_fert[isim],all_param$n_mort[isim]])
    load(file = paste(getwd(),  "\\Results_", dirrates, "\\fert_pp_", paste("sim", isim, sep = "_"),   sep =""))          	   
    all_param$MACB_sim[isim] = weighted.mean(10:49 + 0.5, tapply(fert_pp$event, fert_pp$age, sum))
  } 
  
  # save(all_param, file = paste0('all_param_', growthrate, '.rda'), ascii = FALSE)
} # end of loop by growth rate
} # end of condition re-do
  
  # recover all_param depending on growth rates
  load(file = 'all_param_1.rda')
  all_param_1 = all_param; all_param_1$isim = 1:nrow(all_param_1)
  load(file = 'all_param_3.rda')
  all_param_3 = all_param; all_param_3$isim = 1:nrow(all_param_3)
  
	# Re-organize the mortality and fertility parameters
 	alpha_mort = c(0.2, -0.2, -0.6, -1.0)
	beta_mort =  c(0.7, 1.1)
	n_mort = length(alpha_mort)*length(beta_mort)
	mort_param <- expand.grid(alpha_mort=alpha_mort, beta_mort=beta_mort)
	mort_param$n_mort = 1:n_mort
	
	alpha_fert =  c(-0.5, -0.2,0.1,0.4)
	beta_fert =  c(1.0, 1.15, 1.4, 1.8)
	n_fert = length(alpha_fert)*length(beta_fert)
	 
 	fert_param <- expand.grid(alpha_fert=alpha_fert, beta_fert=beta_fert)
	fert_param$n_fert = 1:n_fert
	
	all_param <- expand.grid(n_fert=1:n_fert, n_mort=1:n_mort)
	all_param = merge(all_param, fert_param, by = "n_fert")
	all_param = merge(all_param, mort_param, by = "n_mort")
	all_param$r = 0.01
	
	all_param = rbind(all_param, all_param)
	all_param$r[c(nrow(all_param)/2+1):nrow(all_param)] = 0.03
	
	# exclude small beta_f when r = 0.01
	all_param = all_param[-which(all_param$r == 0.01 & all_param$beta_fert == 1),]
	# exclude high beta_f when r = 0.03
	all_param = all_param[-which(all_param$r == 0.03 & all_param$beta_fert == 1.8),]
	# number of sim
	n_sim = nrow(all_param)
  all_param$sim = 1:nrow(all_param)
  all_param$isim1 = NA; all_param$isim3 = NA
  for(i in which(all_param$r == 0.01)){
    all_param$isim1[i] = all_param_1$isim[all_param_1$alpha_mort == all_param$alpha_mort[i] &
                                          all_param_1$beta_mort  == all_param$beta_mort[i] &
                                          all_param_1$alpha_fert == all_param$alpha_fert[i] &
                                          all_param_1$beta_fert  == all_param$beta_fert[i]]
  }
  for(i in which(all_param$r == 0.03)){
    all_param$isim3[i] = all_param_3$isim[all_param_3$alpha_mort == all_param$alpha_mort[i] &
                                            all_param_3$beta_mort  == all_param$beta_mort[i] &
                                            all_param_3$alpha_fert == all_param$alpha_fert[i] &
                                            all_param_3$beta_fert  == all_param$beta_fert[i]]
  }

# ___________________________________________________________________________________________________________
# Life tables for all simulations
	F_life_tables <- array(0, dim = c(101,1 + 15,n_fert, n_mort)); dimnames(F_life_tables) = list(c(0:max_age),c("npx", "nqx", "nax", "nx", "tx", "lx", "dx", "Lx", "Tx", "ex", "expaLx", "1ca", "age", "age5", "fx", "nb_naiss"), as.character(1:n_fert), as.character(1:n_mort)) 
	M_life_tables <- array(0, dim = c(101,1 + 15,n_fert, n_mort)); dimnames(M_life_tables) = list(c(0:max_age),c("npx", "nqx", "nax", "nx", "tx", "lx", "dx", "Lx", "Tx", "ex", "expaLx", "1ca", "age", "age5", "fx", "nb_naiss"), as.character(1:n_fert), as.character(1:n_mort)) 

	# Exact ages
	F_life_tables[,13,,] = c(0:(dim(F_life_tables)[1] - 1)) 
	M_life_tables[,13,,] = c(0:(dim(M_life_tables)[1] - 1)) 
		            
	# Survivors
	F_life_tables[1,6,,] = 10000; M_life_tables[1,6,,] = 10000 
	for(i in 2:(nrow(F_life_tables))) {
		for(j in 1:dim(F_life_tables)[4]) {
		F_life_tables[i,6,,j] = (1/(1+exp(2*mort_param[j,1]+2*mort_param[j,2]*as.numeric(logitsBrass[,]))))[i-1]*10000
		M_life_tables[i,6,,j] = (1/(1+exp(2*mort_param[j,1]+2*mort_param[j,2]*as.numeric(logitsBrass[,]))))[i-1]*10000
	}}

	# Deaths
	for(i in 1:(nrow(F_life_tables)-1)) {
			F_life_tables[i,7,,] = F_life_tables[i,6,,] - F_life_tables[i+1,6,,]
			M_life_tables[i,7,,] = M_life_tables[i,6,,] - M_life_tables[i+1,6,,]
	}
					
	F_life_tables[nrow(F_life_tables),7,,] = F_life_tables[nrow(F_life_tables),6,,]
	M_life_tables[nrow(M_life_tables),7,,] = M_life_tables[nrow(M_life_tables),6,,]
		
	# Survival probs
	F_life_tables[,1,,] = (F_life_tables[,6,,] - F_life_tables[,7,,]) / F_life_tables[,6,,]
	M_life_tables[,1,,] = (M_life_tables[,6,,] - M_life_tables[,7,,]) / M_life_tables[,6,,]
		
	# Smoothing at old ages
	for(i in 1:dim(F_life_tables)[3]) {
			for(j in 1:dim(F_life_tables)[4]) {
			F_life_tables[60:101,1,i,j] = spline(c(60:90, 101), c(F_life_tables[60:90,1,i,j], 0), xout = 60: 101)$y
			M_life_tables[60:101,1,i,j] = spline(c(60:90, 101), c(M_life_tables[60:90,1,i,j], 0), xout = 60: 101)$y
	}}
			
	# nqx
	F_life_tables[,2,,] = 1-F_life_tables[,1,,]
	M_life_tables[,2,,] = 1-M_life_tables[,1,,]
		
	M_life_tables[,6,i,j] = M_life_tables[1,6,i,j] * cumprod(c(1,1-M_life_tables[,2,i,j]))[1:101]
	F_life_tables[,6,i,j] = F_life_tables[1,6,i,j] * cumprod(c(1,1-F_life_tables[,2,i,j]))[1:101]

	for(i in 1:(nrow(F_life_tables)-1)) {
			F_life_tables[i,7,,] = F_life_tables[i,6,,] - F_life_tables[i+1,6,,]
			M_life_tables[i,7,,] = M_life_tables[i,6,,] - M_life_tables[i+1,6,,]
	}
			
	F_life_tables[nrow(F_life_tables),7,,] = F_life_tables[nrow(F_life_tables),6,,]
	M_life_tables[nrow(M_life_tables),7,,] = M_life_tables[nrow(M_life_tables),6,,]
		
	# Age groups
	F_life_tables[,14,,] = c(0, rep(1,4), rep(seq(5, 95, by = 5), each = 5), 100)
	M_life_tables[,14,,] = c(0, rep(1,4), rep(seq(5, 95, by = 5), each = 5), 100)
			
	# nAx
	for(i in 2:(nrow(F_life_tables)-1)) {
		F_life_tables[i,3,,] = ((-1/24)*F_life_tables[i-1,7,,] + 0.5*F_life_tables[i,7,,] + (1/24)*F_life_tables[i+1,7,,])/F_life_tables[i,7,,]
		M_life_tables[i,3,,] = ((-1/24)*M_life_tables[i-1,7,,] + 0.5*M_life_tables[i,7,,] + (1/24)*M_life_tables[i+1,7,,])/M_life_tables[i,7,,]
	}
			
	# Intervals
		F_life_tables[,4,,] = 1; M_life_tables[,4,,] = 1
		
	# m-rates
	  F_life_tables[,5,,] = F_life_tables[,2,,] / (F_life_tables[,4,,]-F_life_tables[,2,,]*(F_life_tables[,4,,]-F_life_tables[,3,,]))
	  M_life_tables[,5,,] = M_life_tables[,2,,] / (M_life_tables[,4,,]-M_life_tables[,2,,]*(M_life_tables[,4,,]-M_life_tables[,3,,]))
		
			F_life_tables[,5,,][F_life_tables[,5,,] > 1]  = 1; F_life_tables[101,5,,] = 1 
			M_life_tables[,5,,][M_life_tables[,5,,] > 1]  = 1; M_life_tables[101,5,,] = 1
			
	# Lx
			for(i in 1:(nrow(F_life_tables)-1)) {
			F_life_tables[i,8,,] = F_life_tables[i,4,,]*F_life_tables[i+1,6,,] + F_life_tables[i,3,,]*F_life_tables[i,7,,]
			M_life_tables[i,8,,] = M_life_tables[i,4,,]*M_life_tables[i+1,6,,] + M_life_tables[i,3,,]*M_life_tables[i,7,,]
			}
		
				F_life_tables[nrow(F_life_tables),8,,] = F_life_tables[nrow(F_life_tables),3,,]*F_life_tables[nrow(F_life_tables),7,,]
				M_life_tables[nrow(M_life_tables),8,,] = M_life_tables[nrow(M_life_tables),3,,]*M_life_tables[nrow(M_life_tables),7,,]
		
	#  Tx
	for(i in 1:dim(F_life_tables)[3]) {
		for(j in 1:dim(F_life_tables)[4]) {
			F_life_tables[1,9,i,j] = sum(F_life_tables[,8,i,j]) 
			M_life_tables[1,9,i,j] = sum(M_life_tables[,8,i,j]) 
			}}
			
	for(i in 1:dim(F_life_tables)[3]) {
		for(j in 1:dim(F_life_tables)[4]) {
			F_life_tables[2,9,i,j] = sum(F_life_tables[,8,i,j]) - F_life_tables[1,8,i,j]
			M_life_tables[2,9,i,j] = sum(M_life_tables[,8,i,j]) - M_life_tables[1,8,i,j]
			}}
		
	for(i in 3:(nrow(F_life_tables))) {
		for(j in 1:dim(F_life_tables)[4]) {
			F_life_tables[i,9,,j] = apply(F_life_tables[,8,,j], c(2), sum) - apply(F_life_tables[1:(i-1),8,,j], c(2), sum)			
			M_life_tables[i,9,,j] = apply(M_life_tables[,8,,j], c(2), sum) - apply(M_life_tables[1:(i-1),8,,j], c(2), sum)
			}}
					
	# ex
			F_life_tables[,10,,] = F_life_tables[,9,,]/F_life_tables[,6,,]
			M_life_tables[,10,,] = M_life_tables[,9,,]/M_life_tables[,6,,]
			
			round(range(M_life_tables[1,10,1,]),1); round(mean(M_life_tables[1,10,1,]))
			round(range(M_life_tables[16,10,1,])); round(mean(M_life_tables[16,10,1,]), 1)

			
#____________________________________________________________________________________________________________________________________________________________________________
# ASFRs
	# Estimation de F/TF = exp(-exp(-Y)) avec Y = alpha + beta*Ys
	for(i in 1:dim(F_life_tables)[3]) {
		F_life_tables[12:50,15,i,]  =  exp(-exp(-(fert_param[i,1] + fert_param[i,2] * Boothstandard[,2])))
	}

		F_life_tables[12:49,15,,]  =  F_life_tables[13:50,15,,1] - F_life_tables[12:49,15,,1]
		F_life_tables[50,15,,] = 1 - F_life_tables[50,15,,]
	
	# Mean age at childbearing
	 for(i in 1:n_fert) {		
		 fert_param$MAC[i] = weighted.mean(F_life_tables[,13,i,1]+0.5, F_life_tables[,15,i,1])
	 }
		
	
# Comparison of estimation approaches  
#___________________________________________________________________________       
# 1. Rates calculated directly for the recent period prior to the survey
         SSH = all_param
         # 35q15
         SSH$q35_15input = NA
         SSH$q35_15sim = NA
         SSH$q35_15SSH = NA
         SSH$seq35_15sim = NA
         SSH$seq35_15SSH = NA
         # 30q15
         SSH$q30_15input = NA
         SSH$q30_15sim = NA
         SSH$q30_15SSH = NA
         SSH$seq30_15sim = NA
         SSH$seq30_15SSH = NA
         # nb of respondents and deaths
         SSH$nresp = NA
         SSH$ndeaths35_15sim = NA; SSH$ndeaths30_15sim = NA
         SSH$ndeaths35_15SSH = NA; SSH$ndeaths30_15SSH = NA
         # Age specific nqx
         popnqx = as.list(1)
         
         for (isim in 1:n_sim) {
           # isim = 1
           
           cat(isim); cat("\n")
           # probabilities 35q15 and 30q15 in the life tables
           SSH$q35_15input[isim] =  1-(F_life_tables['50',"lx",all_param$n_fert[isim],all_param$n_mort[isim]]/F_life_tables['15',"lx",all_param$n_fert[isim],all_param$n_mort[isim]])
           SSH$q30_15input[isim] =  1-(F_life_tables['45',"lx",all_param$n_fert[isim],all_param$n_mort[isim]]/F_life_tables['15',"lx",all_param$n_fert[isim],all_param$n_mort[isim]])
           
           if(SSH$r[isim] == 0.01){
           load(file = paste("Results_TZA01\\", "all_",  paste("sim", SSH$isim1[isim], sep = "_"),   sep = ""))
           }
           if(SSH$r[isim] == 0.03){
             load(file = paste("Results_TZA03\\", "all_",  paste("sim", SSH$isim3[isim], sep = "_"),   sep = ""))
           }
           start_sim = max(all$month_dob[all$mother == 0]+1, na.rm = T); end_sim = max(c(all$month_dod, all$month_dob),na.rm = TRUE) + 1; length_sim = end_sim - start_sim
           w = all
           
           bio.mort = w[,-c(grep("month|jour|dod_m|dob_m|bornfromHIV|intrapartum", colnames(w)))]
           bio.mort$fu.time = (bio.mort$dod - bio.mort$dob)
           bio.mort$fu.time[is.na(bio.mort$fu.time)] = (end_sim - bio.mort$dob[is.na(bio.mort$fu.time)])
           bio.mort$fu.time[!is.na(bio.mort$dod) &  (bio.mort$dob < start_sim)] = (bio.mort$dod[!is.na(bio.mort$dod) &  (bio.mort$dob < start_sim)] - start_sim)
           stopifnot(nrow(bio.mort[is.na(bio.mort$fu.time),]) == 0)
           bio.mort$event = as.numeric(!is.na(bio.mort$dod))
           bio.mort_Surv = with(bio.mort, Surv(time = bio.mort$fu.time, event = bio.mort$event, type = "right"))
           bio.mort$startage = NA
           bio.mort$startage[which(bio.mort$dob >= start_sim)] = 0 
           bio.mort$startage[which(bio.mort$dob < start_sim)] = (start_sim - bio.mort$dob[which(bio.mort$dob < start_sim)])/12
           stopifnot(nrow(bio.mort[is.na(bio.mort$startage),]) == 0)
           age <- tcut(bio.mort$startage*12, age_ref, labels=as.character(age_ref[1:(length(age_ref)-1)]))
           time.cut = seq(start_sim, end_sim + 12, by = 12)
           time.cut.labels = as.character(time.cut[1:(length(time.cut)-1)])
           annee <- tcut(bio.mort$dob, time.cut, labels=time.cut.labels)
           annee[bio.mort$dob < start_sim] <- tcut(start_sim, time.cut, labels=time.cut.labels)
           bio.mort_pp <- pyears(bio.mort_Surv ~ annee + age + sex , bio.mort, scale = 12, data.frame=TRUE)
           stopifnot(round(sum(bio.mort$fu.time)/12 - sum(bio.mort_pp$data$pyears), 4) == 0)
           mort_pp <- bio.mort_pp$data; stopifnot(nrow(mort_pp[mort_pp$event > 1 & mort_pp$pyears == 0,]) == 0)
           mort_pp$age <- as.numeric(as.character(mort_pp$age))/12
           mort_pp$age[mort_pp$age > 85] = 85
           mort_pp$agegrp = trunc(mort_pp$age/5)*5
           mort_pp$agegrp[mort_pp$age >= 1 & mort_pp$age < 5 ] = 1
           mort_pp$annee<- (as.numeric(as.character(mort_pp$annee))-start_sim)/12 + 1870
           mort_pp$sim.id = isim
           # last 5 completed years
           mort_pp = mort_pp[mort_pp$annee >= 2015,]
           
           deaths_F = tapply(mort_pp$event, list(mort_pp$agegrp, mort_pp$sex), sum)[,2]
           exp_F = tapply(mort_pp$pyears, list(mort_pp$agegrp, mort_pp$sex), sum)[,2]
           taux_F = deaths_F/exp_F
           q_F = 5*taux_F/(1+(5-2.5)*taux_F) 
           q_F[1] = 1*taux_F[1]/(1+(1-0.35)*taux_F[1])
           q_F[2] = 4*taux_F[2]/(1+(4-1.361)*taux_F[2])
           
           popnqx[[isim]] = as.data.frame(cbind(q_F[names(q_F) %in% 10:49], seq(10, 45, 5), rep(isim, 8)))
           names(popnqx[[isim]]) = c('nqx', 'x', 'isim')                         
                 
           # probability re-estimated from the whole population
           SSH$q35_15sim[isim] = 1-prod(1-q_F[names(q_F) %in% 15:49])
           SSH$q30_15sim[isim] = 1-prod(1-q_F[names(q_F) %in% 15:44])
           SSH$ndeaths35_15sim[isim] = sum(deaths_F[names(deaths_F) %in% 15:49])
           SSH$ndeaths30_15sim[isim] = sum(deaths_F[names(deaths_F) %in% 15:44])
             
           plot(as.numeric(names(q_F)), q_F, log = "y", type = 'l', col = 1, lty = 1, lwd = 1, ylim = c(min(q_F)/2, max(q_F)*2), xlab = 'x', ylab = "nqx")
           nqxsim = c(1-F_life_tables[,"lx",all_param$n_fert[isim],all_param$n_mort[isim]][2]/F_life_tables[,"lx",all_param$n_fert[isim],all_param$n_mort[isim]][1],
                      1-F_life_tables[,"lx",all_param$n_fert[isim],all_param$n_mort[isim]][6]/F_life_tables[,"lx",all_param$n_fert[isim],all_param$n_mort[isim]][2],
                      1-F_life_tables[,"lx",all_param$n_fert[isim],all_param$n_mort[isim]][seq(11, 101, 5)]/F_life_tables[,"lx",all_param$n_fert[isim],all_param$n_mort[isim]][seq(6, 96, 5)])
           points(c(0, 1, seq(5, 95, 5)), nqxsim, type = 'l', col = 1, lwd = 2)
           points(as.numeric(names(q_F)), q_F+1.96*q_F*sqrt((1/deaths_F)*(1-q_F)),  type = 'l', col = 1, lty = 2)
           points(as.numeric(names(q_F)), q_F-1.96*q_F*sqrt((1/deaths_F)*(1-q_F)),  type = 'l', col = 1, lty = 2)                   
           
           # S.E around summary indices: 35q15 and 30q15
           chiang = as.data.frame(cbind(as.numeric(names(deaths_F)), deaths_F, q_F))
           names(chiang) = c('x', 'deaths_F', 'q_F')
           chiang = chiang[chiang$x %in% 15:50,]; stopifnot(chiang$deaths_F > 0)
           chiang$S2q = (1/chiang$deaths_F)*(1-chiang$q_F)*(chiang$q_F^2)
           chiang$np15 = c(1, rep(NA, nrow(chiang)-1))
           for(i in 2:nrow(chiang)){chiang$np15[i]= prod(1-chiang$q_F[1:(i-1)])}
           chiang$sump15S15 = c(0, rep(NA, nrow(chiang)-1))
           for(i in 2:nrow(chiang)){chiang$sump15S15[i]= (1-chiang$q_F[i-1])^(-2)*chiang$S2q[i-1]+chiang$sump15S15[i-1]}
           chiang$S2q15 = c(0, rep(NA, nrow(chiang)-1))
           for(i in 2:nrow(chiang)){chiang$S2q15[i]= chiang$sump15S15[i]*chiang$np15[i]^2}
           chiang$SEq15 = c(0, rep(NA, nrow(chiang)-1))
           for(i in 2:nrow(chiang)){chiang$SEq15[i]= sqrt(chiang$S2q15[i])}
           SSH$seq35_15sim[isim] = chiang[nrow(chiang),]$SEq15
           SSH$seq30_15sim[isim] = chiang[nrow(chiang)-1,]$SEq15
           
           # rates calculated from SSH
           if(SSH$r[isim] == 0.01){
             load(file = paste("Results_TZA01\\", "all_",  paste("sim", SSH$isim1[isim], sep = "_"),   sep = ""))
           }
           if(SSH$r[isim] == 0.03){
             load(file = paste("Results_TZA03\\", "all_",  paste("sim", SSH$isim3[isim], sep = "_"),   sep = ""))
           }
           ego = all[is.na(all$dod) & all$sex == 1,grep('^pid|mother|^dob$', names(all))] 
           end_sim = max(all$month_dob)+1
           start_sim = min(all$month_dob)
           names(ego) = sub("pid", "ego", names(ego)); names(ego) = sub("dob", "dob_ego", names(ego))
           #  All those with a mother identified through a surviving sister 
           sib = all[c(all$mother %in% ego$mother),grep('pid|mother|^dod$|^dob$|^sex$', names(all))] 
           sib = sib[!duplicated(sib$pid),]; names(sib) = sub("pid", "alter", names(sib)); names(sib) = sub("dod", "dod_alter", names(sib)); names(sib) = sub("dob", "dob_alter", names(sib)); names(sib) = sub("sex", "sex_alter", names(sib))
           # merge so that ego is repeated multiple times
           sib = merge(ego, sib, by = c("mother"), all.x = T, all.y = F)
           sib$mother_alter = sib$mother; sib$type = 1
           sib$type[sib$ego == sib$alter] = 0
           sib = sib[, order(names(sib))]
           # remove respondent and keep when respondent is aged 15-49
           data= sib[sib$type == 1 & (end_sim - sib$dob_ego)/12 >= 15 & (end_sim - sib$dob_ego)/12 < 49,]
           SSH$nresp[isim] = length(unique(data$ego))
           data$v005 = 1
           data$v008 = end_sim
           data$death = as.numeric(!is.na(data$dod_alter))
           data$sex_alter = as.factor(data$sex_alter)
           data$tstop <- ifelse(data$death, data$dod_alter, data$v008)
           data$dob <- data$dob_alter
           data$intv <- data$v008
           clusters = ~1
           by = ~sex_alter
           data$weights <- data$v005/mean(data$v005)
           vars <- unique(unlist(lapply(c(by, strata, clusters), all.vars)))
           f <- formula(paste("~", paste(vars, collapse = "+")))
           mf <- model.frame(f, data = data, na.action = na.pass, death = death, 
                             weights = weights, dob = dob, intv = intv, tstop = tstop)
           aggr <- demog_pyears(f, mf, period = NULL, agegr = seq(15,55,5), 
                                tips = c(0,7), event = "(death)", tstart = "(dob)", tstop = "(tstop)", 
                                weights = "(weights)", origin = 1900, scale = 12)$data
           aggr = aggr[aggr$sex_alter == 1,]
           mx = tapply(aggr$event, aggr$agegr, sum)/tapply(aggr$pyears, aggr$agegr, sum)       
           qx = 5*mx/(1+(5-2.5)*mx); names(qx) = seq(15, 50, 5)
           deaths_sib = tapply(aggr$event, aggr$agegr, sum); names(deaths_sib) = seq(15, 50, 5)
           SSH$q35_15SSH[isim] = 1-prod(1-qx[names(qx) %in% 15:49])
           SSH$ndeaths35_15SSH[isim] = sum(deaths_sib[names(deaths_sib) %in% 15:49])
           SSH$ndeaths30_15SSH[isim] = sum(deaths_sib[names(deaths_sib) %in% 15:44])
           SSH$q30_15SSH[isim] = 1-prod(1-qx[names(qx) %in% 15:44])
           
           points(seq(15, 45, 5), qx[names(qx) %in% 15:49], type = 'l', col = 2, lty = 1, lwd = 3)
           points(seq(15, 45, 5), qx[names(qx) %in% 15:49]+1.96*qx[names(qx) %in% 15:49]*sqrt((1/deaths_sib[names(qx) %in% 15:49])*(1-qx[names(qx) %in% 15:49])),  type = 'l', col = 2, lty = 2, lwd = 2)
           points(seq(15, 45, 5), qx[names(qx) %in% 15:49]-1.96*qx[names(qx) %in% 15:49]*sqrt((1/deaths_sib[names(qx) %in% 15:49])*(1-qx[names(qx) %in% 15:49])),  type = 'l', col = 2, lty = 2, lwd = 2)                   
           
           # S.E around summary indices
           chiang = as.data.frame(cbind(as.numeric(names(deaths_sib)), deaths_sib, qx))
           names(chiang) = c('x', 'deaths_F', 'q_F')
           chiang = chiang[chiang$x %in% 15:50,]; stopifnot(chiang$deaths_F > 0)
           chiang$S2q = (1/chiang$deaths_F)*(1-chiang$q_F)*(chiang$q_F^2)
           chiang$np15 = c(1, rep(NA, nrow(chiang)-1))
           for(i in 2:nrow(chiang)){chiang$np15[i]= prod(1-chiang$q_F[1:(i-1)])}
           chiang$sump15S15 = c(0, rep(NA, nrow(chiang)-1))
           for(i in 2:nrow(chiang)){chiang$sump15S15[i]= (1-chiang$q_F[i-1])^(-2)*chiang$S2q[i-1]+chiang$sump15S15[i-1]}
           chiang$S2q15 = c(0, rep(NA, nrow(chiang)-1))
           for(i in 2:nrow(chiang)){chiang$S2q15[i]= chiang$sump15S15[i]*chiang$np15[i]^2}
           chiang$SEq15 = c(0, rep(NA, nrow(chiang)-1))
           for(i in 2:nrow(chiang)){chiang$SEq15[i]= sqrt(chiang$S2q15[i])}
           
           SSH$seq35_15SSH[isim] = chiang[nrow(chiang),]$SEq15
           SSH$seq30_15SSH[isim] = chiang[nrow(chiang)-1,]$SEq15
           
         } # end of loop by sim  
         
#___________________________________________________________________________       
# 2. Original indirect method from age 15, Timaeus et al. 2001 paper
SIM_sib_ind = as.list(1)
SIM_sib_dir = as.list(1)
prop1519 = as.list(1)
         for (isim in 1:n_sim) {
           
           print(isim)
           SIM_sib_dir[[isim]] = matrix(NA, 3, 4)
           
           # direct 6 years
           if(SSH$r[isim] == 0.01){
             load(file = paste("Results_TZA01\\", "all_",  paste("sim", SSH$isim1[isim], sep = "_"),   sep = ""))
           }
           if(SSH$r[isim] == 0.03){
             load(file = paste("Results_TZA03\\", "all_",  paste("sim", SSH$isim3[isim], sep = "_"),   sep = ""))
           }
           ego = all[is.na(all$dod) & all$sex == 1,grep('^pid|mother|^dob$', names(all))] 
           end_sim = max(all$month_dob)+1
           start_sim = min(all$month_dob)
           names(ego) = sub("pid", "ego", names(ego)); names(ego) = sub("dob", "dob_ego", names(ego))
           #  All those with a mother identified through a surviving sister  
           sib = all[c(all$mother %in% ego$mother),grep('pid|mother|^dod$|^dob$|^sex$', names(all))] 
           sib = sib[!duplicated(sib$pid),]; names(sib) = sub("pid", "alter", names(sib)); names(sib) = sub("dod", "dod_alter", names(sib)); names(sib) = sub("dob", "dob_alter", names(sib)); names(sib) = sub("sex", "sex_alter", names(sib))
           # merge so that ego is repeated multiple times
           sib = merge(ego, sib, by = c("mother"), all.x = T, all.y = F)
           sib$mother_alter = sib$mother; sib$type = 1
           sib$type[sib$ego == sib$alter] = 0
           sib = sib[, order(names(sib))]
           # remove respondent and keep when respondent is aged 15-49
           data= sib[sib$type == 1 & (end_sim - sib$dob_ego)/12 >= 15 & (end_sim - sib$dob_ego)/12 < 50,]
           SSH$nresp[isim] = length(unique(data$ego))
           data$v005 = 1
           data$v008 = end_sim
           data$death = as.numeric(!is.na(data$dod_alter))
           data$sex_alter = as.factor(data$sex_alter)
           data$tstop <- ifelse(data$death, data$dod_alter, data$v008)
           data$dob <- data$dob_alter
           data$intv <- data$v008
           clusters = ~1
           by = ~sex_alter
           data$weights <- data$v005/mean(data$v005)
           vars <- unique(unlist(lapply(c(by, strata, clusters), all.vars)))
           f <- formula(paste("~", paste(vars, collapse = "+")))
           mf <- model.frame(f, data = data, na.action = na.pass, death = death, 
                             weights = weights, dob = dob, intv = intv, tstop = tstop)
           aggr <- demog_pyears(f, mf, period = NULL, agegr = seq(15,55,5), 
                                tips = seq(0,18, 6), event = "(death)", tstart = "(dob)", tstop = "(tstop)", 
                                weights = "(weights)", origin = 1900, scale = 12)$data
           aggr = aggr[aggr$sex_alter == 1,]
           for(z in 1:3){
                  thistips = c('0-5', '6-11', '12-17')[z]
                  thisaggr = aggr[aggr$tips == thistips,]
                  mx = tapply(thisaggr$event, thisaggr$agegr, sum)/tapply(thisaggr$pyears, thisaggr$agegr, sum)       
                  qx = 5*mx/(1+(5-2.5)*mx); names(qx) = seq(15, 50, 5)
               deaths_sib = tapply(thisaggr$event, thisaggr$agegr, sum); names(deaths_sib) = seq(15, 50, 5)
               chiang = as.data.frame(cbind(as.numeric(names(deaths_sib)), deaths_sib, qx))
               names(chiang) = c('x', 'deaths_F', 'q_F')
               chiang = chiang[chiang$x %in% 15:50,]; chiang$deaths_F[chiang$deaths_F == 0]=0.0001
               chiang$S2q = (1/chiang$deaths_F)*(1-chiang$q_F)*(chiang$q_F^2)
               chiang$np15 = c(1, rep(NA, nrow(chiang)-1))
               for(i in 2:nrow(chiang)){chiang$np15[i]= prod(1-chiang$q_F[1:(i-1)])}
               chiang$sump15S15 = c(0, rep(NA, nrow(chiang)-1))
               for(i in 2:nrow(chiang)){chiang$sump15S15[i]= (1-chiang$q_F[i-1])^(-2)*chiang$S2q[i-1]+chiang$sump15S15[i-1]}
               chiang$S2q15 = c(0, rep(NA, nrow(chiang)-1))
               for(i in 2:nrow(chiang)){chiang$S2q15[i]= chiang$sump15S15[i]*chiang$np15[i]^2}
               chiang$SEq15 = c(0, rep(NA, nrow(chiang)-1))
               for(i in 2:nrow(chiang)){chiang$SEq15[i]= sqrt(chiang$S2q15[i])}
               SIM_sib_dir[[isim]][z,] = cbind(isim, z, 1-prod(1-qx[names(qx) %in% 15:49]), chiang[nrow(chiang),]$SEq15)
           }
           
           # Indirect
           data= sib[sib$type == 1 & (end_sim - sib$dob_ego)/12 >= 15 & (end_sim - sib$dob_ego)/12 < 49,]
           #  keep only those who survived to age 15 
           data = data[which((is.na(data$dod_alter) & (end_sim - data$dob_alter  ) >= 15*12) | (!is.na(data$dod_alter) & (data$dod_alter - data$dob_alter) >= 15*12)),]
           data$ageego = trunc(trunc(as.numeric(end_sim - data$dob_ego)/12)/5)*5
           data$mm2 = as.numeric(is.na(data$dod_alter))
           data$mm1 = data$sex_alter 
           
           tab_bro = xtabs(~data$ageego + data$mm2 + data$mm1)[,,1]
           tab_bro = as.data.frame(cbind((tab_bro), as.numeric(dimnames(tab_bro)[[1]])))
           colnames(tab_bro)[length(colnames(tab_bro))] = "Age_group"
           tab_bro$Age_group <- factor(tab_bro$Age_group, levels = seq(15, 45, 5), labels = paste(seq(15, 45, by = 5), seq(19, 49, by = 5), sep = "-")) 
           tab_bro$Sex = "Males"
           tab_sis = xtabs(~data$ageego + data$mm2 + data$mm1)[,,2]
           tab_sis = as.data.frame(cbind((tab_sis), as.numeric(dimnames(tab_sis)[[1]])))
           colnames(tab_sis)[length(colnames(tab_sis))] = "Age_group"
           tab_sis$Age_group <- factor(tab_sis$Age_group, levels = seq(15, 45, 5), labels = paste(seq(15, 45, by = 5), seq(19, 49, by = 5), sep = "-")) 
           tab_sis$Sex = "Females"
           tab = rbind(tab_bro, tab_sis)
           tab$v008 = end_sim
           tab$sim.id = isim
           prop1519[[isim]] = prop.table(tab[tab$Age_group == "15-19" & tab$Sex == "Females",c(1,2)])[2]
           
           # Estimates are time-located with Timaeus (2001) coefficients
           temp = merge(tab, location_TZ2001, by = "Age_group", all.x = TRUE)		
           temp$S5n_5 = temp$'1'/(temp$'1' + temp$'0')
           temp$S5n_5_l = temp$S5n_5 - 1.96*sqrt((temp$S5n_5*(1-temp$S5n_5))/(temp[,2]+temp[,3]))
           temp$S5n_5_u = temp$S5n_5 + 1.96*sqrt((temp$S5n_5*(1-temp$S5n_5))/(temp[,2]+temp[,3]))
           
           temp$years_back = temp$c - temp$d * log(temp$S5n_5)
           temp$time = (2020 - temp$years_back)
           
           temp = merge(temp, coeff_TZ2001[, grep("^a$|^b$|Age_group", names(coeff_TZ2001))], by = "Age_group", all.x = TRUE)
           temp$n[temp$Age_group == "15-19"] = 20
           temp = temp[order(temp$n),]
           temp$xp15 = temp$a + temp$b*temp$S5n_5
           temp$xp15_l = temp$a + temp$b*temp$S5n_5_l
           temp$xp15_u = temp$a + temp$b*temp$S5n_5_u
           
           temp$xq15 = 1- temp$xp15
           temp$xq15_l = 1- temp$xp15_u
           temp$xq15_u = 1- temp$xp15_l
           
           SIM_sib_ind[[isim]] = temp
         }

SIM_sib_ind = do.call(rbind, SIM_sib_ind) 
SIM_sib_ind$W30q15 = NA
SIM_sib_ind$W30q15_l = NA
SIM_sib_ind$W30q15_u = NA

SIM_sib_ind$W35q15 = NA
SIM_sib_ind$W35q15_l = NA
SIM_sib_ind$W35q15_u = NA


SIM_sib_ind$x  = SIM_sib_ind$n-15
# remove 15-19
SIM_sib_ind = SIM_sib_ind[-which(SIM_sib_ind$Age_group == "15-19"),]

#30q15
for(i in seq(10, 35, 5)){
  # males - estimate
  SIM_sib_ind$W30q15[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Males"] =
    approx(x = 1-cdmltw(sex = "M")$lx[,as.character(15+i)]/cdmltw(sex = "M")$lx[,'15'],
           y = 1-cdmltw(sex = "M")$lx[,'45']/cdmltw(sex = "M")$lx[,'15'],
           xout = 1-SIM_sib_ind$xp15[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Males"])$y
  SIM_sib_ind$W30q15_l[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Males"] =
    approx(x = 1-cdmltw(sex = "M")$lx[,as.character(15+i)]/cdmltw(sex = "M")$lx[,'15'],
           y = 1-cdmltw(sex = "M")$lx[,'45']/cdmltw(sex = "M")$lx[,'15'],
           xout = 1-SIM_sib_ind$xp15_u[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Males"])$y
  SIM_sib_ind$W30q15_u[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Males"] =
    approx(x = 1-cdmltw(sex = "M")$lx[,as.character(15+i)]/cdmltw(sex = "M")$lx[,'15'],
           y = 1-cdmltw(sex = "M")$lx[,'45']/cdmltw(sex = "M")$lx[,'15'],
           xout = 1-SIM_sib_ind$xp15_l[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Males"])$y
  # females - estimate
  SIM_sib_ind$W30q15[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Females"] =
    approx(x = 1-cdmltw(sex = "F")$lx[,as.character(15+i)]/cdmltw(sex = "F")$lx[,'15'],
           y = 1-cdmltw(sex = "F")$lx[,'45']/cdmltw(sex = "F")$lx[,'15'],
           xout = 1-SIM_sib_ind$xp15[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Females"])$y
  SIM_sib_ind$W30q15_l[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Females"] =
    approx(x = 1-cdmltw(sex = "F")$lx[,as.character(15+i)]/cdmltw(sex = "F")$lx[,'15'],
           y = 1-cdmltw(sex = "F")$lx[,'45']/cdmltw(sex = "F")$lx[,'15'],
           xout = 1-SIM_sib_ind$xp15_u[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Females"])$y
  SIM_sib_ind$W30q15_u[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Females"] =
    approx(x = 1-cdmltw(sex = "F")$lx[,as.character(15+i)]/cdmltw(sex = "F")$lx[,'15'],
           y = 1-cdmltw(sex = "F")$lx[,'45']/cdmltw(sex = "F")$lx[,'15'],
           xout = 1-SIM_sib_ind$xp15_l[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Females"])$y
  
}

#35q15
for(i in seq(10, 35, 5)){
  # males - estimate
  SIM_sib_ind$W35q15[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Males"] =
    approx(x = 1-cdmltw(sex = "M")$lx[,as.character(15+i)]/cdmltw(sex = "M")$lx[,'15'],
           y = 1-cdmltw(sex = "M")$lx[,'50']/cdmltw(sex = "M")$lx[,'15'],
           xout = 1-SIM_sib_ind$xp15[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Males"])$y
  SIM_sib_ind$W35q15_l[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Males"] =
    approx(x = 1-cdmltw(sex = "M")$lx[,as.character(15+i)]/cdmltw(sex = "M")$lx[,'15'],
           y = 1-cdmltw(sex = "M")$lx[,'50']/cdmltw(sex = "M")$lx[,'15'],
           xout = 1-SIM_sib_ind$xp15_u[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Males"])$y
  SIM_sib_ind$W35q15_u[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Males"] =
    approx(x = 1-cdmltw(sex = "M")$lx[,as.character(15+i)]/cdmltw(sex = "M")$lx[,'15'],
           y = 1-cdmltw(sex = "M")$lx[,'50']/cdmltw(sex = "M")$lx[,'15'],
           xout = 1-SIM_sib_ind$xp15_l[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Males"])$y
  # females - estimate
  SIM_sib_ind$W35q15[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Females"] =
    approx(x = 1-cdmltw(sex = "F")$lx[,as.character(15+i)]/cdmltw(sex = "F")$lx[,'15'],
           y = 1-cdmltw(sex = "F")$lx[,'50']/cdmltw(sex = "F")$lx[,'15'],
           xout = 1-SIM_sib_ind$xp15[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Females"])$y
  SIM_sib_ind$W35q15_l[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Females"] =
    approx(x = 1-cdmltw(sex = "F")$lx[,as.character(15+i)]/cdmltw(sex = "F")$lx[,'15'],
           y = 1-cdmltw(sex = "F")$lx[,'50']/cdmltw(sex = "F")$lx[,'15'],
           xout = 1-SIM_sib_ind$xp15_u[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Females"])$y
  SIM_sib_ind$W35q15_u[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Females"] =
    approx(x = 1-cdmltw(sex = "F")$lx[,as.character(15+i)]/cdmltw(sex = "F")$lx[,'15'],
           y = 1-cdmltw(sex = "F")$lx[,'50']/cdmltw(sex = "F")$lx[,'15'],
           xout = 1-SIM_sib_ind$xp15_l[SIM_sib_ind$x == i & SIM_sib_ind$Sex == "Females"])$y
  
}

SIM_sib_ind$inputLT = NA
for(isim in 1:192){
  SIM_sib_ind[SIM_sib_ind$Sex == "Females" & SIM_sib_ind$sim.id == isim,]$inputLT = 
    1- F_life_tables[seq(25, 50, 5)+1,"lx",all_param$n_fert[isim],all_param$n_mort[isim]]/F_life_tables[16,"lx",all_param$n_fert[isim],all_param$n_mort[isim]]
  SIM_sib_ind[SIM_sib_ind$Sex == "Males" & SIM_sib_ind$sim.id == isim,]$inputLT = 
    1- M_life_tables[seq(25, 50, 5)+1,"lx",all_param$n_fert[isim],all_param$n_mort[isim]]/M_life_tables[16,"lx",all_param$n_fert[isim],all_param$n_mort[isim]]
}
SIM_sib_ind$inputLT30q15 = NA
for(isim in 1:192){
  SIM_sib_ind[SIM_sib_ind$Sex == "Females" & SIM_sib_ind$sim.id == isim,]$inputLT30q15 = 
    1- F_life_tables[46,"lx",all_param$n_fert[isim],all_param$n_mort[isim]]/F_life_tables[16,"lx",all_param$n_fert[isim],all_param$n_mort[isim]]
}
SIM_sib_ind$inputLT35q15 = NA
for(isim in 1:192){
  SIM_sib_ind[SIM_sib_ind$Sex == "Females" & SIM_sib_ind$sim.id == isim,]$inputLT35q15 = 
    1- F_life_tables[51,"lx",all_param$n_fert[isim],all_param$n_mort[isim]]/F_life_tables[16,"lx",all_param$n_fert[isim],all_param$n_mort[isim]]
}

#_________________________________________________________________________________________           
# 3. Synthetic cohort method based on the ratio 5pn-5 = 5Sn(t)/5Sn-5(t)           
           
           SIM_new_syn = as.list(NA)
           for (isim in 1:n_sim) {
             print(isim)
             if(SSH$r[isim] == 0.01){
               load(file = paste("Results_TZA01\\", "all_",  paste("sim", SSH$isim1[isim], sep = "_"),   sep = ""))
             }
             if(SSH$r[isim] == 0.03){
               load(file = paste("Results_TZA03\\", "all_",  paste("sim", SSH$isim3[isim], sep = "_"),   sep = ""))
             }
             
             end_sim = max(all$month_dob)+1
             start_sim = min(all$month_dob)
             
             # synthetic cohort
             ego = all[is.na(all$dod) & all$sex == 1,grep('^pid|mother|^dob$', names(all))]
             names(ego) = sub("pid", "ego", names(ego)); names(ego) = sub("dob", "dob_ego", names(ego))
             #  All those with a mother identified through a surviving sister at time t
             sib = all[c(all$mother %in% ego$mother),grep('pid|mother|^dod$|^dob$|^sex$', names(all))]
             sib = sib[!duplicated(sib$pid),]; names(sib) = sub("pid", "alter", names(sib)); names(sib) = sub("dod", "dod_alter", names(sib)); names(sib) = sub("dob", "dob_alter", names(sib)); names(sib) = sub("sex", "sex_alter", names(sib))
             # merge so that ego is repeated multiple times
             sib = merge(ego, sib, by = c("mother"), all.x = T, all.y = F)
             sib$mother_alter = sib$mother; sib$type = 1
             sib$type[sib$ego == sib$alter] = 0
             sib = sib[, order(names(sib))]
             # remove respondent and keep when respondent is aged 15-49 at time t
             data= sib[sib$type == 1 & (end_sim - sib$dob_ego)/12 >= 15 & (end_sim - sib$dob_ego)/12 < 50,]
             data$ageego = trunc(trunc(as.numeric(end_sim - data$dob_ego)/12)/5)*5 # at time t
             data$mm2 = as.numeric(is.na(data$dod_alter)) # at time t
             data$mm2t_5 =  data$mm2
             data$v008 = end_sim
             data$mm2t_5[!is.na(data$dod_alter) & (data$v008 - data$dod_alter) < 5*12] = 1 # resuscitate if deceased in last 5 years
             data$mm1 = data$sex_alter
             # remove those who did not reach age 15 before time t
             data = data[which((is.na(data$dod_alter) & (end_sim - data$dob_alter  ) >= 15*12) | (!is.na(data$dod_alter) & (data$dod_alter - data$dob_alter) >= 15*12)),]
             
             tab_sis = as.data.frame(cbind(prop.table(xtabs(~ data$ageego  + data$mm2t_5 + data$mm1)[,,2], 1)[,2], prop.table(xtabs(~data$ageego  + data$mm2 + data$mm1)[,,2], 1)[,2]))
             names(tab_sis) = c("t_5", "t")
             tab_sis$ratio = tab_sis$t/ tab_sis$t_5
             tab_sis$Sex = "Females"
             tab_sis$Age_group = 1:7
             tab_sis$Age_group <- factor(tab_sis$Age_group, levels = c(1:7), labels = paste(seq(15, 45, by = 5), seq(19, 49, by = 5), sep = "-"))
             tab_sis$n = seq(15, 45, 5) # age group is age group at time t
             tab_sis$sim.id  = isim
             
             SIM_new_syn[[isim]] = tab_sis
           }
           SIM_new_syn = do.call(rbind, SIM_new_syn) 
           
           # Relate to survivorship probabilities
           # Add nqx recalculated in the last 5 years in simulations
           SIM_new_syn$Sim5qn_5 = NA; SIM_new_syn$Sim5qn = NA
           for(isim in 1:192){
             SIM_new_syn[SIM_new_syn$Sex == "Females" & SIM_new_syn$sim.id == isim,]$Sim5qn_5 = 
               popnqx[[isim]]$nqx[1:7]
             SIM_new_syn[SIM_new_syn$Sex == "Females" & SIM_new_syn$sim.id == isim,]$Sim5qn = 
               popnqx[[isim]]$nqx[2:8]
           }
           SIM_new_syn$Sim5pn_5 = 1- SIM_new_syn$Sim5qn_5
           SIM_new_syn$Sim5pn = 1- SIM_new_syn$Sim5qn
           
           # New regression coefficients based on the hypothetical cohort
           coeff_sn_5 = coeff_TZ2001
           coeff_sn_5$n = coeff_sn_5$n-5 
           coeff_sn_5$r = NA
           coeff_sn_5$cv = NA
           for(i in 1:6){
             mod = lm((SIM_new_syn$Sim5pn_5[SIM_new_syn$Sex == "Females" & SIM_new_syn$n == seq(20, 45, 5)[i]])~SIM_new_syn$ratio[SIM_new_syn$Sex == "Females" & SIM_new_syn$n ==  seq(20, 45, 5)[i]])
             coeff_sn_5[i, "a"] = coef(mod)[1]
             coeff_sn_5[i, "b"] = coef(mod)[2]
             coeff_sn_5[i, "r"] = summary(mod)$r.squared
             coeff_sn_5[i, "cv"] = sqrt(mean(((SIM_new_syn$Sim5pn_5[SIM_new_syn$Sex == "Females" & SIM_new_syn$n == seq(20, 45, 5)[i]]) - mod$fitted.values)^2))/mean((SIM_new_syn$Sim5pn_5[SIM_new_syn$Sex == "Females" & SIM_new_syn$n == seq(20, 45, 5)[i]]))
           }
           save(coeff_sn_5, file = "coeff_sn_5.rda")
           xtable(coeff_sn_5, digits = 4)
           
           coeff_sn = coeff_TZ2001[c(1,1:6), ]
           coeff_sn$Age_group[1] = "15-19"
           coeff_sn$n[1] = 20
           coeff_sn$n = coeff_sn$n-5 
           coeff_sn$r = NA
           coeff_sn$cv = NA
           for(i in 1:7){
             mod = lm((SIM_new_syn$Sim5pn[SIM_new_syn$Sex == "Females" & SIM_new_syn$n == seq(15, 45, 5)[i]])~SIM_new_syn$ratio[SIM_new_syn$Sex == "Females" & SIM_new_syn$n ==  seq(15, 45, 5)[i]])
             coeff_sn[i, "a"] = coef(mod)[1]
             coeff_sn[i, "b"] = coef(mod)[2]
             coeff_sn[i, "r"] = summary(mod)$r.squared
             coeff_sn[i, "cv"] = sqrt(mean(((SIM_new_syn$Sim5pn[SIM_new_syn$Sex == "Females" & SIM_new_syn$n == seq(15, 45, 5)[i]]) - mod$fitted.values)^2))/mean((SIM_new_syn$Sim5pn_5[SIM_new_syn$Sex == "Females" & SIM_new_syn$n == seq(15, 45, 5)[i]]))
           }
           save(coeff_sn, file = "coeff_sn.rda")
           xtable(coeff_sn, digits = 4)
           
           
