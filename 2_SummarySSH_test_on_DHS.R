rm(list = ls()) 
set.seed(3)

library(rdhs)
library(demogsurv)
library(gdata)
library(demogR)
library(xtable)
library(survival)
library(survey)
options(survey.lonely.psu="remove")
library(demogR)
library(stargazer)
library(RColorBrewer)
Sys.setenv(rdhs_DATA_TABLE = "TRUE")

# Identify surveys with SSH
mm_surv <- rbind(dhs_surveys(surveyCharacteristicIds = 1),
                 dhs_surveys(surveyIds = "BF2003DHS")) # BF 2003 has SSH
stopifnot("DR2007SPE" %in% mm_surv$SurveyId  == FALSE)
# Remove surveys not in the public domain
mm_surv = mm_surv[mm_surv$CountryName != "Yemen" | mm_surv$SurveyYear != "2013",]
mm_surv = mm_surv[mm_surv$CountryName != "Eritrea" | mm_surv$SurveyYear != "1995",]
mm_surv = mm_surv[mm_surv$CountryName != "Mauritania" | mm_surv$SurveyYear != "2000",]
mm_surv = mm_surv[mm_surv$CountryName != "Yemen" | mm_surv$SurveyYear != "1997",]
# India 1993 or 1999 have no SSH
mm_surv = mm_surv[mm_surv$CountryName != "India" | mm_surv$SurveyYear != 1993,]
mm_surv = mm_surv[mm_surv$CountryName != "India" | mm_surv$SurveyYear != 1999,]
# Peru Peru 2009/2010/2012 have no DOD, DOB
mm_surv = mm_surv[mm_surv$CountryName != "Peru" | mm_surv$SurveyYear < 2009,]
# Some surveys are not in the correct format
mm_surv = mm_surv[mm_surv$CountryName != "Bolivia" | mm_surv$SurveyYear != 1989,]
mm_surv = mm_surv[mm_surv$CountryName != "Egypt" | mm_surv$SurveyYear != 1988,]
mm_surv = mm_surv[mm_surv$CountryName != "Dominican Republic" | mm_surv$SurveyYear != 1996,]
mm_surv = mm_surv[mm_surv$CountryName != "Sudan" | mm_surv$SurveyYear != 1990,]
mm_surv = mm_surv[mm_surv$CountryName != "Pakistan" | mm_surv$SurveyYear != 2006,]
mm_surv = mm_surv[mm_surv$CountryName != "Ghana" | mm_surv$SurveyYear != 1993,]
mm_surv = mm_surv[mm_surv$CountryName != "India" | mm_surv$SurveyYear != "1999",]
mm_surv = mm_surv[mm_surv$CountryName != "India" | mm_surv$SurveyYear != "2006",]

ird <- dhs_datasets(fileType = "IR", fileFormat = "flat")
ird = ird[ird$SurveyId %in% mm_surv$SurveyId,]
ird$path <- unlist(get_datasets(ird$FileName))

mm_surv = mm_surv[!duplicated(mm_surv),]
stopifnot(nrow(mm_surv) == nrow(ird))

ird = ird[order(ird$SurveyId),]
mm_surv = mm_surv[order(mm_surv$SurveyId),]

nrow(mm_surv)
length(unique(mm_surv$CountryName))

# read datasets and keep the variables needed 
ir <- as.list(rep(NA, length(ird$SurveyId))); names(ir) = ird$SurveyId
for(survid in ird$SurveyId){
  print(survid)
  dat <- readRDS(ird[ird$SurveyId == survid,]$path)
  dat <- dat[grep("caseid|^v0|^v1|^b|^mm", names(dat))]
  dat$SurveyId = as.factor(survid)
  dat$CountryName = as.factor(mm_surv[mm_surv$SurveyId == survid,]$CountryName)
  dat$SurveyYear = as.factor(mm_surv[mm_surv$SurveyId == survid,]$SurveyYear)
  dat$FieldworkStart = as.factor(mm_surv[mm_surv$SurveyId == survid,]$FieldworkStart)
  dat$FieldworkEnd = as.factor(mm_surv[mm_surv$SurveyId == survid,]$FieldworkEnd)
  ir[[which(ird$SurveyId==survid)]] <- dat
}

# reshape datasets for SSH
sib = as.list(1:length(ird$SurveyId)); names(sib) = ird$SurveyId
for(survid in ird$SurveyId){
  print(survid)
  data = ir[[survid]]
  for(k in 1:ncol(data)){
    if(class(data[,k])[1] == "haven_labelled"){class(data[,k]) = "labelled"}
  }
  sib[[which(ird$SurveyId==survid)]] = reshape_sib_data(data,
                                                        widevars = c("SurveyId", "CountryName", "SurveyYear", "FieldworkStart", "FieldworkEnd", "v005", "v008", "v021", "v011", "v013", "v023", 'v024', 'v025'))
}

# Direct mortality estimation of 35q15 and 30q15: trends and last 5 years
q30_15 = as.list(1); q35_15 = as.list(1); q30_15_tips4 = as.list(1); q35_15_tips4 = as.list(1)
# age-specific mortality, last 5 years
q5_15_tips4 = as.list(1);q5_20_tips4 = as.list(1);q5_25_tips4 = as.list(1);q5_30_tips4 = as.list(1); q5_35_tips4 = as.list(1); q5_40_tips4 = as.list(1); q5_45_tips4 = as.list(1)

for(survid in ird$SurveyId){
  set.seed(which(ird$SurveyId == survid))
  print(survid)
  sib[[survid]] = sib[[survid]][sib[[survid]]$mm2 %in% c(1,0) & sib[[survid]]$mm1 %in% c(1,2),]
  sib[[survid]]$death <- sib[[survid]]$mm2 == 0
  #sib[[survid]]$death <- !is.na(sib[[survid]]$mm9) & sib[[survid]]$mm9%in% 2:6
  sib[[survid]]$mm8 <- sib[[survid]]$mm8 + 0.5
  sib[[survid]]$mm1 = as.factor(sib[[survid]]$mm1)
  nostrata = FALSE
  
  if(survid %in% c("AF2015DHS", "PE2004DHS", "PE2007DHS", "BD1994DHS", "BD2007DHS", "CO2015DHS", "IA1993DHS", "UA2007DHS", "BR1991DHS")){ nostrata = TRUE}
  if(!nostrata){
    q30_15[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jkn',
                                agegr=seq(15, 45, 5), tips=seq(0, 21,7), dob="mm4", dod="mm8")
    q35_15[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jkn',
                                agegr=seq(15, 50, 5), tips=seq(0, 21,7), dob="mm4", dod="mm8")
    q35_15_tips4[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jkn',
                                      agegr=seq(15, 50, 5), tips=c(0, 5), dob="mm4", dod="mm8")
    q30_15_tips4[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jkn',
                                      agegr=seq(15, 45, 5), tips=c(0, 5), dob="mm4", dod="mm8")
    # by age
    q5_15_tips4[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jkn',
                                     agegr=seq(15, 20, 5), tips=c(0, 5), dob="mm4", dod="mm8")
    q5_20_tips4[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jkn',
                                     agegr=seq(20, 25, 5), tips=c(0, 5), dob="mm4", dod="mm8")
    q5_25_tips4[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jkn',
                                     agegr=seq(25, 30, 5), tips=c(0, 5), dob="mm4", dod="mm8")
    q5_30_tips4[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jkn',
                                     agegr=seq(30, 35, 5), tips=c(0, 5), dob="mm4", dod="mm8")
    q5_35_tips4[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jkn',
                                     agegr=seq(35, 40, 5), tips=c(0, 5), dob="mm4", dod="mm8")
    q5_40_tips4[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jkn',
                                     agegr=seq(40, 45, 5), tips=c(0, 5), dob="mm4", dod="mm8")
    q5_45_tips4[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jkn',
                                     agegr=seq(45, 50, 5), tips=c(0, 5), dob="mm4", dod="mm8")
  }
  
  if(nostrata){
    q30_15[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jk1', strata = NULL,
                                agegr=seq(15, 45, 5), tips=seq(0, 21,7), dob="mm4", dod="mm8")
    q35_15[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jk1', strata = NULL,
                                agegr=seq(15, 50, 5), tips=seq(0, 21,7), dob="mm4", dod="mm8")
    q35_15_tips4[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jk1', strata = NULL,
                                      agegr=seq(15, 50, 5), tips=c(0, 5), dob="mm4", dod="mm8")
    q30_15_tips4[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jk1', strata = NULL,
                                      agegr=seq(15, 45, 5), tips=c(0, 5), dob="mm4", dod="mm8")
    #by age
    q5_15_tips4[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jk1', strata = NULL,
                                     agegr=seq(15, 20, 5), tips=c(0, 5), dob="mm4", dod="mm8")
    q5_20_tips4[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jk1', strata = NULL,
                                     agegr=seq(20, 25, 5), tips=c(0, 5), dob="mm4", dod="mm8")
    q5_25_tips4[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jk1', strata = NULL,
                                     agegr=seq(25, 30, 5), tips=c(0, 5), dob="mm4", dod="mm8")
    q5_30_tips4[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jk1', strata = NULL,
                                     agegr=seq(30, 35, 5), tips=c(0, 5), dob="mm4", dod="mm8")
    q5_35_tips4[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jk1', strata = NULL,
                                     agegr=seq(35, 40, 5), tips=c(0, 5), dob="mm4", dod="mm8")
    q5_40_tips4[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jk1', strata = NULL,
                                     agegr=seq(40, 45, 5), tips=c(0, 5), dob="mm4", dod="mm8")
    q5_45_tips4[[survid]] = calc_nqx(sib[[survid]], by=~SurveyId+CountryName+SurveyYear+mm1, varmethod = 'jk1', strata = NULL,
                                     agegr=seq(45, 50, 5), tips=c(0, 5), dob="mm4", dod="mm8")
    
  }
  q30_15[[survid]]$tips = as.character(q30_15[[survid]]$tips)
  q30_15[[survid]]$start = as.numeric(sapply(strsplit(as.character(q30_15[[survid]]$tips), "-"), function(x) x[1]))
  q30_15[[survid]]$end = as.numeric(sapply(strsplit(as.character(q30_15[[survid]]$tips), "-"), function(x) x[2]))+1
  q30_15[[survid]]$mid = q30_15[[survid]]$start + (q30_15[[survid]]$end-q30_15[[survid]]$start)/2
  q30_15[[survid]]$time = as.numeric(as.Date(sib[[survid]]$FieldworkEnd[1])-as.Date(sib[[survid]]$FieldworkStart[1]))+as.Date(sib[[survid]]$FieldworkStart[1])-q30_15[[survid]]$mid*365.25
  
  q35_15[[survid]]$tips = as.character(q35_15[[survid]]$tips)
  q35_15[[survid]]$start = as.numeric(sapply(strsplit(as.character(q35_15[[survid]]$tips), "-"), function(x) x[1]))
  q35_15[[survid]]$end = as.numeric(sapply(strsplit(as.character(q35_15[[survid]]$tips), "-"), function(x) x[2]))+1
  q35_15[[survid]]$mid = q35_15[[survid]]$start + (q35_15[[survid]]$end-q35_15[[survid]]$start)/2
  q35_15[[survid]]$time = as.numeric(as.Date(sib[[survid]]$FieldworkEnd[1])-as.Date(sib[[survid]]$FieldworkStart[1]))+as.Date(sib[[survid]]$FieldworkStart[1])-q35_15[[survid]]$mid*365.25
  
  q35_15_tips4[[survid]]$tips = as.character(q35_15_tips4[[survid]]$tips)
  q35_15_tips4[[survid]]$start = as.numeric(sapply(strsplit(as.character(q35_15_tips4[[survid]]$tips), "-"), function(x) x[1]))
  q35_15_tips4[[survid]]$end = as.numeric(sapply(strsplit(as.character(q35_15_tips4[[survid]]$tips), "-"), function(x) x[2]))+1
  q35_15_tips4[[survid]]$mid = q35_15_tips4[[survid]]$start + (q35_15_tips4[[survid]]$end-q35_15_tips4[[survid]]$start)/2
  q35_15_tips4[[survid]]$time = as.numeric(as.Date(sib[[survid]]$FieldworkEnd[1])-as.Date(sib[[survid]]$FieldworkStart[1]))+as.Date(sib[[survid]]$FieldworkStart[1])-q35_15_tips4[[survid]]$mid*365.25
  
  q30_15_tips4[[survid]]$tips = as.character(q30_15_tips4[[survid]]$tips)
  q30_15_tips4[[survid]]$start = as.numeric(sapply(strsplit(as.character(q30_15_tips4[[survid]]$tips), "-"), function(x) x[1]))
  q30_15_tips4[[survid]]$end = as.numeric(sapply(strsplit(as.character(q30_15_tips4[[survid]]$tips), "-"), function(x) x[2]))+1
  q30_15_tips4[[survid]]$mid = q30_15_tips4[[survid]]$start + (q30_15_tips4[[survid]]$end-q30_15_tips4[[survid]]$start)/2
  q30_15_tips4[[survid]]$time = as.numeric(as.Date(sib[[survid]]$FieldworkEnd[1])-as.Date(sib[[survid]]$FieldworkStart[1]))+as.Date(sib[[survid]]$FieldworkStart[1])-q30_15_tips4[[survid]]$mid*365.25
  
  # by age
  q5_15_tips4[[survid]]$tips = as.character(q5_15_tips4[[survid]]$tips)
  q5_15_tips4[[survid]]$start = as.numeric(sapply(strsplit(as.character(q5_15_tips4[[survid]]$tips), "-"), function(x) x[1]))
  q5_15_tips4[[survid]]$end = as.numeric(sapply(strsplit(as.character(q5_15_tips4[[survid]]$tips), "-"), function(x) x[2]))+1
  q5_15_tips4[[survid]]$mid = q5_15_tips4[[survid]]$start + (q5_15_tips4[[survid]]$end-q5_15_tips4[[survid]]$start)/2
  q5_15_tips4[[survid]]$time = as.numeric(as.Date(sib[[survid]]$FieldworkEnd[1])-as.Date(sib[[survid]]$FieldworkStart[1]))+as.Date(sib[[survid]]$FieldworkStart[1])-q5_15_tips4[[survid]]$mid*365.25
  
  q5_20_tips4[[survid]]$tips = as.character(q5_20_tips4[[survid]]$tips)
  q5_20_tips4[[survid]]$start = as.numeric(sapply(strsplit(as.character(q5_20_tips4[[survid]]$tips), "-"), function(x) x[1]))
  q5_20_tips4[[survid]]$end = as.numeric(sapply(strsplit(as.character(q5_20_tips4[[survid]]$tips), "-"), function(x) x[2]))+1
  q5_20_tips4[[survid]]$mid = q5_20_tips4[[survid]]$start + (q5_20_tips4[[survid]]$end-q5_20_tips4[[survid]]$start)/2
  q5_20_tips4[[survid]]$time = as.numeric(as.Date(sib[[survid]]$FieldworkEnd[1])-as.Date(sib[[survid]]$FieldworkStart[1]))+as.Date(sib[[survid]]$FieldworkStart[1])-q5_20_tips4[[survid]]$mid*365.25
  
  q5_25_tips4[[survid]]$tips = as.character(q5_25_tips4[[survid]]$tips)
  q5_25_tips4[[survid]]$start = as.numeric(sapply(strsplit(as.character(q5_25_tips4[[survid]]$tips), "-"), function(x) x[1]))
  q5_25_tips4[[survid]]$end = as.numeric(sapply(strsplit(as.character(q5_25_tips4[[survid]]$tips), "-"), function(x) x[2]))+1
  q5_25_tips4[[survid]]$mid = q5_25_tips4[[survid]]$start + (q5_25_tips4[[survid]]$end-q5_25_tips4[[survid]]$start)/2
  q5_25_tips4[[survid]]$time = as.numeric(as.Date(sib[[survid]]$FieldworkEnd[1])-as.Date(sib[[survid]]$FieldworkStart[1]))+as.Date(sib[[survid]]$FieldworkStart[1])-q5_25_tips4[[survid]]$mid*365.25
  
  q5_30_tips4[[survid]]$tips = as.character(q5_30_tips4[[survid]]$tips)
  q5_30_tips4[[survid]]$start = as.numeric(sapply(strsplit(as.character(q5_30_tips4[[survid]]$tips), "-"), function(x) x[1]))
  q5_30_tips4[[survid]]$end = as.numeric(sapply(strsplit(as.character(q5_30_tips4[[survid]]$tips), "-"), function(x) x[2]))+1
  q5_30_tips4[[survid]]$mid = q5_30_tips4[[survid]]$start + (q5_30_tips4[[survid]]$end-q5_30_tips4[[survid]]$start)/2
  q5_30_tips4[[survid]]$time = as.numeric(as.Date(sib[[survid]]$FieldworkEnd[1])-as.Date(sib[[survid]]$FieldworkStart[1]))+as.Date(sib[[survid]]$FieldworkStart[1])-q5_30_tips4[[survid]]$mid*365.25
  
  q5_35_tips4[[survid]]$tips = as.character(q5_35_tips4[[survid]]$tips)
  q5_35_tips4[[survid]]$start = as.numeric(sapply(strsplit(as.character(q5_35_tips4[[survid]]$tips), "-"), function(x) x[1]))
  q5_35_tips4[[survid]]$end = as.numeric(sapply(strsplit(as.character(q5_35_tips4[[survid]]$tips), "-"), function(x) x[2]))+1
  q5_35_tips4[[survid]]$mid = q5_35_tips4[[survid]]$start + (q5_35_tips4[[survid]]$end-q5_35_tips4[[survid]]$start)/2
  q5_35_tips4[[survid]]$time = as.numeric(as.Date(sib[[survid]]$FieldworkEnd[1])-as.Date(sib[[survid]]$FieldworkStart[1]))+as.Date(sib[[survid]]$FieldworkStart[1])-q5_35_tips4[[survid]]$mid*365.25
  
  q5_40_tips4[[survid]]$tips = as.character(q5_40_tips4[[survid]]$tips)
  q5_40_tips4[[survid]]$start = as.numeric(sapply(strsplit(as.character(q5_40_tips4[[survid]]$tips), "-"), function(x) x[1]))
  q5_40_tips4[[survid]]$end = as.numeric(sapply(strsplit(as.character(q5_40_tips4[[survid]]$tips), "-"), function(x) x[2]))+1
  q5_40_tips4[[survid]]$mid = q5_40_tips4[[survid]]$start + (q5_40_tips4[[survid]]$end-q5_40_tips4[[survid]]$start)/2
  q5_40_tips4[[survid]]$time = as.numeric(as.Date(sib[[survid]]$FieldworkEnd[1])-as.Date(sib[[survid]]$FieldworkStart[1]))+as.Date(sib[[survid]]$FieldworkStart[1])-q5_40_tips4[[survid]]$mid*365.25
  
  q5_45_tips4[[survid]]$tips = as.character(q5_45_tips4[[survid]]$tips)
  q5_45_tips4[[survid]]$start = as.numeric(sapply(strsplit(as.character(q5_45_tips4[[survid]]$tips), "-"), function(x) x[1]))
  q5_45_tips4[[survid]]$end = as.numeric(sapply(strsplit(as.character(q5_45_tips4[[survid]]$tips), "-"), function(x) x[2]))+1
  q5_45_tips4[[survid]]$mid = q5_45_tips4[[survid]]$start + (q5_45_tips4[[survid]]$end-q5_45_tips4[[survid]]$start)/2
  q5_45_tips4[[survid]]$time = as.numeric(as.Date(sib[[survid]]$FieldworkEnd[1])-as.Date(sib[[survid]]$FieldworkStart[1]))+as.Date(sib[[survid]]$FieldworkStart[1])-q5_45_tips4[[survid]]$mid*365.25
  
}

q30_15 = do.call(rbind, q30_15[names(q30_15) != ""])
q30_15$Sex = NA
q30_15$Sex[q30_15$mm1 == 1] = "Males"
q30_15$Sex[q30_15$mm1 == 2] = "Females"

q35_15 = do.call(rbind, q35_15[names(q35_15) != ""])
q35_15$Sex = NA
q35_15$Sex[q35_15$mm1 == 1] = "Males"
q35_15$Sex[q35_15$mm1 == 2] = "Females"

q35_15_tips4 = do.call(rbind, q35_15_tips4[names(q35_15_tips4) != ""])
q35_15_tips4$Sex = NA
q35_15_tips4$Sex[q35_15_tips4$mm1 == 1] = "Males"
q35_15_tips4$Sex[q35_15_tips4$mm1 == 2] = "Females"

q30_15_tips4 = do.call(rbind, q30_15_tips4[names(q30_15_tips4) != ""])
q30_15_tips4$Sex = NA
q30_15_tips4$Sex[q30_15_tips4$mm1 == 1] = "Males"
q30_15_tips4$Sex[q30_15_tips4$mm1 == 2] = "Females"

# by age
q5_15_tips4 = do.call(rbind, q5_15_tips4[names(q5_15_tips4) != ""])
q5_15_tips4$Sex = NA
q5_15_tips4$Sex[q5_15_tips4$mm1 == 1] = "Males"
q5_15_tips4$Sex[q5_15_tips4$mm1 == 2] = "Females"

q5_20_tips4 = do.call(rbind, q5_20_tips4[names(q5_20_tips4) != ""])
q5_20_tips4$Sex = NA
q5_20_tips4$Sex[q5_20_tips4$mm1 == 1] = "Males"
q5_20_tips4$Sex[q5_20_tips4$mm1 == 2] = "Females"

q5_25_tips4 = do.call(rbind, q5_25_tips4[names(q5_25_tips4) != ""])
q5_25_tips4$Sex = NA
q5_25_tips4$Sex[q5_25_tips4$mm1 == 1] = "Males"
q5_25_tips4$Sex[q5_25_tips4$mm1 == 2] = "Females"

q5_30_tips4 = do.call(rbind, q5_30_tips4[names(q5_30_tips4) != ""])
q5_30_tips4$Sex = NA
q5_30_tips4$Sex[q5_30_tips4$mm1 == 1] = "Males"
q5_30_tips4$Sex[q5_30_tips4$mm1 == 2] = "Females"

q5_35_tips4 = do.call(rbind, q5_35_tips4[names(q5_35_tips4) != ""])
q5_35_tips4$Sex = NA
q5_35_tips4$Sex[q5_35_tips4$mm1 == 1] = "Males"
q5_35_tips4$Sex[q5_35_tips4$mm1 == 2] = "Females"

q5_40_tips4 = do.call(rbind, q5_40_tips4[names(q5_40_tips4) != ""])
q5_40_tips4$Sex = NA
q5_40_tips4$Sex[q5_40_tips4$mm1 == 1] = "Males"
q5_40_tips4$Sex[q5_40_tips4$mm1 == 2] = "Females"

q5_45_tips4 = do.call(rbind, q5_45_tips4[names(q5_45_tips4) != ""])
q5_45_tips4$Sex = NA
q5_45_tips4$Sex[q5_45_tips4$mm1 == 1] = "Males"
q5_45_tips4$Sex[q5_45_tips4$mm1 == 2] = "Females"

#______________________________________________________________________________
# Indirect mortality estimation using Timaeus 2001 approach (lifetime data)
coeff_TZ2001 = data.frame(Age_group = c('20-24', '25-29', '30-34', '35-39', '40-44', '45-49'),
                          n = seq(25, 50, 5),
                          a = c(-0.0003, -0.1546, -0.1645,-0.1388, -0.1140,-0.1018),
                          b = c(1.0011,1.1560,1.1660,1.1406,1.1168, 1.1066))
location_TZ2001 = data.frame(Age_group = c('20-24', '25-29', '30-34', '35-39', '40-44', '45-49'),
                             n = seq(25, 50, 5),
                             c = c(3.23,5.46,7.52,9.38,11.00, 12.32),
                             d = c(1.12,1.95,2.78,3.62,4.45,5.28))

ind = as.list(1)
for(survid in ird$SurveyId){
  set.seed(which(ird$SurveyId == survid))
  print(survid)  
  w = sib[[survid]]
  w = w[which((is.na(w$mm8) & (w$v008 - w$mm4) >= 15*12) | (!is.na(w$mm8) & (w$mm8 - w$mm4) >= 15*12)),]
  w = w[!is.na(w$mm2) & w$mm2 < 8,]
  w$v005 = (w$v005/1000000)
  
  # no strata info
  if(length(unique(w$v023)) == 1){w$v023 = w$v024*10+w$v025}
  
  # account for the survey design
  bro = w[w$mm1 == 1,]
  bro$surviving = as.numeric(is.na(bro$mm8))
  DHSdesignbro<-svydesign(id=bro$v021, strata=bro$v023, weights=bro$v005, data=bro, nest = TRUE)
  tab_bro = svyby(~surviving, ~v013, DHSdesignbro, svymean, vartype=c("se","ci"))
  tab_bro = as.data.frame(cbind((tab_bro), as.numeric(dimnames(tab_bro)[[1]])))
  colnames(tab_bro)[length(colnames(tab_bro))] = "Age_group"
  tab_bro$Age_group <- factor(tab_bro$Age_group, levels = c(1:14), labels = paste(seq(15, 80, by = 5), seq(19, 85, by = 5), sep = "-")) 
  tab_bro$Sex = "Males"
  
  sis = w[w$mm1 == 2,]
  sis$surviving = as.numeric(is.na(sis$mm8))
  DHSdesignsis<-svydesign(id=sis$v021, strata=sis$v023, weights=sis$v005, data=sis, nest = TRUE)
  tab_sis = svyby(~surviving, ~v013, DHSdesignsis, svymean, vartype=c("se","ci"))
  tab_sis = as.data.frame(cbind((tab_sis), as.numeric(dimnames(tab_sis)[[1]])))
  colnames(tab_sis)[length(colnames(tab_sis))] = "Age_group"
  tab_sis$Age_group <- factor(tab_sis$Age_group, levels = c(1:14), labels = paste(seq(15, 80, by = 5), seq(19, 85, by = 5), sep = "-")) 
  tab_sis$Sex = "Females"

  tab = rbind(tab_bro, tab_sis)
  tab$v008 = as.numeric(as.Date(w$FieldworkEnd[1])-as.Date(w$FieldworkStart[1]))+as.Date(w$FieldworkStart[1])
  tab$SurveyId = as.character(w$SurveyId[1])
  
  # Estimates are time-located with Timaeus (2001) coefficients
  DHS_sib_ind = merge(tab, location_TZ2001, by = "Age_group")		
  DHS_sib_ind$S5n_5 = DHS_sib_ind$surviving
  
  DHS_sib_ind$years_back = DHS_sib_ind$c - DHS_sib_ind$d * log(DHS_sib_ind$S5n_5)
  DHS_sib_ind$time = (DHS_sib_ind$v008 - DHS_sib_ind$years_back*365.25)
  
  DHS_sib_ind = merge(DHS_sib_ind, coeff_TZ2001[, grep("^a$|^b$|Age_group", names(coeff_TZ2001))], by = "Age_group")
  DHS_sib_ind = DHS_sib_ind[order(DHS_sib_ind$n),]
  DHS_sib_ind$xp15 = DHS_sib_ind$a + DHS_sib_ind$b*DHS_sib_ind$S5n_5
  DHS_sib_ind$xp15_l = DHS_sib_ind$a + DHS_sib_ind$b*DHS_sib_ind$ci_l
  DHS_sib_ind$xp15_u = DHS_sib_ind$a + DHS_sib_ind$b*DHS_sib_ind$ci_u
  
  ind[[survid]] = DHS_sib_ind
} 

ind = do.call(rbind, ind[names(ind) != ""])
ind$x = ind$n-15

# convert to 30q15 using the West model
ind$W30q15 = NA; ind$W30q15_l = NA; ind$W30q15_u = NA
for(i in seq(10, 35, 5)){
  # males - estimate
  ind$W30q15[ind$x == i & ind$Sex == "Males"] =
    approx(x = 1-cdmltw(sex = "M")$lx[,as.character(15+i)]/cdmltw(sex = "M")$lx[,'15'],
           y = 1-cdmltw(sex = "M")$lx[,'45']/cdmltw(sex = "M")$lx[,'15'],
           xout = 1-ind$xp15[ind$x == i & ind$Sex == "Males"])$y
  # males - lower
  ind$W30q15_l[ind$x == i & ind$Sex == "Males"] =
    approx(x = 1-cdmltw(sex = "M")$lx[,as.character(15+i)]/cdmltw(sex = "M")$lx[,'15'],
           y = 1-cdmltw(sex = "M")$lx[,'45']/cdmltw(sex = "M")$lx[,'15'],
           xout = 1-ind$xp15_u[ind$x == i & ind$Sex == "Males"])$y
  # males- upper
  ind$W30q15_u[ind$x == i & ind$Sex == "Males"] =
    approx(x = 1-cdmltw(sex = "M")$lx[,as.character(15+i)]/cdmltw(sex = "M")$lx[,'15'],
           y = 1-cdmltw(sex = "M")$lx[,'45']/cdmltw(sex = "M")$lx[,'15'],
           xout = 1-ind$xp15_l[ind$x == i & ind$Sex == "Males"])$y
  
  # females - estimate
  ind$W30q15[ind$x == i & ind$Sex == "Females"] =
    approx(x = 1-cdmltw(sex = "F")$lx[,as.character(15+i)]/cdmltw(sex = "F")$lx[,'15'],
           y = 1-cdmltw(sex = "F")$lx[,'45']/cdmltw(sex = "F")$lx[,'15'],
           xout = 1-ind$xp15[ind$x == i & ind$Sex == "Females"])$y
  # females - lower
  ind$W30q15_l[ind$x == i & ind$Sex == "Females"] =
    approx(x = 1-cdmltw(sex = "F")$lx[,as.character(15+i)]/cdmltw(sex = "F")$lx[,'15'],
           y = 1-cdmltw(sex = "F")$lx[,'45']/cdmltw(sex = "F")$lx[,'15'],
           xout = 1-ind$xp15_u[ind$x == i & ind$Sex == "Females"])$y
  # females - upper
  ind$W30q15_u[ind$x == i & ind$Sex == "Females"] =
    approx(x = 1-cdmltw(sex = "F")$lx[,as.character(15+i)]/cdmltw(sex = "F")$lx[,'15'],
           y = 1-cdmltw(sex = "F")$lx[,'45']/cdmltw(sex = "F")$lx[,'15'],
           xout = 1-ind$xp15_l[ind$x == i & ind$Sex == "Females"])$y
  
}

# convert to 35q15 using the West model
ind$W35q15 = NA; ind$W35q15_l = NA; ind$W35q15_u = NA
for(i in seq(10, 35, 5)){
  # males - estimate
  ind$W35q15[ind$x == i & ind$Sex == "Males"] =
    approx(x = 1-cdmltw(sex = "M")$lx[,as.character(15+i)]/cdmltw(sex = "M")$lx[,'15'],
           y = 1-cdmltw(sex = "M")$lx[,'50']/cdmltw(sex = "M")$lx[,'15'],
           xout = 1-ind$xp15[ind$x == i & ind$Sex == "Males"])$y
  # males - lower
  ind$W35q15_l[ind$x == i & ind$Sex == "Males"] =
    approx(x = 1-cdmltw(sex = "M")$lx[,as.character(15+i)]/cdmltw(sex = "M")$lx[,'15'],
           y = 1-cdmltw(sex = "M")$lx[,'50']/cdmltw(sex = "M")$lx[,'15'],
           xout = 1-ind$xp15_u[ind$x == i & ind$Sex == "Males"])$y
  # males- upper
  ind$W35q15_u[ind$x == i & ind$Sex == "Males"] =
    approx(x = 1-cdmltw(sex = "M")$lx[,as.character(15+i)]/cdmltw(sex = "M")$lx[,'15'],
           y = 1-cdmltw(sex = "M")$lx[,'50']/cdmltw(sex = "M")$lx[,'15'],
           xout = 1-ind$xp15_l[ind$x == i & ind$Sex == "Males"])$y
  
  # females - estimate
  ind$W35q15[ind$x == i & ind$Sex == "Females"] =
    approx(x = 1-cdmltw(sex = "F")$lx[,as.character(15+i)]/cdmltw(sex = "F")$lx[,'15'],
           y = 1-cdmltw(sex = "F")$lx[,'50']/cdmltw(sex = "F")$lx[,'15'],
           xout = 1-ind$xp15[ind$x == i & ind$Sex == "Females"])$y
  # females - lower
  ind$W35q15_l[ind$x == i & ind$Sex == "Females"] =
    approx(x = 1-cdmltw(sex = "F")$lx[,as.character(15+i)]/cdmltw(sex = "F")$lx[,'15'],
           y = 1-cdmltw(sex = "F")$lx[,'50']/cdmltw(sex = "F")$lx[,'15'],
           xout = 1-ind$xp15_u[ind$x == i & ind$Sex == "Females"])$y
  # females - upper
  ind$W35q15_u[ind$x == i & ind$Sex == "Females"] =
    approx(x = 1-cdmltw(sex = "F")$lx[,as.character(15+i)]/cdmltw(sex = "F")$lx[,'15'],
           y = 1-cdmltw(sex = "F")$lx[,'50']/cdmltw(sex = "F")$lx[,'15'],
           xout = 1-ind$xp15_l[ind$x == i & ind$Sex == "Females"])$y
  
}

#___________________________________________________________________________
# Indirect estimates derived from adult lifetime proportions or changes in proportions surviving calculated from a single survey
library(abind)
load("coeff_sn.rda")

syncohort = as.list(1)
for(survid in ird$SurveyId){
  set.seed(which(ird$SurveyId == survid))
  print(survid)  
  w = sib[[survid]]
  if(survid ==  "CO2015DHS"){w = w[w$v013 %in% 1:7,]}
  w = w[!is.na(w$mm2) & w$mm2 < 8,]
  w$v005 = (w$v005/1000000)
  w$mm2t_5 =  w$mm2
  w$mm2t_5[!is.na(w$mm6) & w$mm6 < 5] = 1
  w$mm2t_5[is.na(w$mm6) & !is.na(w$mm8) & w$v008 - w$mm8 < 5*12] = 1
  
  # select those who survived to age 15
  w = w[which((is.na(w$mm8) & (w$v008 - w$mm4) >= 15*12) | (!is.na(w$mm8) & (w$mm8 - w$mm4) >= 15*12)),]
  
  # Sisters - full sample
  tab_sis = as.data.frame(cbind(prop.table(xtabs(w$v005~w$v013 + w$mm2t_5 + w$mm1)[,,2], 1)[,2], prop.table(xtabs(w$v005~w$v013 + w$mm2 + w$mm1)[,,2], 1)[,2]))
  names(tab_sis) = c("t_5", "t")
  tab_sis$ratio = tab_sis$t/ tab_sis$t_5
  tab_sis$pred5pn = coeff_sn$a+coeff_sn$b*tab_sis$ratio
  tab_sis$pred35q15 = 1-prod(tab_sis$pred5pn, na.rm = TRUE)
  tab_sis$Sex = "Females"
  tab_sis$SurveyId = survid
  tab_sis$n = seq(15,45, 5) # age group at time t
  # Jackknife resampling
  sis= w[w$mm1 == 2,]
  temp1 = xtabs(sis$v005~sis$v013 + sis$mm2t_5 + sis$v021)
  temp2 = xtabs(sis$v005~sis$v013 + sis$mm2 + sis$v021)
  temp = abind(temp1, temp2, along  = 2)
  omitone = as.list(1:length(unique(sis$v021)))
  for(i in 1:length(unique(sis$v021))){
    omitone[[i]] = as.data.frame(apply(temp[,,-i], c(1,2), sum))
    omitone[[i]]$t_5 = omitone[[i]][,2]/(omitone[[i]][,2]+omitone[[i]][,1])
    omitone[[i]]$t = omitone[[i]][,4]/(omitone[[i]][,4]+omitone[[i]][,3])
    omitone[[i]]$ratio = omitone[[i]]$t/omitone[[i]]$t_5
    omitone[[i]]$pred5pn = coeff_sn$a+coeff_sn$b*omitone[[i]]$ratio
    omitone[[i]]$pred35q15 = 1-prod(omitone[[i]]$pred5pn, na.rm = TRUE)
    
    omitone[[i]]$pred5pn.fs = tab_sis$pred5pn
    omitone[[i]]$pred35q15.fs = tab_sis$pred35q15
    
    omitone[[i]]$v021 = i
  }
  omitone = do.call(rbind, omitone)
  omitone$age = 1:7
  k = length(unique(sis$v021))
  
  omitone$pred5pn.ri= k*omitone$pred5pn.fs-(k-1)*omitone$pred5pn
  omitone$pred5pn.sqrdiff= (omitone$pred5pn.ri-omitone$pred5pn.fs)^2 
  tab_sis$pred5pn.se = sqrt((1/(k*(k-1)))*tapply(omitone$pred5pn.sqrdiff, omitone$age, sum))
  
  omitone$pred35q15.ri= k*omitone$pred35q15.fs-(k-1)*omitone$pred35q15
  omitone$pred35q15.sqrdiff= (omitone$pred35q15.ri-omitone$pred35q15.fs)^2 
  tab_sis$pred35q15.se = sqrt((1/(k*(k-1)))*tapply(omitone$pred35q15.sqrdiff, omitone$age, sum))
  
  # Brothers - full sample
  tab_bro = as.data.frame(cbind(prop.table(xtabs(w$v005~w$v013 + w$mm2t_5 + w$mm1)[,,1], 1)[,2], prop.table(xtabs(w$v005~w$v013 + w$mm2 + w$mm1)[,,1], 1)[,2]))
  names(tab_bro) = c("t_5", "t")
  tab_bro$ratio = tab_bro$t/ tab_bro$t_5
  tab_bro$pred5pn = coeff_sn$a+coeff_sn$b*tab_bro$ratio
  tab_bro$pred35q15 = 1-prod(tab_bro$pred5pn, na.rm = TRUE)
  tab_bro$Sex = "Males"
  tab_bro$SurveyId = survid
  tab_bro$n = seq(15,45, 5) # age group at time t
  # Jackknife resampling
  bro= w[w$mm1 == 1,]
  temp1 = xtabs(bro$v005~bro$v013 + bro$mm2t_5 + bro$v021)
  temp2 = xtabs(bro$v005~bro$v013 + bro$mm2 + bro$v021)
  temp = abind(temp1, temp2, along  = 2)
  omitone = as.list(1:length(unique(bro$v021)))
  for(i in 1:length(unique(bro$v021))){
    omitone[[i]] = as.data.frame(apply(temp[,,-i], c(1,2), sum))
    omitone[[i]]$t_5 = omitone[[i]][,2]/(omitone[[i]][,2]+omitone[[i]][,1])
    omitone[[i]]$t = omitone[[i]][,4]/(omitone[[i]][,4]+omitone[[i]][,3])
    omitone[[i]]$ratio = omitone[[i]]$t/omitone[[i]]$t_5
    omitone[[i]]$pred5pn = coeff_sn$a+coeff_sn$b*omitone[[i]]$ratio
    omitone[[i]]$pred35q15 = 1-prod(omitone[[i]]$pred5pn, na.rm = TRUE)
    
    omitone[[i]]$pred5pn.fs = tab_bro$pred5pn
    omitone[[i]]$pred35q15.fs = tab_bro$pred35q15
    
    omitone[[i]]$v021 = i
  }
  omitone = do.call(rbind, omitone)
  omitone$age = 1:7
  k = length(unique(bro$v021))
  
  omitone$pred5pn.ri= k*omitone$pred5pn.fs-(k-1)*omitone$pred5pn
  omitone$pred5pn.sqrdiff= (omitone$pred5pn.ri-omitone$pred5pn.fs)^2 
  tab_bro$pred5pn.se = sqrt((1/(k*(k-1)))*tapply(omitone$pred5pn.sqrdiff, omitone$age, sum))
  
  omitone$pred35q15.ri= k*omitone$pred35q15.fs-(k-1)*omitone$pred35q15
  omitone$pred35q15.sqrdiff= (omitone$pred35q15.ri-omitone$pred35q15.fs)^2 
  tab_bro$pred35q15.se = sqrt((1/(k*(k-1)))*tapply(omitone$pred35q15.sqrdiff, omitone$age, sum))
  
  syncohort[[survid]] = rbind(tab_sis, tab_bro)
  
}
syncohort = do.call(rbind, syncohort[names(syncohort) != ""])
syncohort$time = NA
for(survid in ird$SurveyId){
  syncohort$time[syncohort$SurveyId == survid] = 
    as.Date(q35_15_tips4$time[q35_15_tips4$SurveyId == survid][1])
}

# q35_15 contains the direct estimates for "0-5", "6-11", "12-17" y before the survey
# q35_15_tips4 contains the direct estimates for the last 5 years
# ind contains the estimates  from the indirect method from TZA2001 (varying refs)
# syncohort contains the synthetic cohort built for the last 5 years

