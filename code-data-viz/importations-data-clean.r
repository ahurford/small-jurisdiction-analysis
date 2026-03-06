# Title: Data cleaning
# Date: March 6, 2026
# Description: Based on older files from Hurford et al. Pandemic modelling... JTB
#    which is a github repository pandemic-COVID-zero.
#    Produces data files for travel-related cases, close contacts of travellers, and
#.   new community cases for Atlantic Canada and the territories
#.   Nunavut is added relative to the initial analysis in Hurford JTB. Nunavut does not
#    report any travel-related cases or close contacts of travelers during the study period.
#======================

library(dplyr)

## PULLING THE DATA FILES
# This is to pull the data a copy of the PHAC data for new cases
PHAC.data <- read.csv('https://raw.githubusercontent.com/ahurford/covid-nl/master/covid19-download.csv')%>%
  select(date,numtoday, prname)

## These datasets are inidividual-level from the COVID-19 Canada Open data working group. They give travel-related
# cases. The are large files.
CCODWG.2020=read.csv('https://raw.githubusercontent.com/ishaberry/Covid19Canada/master/retired_datasets/individual_level/cases_2020.csv', fill=TRUE)
CCODWG.2021=read.csv('https://raw.githubusercontent.com/ishaberry/Covid19Canada/master/retired_datasets/individual_level/cases_2021_1.csv', fill = TRUE)
CCODWG.2021b=read.csv('https://raw.githubusercontent.com/ishaberry//Covid19Canada/master/retired_datasets/individual_level/cases_2021_2.csv', fill = TRUE)
date_report <- data.frame(date_report=as.Date(seq(from=as.Date("2020-06-15"), to =as.Date("2021-05-31"), by = "day")))


## TRAVEL-RELATED CASES FROM CCODWG
# Correcting some syntax inconsistencies
CCODWG.2020$locally_acquired[CCODWG.2020$locally_acquired =="Close contact"] = "Close Contact"
CCODWG.2020$locally_acquired[CCODWG.2020$locally_acquired =="close contact"] = "Close Contact"
CCODWG.2021$locally_acquired[CCODWG.2021$locally_acquired =="Close contact"] = "Close Contact"
CCODWG.2021$locally_acquired[CCODWG.2021$locally_acquired =="close contact"] = "Close Contact"
CCODWG.2021$locally_acquired[CCODWG.2021$locally_acquired =="Close Contact "] = "Close Contact"
CCODWG.2021b$locally_acquired[CCODWG.2021b$locally_acquired =="Close contact"] = "Close Contact"
CCODWG.2021b$locally_acquired[CCODWG.2021b$locally_acquired =="close contact"] = "Close Contact"
CCODWG.2020$travel_history_country[CCODWG.2020$travel_history_country =="Not Reported "] = "Not Reported"
CCODWG.2021$travel_history_country[CCODWG.2021$travel_history_country =="Not Repoted"] = "Not Reported"
CCODWG.2021$travel_history_country[CCODWG.2021$travel_history_country =="Not reported"] = "Not Reported"
i = which(CCODWG.2021$travel_history_country =="Close contact")
CCODWG.2021 = CCODWG.2021[-i,]

importations=function(province){
  travel.data.2020 <- CCODWG.2020[CCODWG.2020$province==province,]
  travel.data.2021 <- CCODWG.2021[CCODWG.2021$province==province,]
  travel.data.2021b <- CCODWG.2021b[CCODWG.2021b$province==province,]
  travel.data = rbind(travel.data.2020,travel.data.2021, travel.data.2021b)
  travel.data$date_report=format(as.Date(travel.data$date_report, format = "%d-%m-%Y"),"%Y-%m-%d")
  # Only travel-related
  travel = travel.data[travel.data$travel_yn==1 & travel.data$locally_acquired!="Close Contact",]
  travel = dplyr::select(travel,date_report,travel_yn,locally_acquired)%>%
    group_by(date_report)%>%
    add_count()%>%
    distinct()%>%
    rename(travel = n)%>%
    dplyr::select(date_report, travel)%>%
    as.data.frame()

  contacts = travel.data[travel.data$travel_yn!=1 & travel.data$locally_acquired=="Close Contact",]
  contacts = dplyr::select(contacts,date_report,travel_yn,locally_acquired)%>%
    group_by(date_report)%>%
    add_count()%>%
    distinct()%>%
    rename(contacts = n)%>%
    dplyr::select(date_report, contacts)%>%
    as.data.frame()
  
  if(province=="NL"){
    province="Newfoundland and Labrador"
  }
  if(province=="PEI"){
    province="Prince Edward Island"
  }
  if(province=="NWT"){
    province="Northwest Territories"
  }
  data = filter(PHAC.data, prname == province)%>%
    select(date,numtoday)%>%
    rename(date_report = date)%>%
    as.data.frame()
  data$date_report = as.Date(data$date_report)
  travel$date_report=as.Date(travel$date_report)
  contacts$date_report=as.Date(contacts$date_report)
  
  data=left_join(data,travel)%>%left_join(contacts)
  data[is.na(data)]=0
  return(data)
}

NL.travel = importations("NL")%>%
  rename("NL_import"=travel)%>%
  rename("NL_contacts"=contacts)%>%
  rename("NL_new_cases"=numtoday)
NS.travel = importations("Nova Scotia")%>%
  rename("NS_import"=travel)%>%
  rename("NS_contacts"=contacts)%>%
  rename("NS_new_cases"=numtoday)
YT.travel = importations("Yukon")%>%
  rename("YT_import"=travel)%>%
  rename("YT_contacts"=contacts)%>%
  rename("YT_new_cases"=numtoday)
NB.travel = importations("New Brunswick")%>%
  rename("NB_import"=travel)%>%
  rename("NB_contacts"=contacts)%>%
  rename("NB_new_cases"=numtoday)
PEI.travel = importations("PEI")%>%
  rename("PEI_import"=travel)%>%
  rename("PEI_contacts"=contacts)%>%
  rename("PEI_new_cases"=numtoday)
NWT.travel = importations("NWT")%>%
  rename("NWT_import"=travel)%>%
  rename("NWT_contacts"=contacts)%>%
  rename("NWT_new_cases"=numtoday)
NU.travel = importations("Nunavut")%>%
  rename("NU_import"=travel)%>%
  rename("NU_contacts"=contacts)%>%
  rename("NU_new_cases"=numtoday)

travel.day = full_join(NL.travel,NS.travel)%>%
  full_join(NB.travel)%>%
  full_join(PEI.travel)%>%
  full_join(YT.travel)%>%
  full_join(NWT.travel)%>%
  full_join(NU.travel)%>%
  arrange(date_report)%>%
  filter(date_report>"2020-07-01"&date_report<"2021-06-01")%>%
  as.data.frame()


# (1) Data to make the graph of time series of travel-related cases, close contacts, and new cases, aggregated by day
# for Atlantic Canada and the territories
write.csv(travel.day, "data/travel_day.csv")

