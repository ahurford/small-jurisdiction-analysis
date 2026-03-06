library(ggplot2)
library(scales)
library(dplyr)
library(patchwork)

# Population sizes from
#https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1710000901
#Q1 2021


## Color blind palette
cb <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#999999","#D55E00", "#C779A7")
travel_day <- read.csv('data/travel_day.csv')[,-1]

g.NL = ggplot(data = travel_day, aes(x=as.Date(date_report), group =1))+
  geom_ribbon(aes(ymin = 0, ymax = NL_new_cases - NL_contacts-NL_import), fill="#E69F00")+
  geom_line(aes(y=NL_contacts+NL_import), color = "black", size=.2)+
  scale_x_datetime(labels = date_format("%b"),date_breaks = "1 month")+
  scale_y_log10()+
  ylab("incidence")+
  xlab("")+
  ggtitle("Newfoundland and Labrador (525,895)")+theme_bw()

g.NS = ggplot(data = travel_day, aes(x=as.Date(date_report), group =1))+
  geom_ribbon(aes(ymin = 0, ymax = NS_new_cases-NS_contacts-NS_import), fill="#56B4E9")+
  geom_line(aes(y=NS_contacts+NS_import), color = "black", size=0.2)+
  scale_x_datetime(labels = date_format("%b"),date_breaks = "1 month")+
  scale_y_log10()+
  ylab("incidence")+
  xlab("")+
  ggtitle("Nova Scotia (990,025)")+theme_bw()

g.NB = ggplot(data = travel_day, aes(x=as.Date(date_report), group =1))+
  geom_ribbon(aes(ymin = 0, ymax =  NB_new_cases-NB_contacts-NB_import), fill="#999999", alpha = 0.8)+
  geom_line(aes(y=NB_contacts+NB_import), color = "black", size=0.2)+
  scale_x_datetime(labels = date_format("%b"),date_breaks = "1 month")+
  scale_y_log10()+
  ylab("incidence")+
  xlab("")+
  ggtitle("New Brunswick (784,950)")+theme_bw()

g.PEI = ggplot(data = travel_day, aes(x=as.Date(date_report), group =1))+
  geom_ribbon(aes(ymin = 0, ymax =  PEI_new_cases-PEI_contacts-PEI_import), fill="#F0E442")+
  geom_line(aes(y=PEI_contacts+PEI_import), color = "black", size=0.2)+
  scale_x_datetime(labels = date_format("%b"),date_breaks = "1 month")+
  scale_y_log10()+
  ylab("incidence")+
  xlab("")+
  ggtitle("Prince Edward Island (159,240)")+theme_bw()

g.NWT = ggplot(data = travel_day, aes(x=as.Date(date_report), group =1))+
  geom_ribbon(aes(ymin = 0, ymax = NWT_new_cases-NWT_contacts-NWT_import), fill="#009E73")+
  geom_line(aes(y=NWT_contacts+NWT_import), color = "black", size=0.2)+
  scale_x_datetime(labels = date_format("%b"),date_breaks = "1 month")+
  scale_y_log10()+
  ylab("incidence")+
  xlab("")+
  ggtitle("Northwest territories (44,432)")+theme_bw()

g.YT = ggplot(data = travel_day, aes(x=as.Date(date_report), group =1))+
  geom_ribbon(aes(ymin = 0, ymax =  YT_new_cases-YT_contacts-YT_import), fill="#0072B2")+
  geom_line(aes(y=YT_contacts+YT_import), color = "black",size=0.2)+
  scale_x_datetime(labels = date_format("%b"),date_breaks = "1 month")+
  scale_y_log10()+
  ylab("incidence")+
  xlab("")+
  ggtitle("Yukon (42,157)")+theme_bw()

g.NU = ggplot(data = travel_day, aes(x=as.Date(date_report), group =1))+
  geom_ribbon(aes(ymin = 0, ymax =  NU_new_cases-NU_contacts-NU_import), fill="#C779A7")+
  geom_line(aes(y=NU_contacts+NU_import), color = "black", size=0.2)+
  scale_x_datetime(labels = date_format("%b"),date_breaks = "1 month")+
  scale_y_log10()+
  ylab("incidence")+
  xlab("")+
  ggtitle("Nunavut (39,750)")+theme_bw()

g1=(g.NS+g.NB)/(g.NL+g.PEI)/(g.YT+g.NWT)/(g.NU+plot_spacer())

ggsave('code-data-viz/travel-related-vs-new.png', width=10,height=10)

PEI_daily_rolling28 = data.frame(date_report=tail(travel_day$date_report,-27), community = rollmean(travel_day$PEI_new_cases-travel_day$PEI_contacts-travel_day$PEI_import,28), travel_contact = rollmean(travel_day$PEI_contacts+travel_day$PEI_import,28))
YT_daily_rolling28 = data.frame(date_report=tail(travel_day$date_report,-27), community = rollmean(travel_day$YT_new_cases-travel_day$YT_contacts-travel_day$YT_import,28), travel_contact = rollmean(travel_day$YT_contacts+travel_day$YT_import,28))
NWT_daily_rolling28 = data.frame(date_report=tail(travel_day$date_report,-27), community = rollmean(travel_day$NWT_new_cases-travel_day$NWT_contacts-travel_day$NWT_import, 28), travel_contact = rollmean(travel_day$NWT_contacts+travel_day$NWT_import,28))

