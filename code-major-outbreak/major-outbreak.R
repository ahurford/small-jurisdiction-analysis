# Figure 4: Expected rate of infectious travel-related cases in the NL community

library(ggplot2)
library(scales)
library(zoo)
library(bbmle)
library(dplyr)
library(tidyr)
library(imputeTS)
library(patchwork)

## Color map
cb = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")

##---------- 1. CASES DATA
n1 <- read.csv('data/travel-day.csv')[-1]%>%
  rename(date=date_report)
n1$date <- as.Date(n1$date)
dates<-data.frame(date=n1$date)

# Newfoundland & Labrador
cases <- select(n1,date, NL_new_cases, NL_import, NL_contacts)%>%
  mutate(community = NL_new_cases-NL_contacts-NL_import)%>%
  select(-NL_new_cases)%>%
  rename(contacts=NL_contacts, import = NL_import)%>%
  pivot_longer(
    cols = `import`:`community`, 
    names_to = "cases",
    values_to = "value"
  )%>%
  mutate(jurisdiction = "NL")

# Nova Scotia
n2 <- select(n1,date, NS_new_cases, NS_import, NS_contacts)%>%
  mutate(community = NS_new_cases-NS_contacts-NS_import)%>%
  select(-NS_new_cases)%>%
  rename(contacts=NS_contacts, import = NS_import)%>%
  pivot_longer(
    cols = `import`:`community`, 
    names_to = "cases",
    values_to = "value"
  )%>%
  mutate(jurisdiction = "NS")
cases <- rbind(cases,n2)

# New Brunswick
n2 <- select(n1,date, NB_new_cases, NB_import, NB_contacts)%>%
  mutate(community = NB_new_cases-NB_contacts-NB_import)%>%
  select(-NB_new_cases)%>%
  rename(contacts=NB_contacts, import = NB_import)%>%
  pivot_longer(
    cols = `import`:`community`, 
    names_to = "cases",
    values_to = "value"
  )%>%
  mutate(jurisdiction = "NB")
cases <- rbind(cases,n2)

# PEI
n2 <- select(n1,date, PEI_new_cases, PEI_import, PEI_contacts)%>%
  mutate(community = PEI_new_cases-PEI_contacts-PEI_import)%>%
  select(-PEI_new_cases)%>%
  rename(contacts=PEI_contacts, import = PEI_import)%>%
  pivot_longer(
    cols = `import`:`community`, 
    names_to = "cases",
    values_to = "value"
  )%>%
  mutate(jurisdiction = "PEI")
cases <- rbind(cases,n2)

#YT
n2 <- select(n1,date, YT_new_cases, YT_import, YT_contacts)%>%
  mutate(community = YT_new_cases-YT_contacts-YT_import)%>%
  select(-YT_new_cases)%>%
  rename(contacts=YT_contacts, import = YT_import)%>%
  pivot_longer(
    cols = `import`:`community`, 
    names_to = "cases",
    values_to = "value"
  )%>%
  mutate(jurisdiction = "YT")
cases <- rbind(cases,n2)

#YWT
n2 <- select(n1,date, NWT_new_cases, NWT_import, NWT_contacts)%>%
  mutate(community = NWT_new_cases-NWT_contacts-NWT_import)%>%
  select(-NWT_new_cases)%>%
  rename(contacts=NWT_contacts, import = NWT_import)%>%
  pivot_longer(
    cols = `import`:`community`, 
    names_to = "cases",
    values_to = "value"
  )%>%
  mutate(jurisdiction = "NWT")
cases <- rbind(cases,n2)

#NU
n2 <- select(n1,date, NU_new_cases, NU_import, NU_contacts)%>%
  mutate(community = NU_new_cases-NU_contacts-NU_import)%>%
  select(-NU_new_cases)%>%
  rename(contacts=NU_contacts, import = NU_import)%>%
  pivot_longer(
    cols = `import`:`community`, 
    names_to = "cases",
    values_to = "value"
  )%>%
  mutate(jurisdiction = "NU")
cases <- rbind(cases,n2)


## No substantial amount of vaccination occurred during this time.

## 2. VARIANT DATA--------------------
# Variant data (source: PHAC)
variant <- read.csv('data/covid19-epiSummary-variants.csv')[,-1]%>%
  rename(date = "Collection..week.", fraction = "X.CT.Count.of.Sample..")
variant$date <- as.Date(variant$date, format = "%Y-%m-%d")

variant.clean  = function(var){
  variant%>%filter(X_Identifier==var)%>%dplyr::select(date,fraction)%>%
  group_by(date)%>%
  add_tally(fraction)%>%
  dplyr::select(date,n)%>%
  distinct()%>%
  arrange(date)%>%
    rename(freq=n)%>%
  as.data.frame()
}

alpha = variant.clean("Alpha")
delta = variant.clean("Delta")

# Join the alpha variant
data <- left_join(dates,alpha)

var.interp = function(){
  data$freq[1]=0
  data$freq = na_interpolation(data$freq)
  return(data)
}

variant = var.interp()%>%mutate(variant="alpha")
# Join the delta variant
data <- left_join(dates,delta)
delta1 <- var.interp()%>%mutate(variant="delta")
original <- data.frame(date=dates$date, freq=1-variant$freq-delta1$freq, variant="original")

variant <- rbind(variant,delta1,original)

## 3. INTERVENTIONS DATA -----------------
### intervention data: https://www.bankofcanada.ca/markets/market-operations-liquidity-provision/covid-19-actions-support-economy-financial-system/covid-19-stringency-index/
interventions <- read.csv("data/COVID-19_STRINGENCY_INDEX.csv", skip = 26)%>%filter(date>=dates$date[1]&date<=tail(dates$date,1))
# This column corresponds to NL
NPIs <- data.frame(date = interventions$date, stringency=interventions$SAN_LYOJ20210218_C2_S7)%>%
  mutate(jurisdiction="NL")

NPIs.NS <- data.frame(date = interventions$date, stringency=interventions$SAN_LYOJ20210218_C2_S9)%>%
  mutate(jurisdiction="NS")
NPIs <- rbind(NPIs,NPIs.NS)

NPIs.NB <- data.frame(date = interventions$date, stringency=interventions$SAN_LYOJ20210218_C2_S10)%>%
  mutate(jurisdiction="NB")
  NPIs <- rbind(NPIs,NPIs.NB)
  
NPIs.PEI <- data.frame(date = interventions$date, stringency=interventions$SAN_LYOJ20210218_C2_S8)%>%
    mutate(jurisdiction="PEI")
  NPIs <- rbind(NPIs,NPIs.PEI)

# Stringency index unavailable for the territories - assume it to be the same as Atlantic Canada
  NPIs.NU <- data.frame(date = interventions$date, stringency=interventions$SAN_LYOJ20210218_C2_S2)%>%
    mutate(jurisdiction="NU")
  NPIs <- rbind(NPIs,NPIs.NU)
  NPIs.YT <- data.frame(date = interventions$date, stringency=interventions$SAN_LYOJ20210218_C2_S2)%>%
    mutate(jurisdiction="YT")
  NPIs <- rbind(NPIs,NPIs.YT)
  NPIs.NWT <- data.frame(date = interventions$date, stringency=interventions$SAN_LYOJ20210218_C2_S2)%>%
    mutate(jurisdiction="NWT")
  NPIs <- rbind(NPIs,NPIs.NWT)
  
## TRANSMISSION
delta.trans <- 1 # reference value
alpha.trans <- delta.trans/1.97
original.trans <- alpha.trans/1.77 #correction for variant (these are justified in Hurford et al. JTB)

variant<-variant%>%as.data.frame()

alpha <- variant%>%filter(variant=="alpha")
alpha <- alpha$freq
delta <- variant%>%filter(variant=="delta")
delta<-delta$freq
original <- variant%>%filter(variant=="original")
original <-original$freq

NPIs<-NPIs%>%as.data.frame()

### 4. CALCULATE QUANTITIES OF INTEREST
P.major <- function(){
  NPIs <- filter(NPIs, jurisdiction==j)
  NPIs <- NPIs$stringency
  cases <- filter(cases, jurisdiction==j)
  imports <- filter(cases, cases=="import")
  imports <- imports$value
  contacts <- filter(cases, cases=="contacts")
  contacts <- contacts$value
  community <- filter(cases, cases=="community")
  community <- community$value
  community[community<0]=0
  community.max = max(community)
  Rt <- local.trans*NPIs*(alpha*alpha.trans+original*original.trans+delta*delta.trans)/100
  n <- (contacts+imports)
  n <- c(head(n,3),rollsum(n,4))
  P.major <- traveller.trans*n*(1-1/Rt) # assume more than one spillover negligible, and multiply by the probability of one spillover
  P.major[n==0]<-0
  P.major[Rt<=1]<-0
  P.major[P.major>1]<-1
  output <- data.frame(date=dates$date, n=n, Rt=Rt, P.major=P.major, community=community/community.max)
  return(output)
}

j<-"NL"
traveller.trans <- 1/100
local.trans <- 5 
NL <-P.major()%>% mutate(jurisdiction="NL")

j<-"NS"
traveller.trans <- 1/50
local.trans <- 6 
NS<-P.major()%>%mutate(jurisdiction="NS")

j<-"NB"
traveller.trans <- 1/50
local.trans <- 4 
NB<-P.major()

j<-"PEI"
traveller.trans <- 1/100
local.trans <- 4 
PEI<-P.major()

j<-"YT"
traveller.trans <- 1/100
local.trans <- 4 
YT<-P.major()

j<-"NWT"
traveller.trans <- 1/100
local.trans <- 4 
NWT<-P.major()

j<-"NU"
traveller.trans <- 1/100
local.trans <- 4 
NU<-P.major()

make.figure = function(data, color, title){
g1 = ggplot(data=data,aes(x=date))+
  geom_hline(yintercept = 1,lty=2)+
  geom_line(aes(y=Rt),color=color)+
  scale_x_datetime(labels = date_format("%b"),date_breaks = "1 month", guide = guide_axis(angle = 90))+
  xlab("")+
  ylab("effective reproduction number")+
  theme_classic()

g2  = ggplot(data=data,aes(x=date))+
  geom_line(aes(y=n), color=color)+
  scale_x_datetime(labels = date_format("%b"),date_breaks = "1 month", guide = guide_axis(angle = 90))+
  xlab("")+
  ylab("infectious travelers and contacts")+
  theme_classic()

g3  = ggplot(data=data,aes(x=date))+
  geom_ribbon(aes(ymax=community*.25, ymin=0), fill="darkgrey")+
  geom_line(aes(y=P.major), color=color)+
  scale_x_datetime(labels = date_format("%b"),date_breaks = "1 month", guide = guide_axis(angle = 90))+
  xlab("")+
  ylim(c(0,.25))+
  ylab("major outbreak probability")+
  theme_classic()

g=(g1+g2+g3+plot_annotation(title=title, theme = theme(plot.title = element_text(hjust = 0.5))))
return(g)
}

gNL = make.figure(NL,"#E69F00", "Newfoundland and Labrador")
ggsave('code-major-outbreak/NL.png', width=8,height=3)

gNS = make.figure(NS, "#56B4E9","Nova Scotia")
ggsave('code-major-outbreak/NS.png', width=8,height=3)

gNB = make.figure(NB, "black","New Brunswick")
ggsave('code-major-outbreak/NB.png', width=8,height=3)

gPEI = make.figure(PEI, "#F0E442","Prince Edward Island")
ggsave('code-major-outbreak/PEI.png', width=8,height=3)

gYT = make.figure(YT,"#0072B2", "Yukon")
ggsave('code-major-outbreak/YT.png', width=8,height=3)

gNWT = make.figure(NWT, "#009E73","Northwest territories")
ggsave('code-major-outbreak/NWT.png', width=8,height=3)

gNU = make.figure(NU, "#C779A7","Nunavut")
ggsave('code-major-outbreak/NU.png', width=8,height=3)

