# Figure 4: Expected rate of infectious travel-related cases in the NL community

library(ggplot2)
library(scales)
library(zoo)
library(bbmle)
library(dplyr)
library(imputeTS)
library(patchwork)

## Color map
cb = c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2")

##---------- 1. LOAD AND CLEAN DATA
## Load the data & clean to minimum level with date and the time variable
# Travel-related cases arriving in NL, n (source: NLCHI data)
n <- read.csv('~/Desktop/Work/Research/Research_Projects/2022/reopening/pandemic-COVID-zero/data/NL-travel.csv')[,-1]%>%
  rename(date=REPORTED_DATE)%>%filter(date<"2021-12-25")
n <- data.frame(date=n$date, n = n$TRAVEL, c=n$CLOSE_CONTACT)
n$date <- as.Date(n$date)

# Load community cases
community <- read.csv('~/Desktop/Work/Research/Research_Projects/2022/reopening/pandemic-COVID-zero/data/NLCHI_cases.csv')[,-1]%>%
  rename(date=REPORTED_DATE)%>%filter(date<"2021-12-25")%>%
  dplyr::select(date,COMMUNITY)
community$date = as.Date(community$date)

# Vaccination data (source: PHAC)
vaccination <- read.csv('~/Desktop/Work/Research/Research_Projects/2022/reopening/pandemic-COVID-zero/data/vaccination-coverage-map.csv')%>%
  rename(date = week_end)
vaccination$date = as.Date(vaccination$date)
vaccination$proptotal_additional[is.na(vaccination$proptotal_additional)]=0
vaccination$proptotal_fully = as.numeric(vaccination$proptotal_fully)
vaccination <- mutate(vaccination, proptotal_fully = proptotal_fully/100)%>%
  mutate(proptotal_partially = proptotal_partially/100)%>%
  mutate(proptotal_additional = proptotal_additional/100)

vaccination.Canada <- filter(vaccination,prename=="Canada")%>%
  dplyr::select(date, proptotal_partially, proptotal_fully, proptotal_additional)
vaccination.NL <- filter(vaccination,prename=="Newfoundland and Labrador")%>%
  dplyr::select(date, proptotal_partially, proptotal_fully, proptotal_additional)

# Variant data (source: PHAC)
variant <- read.csv('~/Desktop/Work/Research/Research_Projects/2022/reopening/pandemic-COVID-zero/data/covid19-epiSummary-variants.csv')[,-1]%>%
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
omicron = variant.clean("BA.1")

#-------------
# Align all the data sources, fill NAs with 0s or interpolate missing values
start.date = as.Date(min(n$date))
data = data.frame(date = seq(from=start.date, to = as.Date("2021-12-24"), by = "days"))
data = left_join(data, n)

# replace any NAs created by missing dates with 0 travel-related cases
data[is.na(data)]=0

# Add the community cases
data <- left_join(data, community)
data[is.na(data)]=0

# Join the Canadian vaccination data
data <- left_join(data, vaccination.Canada)

vacc.interp = function(){
# Linear interpolation for the vaccination data where NAs were created during the left_join
data$proptotal_partially = na_interpolation(data$proptotal_partially)
data$proptotal_fully = na_interpolation(data$proptotal_fully)
data$proptotal_additional = na_interpolation(data$proptotal_additional)
# Exclude additional from full
data$proptotal_fully = data$proptotal_fully- data$proptotal_additional
return(data)
}

# Interpolate and clean-up the Canadian vaccination data
data = vacc.interp()
# Rename the Canadian vaccination data columns
data = data %>% rename(CAN.partial = proptotal_partially, CAN.full = proptotal_fully, CAN.additional = proptotal_additional)%>%
  mutate(CAN.unvax = 1 - CAN.full - CAN.partial - CAN.additional)

# Join the NL vaccination data:
data <- left_join(data, vaccination.NL)
# Perform the linear interpolation for the NL vaccination data:
data = vacc.interp()%>% rename(NL.partial = proptotal_partially, NL.full = proptotal_fully, NL.additional = proptotal_additional)%>%
  mutate(NL.unvax = 1 - NL.full - NL.partial - NL.additional)

# Join the alpha variant
data <- left_join(data,alpha)

var.interp = function(){
  data$freq[1]=0
  data$freq = na_interpolation(data$freq)
  return(data)
}

data = var.interp()%>%rename(alpha = freq)
# Join the delta variant
data <- left_join(data,delta)
data <- var.interp()%>%rename(delta = freq)
# Join the omicron variant
data <- left_join(data,omicron)
data <- var.interp()%>%rename(omicron = freq)%>%
  mutate(original = round(1- omicron - delta-alpha,2))
# rounding is to remove -ve numbers
# contacts from travellers:
contacts.data = inner_join(n, data.frame(date=data$date,original = data$original,alpha=data$alpha,delta=data$delta,omicron=data$omicron))%>%filter(n>0)
c.original = filter(contacts.data, original>0.01) %>% mutate(contacts.per = c/n)
c.original = mean(c.original$contacts.per)
c.alpha = filter(contacts.data, alpha>0.01) %>% mutate(contacts.per = c/n)
c.alpha = mean(c.alpha$contacts.per)
c.delta = filter(contacts.data, delta>0.01) %>% mutate(contacts.per = c/n)
c.delta = mean(c.delta$contacts.per)
c.omicron = filter(contacts.data, omicron>0.01) %>% mutate(contacts.per = c/n)
c.omicron = mean(c.omicron$contacts.per)

### intervention data: https://www.bankofcanada.ca/markets/market-operations-liquidity-provision/covid-19-actions-support-economy-financial-system/covid-19-stringency-index/
interventions = read.csv("https://raw.githubusercontent.com/ahurford/pandemic-COVID-zero/main/data/COVID-19_STRINGENCY_INDEX.csv", skip = 26)
# This column corresponds to NL
NPIs = data.frame(date = interventions$date, stringency=interventions$SAN_LYOJ20210218_C2_S7)
# Extract the dates corresponding to data
NPIs = NPIs%>%filter(date>=data$date[1]&date<=tail(data$date,1))


# Clean data from actual community outbreaks
community.outbreaks = data.frame(date= data$date, community = data$COMMUNITY)
ymax=-.05*7
ymin=-.075*7
community.outbreaks$community[which(community.outbreaks$community<=5)]=ymin
community.outbreaks$community[which(community.outbreaks$community>5)]=ymax
community.alpha = community.outbreaks%>%
  filter(date<"2021-05-01")%>%
  rename(alpha=community)
community.delta= community.outbreaks%>%
  filter(date>="2021-05-01"&date<="2021-12-01")%>%
  rename(delta=community)
community.omicron= community.outbreaks%>%
  filter(date>="2021-12-01")%>%
  rename(omicron=community)
community.outbreaks=left_join(community.outbreaks, community.alpha)%>%
  left_join(community.delta)%>%
  left_join(community.omicron)
community.outbreaks[is.na(community.outbreaks)]=ymin

###------- POST-ARRIVAL TRAVEL RESTRICTIONS
# Vaccine efficacies and calculating, T_j,k, vaccination status and variant given infection
# These are the Pfizer values
# https://www.nejm.org/doi/full/10.1056/NEJMoa2119451 (Omicron values, 2-does 25+ weeks)
# rows are vaccination status: unvax, 1-dose, 2-doses and 3-doses.
VE = data.frame(original = c(0, .49, .93, .93), alpha = c(0, .49, .93, .93), delta = c(0, .33, 0.88, 0.88), omicron = c(0, 0, 0.09, 0.67))
T.original1 = data.frame(unvax = data$CAN.unvax*VE$original[1], partial = data$CAN.partial*VE$original[2], full = data$CAN.full*VE$original[3], additional = data$CAN.additional*VE$original[4])*data$original
T.alpha1 = data.frame(unvax = data$CAN.unvax*VE$alpha[1], partial = data$CAN.partial*VE$alpha[2], full = data$CAN.full*VE$alpha[3], additional = data$CAN.additional*VE$alpha[4])*data$alpha
T.delta1 = data.frame(unvax = data$CAN.unvax*VE$delta[1], partial = data$CAN.partial*VE$delta[2], full = data$CAN.full*VE$delta[3], additional = data$CAN.additional*VE$delta[4])*data$delta
T.omicron1 = data.frame(unvax = data$CAN.unvax*VE$omicron[1], partial = data$CAN.partial*VE$omicron[2], full = data$CAN.full*VE$omicron[3], additional = data$CAN.additional*VE$omicron[4])*data$omicron

# Normalize by dividing by the sum of the row sums:
T.sum = rowSums(T.original1+T.alpha1+T.delta1+T.omicron1)
T.original = T.original1/T.sum
T.alpha = T.alpha1/T.sum
T.delta = T.delta1/T.sum
T.omicron = T.omicron1/T.sum

test1=rowSums(T.original+T.alpha+T.delta+T.omicron)

# Self-isolation and testing
# Infectivity on each day based on Ferreti
#dweibull(seq(0,10), 2.83, scale = 5.67, log = FALSE)

# where i is the number of days exposure was prior to arrival
# t is the day of a test.
# The probability of a true positive when days since exposure is
# uniformly distribution from 0 to 10 is the mean of $t$, where $t_i$ is the probability of testing positive
# given infection $i$ days ago:

t.sens <- c(0, 0.05, .1, .55, .78,.77,.726, .682, .638, .594, .55, .49, .43, .37, .31, .25, .22, .19, .16, .13, .1, .09, .08, .07, 0.06, .05)

inf.fun = function(x,t){
  incr = 0.1
  L = 40
  compl=0.7
  exposure <- seq(0,L, incr)
  ysum=NULL
  for(i in seq(0,10)){
    seq1 <- rep(0,length(exposure))
    seq2 <-rep(0,length(exposure))
    # if compliant with self-isolation
  j = min(which(exposure>=(i+x)))
  seq1[j:length(seq1)] = compl*dweibull(exposure[j:length(seq1)], 2.83, scale = 5.67, log = FALSE)
  if((i+t+1)<(length(t.sens)-1)){
  #   # A false negative when the test is administered on day i+t of the infection
    seq1 = seq1*(1-t.sens[(i+t+1)])
  }
  ysum[i+1]=sum(seq1[which(exposure>=i)])*incr
  j1 = min(which(exposure>=(i)))
  # Non-compliance
  seq2[j1:length(seq2)] = (1-compl)*dweibull(exposure[j1:length(seq2)], 2.83, scale = 5.67, log = FALSE)
  ysum[i+1]=sum(seq2[which(exposure>=i)])*incr+ysum[i+1]
  }
  ymean = mean(ysum)
}

L=length(data$date)

## Restrictions for unvaccinated travellers
m.unvax = rep(inf.fun(14,100),L)
# on August 1, test on day 8 + isolation to negative for unvaccinated
i = which(data$date=="2021-08-01")
m.unvax[i:L] <- inf.fun(8,8)


## Restrictions for travellers with 1 dose
# Same restrictions prior re-opening
m.1 = rep(inf.fun(14,100),L)

# July 1 - negative test step 1
i = which(data$date=="2021-07-01")
m.1[i:L] <- inf.fun(0,0)

# August 1 - no measures
i = which(data$date=="2021-08-01")
m.1[i:L] <- inf.fun(0,100)
# Sept 30 - same as unvaccinated
i = which(data$date=="2021-09-30")
m.1[i:L] <- m.unvax[i:L]

## Two dose travellers (same as 3-dose also)
# Same restrictions as other travellers prior to reopening
m.2 <- rep(inf.fun(14,100),L)

# No requirements after July 1.
i = which(data$date=="2021-07-01")
m.2[i:L] <- inf.fun(0,100)

# Dec 21: 5 RAT and isolate for 5 days.
# Since the study ends on Dec 24, the only impact
# is due to self-isolation
i = which(data$date=="2021-12-21")
m.2[i:L] = inf.fun(5,100)*.1
m.unvax[i:L] = min(inf.fun(5,100)*.1, m.unvax[i:L])
m.1[i:L] = min(inf.fun(5,100)*.1, m.1[i:L])

traveller.measures = data.frame(date = data$date,m.unvax, m.1, m.2, m.3=m.2)

####### PLOTS OF POST-ARRIVAL MEASURES
g.travel.measures =ggplot(traveller.measures,aes(as.Date(date),group=1)) +
  geom_line(aes(y = 1-m.unvax), col="grey")+
  #geom_ribbon(aes(ymax = m.unvax, ymin=0), fill="grey", alpha = 0.3)+
  geom_line(aes(y = 1-m.1), col="darkorchid")+
  #geom_ribbon(aes(ymax = m.1, ymin=m.unvax), fill="darkorchid", alpha = 0.3)+
  geom_line(aes(y = 1-m.3), col="dodgerblue")+
  geom_line(aes(y = 1-m.2), col=palette.colors(2)[2])+
  #geom_ribbon(aes(ymax = m.2, ymin=m.1), fill=palette.colors(2)[2], alpha = 0.3)+
  geom_line(data = data.frame(x = c(as.Date("2021-07-01"), as.Date("2021-07-01")), y = c(0, 1)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  geom_line(data = data.frame(x = c(as.Date("2021-08-01"), as.Date("2021-08-01")), y = c(0, 1)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  geom_line(data = data.frame(x = c(as.Date("2021-09-30"), as.Date("2021-09-30")), y = c(0, 1)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  geom_line(data = data.frame(x = c(as.Date("2021-12-21"), as.Date("2021-12-21")), y = c(0, 1)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  scale_x_date(breaks = date_breaks("1 month"),
               labels = date_format("%b %Y"), limits = c(as.Date("2021-04-01"), as.Date("2021-12-24")))+
  xlab("") +
  ylab("Stringency")+
  ggtitle("NL post-arrival travel restrictions")+
  #coord_cartesian(ylim=c(0, 25))+
  annotate("text", x = as.Date("2021-09-01"), y = 0.95, label = "0 doses", col = "darkgrey", angle=0, fontface=2)+
  annotate("text", x = as.Date("2021-11-01"), y = 0.95, label = "1 dose", col = "darkorchid", angle=0, fontface=2)+
  annotate("text", x = as.Date("2021-11-01"), y = 0.6, label = "2+ doses", col = palette.colors(2)[2], angle=0, fontface=2)+
  annotate("text", x = as.Date("2021-06-23"), y = 0.2, label = "Step 1", angle=90, col  ="darkgrey", fontface=2)+
  annotate("text", x = as.Date("2021-07-23"), y = 0.2, label = "Step 2", angle=90, col  ="darkgrey", fontface=2)+
  annotate("text", x = as.Date("2021-09-23"), y = 0.2, label = "Step 2a", angle=90, col  ="darkgrey", fontface=2)+
  annotate("text", x = as.Date("2021-12-12"), y = 0.2, label = "Step 2b", angle=90, col  ="darkgrey", fontface=2)+
  #geom_line(data = data.frame(x = c(as.Date("2021-12-15"), as.Date("2021-12-15")), y = c(0, 1)), aes(x = x, y = y),lty=2, col = palette.colors(7)[7])+
  #geom_line(data = data.frame(x = c(as.Date("2021-07-01"), as.Date("2021-07-01")), y = c(0, 1)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  #annotate("text", x = as.Date("2021-12-24"), y = .75, label = "+1", col = "black", fontface=2)+
  theme_classic() + theme(axis.text.x = element_text(angle = 90, size=rel(1)), legend.title = element_blank(),legend.text=element_text(size=rel(1.2)),plot.title=element_text(size=rel(1.2), face="bold"),axis.title = element_text(size=rel(1.2)),axis.text.y = element_text(size=rel(1.2)))

travel.unvax= m.unvax*(T.original$unvax+T.alpha$unvax+T.delta$unvax+T.omicron$unvax)*data$CAN.unvax
travel.partial = m.1*(T.original$partial+T.alpha$partial+T.delta$partial+T.omicron$partial)*data$CAN.partial
travel.full = m.2*(T.original$full+T.alpha$full+T.delta$full+T.omicron$full)*data$CAN.full
travel.additional = m.2*(T.original$additional+T.alpha$additional+T.delta$additional+T.omicron$additional)*data$CAN.additional
traveller.measures2 = data.frame(date = data$date,unvax = travel.unvax ,partial=travel.partial,full=travel.full, additional = travel.additional)

######### NL COMMUNITY MEASURES
# Vulnerability of the NL community to different variants
NL.original = data$NL.unvax*(1-VE$original[1]) + data$NL.partial*(1-VE$original[2]) + data$NL.full*(1-VE$original[3]) + data$NL.additional*(1-VE$original[4])
NL.alpha = data$NL.unvax*(1-VE$alpha[1]) + data$NL.partial*(1-VE$alpha[2]) + data$NL.full*(1-VE$alpha[3]) + data$NL.additional*(1-VE$alpha[4])
NL.delta = data$NL.unvax*(1-VE$delta[1]) + data$NL.partial*(1-VE$delta[2]) + data$NL.full*(1-VE$delta[3]) + data$NL.additional*(1-VE$delta[4])
NL.omicron = data$NL.unvax*(1-VE$omicron[1]) + data$NL.partial*(1-VE$omicron[2]) + data$NL.full*(1-VE$omicron[3]) + data$NL.additional*(1-VE$omicron[4])

PIs.NPIs = data.frame(date = data$date, NPIs=NPIs$stringency/100, NL.original, NL.alpha, NL.delta, NL.omicron)

n.original = data$n*(T.original$unvax+T.original$partial+T.original$full+T.original$additional)
n.alpha = data$n*(T.alpha$unvax+T.alpha$partial+T.alpha$full+T.alpha$additional)
n.delta = data$n*(T.delta$unvax+T.delta$partial+T.delta$full+T.delta$additional)
n.omicron = data$n*(T.omicron$unvax+T.omicron$partial+T.omicron$full+T.omicron$additional)

n.variant = data.frame(date = data$date, original = n.original, alpha = n.alpha, delta = n.delta, omicron = n.omicron)
n.variant$original = c(n.variant$original[1:6], rollmean(n.variant$original, 7))
n.variant$alpha = c(n.variant$alpha[1:6], rollmean(n.variant$alpha, 7))
n.variant$delta = c(n.variant$delta[1:6], rollmean(n.variant$delta, 7))
n.variant$omicron = c(n.variant$omicron[1:6], rollmean(n.variant$omicron, 7))

## TRANSMISSION
base.trans = .2
delta.trans = 1 # reference value
omicron.trans = delta.trans*1.51
alpha.trans =delta.trans/1.97
original.trans = alpha.trans/1.77 #correction for variant

#Expected spillovers
E.original = base.trans*original.trans*(1-NPIs$stringency/100)*(PIs.NPIs$NL.original)*(1+c.original)*(data$n)*(T.original$unvax*m.unvax+T.original$partial*m.1+T.original$full*m.2+T.original$additional*m.2)
E.alpha = base.trans*alpha.trans*(1-NPIs$stringency/100)*(PIs.NPIs$NL.alpha)*(data$n)*(1+c.alpha)*(T.alpha$unvax*m.unvax+T.alpha$partial*m.1+T.alpha$full*m.2+T.alpha$additional*m.2)
E.delta = base.trans*delta.trans*(1-NPIs$stringency/100)*(PIs.NPIs$NL.delta)*(data$n)*(1+c.delta)*(T.delta$unvax*m.unvax+T.delta$partial*m.1+T.delta$full*m.2+T.delta$additional*m.2)
E.omicron = base.trans*omicron.trans*(1-NPIs$stringency/100)*(PIs.NPIs$NL.omicron)*(data$n)*(1+c.omicron)*(T.omicron$unvax*m.unvax+T.omicron$partial*m.1+T.omicron$full*m.2+T.omicron$additional*m.2)

# Take rolling mean
E.original = c(E.original[1:6], rollmean(E.original, 7))
E.alpha = c(E.alpha[1:6], rollmean(E.alpha, 7))
E.delta = c(E.delta[1:6], rollmean(E.delta, 7))
E.omicron = c(E.omicron[1:6], rollmean(E.omicron, 7))

E.spillovers = data.frame(date = data$date, original = E.original, alpha = E.alpha, delta = E.delta, omicron = E.omicron, total = E.original+E.alpha+E.delta+E.omicron)



### PLOTS
g.var =ggplot(data,aes(as.Date(date),group=1)) +
  geom_ribbon(aes(ymax=1, ymin=alpha+delta+omicron), fill="grey", alpha=0.3)+
  geom_ribbon(aes(ymax = alpha+delta+omicron, ymin=alpha+delta), fill=palette.colors(7)[7], alpha=.8)+
  geom_ribbon(aes(ymax = alpha+delta, ymin=alpha), fill=palette.colors(3)[3], alpha=.8)+
  geom_ribbon(aes(ymax = alpha, ymin=0), fill=palette.colors(4)[4], alpha=.8)+
  scale_x_date(breaks = date_breaks("1 month"),
               labels = date_format("%b %Y"),limits = c(as.Date("2021-01-01"), as.Date("2021-12-24")))+
  xlab("") +
  ylab("proportion")+
  ggtitle("Variants")+
  #coord_cartesian(ylim=c(0, 25))+
  annotate("text", x = as.Date("2021-12-10"), y = .9, label = "BA.1", fontface=2)+
  annotate("text", x = as.Date("2021-09-01"), y = .6, label = "Delta", col = "black", fontface=2)+
  annotate("text", x = as.Date("2021-04-07"), y = .25, label = "Alpha", col = "black", fontface=2)+
  annotate("text", x = as.Date("2021-03-01"), y = .9, label = "Original", col = "black", fontface=2)+
  theme_classic() + theme(axis.text.x = element_text(angle = 90, size=rel(1)), legend.title = element_blank(),legend.text=element_text(size=rel(1.2)),plot.title=element_text(size=rel(1.2)),axis.title = element_text(size=rel(1)))

g.vax =ggplot(data,aes(as.Date(date),group=1)) +
  geom_ribbon(aes(ymax=1, ymin=1-CAN.unvax), fill="grey", alpha=0.3)+
  geom_ribbon(aes(ymax = CAN.full+CAN.partial+CAN.additional, ymin=CAN.partial+CAN.full), fill="dodgerblue", alpha=.5)+
  geom_ribbon(aes(ymax = CAN.partial+CAN.full, ymin=CAN.partial), fill=palette.colors(2)[2], alpha=.5)+
  geom_ribbon(aes(ymax = CAN.partial, ymin=0), fill="darkorchid", alpha=.5)+
  geom_line(aes(y=NL.full+NL.partial+NL.additional), col="green")+
  geom_line(aes(y=NL.partial+NL.full), col=palette.colors(2)[2])+
  geom_line(aes(y=NL.partial), col="darkorchid")+
  scale_x_date(breaks = date_breaks("1 month"),
               labels = date_format("%b %Y"),limits = c(as.Date("2021-01-01"), as.Date("2021-12-24")))+
  xlab("") +
  ylab("proportion")+
  ggtitle("Vaccination")+
  #coord_cartesian(ylim=c(0, 25))+
  annotate("text", x = as.Date("2021-11-01"), y = .6, label = "2 doses", fontface=2)+
  annotate("text", x = as.Date("2021-05-25"), y = .1, label = "1 dose", fontface=2)+
  annotate("text", x = as.Date("2021-03-01"), y = .9, label = "0 doses", col = "black", fontface=2)+
  annotate("text", x = as.Date("2021-12-24"), y = .75, label = "3", col = "black", fontface=2)+
  theme_classic() + theme(axis.text.x = element_text(angle = 90, size=rel(1)), legend.title = element_blank(),legend.text=element_text(size=rel(1.2)),plot.title=element_text(size=rel(1.2)),axis.title = element_text(size=rel(1)))

g.NPIs =ggplot(PIs.NPIs,aes(as.Date(date),group=1)) +
  #geom_ribbon(aes(ymax = NPIs, ymin = 0), fill="black", alpha =.1)+
  geom_line(data = data.frame(x = c(as.Date("2021-07-01"), as.Date("2021-07-01")), y = c(0, 1)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  geom_line(data = data.frame(x = c(as.Date("2021-08-01"), as.Date("2021-08-01")), y = c(0, 1)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  geom_line(data = data.frame(x = c(as.Date("2021-09-30"), as.Date("2021-09-30")), y = c(0, 1)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  geom_line(data = data.frame(x = c(as.Date("2021-12-21"), as.Date("2021-12-21")), y = c(0, 1)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  annotate("text", x = as.Date("2021-06-23"), y = 0.78, label = "Step 1", angle=90, col  ="darkgrey", fontface=2)+
  annotate("text", x = as.Date("2021-07-23"), y = 0.78, label = "Step 2", angle=90, col  ="darkgrey", fontface=2)+
  annotate("text", x = as.Date("2021-09-23"), y = 0.8, label = "Step 2a", angle=90, col  ="darkgrey", fontface=2)+
  annotate("text", x = as.Date("2021-12-08"), y = 0.78, label = "Step 2b", angle=90, col  ="darkgrey", fontface=2)+
  geom_line(aes(y = NPIs), col="black")+
  scale_x_date(breaks = date_breaks("1 month"),
               labels = date_format("%b %Y"), limits = c(as.Date("2021-01-01"), as.Date("2021-12-24")))+
  xlab("") +
  ylab("Stringency")+
  ggtitle("NL community measures: non-pharmaceutical")+
  ylim(c(0,1))+
  theme_classic() + theme(axis.text.x = element_text(angle = 90, size=rel(1)), axis.text.y = element_text(size=rel(1.2)), legend.title = element_blank(),legend.text=element_text(size=rel(1.2)),plot.title=element_text(size=rel(1.2),face="bold"),axis.title = element_text(size=rel(1.2)))

g.PIs =ggplot(PIs.NPIs,aes(as.Date(date),group=1)) +
  #geom_ribbon(aes(ymax = NL.original, ymin = 0), fill="grey", alpha = 0.3)+
  #geom_ribbon(aes(ymax = NL.alpha, ymin = 0), fill=palette.colors(4)[4], alpha = 0.3)+
  #geom_ribbon(aes(ymax = NL.delta, ymin = NL.alpha), fill=palette.colors(3)[3],alpha=0.3)+
  #geom_ribbon(aes(ymax = NL.omicron, ymin = NL.delta), fill=palette.colors(7)[7], alpha = 0.3)+
  geom_line(data = data.frame(x = c(as.Date("2021-07-01"), as.Date("2021-07-01")), y = c(0, 1)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  geom_line(data = data.frame(x = c(as.Date("2021-08-01"), as.Date("2021-08-01")), y = c(0, 1)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  geom_line(data = data.frame(x = c(as.Date("2021-09-30"), as.Date("2021-09-30")), y = c(0, 1)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  geom_line(data = data.frame(x = c(as.Date("2021-12-21"), as.Date("2021-12-21")), y = c(0, 1)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  annotate("text", x = as.Date("2021-06-23"), y = 0.25, label = "Step 1", angle=90, col  ="darkgrey", fontface=2)+
  annotate("text", x = as.Date("2021-07-23"), y = 0.25, label = "Step 2", angle=90, col  ="darkgrey", fontface=2)+
  annotate("text", x = as.Date("2021-09-23"), y = 0.65, label = "Step 2a", angle=90, col  ="darkgrey", fontface=2)+
  annotate("text", x = as.Date("2021-12-08"), y = 0.55, label = "Step 2b", angle=90, col  ="darkgrey", fontface=2)+
  geom_line(aes(y = NL.original), col="grey")+
  geom_line(aes(y = NL.alpha), col=palette.colors(4)[4])+
  geom_line(aes(y = NL.delta), col=palette.colors(3)[3])+
  geom_line(aes(y = NL.omicron), col=palette.colors(7)[7])+
  scale_x_date(breaks = date_breaks("1 month"),
               labels = date_format("%b %Y"), limits = c(as.Date("2021-01-01"), as.Date("2021-12-24")))+
  ylab("prob. of symptom. infection")+
  xlab("")+
  ggtitle("NL community measures: pharmaceutical")+
  annotate("text", x = as.Date("2021-12-1"), y = .85, label = "BA.1", col=palette.colors(7)[7], fontface=2)+
  annotate("text", x = as.Date("2021-08-25"), y = .5, label = "Delta", col=palette.colors(3)[3], fontface=2)+
  annotate("text", x = as.Date("2021-03-20"), y = .85, label = "Alpha", col=palette.colors(4)[4], fontface=2)+
  #annotate("text", x = as.Date("2021-03-01"), y = .9, label = "Original", col = "grey")+
  ylim(c(0,1))+
  theme_classic() + theme(axis.text.x = element_text(angle = 90, size=rel(1)), legend.title = element_blank(),legend.text=element_text(size=rel(1.2)),plot.title=element_text(size=rel(1.2), face="bold"),axis.title = element_text(size=rel(1.2)), axis.text.y = element_text(size=rel(1.2)))

g.NPIs + g.PIs +plot_annotation(tag_levels = 'A', title = "NL community measures",theme = theme(plot.title = element_text(hjust = 0.5,face = "bold")))
ggsave("community.png", width = 8, height=4)

# For plotting the shading under the plot
val = which(E.spillovers$omicron-E.spillovers$delta>0)
min.omicron=rep(0,length(E.spillovers[,1]))
max.omicron = min.omicron
min.omicron[val] = E.spillovers$delta[val]
max.omicron[val] = E.spillovers$omicron[val]
val = which(E.spillovers$delta-E.spillovers$alpha>0)
min.delta=rep(0,length(E.spillovers[,1]))
max.delta = min.delta
min.delta[val] = E.spillovers$alpha[val]
max.delta[val] = E.spillovers$delta[val]
val = which(E.spillovers$alpha-E.spillovers$original>0)
min.alpha=rep(0,length(E.spillovers[,1]))
max.alpha = min.alpha
min.alpha[val] = E.spillovers$original[val]
max.alpha[val] = E.spillovers$alpha[val]

# 7 day rolling mean (daily) x 7 days/week
g1 = ggplot(E.spillovers,aes(as.Date(date),group=1))+
  geom_ribbon(aes(ymax = 7*max.omicron, ymin = 7*min.omicron), fill=palette.colors(7)[7], alpha = 0.2)+
  geom_ribbon(aes(ymax = 7*max.delta, ymin = 7*min.delta), fill=palette.colors(3)[3], alpha = 0.2)+
  geom_ribbon(aes(ymax = 7*max.alpha, ymin = 7*min.alpha), fill=palette.colors(4)[4], alpha = 0.2)+
  geom_ribbon(aes(ymax = 7*original, ymin = 0), fill="grey", alpha = 0.2)+
  geom_line(data = data.frame(x = c(as.Date("2021-07-01"), as.Date("2021-07-01")), y = c(0, 1.8)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  geom_line(data = data.frame(x = c(as.Date("2021-08-01"), as.Date("2021-08-01")), y = c(0, 1.8)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  geom_line(data = data.frame(x = c(as.Date("2021-09-30"), as.Date("2021-09-30")), y = c(0, 1.8)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  geom_line(data = data.frame(x = c(as.Date("2021-12-21"), as.Date("2021-12-21")), y = c(0, 1.8)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  geom_line(aes(y=7*total), col = "black")+
  geom_line(aes(y=7*omicron), col=palette.colors(7)[7])+
  geom_line(aes(y=7*delta), col=palette.colors(3)[3])+
  geom_line(aes(y = 7*alpha), col=palette.colors(4)[4])+
  geom_line(aes(y = 7*original), col="grey")+
  geom_ribbon(aes(ymax=community.outbreaks$alpha, ymin = ymin), fill =palette.colors(4)[4])+
  geom_ribbon(aes(ymax=community.outbreaks$delta, ymin = ymin), fill =palette.colors(3)[3])+
  geom_ribbon(aes(ymax=community.outbreaks$omicron, ymin = ymin), fill =palette.colors(7)[7])+
  scale_x_date(breaks = date_breaks("1 month"),
               labels = date_format("%b %Y"),limits = c(as.Date("2021-01-01"), as.Date("2021-12-24")))+
  xlab("") +
  scale_y_continuous(breaks=c(0,.5, 1,1.5,2))+
  ylab("")+
  ggtitle("Expected weekly spillovers to NL community")+
  #geom_line(data = data.frame(x = c(as.Date("2021-12-15"), as.Date("2021-12-15")), y = c(0, .2)), aes(x = x, y = y),lty=2, col = palette.colors(7)[7])+
  annotate("text", x = as.Date("2021-12-10"), y = .1*7, label = "BA.1",col = palette.colors(7)[7], fontface=2)+
  annotate("text", x = as.Date("2021-09-01"), y = .45, label = "Delta", col = palette.colors(3)[3], fontface=2)+
  annotate("text", x = as.Date("2021-05-05"), y = .35, label = "Alpha", col = palette.colors(4)[4], fontface=2)+
  annotate("text", x = as.Date("2021-07-09"), y = .08*7, label = "All", col = "black", fontface=2)+
  annotate("text", x = as.Date("2021-01-10"), y = -.03*7, label = "Community\noutbreaks", col = "black", size=3, fontface=2)+
  annotate("text", x = as.Date("2021-02-20"), y = -.03*7, label = "Mt. Pearl", col = palette.colors(4)[4], size=3, fontface=2)+
  annotate("text", x = as.Date("2021-06-01"), y = -.03*7, label = "Lewisporte", col = palette.colors(3)[3], size=3, fontface=2)+
  annotate("text", x = as.Date("2021-09-29"), y = -.03*7, label = "Baie Verte", col = palette.colors(3)[3], size=3, fontface=2)+
  annotate("text", x = as.Date("2021-09-05"), y = -.03*7, label = "Labrador", col = palette.colors(3)[3], size=3, fontface=2)+
  annotate("text", x = as.Date("2021-10-26"), y = -.03*7, label = "Burin Pen.", col = palette.colors(3)[3], size=3, fontface=2)+
  annotate("text", x = as.Date("2021-12-18"), y = -.03*7, label = "St. John's", col = palette.colors(7)[7], size=3, fontface=2)+
  #annotate("text", x = as.Date("2021-12-10"), y = 0.15, label = "1st BA.1", col = palette.colors(7)[7], angle=90)+
  annotate("text", x = as.Date("2021-06-26"), y = 1.4, label = "Step 1", angle=90, col  ="darkgrey", fontface=2)+
  annotate("text", x = as.Date("2021-07-26"), y = 1.4, label = "Step 2", angle=90, col  ="darkgrey", fontface=2)+
  annotate("text", x = as.Date("2021-09-25"), y = 1.4, label = "Step 2a", angle=90, col  ="darkgrey", fontface=2)+
  annotate("text", x = as.Date("2021-12-24"), y = 0.8, label = "Step 2b", angle=90, col  ="darkgrey", fontface=2)+
  theme_classic() + theme(axis.text.x = element_text(angle = 90, size=rel(1.2)),axis.text.y = element_text(size=rel(1.2)), legend.title = element_blank(),legend.text=element_text(size=rel(1.2)),plot.title=element_text(size=rel(1.5), face="bold"),axis.title = element_text(size=rel(1.2)))

g.n =ggplot(n.variant,aes(as.Date(date),group=1)) +
  geom_line(data = data.frame(x = c(as.Date("2021-12-15"), as.Date("2021-12-15")), y = c(0, 10)), aes(x = x, y = y),lty=2, col = palette.colors(7)[7])+
  geom_line(data = data.frame(x = c(as.Date("2021-07-01"), as.Date("2021-07-01")), y = c(0, 10)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  geom_line(data = data.frame(x = c(as.Date("2021-08-01"), as.Date("2021-08-01")), y = c(0, 10)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  geom_line(data = data.frame(x = c(as.Date("2021-09-30"), as.Date("2021-09-30")), y = c(0, 10)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  geom_line(data = data.frame(x = c(as.Date("2021-12-21"), as.Date("2021-12-21")), y = c(0, 10)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  geom_line(aes(y = alpha+delta+omicron+original), col=palette.colors(7)[7])+
  geom_line(aes(y = alpha+delta+original), col=palette.colors(3)[3])+
  geom_line(aes(y = alpha+original), col=palette.colors(4)[4])+
  geom_line(aes(y = original), col="grey")+
  geom_ribbon(aes(ymax = alpha+delta+omicron+original, ymin=alpha+delta+original), fill=palette.colors(7)[7], alpha=.3)+
  geom_ribbon(aes(ymax = alpha+delta+original, ymin=alpha+original), fill=palette.colors(3)[3], alpha=.3)+
  geom_ribbon(aes(ymax = alpha+original, ymin=original), fill=palette.colors(4)[4], alpha=.3)+
  geom_ribbon(aes(ymax = original, ymin=0), fill="grey", alpha=.3)+
  scale_x_date(breaks = date_breaks("1 month"),
               labels = date_format("%b %Y"),limits = c(as.Date("2021-01-01"), as.Date("2021-12-24")))+
  #scale_y_continuous(trans='log2')+
  xlab("") +
  annotate("text", x = as.Date("2021-06-23"), y = 7.5, label = "Step 1", angle=90, col  ="darkgrey", fontface=2)+
  annotate("text", x = as.Date("2021-07-23"), y = 7.5, label = "Step 2", angle=90, col  ="darkgrey", fontface=2)+
  annotate("text", x = as.Date("2021-09-23"), y = 7.5, label = "Step 2a", angle=90, col  ="darkgrey", fontface=2)+
  annotate("text", x = as.Date("2021-12-09"), y = 5, label = "1st BA.1", col = palette.colors(7)[7], angle=90, fontface=2)+
  #annotate("text", x = as.Date("2021-06-23"), y = 5, label = "Reopening", angle=90, col  ="darkgrey")+
  #geom_line(data = data.frame(x = c(as.Date("2021-07-01"), as.Date("2021-07-01")), y = c(0, 10)), aes(x = x, y = y),lty=2, col = "darkgrey")+
  ylab("7-day rolling mean, daily")+
  ggtitle("Importations to NL")+
  #coord_cartesian(ylim=c(0, .5))+
  annotate("text", x = as.Date("2021-12-07"), y = 9, label = "BA.1", fontface=2, col = palette.colors(7)[7])+
  annotate("text", x = as.Date("2021-10-01"), y = 2.5, label = "Delta", col =palette.colors(3)[3], fontface=2)+
  annotate("text", x = as.Date("2021-04-21"), y = 7, label = "Alpha", fontface=2, col=palette.colors(4)[4])+
  annotate("text", x = as.Date("2021-02-01"), y = 1.4, label = "Original", col = "grey", fontface=2)+
  theme_classic() + theme(axis.text.x = element_text(angle = 90, size=rel(1)), legend.title = element_blank(),legend.text=element_text(size=rel(1.2)),plot.title=element_text(size=rel(1.2), face="bold"),axis.title = element_text(size=rel(1.2)),axis.text.y = element_text(size=rel(1.2)))

gout1 = (g.n+g.travel.measures+g.PIs+g.NPIs)/g1 + plot_annotation(tag_levels = 'A', theme = theme(plot.title = element_text(hjust = 0.5,face = "bold")))
#+plot_layout(height = c(1, 1, 2))
ggsave("~/Desktop/community-outbreak.png", width = 12, height = 12)

(g.var + g.vax)+ plot_annotation(tag_levels = 'A')
ggsave("~/Desktop/vax_var.png", width = 10)

