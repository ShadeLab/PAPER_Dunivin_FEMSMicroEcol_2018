setwd("/Users/dunivint/Documents/ShadeLab/Meta_centralia/")

#read data
data=read.delim("cen_18s.txt")

#calculate
data$norm=data$X18S/data$Gbases
data$oddsratio=data$norm/(mean(data$norm))
data$nonodds=data$X18S/(mean(data$X18S))

#plot
ggplot(data, aes(x=Temp, y=oddsratio)) +
  geom_point()

ggplot(data, aes(x=Temp, y=norm)) +
  geom_point(aes(shape=History))
ggplot(data, aes(x=Temp, y=norm)) +
  geom_point(aes(shape=History)) +
  ylab("Assembled 18S Count (JGI)") +
  xlab("Temperature")


##############
MicrobeCensus
##############

census=read.delim("microbe_census.txt")

census=inner_join(census, data)

ggplot(census, aes(x=Temp, y=AGS)) +
  geom_point(aes(shape=History, size=3)) +
  ylab("Average genome size") +
  theme_gray(base_size = 12)





