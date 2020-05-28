
############################################# Scale on WE

rm(list = ls(all.names = TRUE)) 

require(ggplot2)

## Create a data vector of the frequency of cell sizes (a,b,c,d,) used in the AOO WE literature 
a<-rep(1, 2)
b<-rep(2, 6)
c<-rep(3, 6)
d<-rep(4, 9)
data<-c(a,b,c,d)

## Plotting the frequency bar chart 
p<-ggplot(data.frame(data), aes(x=factor(data))) +
  geom_bar(width = 0.4)+
  theme_bw()+
  scale_x_discrete(labels=c("0.0001 - 99", "100 - 999", "1000 - 9999", "10 000 +"))+
  xlab("Grid cell size (km2)") + 
  ylab("Frequency") +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank(), axis.ticks.x = element_blank(),
        axis.text=element_text(size=14), axis.title=element_text(size=15))+
  scale_y_continuous(breaks=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9), expand = expansion(mult = c(0, 0.1)))
p
  




