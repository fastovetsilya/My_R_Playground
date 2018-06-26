library('readxl')
library('ggplot2')
library('reshape2')
library('RColorBrewer')

D <- read_excel('Vytyazhki.xlsx', sheet = 1)
D_reshaped <- melt(D, id.vars=1 )
colnames(D_reshaped) <- c('Group', 'Extractant', 'Concentration')
#D_reshaped <- as.matrix(D_reshaped)

ggplot(D_reshaped) + 
  aes(fill=Group, y=Concentration, x=Extractant) +
  scale_fill_brewer(palette='BuPu')+
  geom_bar(position="dodge", stat="identity") +
  #geom_errorbar(aes(ymin=Plot_LCI[,3], ymax=Plot_UCI[,3]), position = 'dodge') +
  xlab('Extractant') + 
  ylab('Concentration of La, mg/kg')

  