# This is the script for creating good mean plots with 
#   95% confidence interval bars. The confidence intervals
#   may be either t-parametric or bca bootstrap intervals.
# The algorithm was developed for cell analysis mean plots
#   of mitotic index (MI) and frequency of abberant cells (FAC).
#
# Copyright: Ilia Fastovets (2018)
# Contact: fastovetsilya@yandex.ru


library('readxl')
library('ggplot2')
#library('extrafont')
library('boot')

Data <- read_excel('Correlations_La_Ce_rastv_pochv.xlsx')
Datamat <- as.matrix(Data)

Element <- 'Ce'
Medium <- 'pochv'
Feature <- 'RootLength'

Index_Group <- paste(Element, Medium, sep = '_')
Index_Feature <- paste(Element, Medium, Feature, sep = '_')
D_r <- na.exclude(Datamat[,c(Index_Group, Index_Feature)])

create_summary <- function(D_r, method = 't') # 't' or 'boot'
{
  
  if(method == 't')
  {
    conf_interval_low <- function(x)
    {
      return(t.test(x)$conf.int[1])
    }
    conf_interval_high <- function(x)
    {
      return(t.test(x)$conf.int[2])
    }
    
    D_mean <- aggregate(x = D_r[,2], by = list(D_r[,1]), mean)
    D_low <- aggregate(x = D_r[,2], by = list(D_r[,1]), 
                       conf_interval_low)
    D_high <- aggregate(x = D_r[,2], by = list(D_r[,1]), 
                        conf_interval_high)
    D_agg <- cbind(D_mean, D_low[2], D_high[2])
    colnames(D_agg) <- c('Group', 'Mean', 'Low', 'High')
    D_agg[, 'Low'][D_agg[, 'Low'] < 0] <- 0
    return(as.data.frame(D_agg))
  }
  
  if(method == 'boot')
  {
    
    conf_interval_low <- function(x, R = 10^5)
    {
      bootfun <- function(data, index)
      {
        d <- data[index]
        return(mean(d))
      }
      bootresult <- boot(na.omit(x), bootfun, R = R)
      bootci <- boot.ci(bootresult, type = 'bca')
      return(bootci$bca[4])
    }
    
    conf_interval_high <- function(x, R = 10^5)
    {
      bootfun <- function(data, index)
      {
        d <- data[index]
        return(mean(d))
      }
      bootresult <- boot(na.omit(x), bootfun, R = R)
      bootci <- boot.ci(bootresult, type = 'bca')
      return(bootci$bca[5])
    }
    
    D_mean <- aggregate(x = D_r[,2], by = list(D_r[,1]), mean)
    D_low <- aggregate(x = D_r[,2], by = list(D_r[,1]), 
                       conf_interval_low)
    D_high <- aggregate(x = D_r[,2], by = list(D_r[,1]), 
                        conf_interval_high)
    D_agg <- cbind(D_mean, D_low[2], D_high[2])
    colnames(D_agg) <- c('Group', 'Mean', 'Low', 'High')
    D_agg[, 'Low'][D_agg[, 'Low'] < 0] <- 0
    return(as.data.frame(D_agg))
  }
  
}

D_agg <- create_summary(D_r, 't')


#loadfonts(device = "win")
#font_install('TT Arial')
good_plot <- function(D_agg, title, xlab, ylab)
{
  ggplot(D_agg, aes(x=Group, y=Mean, vjust = 1.5)) +
    geom_line(aes(group=1)) +
    geom_point(size=2) +
    geom_errorbar(aes(ymin=Low, ymax=High), 
                  width=2) +
    xlab(paste0(xlab)) + 
    ylab(paste0(ylab)) +
    ggtitle(title) +
    scale_x_discrete(limits = D_agg$Group, expand = c(0.05,0.05)) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text = element_text(face='plain',size=24, colour='black'),
          text = element_text(face="bold",size=24, colour='black'),
          axis.title.x = element_text(margin = margin(t = 20)),
          axis.title.y = element_text(margin = margin(r = 20))) +
    ylim(0,max(D_agg$High)+1)
         
}

# good_plot(D_agg, '', 
#           paste0(Element,' concentration', ', mg/l'), 
#           paste0(Feature, ', %'))

# good_plot(D_agg, '',
#           paste0(Element,' extraneous concentration', ', mg/kg'),
#           paste0(Feature, ', %'))

# good_plot(D_agg, '',
#           paste0(Element,' concentration', ', mg/l'),
#           'Root length, mm')

good_plot(D_agg, '',
          paste0(Element,' extraneous concentration', ', mg/kg'),
          'Root length, mm')
