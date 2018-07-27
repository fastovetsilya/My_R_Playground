library('kSamples')
library('readxl')
library('nparcomp')
library('ggplot2')

# Initialize 
init_dir <- getwd()
Data <- read_excel('Data.xlsx')
N_sim <- 10^6
confidence_intervals <- 't' # t or boot

# Perform Kruskal-Wallis and Gao tests, save .txt files
tests_dir <- paste0(getwd(), '/Tests')
dir.create(tests_dir, showWarnings = FALSE)
setwd(tests_dir)
cat('Kruskal-Wallis tests for all data', file = 'Kruskal.txt')
for (i in colnames(Data)[-1])
{
formula <- as.formula(paste(i, '~ Group')) 
output <- qn.test(formula, data = Data, test = 'KW', 
                  method = 'simulated', Nsim = N_sim)
cat('\n\n', file = 'Kruskal.txt', append = TRUE)
cat(c(as.character(formula)[[2]], as.character(formula)[[3]]), 
    file = 'Kruskal.txt', append = TRUE)
capture.output(output, file = 'Kruskal.txt', append = TRUE)

if(as.numeric(output$qn[3]) <= 0.05)
{
  print('SIGNIFICANT DIFFERENCE FOUND!')
  gaocs_output <- gao_cs(formula, data = Data)
  gao_output <- gao(formula, data = Data)
  
  cat('Gao_cs tests for significant differences', file = 'Gao_cs.txt',
      append = TRUE)
  cat(c(as.character(formula)[[2]], as.character(formula)[[3]]), 
      file = 'Gao_cs.txt', append = TRUE)
  capture.output(gaocs_output, file = 'Gao_cs.txt', append = TRUE)
  
  cat('Gao tests for significant differences', file = 'Gao.txt',
      append = TRUE)
  cat(c(as.character(formula)[[2]], as.character(formula)[[3]]), 
      file = 'Gao.txt', append = TRUE)
  capture.output(gao_output, file = 'Gao.txt', append = TRUE)
  
 
}

}
setwd(init_dir)

# Plot good meanplots
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
good_plot <- function(D_agg, title, xlab, ylab)
{
  errorbar_width = (max(D_agg[1]) - min(D_agg[1])) / 100
  plot <- ggplot(D_agg, aes(x=Group, y=Mean, vjust = 1.5)) +
    geom_line(aes(group=1)) +
    geom_point(size=2) +
    geom_errorbar(aes(ymin=Low, ymax=High), 
                  width=errorbar_width) +
    xlab(paste0(xlab)) + 
    ylab(paste0(ylab)) +
    ggtitle(title) +
    scale_x_discrete(limits = D_agg$Group, expand = c(0.05,0.05)) +
    theme(plot.title = element_text(hjust = 0.5))
          #axis.text = element_text(face='plain',size=24, colour='black'),
          #text = element_text(face="bold",size=36, colour='black'),
          #axis.title.x = element_text(margin = margin(t = 20)),
          #axis.title.y = element_text(margin = margin(r = 20)))
    ylim(min(D_agg$Low),max(D_agg$High))
  return(plot)
}

plot_dir <- paste0(getwd(), '/Plots')
dir.create(plot_dir, showWarnings = FALSE)
setwd(plot_dir)
for (i in colnames(Data)[-1])
{
D_r <- as.matrix(na.exclude(Data[,c('Group', i)]))
D_agg <- try(create_summary(D_r))
goodplot <- try(good_plot(D_agg, i,'',''))
try(ggsave(paste(i, '.png'), plot = goodplot))
}
setwd(init_dir)
