# This algorithm computes KW tests and creates pictures of them
# The input is Excel file with columns as variables in mainwd dir
# The grouping variable must be named 'Group' 
# The mainwd dir may be changed
# The algorithm by default ignores columns with all missing values
#
# Copyright: Ilia Fastovets (2017)
# Contact: fastovetsilya@yandex.ru


library('readxl')
library('coin')

# Create directories and specify what de do
mainwd <- getwd()
#mainwd <- setwd("~/Dropbox/Общая/Наука/Институт/Данные/Церий/Церий в почве с растениями/ICP+Лук")
what <- 'La_Leaves'
elements <- c('Ca', 'K', 'Na', 'Mg', 'P', 'Fe', 'Mn', 'Cu', 'Zn') # Choose specific elements or comment this line
dir.create(paste0(mainwd,'/ANOVAresults_', what))
nperm <- 10^7
p_adj_method = 'BH'
sheet <- 1

# Load the data 
main_data <- read_excel('Analysis_La.xlsx', sheet = sheet)
main_data <- main_data[c('Group', elements)] # Valid only of elements are specified
# Now remove NA columns from dataframe
main_data <- main_data[colSums(!is.na(main_data))!=0]
clist <- colnames(main_data)[2:length(colnames(main_data))]
fact_levels <- levels(factor(main_data[['Group']]))

# Now compute descriptive statistics for plotting
# Create new data frame with main values (zeros)
Means <- as.data.frame(matrix(0, length(fact_levels), 
                              ncol(main_data)))
colnames(Means) <- colnames(main_data)
Means['Group'] <- as.numeric(fact_levels)
# ... and with lower confidence bounds
Lower <- as.data.frame(matrix(0, length(fact_levels), 
                              ncol(main_data)))
colnames(Lower) <- colnames(main_data)
Lower['Group'] <- as.numeric(fact_levels)
# ... and with upper confidence bounds
Upper <- as.data.frame(matrix(0, length(fact_levels), 
                              ncol(main_data)))
colnames(Upper) <- colnames(main_data)
Upper['Group'] <- as.numeric(fact_levels)

# Construct loop that computes means and CI and puts them to matrices
for (i in clist){
  
  for (k in fact_levels){
    mean <- mean(na.omit(main_data[main_data['Group']==k, i][[1]]))
    lower_ci <- t.test(na.omit(main_data[main_data['Group']==k, i][[1]]))$conf.int[1]
    upper_ci <- t.test(na.omit(main_data[main_data['Group']==k, i][[1]]))$conf.int[2]
    Means[Means['Group']==k, i] <- mean
    Lower[Lower['Group']==k, i] <- lower_ci
    Upper[Upper['Group']==k, i] <- upper_ci
  }
}

# Create table for printing

PrintTable <- function(main_data){
output <- main_data
output.mean <- as.data.frame(matrix(0,length(output[['Group']])))
output.lower <- as.data.frame(matrix(0,length(output[['Group']])))
output.upper <- as.data.frame(matrix(0,length(output[['Group']])))
output <- cbind(output, output.lower)

for (i in clist){
  colnames(output.mean) <- paste0(i, '_Mean')  
  colnames(output.lower) <- paste0(i, '_LowerCI')
  colnames(output.upper) <- paste0(i, '_UpperCI')

  target <- which(names(output) == i)[1]
  output <- cbind(output[,1:target,drop=F], output.upper, 
                output[,(target+1):length(output),drop=F])
  output <- cbind(output[,1:target,drop=F], output.lower, 
                output[,(target+1):length(output),drop=F])
  output <- cbind(output[,1:target,drop=F], output.mean, 
                output[,(target+1):length(output),drop=F])
}
output <- output[-length(output)]
return(output)
}
output <- PrintTable(main_data)

# Fill this table with means and confidence intervals

for (i in clist){  
  for (k in fact_levels){
    output.index <- max(which(main_data['Group'] == k))
    name_mean <- paste0(i, '_Mean')
    name_lower <- paste0(i, '_LowerCI')
    name_upper <- paste0(i, '_UpperCI')
    output[output.index, name_mean] <- Means[Means['Group']==k, i]
    output[output.index, name_lower] <- Lower[Lower['Group']==k, i]
    output[output.index, name_upper] <- Upper[Upper['Group']==k, i]
  }
}

# Now remove extra zeros from the output
for (i in clist){
  name_mean <- paste0(i, '_Mean')
  name_lower <- paste0(i, '_LowerCI')
  name_upper <- paste0(i, '_UpperCI')
  output[which(output[name_mean]==0), name_mean] <- NA
  output[which(output[name_lower]==0), name_lower] <- NA
  output[which(output[name_upper]==0), name_upper] <- NA
}

# Now save the file with the table
write.csv(output, file = paste0(mainwd,'/ANOVAresults_', what, 
                                '/Descriptives output_', what, '.csv'), row.names = F, na = '')

#Create loop that performs KW test and adjust p-values
pval_raw <- numeric()
pval_adj <- numeric()
for (i in clist){
  tst <- kruskal_test(main_data[[i]]~factor(main_data[['Group']]),
                      distribution = approximate(B=nperm))
  pval <- tst@distribution@pvalue(tst@statistic@teststatistic)[1]
  pval_raw[i] <- pval
  pval_adj[i] <- pval 
}
pval_adj <- p.adjust(pval_adj, method = p_adj_method)

#creates pic with KW p-val
for (i in clist){
  CI_dn <- Lower[[i]]
  CI_up <- Upper[[i]]
  filename <- paste0(mainwd,'/ANOVAresults_', what, '/Plots_', i, '.png')
  png()
  png(filename=filename)
  plot(Means[[i]]~Means[['Group']],
       main = i, xlab = 'Groups', ylab = paste0('Concentration of ', i, ' in ',what),
       type='b', ylim = c(min(CI_dn), max(CI_up)))
  arrows(Means[['Group']], CI_dn, Means[['Group']], CI_up,
         code = 3, length = 0.05, angle = 90)
  mtext(paste0('P values adjusted using: ', p_adj_method), 3, line = -2, col = 'green')
  mtext(paste0('P raw: ', pval_raw[[i]]), 3, line = -4, col = 'black')

    if (pval_adj[[i]] < 0.05){
    mtext(paste0('P-value=', pval_adj[[i]]), 3, line = -3, col = 'red')
    dev.off()
    }
  else{
    mtext(paste0('P-value=', pval_adj[[i]]), 3, line = -3, col='black')
    dev.off()
  }
  dev.off()
}






