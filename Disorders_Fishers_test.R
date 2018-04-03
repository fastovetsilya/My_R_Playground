# Statistical analysis of mitotic index and disorders in cells based on
#   Fisher's exact test an Clopper-Pearson exact confidence intervals.
# The input is an .xlsx data matrix with any number of raws/columns.
#   However, the first 3 columns must be 'Group', 'Total', 'Dividing'.
#   The data is raw numbers of cells by replications. The data is 
#   grouped automatically by the algorithm.
# The output are plots of MI and disorders (individual and grouped)
#   and tables of proportions, confidence intervals and p-values. The
#   algorithms outputs raw p-values as well as adjusted p-values. By
#   default, Benjamini-Hochberg & Yekutielli p-value adjustment is used.
# This program is developed for citotoxicity/genotoxicity analysis of
#   allium root tip cells (allium-test).
#
# Copyright: Ilia Fastovets (2018)
# Contact: fastovetsilya@yandex.ru


library('readxl')
library('ggplot2')
library('reshape2')
library('binom')

Data <- read_excel('Overall_data.xlsx', sheet = 5)
D_r <- as.matrix(aggregate(. ~ Group, data = Data, sum))

#///DEFINE FUNCTIONS///
create_MI <- function(D_r)
{
  D_rr <- D_r[, c(1,2,3)]
  MI <- binom.confint(D_rr[,3], D_rr[,2], method='exact')[, -c(1,2,3)]
  MI <- cbind(D_r[,1], MI)
  colnames(MI) <- c('Group', 'Mean', 'Low', 'High')
  return(MI)
}

MI_by_group_p <- function(D_r)
{
  D_rr <- D_r[, c(1,2,3)]
  mat <- matrix(nrow=1, ncol=(nrow(D_rr)-1))
  for(i in 1:(nrow(D_r)-1))
  {
    comp_table <- rbind(D_rr[1,c(2,3)], D_rr[(i+1),c(2,3)])
    comp_table[,1] <- comp_table[,1] - comp_table[,2]
    mat[1,i] <- fisher.test(comp_table)$p.value
    colnames(mat) <- D_rr[,1][-1]
  }
  return(mat)
}

create_MI_comb_p <- function(D_r)
{
  D_rr <- D_r[, c(1,2,3)]
  comb <- combn(D_rr[,'Group'], 2)
  mat <- matrix(nrow=1, ncol=ncol(comb))
  for(i in 1:ncol(mat))
  {
    comp_table <- rbind(D_rr[D_rr[,1] == comb[1,i],], D_rr[D_r[,1] == comb[2,i],])[,-1]
    comp_table[,1] <- comp_table[,1] - comp_table[,2]
    mat[1,i] <- fisher.test(comp_table)$p.value
  }
  mat <- rbind(comb, mat)
  mat[c(1,2),] <- as.character(mat[c(1,2),])
  return(mat)
}

adjust_comb_p <- function(comb_p, method = 'none')
{
  comb_p_adj <- comb_p
  comb_p_adj[3,] <- p.adjust(comb_p[3,], method=method)
  rownames(comb_p_adj) <- c('Group1', 'Group2', 'adj p-val')
  return(comb_p_adj)
}

create_proportions <- function(D_r)
{
  Proportions <- matrix(0, nrow = nrow(D_r), ncol = ncol(D_r)-3)
  for(i in seq(4, ncol(D_r)))
  {
    Proportions[, i-3] <- D_r[, i] / D_r[, 'Dividing']
  }
  Proportions <- cbind(D_r[, 'Group'], Proportions)
  colnames(Proportions) <- colnames(D_r)[c(1, 4:ncol(D_r))]
  return(Proportions)
}

create_Lower_CI <- function(D_r)
{
  Lower_CI <- matrix(0, nrow = nrow(D_r), ncol = ncol(D_r)-3)
  for(i in seq(4, ncol(D_r)))
  {
    Lower_CI[, i-3] <- binom.confint(D_r[, i], D_r[, 'Dividing'], 
                                     method = 'exact')$lower
  }
  Lower_CI <- cbind(D_r[, 'Group'], Lower_CI)
  colnames(Lower_CI) <- colnames(D_r)[c(1, 4:ncol(D_r))]
  Lower_CI[Lower_CI < 0]
  return(Lower_CI)
}

create_Upper_CI <- function(D_r)
{
  Upper_CI <- matrix(0, nrow = nrow(D_r), ncol = ncol(D_r)-3)
  for(i in seq(4, ncol(D_r)))
  {
    Upper_CI[, i-3] <- binom.confint(D_r[, i], D_r[, 'Dividing'], 
                                     method = 'exact')$upper
  }
  Upper_CI <- cbind(D_r[, 'Group'], Upper_CI)
  colnames(Upper_CI) <- colnames(D_r)[c(1, 4:ncol(D_r))]
  return(Upper_CI)
}

create_by_d_prop <- function(D_r)
{
  by_g_prop <- binom.confint(colSums(D_r[, 4:ncol(D_r)]), sum(D_r[,3]),
                             method = 'exact')[-c(1, 2, 3)]
  rownames(by_g_prop) <- colnames(D_r)[-c(1, 2, 3)]
  by_g_prop <- cbind(rownames(by_g_prop), by_g_prop)
  colnames(by_g_prop) <- c('Group', 'Mean', 'Low', 'High')
  return(by_g_prop)
}

create_by_g_prop <- function(D_r)
{
  D_r_t <- t(D_r)
  colnames(D_r_t) <- c(D_r_t['Group',])
  
  by_g_prop <- binom.confint(colSums(D_r_t[-c(1,2,3),]), D_r_t[3,],
                             method = 'exact')[-c(1, 2, 3)]
  
  rownames(by_g_prop) <-  D_r[, "Group"]
  by_g_prop <- cbind(rownames(by_g_prop), by_g_prop)
  colnames(by_g_prop) <- c('Group', 'Mean', 'Low', 'High')
  return(by_g_prop)
}

adjust_by_disorder_p <- function(by_disorder_p, method = 'none')
{
  m <- matrix(0, ncol = ncol(by_disorder_p), nrow = 1)
  by_group_p_adj <- p.adjust(by_disorder_p, method = method)
  m[1, ] <- by_group_p_adj
  colnames(m) <- colnames(by_disorder_p)
  return(m)
}

create_by_group_p <- function(D_r)
{
  num <- cbind(D_r[,c(1,2,3)], rowSums(D_r[,-c(1,2,3)]))
  colnames(num) <- c(colnames(D_r[,c(1,2,3)]), 'No Disorders')
  mat <- matrix(nrow=1, ncol=(nrow(D_r)-1))
  colnames(mat) <- D_r[,'Group'][-1]
  
  for(i in 1:ncol(mat))
  {
    sample <- num[c(1,(i+1)),-c(1, 2)]
    sample[,1] <- sample[,1] - sample[,2]
    mat[1,i] <- fisher.test(sample)$p.value
  }
  return(mat)
}

create_combinations_p <- function(D_r)
{
  comb <- combn(D_r[,'Group'], 2)
  for(i in 1:ncol(comb))
  {
    comb <- rbind(comb, matrix(nrow=1, ncol=ncol(comb)))
    comparemat <- rbind(D_r[D_r[, 'Group'] == comb[1,i]], 
                        D_r[D_r[, 'Group'] == comb[2,i]])
    colnames(comparemat) <- colnames(D_r)
    summat <- cbind(comparemat[,c(1,2,3)], rowSums(comparemat[, -c(1, 2, 3)]))
    colnames(summat) <- c(colnames(comparemat[,c(1,2,3)]), 'No_disorders')
    comp_table <- summat[, c(3, 4)]
    comp_table[,1] <- comp_table[,1] - comp_table[,2]
    comb[3,i] <- fisher.test(comp_table)$p.value
  }
  comb[c(1,2),] <- as.character(comb[c(1,2),])
  comb <- na.exclude(comb)
  rownames(comb) <- c('Group1', 'Group2', 'raw p-val')
  return(comb)
}

create_by_dis_combinations_p <- function(D_r)
{
  comb <- combn(colnames(D_r)[-c(1,2,3)], 2)
  for(i in 1:ncol(comb))
  {
    comb <- rbind(comb, matrix(nrow=1, ncol=ncol(comb)))
    comparemat <- rbind(D_r[,comb[1,i]], 
                        D_r[,comb[2,i]])
    
    comp_table <- cbind(sum(D_r[, 3]) - rowSums(comparemat), rowSums(comparemat))
    comb[3,i] <- fisher.test(comp_table)$p.value
  }
  comb[c(1,2),] <- as.character(comb[c(1,2),])
  comb <- na.exclude(comb)
  rownames(comb) <- c('Group1', 'Group2', 'raw p-val')
  return(comb)
}

create_individ_p_values <- function(Proportions)
{
  p_values <- matrix(0, nrow = nrow(Proportions)-1, ncol = ncol(Proportions)-1)
  rownames(p_values) <- D_r[,'Group'][-1]
  colnames(p_values) <- colnames(D_r)[4:ncol(D_r)]
  for(i in as.character(D_r[,'Group'][-1]))
  {
    for(k in 1:ncol(p_values))
    {
      table <- rbind(cbind(D_r[,'Dividing'][D_r[, 'Group'] == '0'], D_r[, k+3][D_r[, 'Group'] == '0']),
                     cbind(D_r[,'Dividing'][D_r[, 'Group'] == i], D_r[, k+3][D_r[, 'Group'] == i]))
      table[,1] <- table[,1] - table[,2]
      table <- as.table(table)
      p <- fisher.test(table)$p.value
      p_values[i, k] <- p
    }
  }
  return(p_values)
}

adjust_p_values <- function(p_values, method = 'none')
{
  p_val_vect <- melt(p_values, id.vars = 'Disorder')[3][[1]]
  p_val_adj_vect <- p.adjust(p_val_vect, method = method)
  p_val_adj <- matrix(p_val_adj_vect, nrow = nrow(p_values), ncol = ncol(p_values))
  rownames(p_val_adj) <- rownames(p_values)
  colnames(p_val_adj) <- colnames(p_values)
  return(p_val_adj)
}

good_plot <- function(D_agg, title, xlab, ylab)
{
  ggplot(D_agg, aes(x=Group, y=Mean)) +
    geom_line(aes(group=1)) +
    geom_point(size=2) +
    geom_errorbar(aes(ymin=Low, ymax=High), position = 'dodge') +
    xlab(xlab) + 
    ylab(ylab) +
    ggtitle(title) +
    scale_x_discrete(limits = D_agg$Group, expand = c(0.05,0.05)) +
    theme(plot.title = element_text(hjust = 0.5), 
          axis.text = element_text(face='plain',size=14, colour='black'),
          text = element_text(face="bold",size=14, colour='black'))
}

write_output <- function()
{
  write.table('//////////FISHER EXACT TESTS & CLOPPER_PEARSON INTERVALS//////////', 
              file = 'Output.csv', row.names = FALSE, col.names = FALSE)
  write.table('//////////MITOTIC INDEX//////////', 
              file = 'Output.csv', row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(matrix(nrow = 3, ncol = 5), 
              file = 'Output.csv', na = '', col.names = FALSE, row.names = FALSE, append = TRUE)
  
  write.table('MI(+- CLOPPER PEARSON CONFIDENCE INTERVAL)', file = 'Output.csv', row.names = FALSE, col.names = FALSE, append=TRUE)
  write.table(rbind(MI, matrix(nrow = 3, ncol = ncol(MI))),
              file = 'Output.csv', na = '', row.names = FALSE, append = TRUE)
  
  
  write.table('RAW P-VAL MI COMPARISONS WITH CONTROL', file = 'Output.csv', row.names = FALSE, append=TRUE)
  write.table(rbind(MI_by_g_p, matrix(nrow = 3, ncol = ncol(MI_by_g_p))),
              file = 'Output.csv', na = '', row.names = FALSE, append = TRUE)
  
  write.table('RAW P-VAL MI PAIRWISE COMPARISONS', file = 'Output.csv', row.names = FALSE, col.names = FALSE, append=TRUE)
  write.table(rbind(MI_comb_p, matrix(nrow = 3, ncol = ncol(MI_comb_p))),
              file = 'Output.csv', na = '', row.names = FALSE, col.names = FALSE, append = TRUE)
  
  write.table('ADJ P-VAL MI COMPARISONS WITH CONTROL, BENJAMINI-HOCHBERG-YEKUTIELI ADJUSTMENT',
              file = 'Output.csv', row.names = FALSE, col.names = FALSE, append=TRUE)
  write.table(rbind(MI_by_g_p_adj, matrix(nrow = 3, ncol = ncol(MI_by_g_p_adj))),
              file = 'Output.csv', na = '', row.names = FALSE, append = TRUE)
  
  write.table('ADJ P-VAL MI PAIRWISE COMPARISONS, BENJAMINI-HOCHBERG-YEKUTIELI ADJUSTMENT', 
              file = 'Output.csv', row.names = FALSE, col.names = FALSE, append=TRUE)
  write.table(rbind(MI_comb_p_adj, matrix(nrow = 3, ncol = ncol(MI_comb_p_adj))),
              file = 'Output.csv', na = '', row.names = FALSE, col.names = FALSE, append = TRUE)
  
  write.table(matrix(nrow = 3, ncol = 5), 
              file = 'Output.csv', na = '', col.names = FALSE, row.names = FALSE, append = TRUE)
  write.table('//////////FAC//////////', 
              file = 'Output.csv', row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(matrix(nrow = 3, ncol = 5), 
              file = 'Output.csv', na = '',col.names = FALSE, row.names = FALSE, append = TRUE)
  
  write.table('DATA TABLE', file = 'Output.csv', row.names = FALSE, col.names = FALSE, append=TRUE)
  write.table(rbind(D_r, matrix(nrow = 3, ncol = ncol(D_r))),
              file = 'Output.csv', na = '', row.names = FALSE, append = TRUE)  
  
  write.table('INDIVIDUAL FAC', file = 'Output.csv', row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(rbind(Proportions, matrix(nrow = 3, ncol = ncol(Proportions))),
              file = 'Output.csv', na = '', row.names = FALSE, append = TRUE)
  
  write.table(matrix(nrow = 3, ncol = 5), 
              file = 'Output.csv', na = '', col.names = FALSE, row.names = FALSE, append = TRUE)
  write.table('CONFIDENCE INTERVALS', 
              file = 'Output.csv', row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(matrix(nrow = 3, ncol = 5), 
              file = 'Output.csv', na = '',col.names = FALSE, row.names = FALSE, append = TRUE)
  
  
  write.table('OVERALL CLOPPER_PEARSON 95% CI', 
              file = 'Output.csv', row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(rbind(Plot_overall, matrix(nrow = 3, ncol = ncol(Plot_overall))), 
              file = 'Output.csv', na = '', row.names = FALSE, append = TRUE)
  
  write.table('BY-DISORDER CLOPPER_PEARSON 95% CI', 
              file = 'Output.csv', row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(rbind(by_d_prop, matrix(nrow = 3, ncol = ncol(by_d_prop))), 
              file = 'Output.csv', na = '', row.names = FALSE, append = TRUE)
  
  write.table('BY-CONCENTRATION CLOPPER_PEARSON 95% CI', 
              file = 'Output.csv', row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(rbind(by_g_prop, matrix(nrow = 3, ncol = ncol(by_g_prop))), 
              file = 'Output.csv', na = '', row.names = FALSE, append = TRUE)
  
  
  write.table(matrix(nrow = 3, ncol = 5), 
              file = 'Output.csv', na = '', col.names = FALSE, row.names = FALSE, append = TRUE)
  write.table('COMPARISONS', 
              file = 'Output.csv', row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(matrix(nrow = 3, ncol = 5), 
              file = 'Output.csv', na = '',col.names = FALSE, row.names = FALSE, append = TRUE)
  
  
  write.table('RAW INDIVIDUAL P-VALUES (COMPARISONS TO CONTROL GROUP)', 
              file = 'Output.csv', row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(rbind(p_values, matrix(nrow = 3, ncol = ncol(p_values))), 
              file = 'Output.csv', na = '', append = TRUE)


  write.table('RAW BY-GROUP P_VAL (COMPARISONS TO CONTROL GROUP)', file = 'Output.csv', row.names = FALSE, 
              col.names = FALSE, append = TRUE)
  write.table(rbind(by_g_p, matrix(nrow = 3, ncol = ncol(by_g_p))),
              file = 'Output.csv', na = '', row.names = FALSE, append = TRUE)
  
  write.table('RAW PAIRWISE BETWEEN-DISORDERS P-VALUES', 
              file = 'Output.csv', row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(rbind(by_dis_combinations_p, matrix(nrow = 3, ncol = ncol(by_dis_combinations_p))), 
              file = 'Output.csv', na = '', col.names = FALSE, row.names = FALSE, append = TRUE)
  
  write.table('RAW PAIRWISE BETWEEN-CONCENTRATIONS P-VALUES', 
              file = 'Output.csv', row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(rbind(comb_p, matrix(nrow = 3, ncol = ncol(comb_p))), 
              file = 'Output.csv', na = '', col.names = FALSE, row.names = FALSE, append = TRUE)
  
  write.table('ADJ INDIVIDUAL P-VALUES (COMPARISONS TO CONTROL GROUP): 
              BENJAMINI-HOCHBERG-YEKUTIELI ADJUSTMENT', 
              file = 'Output.csv', row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(rbind(p_val_adj, matrix(nrow = 3, ncol = ncol(p_val_adj))), 
              file = 'Output.csv', na = '', append = TRUE)
  
  write.table('ADJ BY-GROUP P-VALUES (COMPARISONS TO CONTROL GROUP): 
            BENJAMINI-HOCHBERG-YEKUTIELI ADJUSTMENT', file = 'Output.csv', row.names = FALSE, 
              col.names = FALSE, append = TRUE)
  write.table(rbind(by_g_p_adj, matrix(nrow = 3, ncol = ncol(by_g_p_adj))),
              file = 'Output.csv', na = '', row.names = FALSE, append = TRUE)
  
  write.table('ADJ PAIRWISE BETWEEN-CONCENTRATIONS ADJ P-VALUES, 
            BENJAMINI-HOCHBERG-YEKUTIELI ADJUSTMENT', 
              file = 'Output.csv', row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(rbind(comb_p_adj, matrix(nrow = 3, ncol = ncol(comb_p_adj))), 
              file = 'Output.csv', na = '', col.names = FALSE, row.names = FALSE, append = TRUE)
  
  write.table('ADJ PAIRWISE BETWEEN-DISORDERS P-VALUES, 
            BENJAMINI-HOCHBERG-YEKUTIELI ADJUSTMENT', 
              file = 'Output.csv', row.names = FALSE, col.names = FALSE, append = TRUE)
  write.table(rbind(by_dis_combinations_p_adj, matrix(nrow = 3, ncol = ncol(by_dis_combinations_p_adj))), 
              file = 'Output.csv', na = '', col.names = FALSE, row.names = FALSE, append = TRUE)
}


#///MI///


MI <- create_MI(D_r)

MI_by_g_p <- MI_by_group_p(D_r)
MI_by_g_p_adj <- MI_by_g_p
MI_by_g_p_adj[1,] <- p.adjust(MI_by_g_p, method = 'BY')

MI_comb_p <- create_MI_comb_p(D_r)

MI_comb_p_adj <- adjust_comb_p(MI_comb_p, method='BY')


#///FAC///
# Creating proportion matrices and tables for plotting/writing

Proportions <- create_proportions(D_r)
Lower_CI <- create_Lower_CI(D_r)
Upper_CI <- create_Upper_CI(D_r)

Proportions <- as.data.frame(Proportions)
Plot_table <- melt(Proportions, id.vars='Group')
colnames(Plot_table) <- c('Group', 'Disorder', 'Proportion')

Lower_CI <- as.data.frame(Lower_CI)
Plot_LCI <- melt(Lower_CI, id.vars='Group')
colnames(Plot_LCI) <- c('Group', 'Disorder', 'Proportion')

Upper_CI <- as.data.frame(Upper_CI)
Plot_UCI <- melt(Upper_CI, id.vars='Group')
colnames(Plot_UCI) <- c('Group', 'Disorder', 'Proportion')

Plot_overall <- as.matrix(cbind(Plot_table, Plot_LCI[,3], Plot_UCI[, 3]))
colnames(Plot_overall) <- c('Group', 'Disorder', 'Proportion', 
                            'Lower 95% CI', 'Upper 95% CI')

by_d_prop <- create_by_d_prop(D_r)
by_g_prop <- create_by_g_prop(D_r)

# Create tables with p-values

by_g_p <- create_by_group_p(D_r)
by_g_p_adj <- adjust_by_disorder_p(by_g_p, method = 'BY')

comb_p <- create_combinations_p(D_r)
comb_p_adj <- adjust_comb_p(comb_p, method = 'BY')

by_dis_combinations_p <- create_by_dis_combinations_p(D_r)
by_dis_combinations_p_adj <- adjust_comb_p(by_dis_combinations_p, method = 'BY')


p_values <- create_individ_p_values(Proportions)
p_val_adj <- adjust_p_values(p_values, method = 'BY') 

# Plotting data
good_plot(MI, 'Mitotic index (+- Clopper-Pearson exact confidence intervals)', 
          'Groups of concentrations', 'MI')
ggplot(Plot_table, aes(fill=Disorder, y=Proportion, x=factor(Group))) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin=Plot_LCI[,3], ymax=Plot_UCI[,3]), position = 'dodge') +
  xlab('Group') + 
  ylab('Proportion') 

good_plot(by_d_prop, 'By disorder plot (+- Clopper-Pearson exact confidence intervals)', 
          'Groups of disorders', 'Proportions of disorders')
good_plot(by_g_prop, 'By groups plot (+- Clopper-Pearson exact confidence intervals)', 
          'Groups of concentrations', 'Proportions of disorders')

# Writing tables
MI <- as.matrix(MI)
Proportions <- as.matrix(Proportions)
by_d_prop <- as.matrix(by_d_prop)
by_g_prop <- as.matrix(by_g_prop)

write_output()

 

