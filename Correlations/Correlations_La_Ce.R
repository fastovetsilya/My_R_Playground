# This program performs correlation analysis of specified features.
# It was developed for performing correlation analysis for Mitotic
#   index (MI) and frequency of abberant cells (FAC) of onion root
#   tip cells (Allium-test) and length of onion roots. 
# For nonlinear relationships, this algorithm computes permutation
#   test for Spearman rank-based correlations.
# The input is a file with columns named with structure 
#   [Element]_[Medium]_[Feature] as in the example. 
# The output is data table with adjusted p-values (BHY adjustment)
#   and the plots of correlations with numerical summaries and
#   LOWESS smoothing curve for assessing linearity of the 
#   relationship.
#
# Copyright: Ilia Fastovets (2018)
# Contact: fastovetsilya@yandex.ru


library('readxl')
library('sm')

Data <- read_excel('Correlations_La_Ce_rastv_pochv.xlsx', sheet = 1)
Datamat <- as.matrix(Data)

# Permutation test for Spearman correlation function
Permfun <- function(x_perm, y_perm, N_perm = 10^5)
{
  N_perm <- N_perm
  n_sample = length(y_perm)
  observed <- cor(x_perm, y_perm, method = 'spearman')
  perm_result <- numeric()
  for (i in 1:N_perm)
  {
    index <- sample(n_sample, replace=FALSE)
    Short.permuted <- y_perm[index]
    perm_result[i] <- cor(Short.permuted, x_perm)
  }
  #hist(perm_result)
  #abline(v = observed, col = 'blue')
  if (observed < 0)
  {
    p_val <- (sum(observed >= perm_result) + 1) / (N_perm + 1) * 2
    return(p_val)
  }
  
  if (observed > 0)
  {
    p_val <- (sum(observed <= perm_result) + 1) / (N_perm + 1) * 2
    return(p_val)
  }
  
  if (observed == 0)
  {
    print('Correlation is 0')
  }
}

PlotData_2 <- function(Datamat, output_mat, save = FALSE)
{ 
  for(i in 1:nrow(output_mat))
{
    Index_Abscissa <- as.character(output_mat[i,1])
    Index_Feature <- as.character(output_mat[i,2])
    spearman_cor <- as.numeric(output_mat[i,3])
    spearman_p <- as.numeric(output_mat[i,4])
    D_Ce_r <- na.exclude(Datamat[,c(Index_Abscissa, Index_Feature)])
    pearson_cor <- as.numeric(output_mat[i,5])
    pearson_rsq <- as.numeric(output_mat[i,6])
    ftest_p <- as.numeric(output_mat[i,7])
    
    if(save == TRUE){
    png(filename = paste0(i,Index_Feature, '_', Index_Abscissa, '.png'))
}
    
    sm.regression(D_Ce_r[,Index_Abscissa], D_Ce_r[,Index_Feature], col = 'red',
                  xlab = Index_Abscissa, ylab = Index_Feature, method = 'cv', lty = 2)
    
    abline(lm(D_Ce_r[,Index_Feature] ~ D_Ce_r[,Index_Abscissa]), col = 'blue')
    
    mtext(paste0('Spearman S: ', spearman_cor), line = -1, col = 'blue', adj = 0)
    mtext(paste0('Pearson R: ', pearson_cor), line = -3, col = 'blue', adj = 0)
    mtext(paste0('R-squared: ', pearson_rsq), line = -4, col = 'blue', adj = 0)
    
    if(ftest_p <= 0.05)
    {
      mtext(paste0('adj p-val (F-test, linear!): ', ftest_p), 
            line = -5, col = 'red', adj = 0)
    }
    
    if(ftest_p > 0.05)
    {
      mtext(paste0('adj p-val (F-test, linear!): ', ftest_p), 
            line = -5, col = 'blue', adj = 0)
    }
    
    if(spearman_p <= 0.05)
    {
      mtext(paste0('adj p-val (Spearman): ', signif(spearman_p, digits = 4)),
            line = -2, col = 'red', adj = 0)
    }
    
    if(spearman_p > 0.05)
    {
      mtext(paste0('adj p-val (Spearman): ', signif(spearman_p, digits = 4)),
            line = -2, col = 'blue', adj = 0)
    }
    legend('topright', c('Linear', 'LOESS'), col = c('blue', 'red'), lty = c(1,2))
    
    if(save == TRUE){
      dev.off()
}
    
    
}
}

output_mat <- matrix(nrow = nrow(expand.grid(colnames(Datamat), colnames(Datamat))),
              ncol = 7)
i <- 1
for(Element in c('La', 'Ce')) 
{ 
for(Medium in c('rastv', 'pochv')) 
{
for(Abscissa in c('RootLength', 'FAC'))
{

for(Feature in c('MI', 'FAC')) 
{

colnames(output_mat) <- c('X', 'Y', 'Spearman S', 'adj sim. p-value (S)', 'Pearson R', 'R-squared', 'adj p-val (R)')
Index_Abscissa <- paste(Element, Medium, Abscissa, sep = '_')
Index_Feature <- paste(Element, Medium, Feature, sep = '_')

D_Ce_r <- na.exclude(Datamat[, c(Index_Abscissa, Index_Feature)])
output_mat[i, 1] <-  colnames(D_Ce_r)[1]
output_mat[i, 2] <-  colnames(D_Ce_r)[2]

spearman_cor <- signif(cor(cbind(D_Ce_r[,Index_Feature], D_Ce_r[,Index_Abscissa]),
                           method = 'spearman')[1,2], digits = 3)
output_mat[i, 3] <- spearman_cor
spearman_p <- Permfun(D_Ce_r[,Index_Abscissa], D_Ce_r[,Index_Feature])
output_mat[i, 4] <- spearman_p
pearson_cor <- signif(cor(cbind(D_Ce_r[,Index_Feature], D_Ce_r[,Index_Abscissa]),
                           method = 'pearson')[1,2], digits = 3)
output_mat[i, 5] <- pearson_cor
output_mat[i, 6] <- signif(pearson_cor^2, digits = 3)
ftest_p <- summary(lm(D_Ce_r[,Index_Feature] ~ D_Ce_r[,Index_Abscissa]))[[4]][8]
output_mat[i, 7] <-ftest_p

i <- i + 1
}
}
}
}
output_mat <- na.exclude(output_mat)
output_mat <- output_mat[output_mat[,6] != 1, ]
output_mat[,4] <- signif(p.adjust(output_mat[,4], method = 'BY'), 4)
output_mat[,7] <- signif(p.adjust(output_mat[,7], method = 'BY'), 4)

# spear_cor_mat <- cor(Data, method = 'spearman', use = 'pairwise.complete.obs')
write.table('BENJAMINI, HOCHBERG & YEKUTIELLI P-VALUE ADJUSTMENT WAS USED',
            file = 'Spearman correlations output.csv', row.names = FALSE, col.names = FALSE)
write.table(output_mat, file = 'Spearman correlations output.csv', row.names = FALSE, append = TRUE) 


PlotData_2(Datamat, output_mat, TRUE)



