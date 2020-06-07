# Load library
library('readxl')

### BEGIN INPUT
r.crit <- 0.3 # Set critical r
r.method <- 'pearson' # Method of correlation. "pearson", "kendall", "spearman"
compute_ranks <- 'simple' # Rank matrices to be computed. 'simple' or 'weighted'
#setwd("~/Dropbox/Общая/Наука/Институт/Формулы сбалансированности") # Set working directory 
Data <- read_excel('SampleData.xlsx', sheet = 1) # Data file to read
dep_var <- 'us' # Dependent variable to correlate with
cmb_data <- Data[4:11] # Set independent variables
### END OF INPUT

cmb <- combn(ncol(cmb_data), 2)
colnames(cmb) <- paste0(names(cmb_data)[cmb[1,]], 
                        '/', names(cmb_data)[cmb[2,]])
ratios <- apply(cmb, 2, function(j) as.matrix(cmb_data[, j[1]])/
              as.matrix(cmb_data[, j[2]]))
cormat <- cor(ratios, Data[,dep_var], method = r.method)
rowrank_mat <- matrix(0, nrow = 2, ncol = length(cmb_data))
colnames(rowrank_mat) <- names(cmb_data)
cormat2 <- cbind(matrix(unlist(strsplit(rownames(cormat), '/')),
       byrow = T, ncol = 2), cormat)
colnames(cormat2) <- c('num', 'den', 'cor')

if (compute_ranks == 'simple')
{
  
for(i in 1:length(cormat))
{
  sign_temp <- sign(as.numeric(cormat2[i,3]))
  num_temp <- cormat2[i,1]
  den_temp <- cormat2[i,2]
  
  if (sign_temp < 0 & abs(as.numeric(cormat2[i,3])) >= r.crit) #crit rho
{
    rowrank_mat[2,num_temp] <- rowrank_mat[2,num_temp] + 1
    rowrank_mat[1,den_temp] <- rowrank_mat[1,den_temp] + 1
}
  if (sign_temp > 0 & abs(as.numeric(cormat2[i,3])) >= r.crit) #crit rho
{
    rowrank_mat[1,num_temp] <- rowrank_mat[1,num_temp] + 1
    rowrank_mat[2,den_temp] <- rowrank_mat[2,den_temp] + 1
}
}
  
}

if (compute_ranks == 'weighted')
{
  
for(i in 1:length(cormat))
{
  cor_temp <- as.numeric(cormat2[i,3])
  sign_temp <- sign(as.numeric(cormat2[i,3]))
  num_temp <- cormat2[i,1]
  den_temp <- cormat2[i,2]
    
  if (sign_temp < 0 & abs(as.numeric(cormat2[i,3])) >= r.crit) #crit rho
{
    rowrank_mat[2,num_temp] <- rowrank_mat[2,num_temp] + abs(cor_temp^2)           
    rowrank_mat[1,den_temp] <- rowrank_mat[1,den_temp] + abs(cor_temp^2)
}
    if (sign_temp > 0 & abs(as.numeric(cormat2[i,3])) >= r.crit) #crit rho
{
    rowrank_mat[1,num_temp] <- rowrank_mat[1,num_temp] + abs(cor_temp^2)
    rowrank_mat[2,den_temp] <- rowrank_mat[2,den_temp] + abs(cor_temp^2)
}
} 

}

def_rank_mat <- rowrank_mat
def_rank_mat[1,] <- -def_rank_mat[1,] + length(cmb_data)
def_rank_mat[2,] <- -def_rank_mat[2,] + length(cmb_data)
rel_rank_mat <- def_rank_mat[1,] / def_rank_mat[2,]
print(c(compute_ranks, r.method, r.crit))
print('The table of frequencies is:', quote = F)
print(rowrank_mat)
print('The table of ranks is:', quote = F)
print(def_rank_mat)
print('The relative ranks (balance formula) are:', quote = F)
print(sort(rel_rank_mat))
