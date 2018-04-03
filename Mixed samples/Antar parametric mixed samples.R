library('readxl')
# Set data
Data <- read_excel('Data_Antar.xlsx', sheet = 2)
Data_mixed <- Data[1:36,]
Data_ind <- Data [37:42,]
Ftr <- 'Car'
alpha <- 0.05
# Process individual samples
std <- sd(as.matrix(Data_ind[Ftr]))
n <- length(as.matrix(Data_ind[Ftr]))
std_mc <- sqrt((n-1)*std^2/qchisq(1-0.05,n-1,lower.tail = FALSE))
# Process comparisons of mixed samples according 
# to z-test with sminimum credence principle
n_mix <- 2
cmb <- combn(nrow(Data_mixed), 2)
colnames(cmb) <- paste0(Data_mixed$Index[cmb[1,]], 
                        '--', Data_mixed$Index[cmb[2,]])
z_test <- abs(as.matrix(apply(cmb, 2, function(j) as.matrix(Data_mixed[j[1], Ftr]) -
             as.matrix(Data_mixed[j[2], Ftr]))))/std_mc/sqrt(2/n_mix)
p_val_row <- pnorm(z_test, mean = 0, sd = 1, lower.tail = FALSE)
p_val_adj <- p.adjust(p_val_row, method = 'BY')
alpha <- matrix(alpha, length(cmb[1,]), 1)
signif <- matrix('Not significant',length(cmb[1,]),1)
signif[p_val_adj <= alpha] <- 'Significant'
z_test <- cbind(z_test, p_val_row, p_val_adj, signif)
colnames(z_test) <- c('z-val', 'p-val_row', 'p-val-adj (BY)', 'p<=0.05')
print('Parametric Minimum Credence (MC) comparisons of mixed samples', quote = F)
print('MC estimate of sigma one-tailed confidence bound was found with chi-suared distibution from the reference sample (n = 6)', quote = F)
print('Normal z-test was then used to compare samples under MC sigma', quote = F)
print('Benjamini, Hochberg, and Yekutieli p-value adjustment was used', quote = F)
print('The number of Significant results is:', quote = F)
print(length(signif [signif == 'Significant']))
#write.table(z_test, file = 'Car.csv')
print('File with the results was created in the home directory')
print('MC standard deviation estimate is (mg/100g):', quote = F)
print(std_mc)



