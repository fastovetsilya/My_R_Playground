library('readxl')
library('ggplot2')
library('RColorBrewer')

Data <- read_excel('Data.xlsx', sheet = 5)
adjust_method <- 'BH'

# Plotting
cols <- brewer.pal(8, 'Set2')

ggplot(Data, aes(x=factor(Fraction), y=U))+
  geom_boxplot(fill = cols)+
  xlab('Fraction')+
  ylab('U')


# Comparisons
fr1 <- c(1,2,5,6,1,2,3,4)
fr2 <- c(3,4,7,8,5,6,7,8)
mat <- matrix(0, length(fr1), 9)
colnames(mat) <- c('FF_p_raw', 'Rdn_p_raw', 'U_p_raw',
                   'FF_p_adj', 'Rdn_p_adj', 'U_p_adj',
                   'FF_signif', 'Rdn_signif', 'U_signif')
rownames(mat) <- paste0(fr1, '-',  fr2)


for (k in 1:3)
{
for (i in 1:length(fr1))
{
D_r <- rbind(Data[Data[,1] == fr1[i], c(1,(k+1))],
             Data[Data[,1] == fr2[i], c(1,(k+1))])
colnames(D_r)[2] <- 'Param'
mat[i,k] <- wilcox.test(Param~Fraction, D_r)$p.val
} 
}

for (i in 1:3)
{
mat[,(i+3)] <- p.adjust(mat[,i], method = adjust_method)
}

for (i in 1:3)
{
mat[,(i+6)] <- as.numeric(mat[,(i+3)] <= 0.05)
mat[,(i+6)][mat[,(i+6)] == 1] <- 'Yes'
mat[,(i+6)][mat[,(i+6)] == 0] <- 'No'
}

write.table(mat, file = 'output.csv')






