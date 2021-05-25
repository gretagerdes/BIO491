library("readxl")
my_data <- read_excel("coronavirus data.xlsx", sheet = 2)
my_data<- data.frame(my_data)
loci <- list()
replicates = 4 # 4 variant strains 
indices = seq(1, 825, by=5) # 825 rows in the excel file, 5 source populations
p_val = 0 # indexing variable 
loci_i = 0 # indexing variable 
loci_i_tot = 0 # indexing variable 
loci_tot = list()
p = list()
p_total = list()
counter =0
coeff = list()

# linear regressions
for (i in 1:length(indices)){
  counter =counter +1
    lr = lm(my_data$Phenotype[indices[i]:(indices[i]+replicates)]~my_data$Allele[indices[i]:(indices[i]+replicates)], data = my_data)
    p_val = p_val +1 
    p_total[p_val] = summary(lr)$coefficients[-1,4]
    loci_i_tot = loci_i_tot + 1
    loci_tot[loci_i_tot] = my_data$Locus[indices[i]]
}
# adjust p-values 
p_total <- p.adjust(p_total, method = "BH", n = length(p_total))

# determine the number of P-values < 0.1 
p_sig = list()
p_s = 0
min_p = min(p_total)
for (i in 1:length(p_total)) {
  if (p_total[i] <= min_p){
    p_s = p_s +1 
    p_sig[p_s] = p_total[i]
  }
}
 
# reformat data frame for manhattan plot 
library("qqman")
SNP <- list()
for (i in 1:length(loci_tot)){
  SNP[i] = paste("rs", as.character(loci_tot[i]), sep = "")
}
Chr <- list()
for (i in 1:length(loci_tot)){
  Chr[i] = 1
}
data_gwas = data.frame(unlist(SNP), unlist(Chr), unlist(loci_tot), unlist(p_total))
names(data_gwas) <- c("SNP", "CHR", "BP", "P")

# create manhattan plot 
manhattan(data_gwas)
