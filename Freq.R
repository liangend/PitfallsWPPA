library(rrBLUP)

args = commandArgs(trailingOnly=TRUE)
nqtl_index = args[1]  # number of QTL index = 1, 2, 3
replicate = args[2] # replicate = 1 to 30

nqtl_index = as.numeric(nqtl_index)
replicate = as.numeric(replicate)

nqtl_value = c(30, 90, 300)
nqtl = nqtl_value[nqtl_index]

datapath = paste0("nqtl_",nqtl,"_replicate_", replicate)

## Import genotype, phenotype and geno map
Genotype = read.csv('controlled_genotype.csv')
Map = read.csv('controlled_map.csv')
Phenotype = read.csv(paste0(datapath, '/Pheno_data.csv'))

## Window-based GWAS using EMMAX
y = Phenotype[,1]
Z = as.matrix(Genotype[,-1])
result = mixed.solve(y = y, Z = Z)
g = result$u
sigmag2 = result$Vu
sigmae2 = result$Ve
X = rep(1, nrow(Phenotype))
A = t(X) %*% X / sigmae2
B = t(X) %*% Z / sigmae2
C = t(Z) %*% X / sigmae2
D = t(Z) %*% Z / sigmae2 + diag(ncol(Z)) / sigmag2
Cgg = solve(D - C %*% solve(A) %*% B) 
# Note that Cgg matrix is very large (> 10 GB), takes very long time to calculate.

savepath = paste0("/Freq/nqtl_",nqtl,"_replicate_", replicate)
saveRDS(Cgg, file = paste0(savepath, "/Cgg_matrix"))


## Organize the genomic data based on window size, 
## works for window size of 1 Mb (window_size = 1000000) and 100 (window_size = 100))
window_dat_org = function (window_size){
  win = c()
  winrow = c()
  if (window_size == 100 | window_size == 1) {
    Map$num = 0:(nrow(Map) - 1)
    class = table(floor(Map$num / window_size))
    clapos = 1
    for (n in 1:length(class)) {
      win = append(win, c(Map$pos[clapos], Map$pos[clapos + class[n] - 1]))
      clapos = clapos + class[n]
    }
    winrow = append(winrow, class)
  } else {
    for (i in 1:12) {
      Map_sub = Map[Map$chrom == i,]
      class = table(floor(Map_sub$pos / window_size))
      clapos = 1
      for (n in 1:length(class)) {
        win = append(win, c(Map_sub$pos[clapos], Map_sub$pos[clapos + class[n] - 1]))
        clapos = clapos + class[n]
      }
      winrow = append(winrow, class)
    }
  }
  
  chrom = c()
  ch = 1
  if (window_size == 100) {
    for (i in 1:(length(winrow) - 1)) {
      if (win[2*i - 1] < win[2*i]){
        chrom[i] = ch
      } else {
        chrom[i] = ch
        ch = ch + 1
      }
    }
    chrom[length(chrom) + 1] = 12
  } else if (window_size == 1) {
    chrom = Map$chrom
  } else {
    for (i in 1:(length(winrow) - 1)) {
      if (floor(abs(win[2*i + 1] - win[2*i])/window_size) < 10){
        chrom[i] = ch
      } else {
        chrom[i] = ch
        ch = ch + 1
      }
    }
    chrom[length(chrom) + 1] = 12
  }
  
  window_dat = data.frame(chr = chrom, start_SNP = win[2 * (1:length(winrow)) - 1],
                          end_SNP = win[2 * (1:length(winrow))], window = 1:length(winrow), 
                          numSNP = winrow)
  
  return(window_dat)
}

## Calculate the chi square test statistic for each window
chi2_cal = function (win_dat, Cgg, g, sigmag2) {
  chi2k = c()
  i = 1
  winpos = 1
  winrow = win_dat$numSNP
  for (k in 1:length(winrow)) {
    gk = g[winpos : (winpos + winrow[k] - 1)]
    Cggk = Cgg[winpos : (winpos + winrow[k] - 1), winpos : (winpos + winrow[k] - 1)]
    winpos = winpos + winrow[k]
    if (length(Cggk) == 1) {
      a = 1
    } else {a = ncol(Cggk)}
    chi2k[i] = t(gk) %*% solve(diag(a) * sigmag2 - Cggk) %*% gk
    i = i + 1
  }
  return(list(chi2k = chi2k, df = winrow))
}

win_1mb = window_dat_org(1000000)
chi2k_1mb = chi2_cal(win_1mb, Cgg, g, sigmag2)
P_value_1mb = -pchisq(chi2k_1mb[[1]], chi2k_1mb[[2]], lower.tail = F, log.p = T)
win_1mb$neg_logP = P_value_1mb

write.csv(win_1mb, paste0(savepath, '/win_1mb.csv'))








