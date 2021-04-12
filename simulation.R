setwd("~/Desktop/MiCM/")

library(snpStats)
library(rms)
genotype <- read.plink("EUR_10kSNP_chr22.bed",
                       "EUR_10kSNP_chr22.bim",
                       "EUR_10kSNP_chr22.fam")

set.seed(1)

### Randomly sample 100 / 9973 SNPs to have non-zero effects
causalSNP <- sort(sample(1:ncol(genotype$genotypes),100,replace = F)) 
trueBeta <- rnorm(100)

### Generate response Y = Xb + E + e
### SNP heritability (Var(Xb) / Var(Y)) = 0.3; Environmental risk factor PVE (Var(E) / Var(Y)) = 0.1
X <- read.table("EUR_10kSNP_chr22.raw",header = T)
rownames(X) <- X$FID
X <- X[,-c(1:6)]
SNPeffect <- as.matrix(X[,causalSNP]) %*% trueBeta
SNPeffect <- scale(SNPeffect)
E <- rnorm(15000)
E <- scale(E)
noise <- rnorm(15000)
noise <- scale(noise)
Y <- sqrt(0.3) * SNPeffect + sqrt(0.1) * E + sqrt(0.6) * noise
var(Y)
summary(lm(Y ~ SNPeffect + E))

training_genotype <- X[1:5000,]
modelselection_genotype <- X[5001:10000,]
test_genotype <- X[10001:15000,]
save(training_genotype, modelselection_genotype, test_genotype, file = "EUR_10kSNP_chr22.raw.RData")

training_Y <- Y[1:5000]
modelselection_Y <- Y[5001:10000]
test_Y <- Y[10001:15000]
write.table(training_Y, "training_Y.txt", quote = F,row.names = F,col.names = F)
write.table(modelselection_Y, "modelselection_Y.txt", quote = F,row.names = F,col.names = F)
write.table(test_Y, "test_Y.txt", quote = F,row.names = F,col.names = F)

training_E <- E[1:5000]
modelselection_E <- E[5001:10000]
test_E <- E[10001:15000]
write.table(training_E, "training_E.txt", quote = F,row.names = F,col.names = F)
write.table(modelselection_E, "modelselection_E.txt", quote = F,row.names = F,col.names = F)
write.table(test_E, "test_E.txt", quote = F,row.names = F,col.names = F)

### Binary outcome with baseline prevalence of 0.2
### logit(P) = Y 
ilogit <- function(x) {
  exp(x) / (1 + exp(x))
}

set.seed(2021)
Z <- sapply(1:15000, function(x) rbinom(1,1,ilogit(log(0.2 / (1-0.2)) + Y[x])))
summary(glm(Z ~ SNPeffect + E, family = binomial))

training_Z <- Z[1:5000]
modelselection_Z <- Z[5001:10000]
test_Z <- Z[10001:15000]
write.table(training_Z, "training_Z.txt", quote = F,row.names = F,col.names = F)
write.table(modelselection_Z, "modelselection_Z.txt", quote = F,row.names = F,col.names = F)
write.table(test_Z, "test_Z.txt", quote = F,row.names = F,col.names = F)

### Time-to-event
set.seed(2020)
time_to_event <- data.frame(Event = Z, Age = rnorm(15000, mean = 60, sd = 5))
time_to_event$Age[time_to_event$Event==1] <- time_to_event$Age[time_to_event$Event==1] - Y[time_to_event$Event==1] * runif(sum(time_to_event$Event==1),min = 1,max = 5)

cphmod <- coxph(Surv(time_to_event$Age, time_to_event$Event) ~ as.numeric(SNPeffect) + as.numeric(E), x = T, y = T)
summary(cphmod)

write.table(time_to_event[1:5000,], "training_timetoevent.txt", quote = F,row.names = F,col.names = T,sep = "\t")
write.table(time_to_event[5001:10000,], "modelselection_timetoevent.txt", quote = F,row.names = F,col.names = T,sep = "\t")
write.table(time_to_event[10001:15000,], "test_timetoevent.txt", quote = F,row.names = F,col.names = T,sep = "\t")

### Association test based on training dataset
### Adjusted for the effect of environmental risk factor
sumstat <- data.frame(SNP = rownames(genotype$map),
                      chromosome = genotype$map$chromosome,
                      position = genotype$map$position,
                      effect_allele = genotype$map$allele.1,
                      other_allele = genotype$map$allele.2,
                      beta = NA,
                      StdErr = NA,
                      P = NA)
for (i in 1:9973) {
  cat("processing SNP",i,"/ 9973",collapse="\n")
  sumstat[i,6:8] <- summary(lm(training_Y ~ training_genotype[,i] + training_E))$coefficients[2,c(1,2,4)]
}

### Visualize summary statistics
plotdat <- data.frame(position = sumstat$position,
                      p = sumstat$P,
                      causal = 0,
                      beta = 0)
plotdat$causal[causalSNP] <- 1
plotdat$beta[causalSNP] <- trueBeta

ggplot(plotdat,aes(x = position, y = -log10(p), 
                   shape = as.factor(causal),
                   size = as.factor(causal),
                   color = as.factor(causal))) + 
  geom_point() + theme_classic() + xlab("Position") +
  ylab(expression(-log[10](P-value))) +
  scale_color_manual(values = c("darkblue","orange")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13)) +
  geom_hline(yintercept = -log10(5e-8), col = "darkred", lty = 2)

ggplot(plotdat,aes(x = position, y = beta, 
                   shape = as.factor(causal),
                   size = as.factor(causal),
                   color = as.factor(causal))) + 
  geom_point() + theme_classic() + xlab("Position") +
  ylab("True effect") +
  scale_color_manual(values = c("darkblue","orange")) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 13))

sumstat$true_beta <- plotdat$beta
write.table(sumstat, "training_GWAS.sumstat.txt", quote = F,row.names = F,col.names = T,sep = "\t")

