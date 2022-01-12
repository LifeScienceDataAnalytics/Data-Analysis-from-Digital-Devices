phewas_catalog <- read_csv("Downloads/phewas-catalog.csv")

ad_snp_catalog = phewas_catalog[phewas_catalog$`phewas phenotype` == "Alzheimer's disease",]
C0002395_disease_vda_summary <- read_excel("Downloads/C0002395_disease_vda_summary.xlsx")

RiskSnps_genotypes <- read.csv("RiskSnps_genotypes.csv",  row.names = 1, header= TRUE)
RiskSnps_genotypes_t <- t(RiskSnps_genotypes)
colnames(RiskSnps_genotypes_t) = RiskSnps_genotypes_t[1,]
write.csv(RiskSnps_genotypes_t, "RiskSnps_genotypes_t.csv")
