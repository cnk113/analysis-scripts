g <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/GBM10/star/Solo.out/Velocyto/spliced/')
g2 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/GBM10/star/Solo.out/Velocyto/unspliced/')
g3 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/GBM10/star/Solo.out/Velocyto/ambiguous/')
g <- g + g3
gbm <- CreateSeuratObject(g+g2+g3)
gbm[['spliced']] <- CreateAssayObject(g)
gbm[['unspliced']] <- CreateAssayObject(g2)
gbm$orig.ident <- 'GBM10'

g <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/GBM9/star/Solo.out/Velocyto/spliced/')
g2 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/GBM9/star/Solo.out/Velocyto/unspliced/')
g3 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/GBM9/star/Solo.out/Velocyto/ambiguous/')
g <- g + g3
gbm2 <- CreateSeuratObject(g+g2+g3)
gbm2[['spliced']] <- CreateAssayObject(g)
gbm2[['unspliced']] <- CreateAssayObject(g2)
gbm2$orig.ident <- 'GBM9'

g <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11159_1/star/Solo.out/Velocyto/spliced/')
g2 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11159_1/star/Solo.out/Velocyto/unspliced/')
g3 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11159_1/star/Solo.out/Velocyto/ambiguous/')
g <- g + g3
gbm3 <- CreateSeuratObject(g+g2+g3)
gbm3[['spliced']] <- CreateAssayObject(g)
gbm3[['unspliced']] <- CreateAssayObject(g2)
gbm3$orig.ident <- 'SF11159_1'

g <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11159_2/star/Solo.out/Velocyto/spliced/')
g2 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11159_2/star/Solo.out/Velocyto/unspliced/')
g3 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11159_2/star/Solo.out/Velocyto/ambiguous/')
g <- g + g3
gbm4 <- CreateSeuratObject(g+g2+g3)
gbm4[['spliced']] <- CreateAssayObject(g)
gbm4[['unspliced']] <- CreateAssayObject(g2)
gbm4$orig.ident <- 'SF11159_2'

g <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11209_1/star/Solo.out/Velocyto/spliced/')
g2 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11209_1/star/Solo.out/Velocyto/unspliced/')
g3 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11209_1/star/Solo.out/Velocyto/ambiguous/')
g <- g + g3
gbm5 <- CreateSeuratObject(g+g2+g3)
gbm5[['spliced']] <- CreateAssayObject(g)
gbm5[['unspliced']] <- CreateAssayObject(g2)
gbm5$orig.ident <- 'SF11209_1'

g <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11209_2/star/Solo.out/Velocyto/spliced/')
g2 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11209_2/star/Solo.out/Velocyto/unspliced/')
g3 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11209_2/star/Solo.out/Velocyto/ambiguous/')
g <- g + g3
gbm6 <- CreateSeuratObject(g+g2+g3)
gbm6[['spliced']] <- CreateAssayObject(g)
gbm6[['unspliced']] <- CreateAssayObject(g2)
gbm6$orig.ident <- 'SF11209_2'

g <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11215_1/star/Solo.out/Velocyto/spliced/')
g2 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11215_1/star/Solo.out/Velocyto/unspliced/')
g3 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11215_1/star/Solo.out/Velocyto/ambiguous/')
g <- g + g3
gbm7 <- CreateSeuratObject(g+g2+g3)
gbm7[['spliced']] <- CreateAssayObject(g)
gbm7[['unspliced']] <- CreateAssayObject(g2)
gbm7$orig.ident <- 'SF11215_1'

g <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11215_2/star/Solo.out/Velocyto/spliced/')
g2 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11215_2/star/Solo.out/Velocyto/unspliced/')
g3 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11215_2/star/Solo.out/Velocyto/ambiguous/')
g <- g + g3
gbm8 <- CreateSeuratObject(g+g2+g3)
gbm8[['spliced']] <- CreateAssayObject(g)
gbm8[['unspliced']] <- CreateAssayObject(g2)
gbm8$orig.ident <- 'SF11215_2'

g <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11232_1/star/Solo.out/Velocyto/spliced/')
g2 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11232_1/star/Solo.out/Velocyto/unspliced/')
g3 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11232_1/star/Solo.out/Velocyto/ambiguous/')
g <- g + g3
gbm9 <- CreateSeuratObject(g+g2+g3)
gbm9[['spliced']] <- CreateAssayObject(g)
gbm9[['unspliced']] <- CreateAssayObject(g2)
gbm9$orig.ident <- 'SF11232_1'

g <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11232_2/star/Solo.out/Velocyto/spliced/')
g2 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11232_2/star/Solo.out/Velocyto/unspliced/')
g3 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11232_2/star/Solo.out/Velocyto/ambiguous/')
g <- g + g3
gbm10 <- CreateSeuratObject(g+g2+g3)
gbm10[['spliced']] <- CreateAssayObject(g)
gbm10[['unspliced']] <- CreateAssayObject(g2)
gbm10$orig.ident <- 'SF11232_2'

g <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11247_1/star/Solo.out/Velocyto/spliced/')
g2 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11247_1/star/Solo.out/Velocyto/unspliced/')
g3 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11247_1/star/Solo.out/Velocyto/ambiguous/')
g <- g + g3
gbm11 <- CreateSeuratObject(g+g2+g3)
gbm11[['spliced']] <- CreateAssayObject(g)
gbm11[['unspliced']] <- CreateAssayObject(g2)
gbm11$orig.ident <- 'SF11247_1'

g <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11247_2/star/Solo.out/Velocyto/spliced/')
g2 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11247_2/star/Solo.out/Velocyto/unspliced/')
g3 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11247_2/star/Solo.out/Velocyto/ambiguous/')
g <- g + g3
gbm12 <- CreateSeuratObject(g+g2+g3)
gbm12[['spliced']] <- CreateAssayObject(g)
gbm12[['unspliced']] <- CreateAssayObject(g2)
gbm12$orig.ident <- 'SF11247_2'

g <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11285_1/star/Solo.out/Velocyto/spliced/')
g2 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11285_1/star/Solo.out/Velocyto/unspliced/')
g3 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11285_1/star/Solo.out/Velocyto/ambiguous/')
g <- g + g3
gbm13 <- CreateSeuratObject(g+g2+g3)
gbm13[['spliced']] <- CreateAssayObject(g)
gbm13[['unspliced']] <- CreateAssayObject(g2)
gbm13$orig.ident <- 'SF11285_1'

g <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11285_2/star/Solo.out/Velocyto/spliced/')
g2 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11285_2/star/Solo.out/Velocyto/unspliced/')
g3 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/SF11285_2/star/Solo.out/Velocyto/ambiguous/')
g <- g + g3
gbm14 <- CreateSeuratObject(g+g2+g3)
gbm14[['spliced']] <- CreateAssayObject(g)
gbm14[['unspliced']] <- CreateAssayObject(g2)
gbm14$orig.ident <- 'SF11285_2'

g <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/Tumor2/star/Solo.out/Velocyto/spliced/')
g2 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/Tumor2/star/Solo.out/Velocyto/unspliced/')
g3 <- Read10X('/media/chang/HDD-6/chang/mg/gbm_aparna/Tumor2/star/Solo.out/Velocyto/ambiguous/')
g <- g + g3
gbm15 <- CreateSeuratObject(g+g2+g3)
gbm15[['spliced']] <- CreateAssayObject(g)
gbm15[['unspliced']] <- CreateAssayObject(g2)
gbm15$orig.ident <- 'Tumor2'

gbm1 <- merge(gbm, c(gbm2,gbm3,gbm4,gbm5,gbm6,gbm7,gbm8,gbm9,gbm10,gbm11,gbm12,gbm13,gbm14,gbm15))
gbm <- gbm1
rm(gbm1,gbm2,gbm3,gbm4,gbm5,gbm6,gbm7,gbm8,gbm9,gbm10,gbm11,gbm12,gbm13,gbm14,gbm15)
gbm1[["percent.mt"]] <- PercentageFeatureSet(gbm1, pattern = "^MT-", assay='RNA')