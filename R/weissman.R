
filtered_allele_table = filtered_allele_table[filtered_allele_table$cellBC %in% 
                                                intersect(filtered_allele_table$cellBC, colnames(mice)),]


r1 = table(filtered_allele_table$r1)
r1 = r1[order(r1, decreasing = T)][1]
filtered_allele_table2$r1[filtered_allele_table2$r1 == 'CCGAA[None]AAATG']  <- 0
r2 = table(filtered_allele_table$r2)
r2 = r2[order(r2, decreasing = T)][1]
filtered_allele_table2$r2[filtered_allele_table2$r2 == 'GATAT[None]CTCTG']  <- 0
r3 = table(filtered_allele_table$r3)
r3 = r3[order(r3, decreasing = T)][1]
filtered_allele_table2$r3[filtered_allele_table2$r3 == 'ATTCG[None]CGGAG']  <- 0

cells = table(filtered_allele_table2$cellBC)
cells = cells[cells==1]

filtered_allele_table3 = filtered_allele_table2[filtered_allele_table2$cellBC %in% names(cells),]
filtered_allele_table3 = data.frame(filtered_allele_table3)
rownames(filtered_allele_table3) = filtered_allele_table3$cellBC
filtered_allele_table3$type = mice3$SCT_snn_res.0.4[filtered_allele_table3$cellBC]
matrix = as.matrix(filtered_allele_table3[,c('r1','r2','r3')])

filtered_allele_table4 = filtered_allele_table3[filtered_allele_table3$type == 2,]
matrix = as.matrix(filtered_allele_table4[,c('r1','r2','r3')])
filtered_allele_table4$type = mice4$SCT_snn_res.1[filtered_allele_table4$cellBC]



library(furrr)
plan(multisession, workers = 8)

tr = phylotime(matrix, t_total = 26)

plot_barcodes(matrix, tr, show_column_names = F)

plot_tr(tr,
        end_alpha_terminal = 0.05,
        edge_width_terminal = 0.02)

sc_celltypes = filtered_allele_table3$type
names(sc_celltypes) = filtered_allele_table3$cell
print(sc_celltypes[1:10])


plot_barcodes(matrix, tr, tip_celltype = sc_celltypes, show_column_names = F)