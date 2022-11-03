import scanpy as sc
import anndata
import sys
import numpy as np
import dropkick as dk
import pandas as pd
import scanpy as sc
import doubletdetection
from scar import model
from scipy import sparse

def buildAnndataFromStar(path, barcode_path=None, batch=None):
    """Generate an anndata object from the STAR aligner output folder"""
    # Load Read Counts
    X = sc.read_mtx(path+'GeneFull_Ex50pAS/raw/matrix.mtx')

    # Transpose counts matrix to have Cells as rows and Genes as cols as expected by AnnData objects
    X = X.X.transpose()

    # Load the 3 matrices containing Spliced, Unspliced and Ambigous reads
    mtxU = np.loadtxt(path+'Velocyto/raw/unspliced.mtx', skiprows=3, delimiter=' ')
    mtxS = np.loadtxt(path+'Velocyto/raw/spliced.mtx', skiprows=3, delimiter=' ')
    mtxA = np.loadtxt(path+'Velocyto/raw/ambiguous.mtx', skiprows=3, delimiter=' ')

    # Extract sparse matrix shape informations from the third row
    shapeU = np.loadtxt(path+'Velocyto/raw/unspliced.mtx', skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)
    shapeS = np.loadtxt(path+'Velocyto/raw/spliced.mtx', skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)
    shapeA = np.loadtxt(path+'Velocyto/raw/ambiguous.mtx', skiprows=2, max_rows = 1 ,delimiter=' ')[0:2].astype(int)

    # Read the sparse matrix with csr_matrix((data, (row_ind, col_ind)), shape=(M, N))
    # Subract -1 to rows and cols index because csr_matrix expects a 0 based index
    # Traspose counts matrix to have Cells as rows and Genes as cols as expected by AnnData objects

    spliced = sparse.csr_matrix((mtxS[:,2], (mtxS[:,0]-1, mtxS[:,1]-1)), shape = shapeS).transpose()
    unspliced = sparse.csr_matrix((mtxU[:,2], (mtxU[:,0]-1, mtxU[:,1]-1)), shape = shapeU).transpose()
    ambiguous = sparse.csr_matrix((mtxA[:,2], (mtxA[:,0]-1, mtxA[:,1]-1)), shape = shapeA).transpose()

    # Load Genes and Cells identifiers
    obs = pd.read_csv(path+'GeneFull_Ex50pAS/raw/barcodes.tsv',
                  header = None, index_col = 0)

    # Remove index column name to make it compliant with the anndata format
    obs.index.name = None

    var = pd.read_csv(path+'GeneFull_Ex50pAS/raw/features.tsv', sep='\t',
                                    names = ('gene_ids', 'feature_types'), index_col = 1)

    # Build AnnData object to be used with ScanPy and ScVelo
    adata = anndata.AnnData(X = X, obs = obs, var = var,
                            layers = {'spliced': spliced, 'unspliced': unspliced, 'ambiguous': ambiguous})

    adata.var_names_make_unique()

    sc.pp.filter_genes(adata, min_cells=1)

    # Run dropkick
    adata_model = dk.dropkick(adata, n_jobs=5)

    # Run scAR
    cell_free = adata[adata.obs["dropkick_label"] == "False"]
    # average and normalize the transcript in cell-free droplets.
    ambient_profile = pd.DataFrame((cell_free.X.sum(axis=0)/cell_free.X.sum()).A1, index=adata.var.index)

    adata = adata[adata.obs["dropkick_label"] == "True"]
    raw_counts = adata.to_df()

    scar = model(raw_count = raw_counts,
                 ambient_profile = ambient_profile,
                 feature_type='mRNA')
    scar.train(epochs=400,
               batch_size=64,
               verbose=True)
    scar.inference()

    adata.layers['scar'] = sparse.csr_matrix(pd.DataFrame(scar.native_counts,
                                                          index=raw_counts.index,
                                                          columns=raw_counts.columns))
    
    # Doublet detection
    clf = doubletdetection.BoostClassifier(
        n_iters=25,
        clustering_algorithm="leiden",
        standard_scaling=True,
        pseudocount=0.1,
        n_jobs=-1)
    adata.obs["doublet"] = clf.fit(adata.layers['scar']).predict()
    adata.obs["doublet_score"] = clf.doublet_score()

    if barcode_path:
        bc_tsv = pd.read_csv(barcode_path, sep='\t', header=0, index_col=1)
        df = pd.DataFrame(0, index=bc_tsv.index.unique(),columns=bc_tsv.CBC.unique())
        for index, row in bc_tsv.iterrows():
            df.loc[index,row.CBC] = row.UMI_Count

        raw_counts = df[df.columns.intersection(adata.obs.index)]
        cell_free = df.drop(df.columns.intersection(raw_counts.columns), axis=1)
        ambient_profile = cell_free.T.sum()/cell_free.T.sum().sum()
        ambient_profile = ambient_profile.to_frame("ambient profile")

        sgRNAs = model(raw_count = raw_counts.T,
               ambient_profile = ambient_profile,
               feature_type = 'sgRNAs')

        sgRNAs.train(epochs=200,
             batch_size=64,
             verbose=True)

        sgRNAs.inference(cutoff=3)
        adata.obs["X_scar_assignment_3"] = "NA"
        adata.obs["X_scar_assignment_3_n_BC"] = "NA"
        adata.obs.loc[raw_counts.columns,"X_scar_assignment_3"] = sgRNAs.feature_assignment.sgRNAs.astype(str)
        adata.obs.loc[raw_counts.columns,"X_scar_assignment_3_n_BC"] = sgRNAs.feature_assignment.n_sgRNAs.astype(str)

        sgRNAs.inference(cutoff=10)
        adata.obs["X_scar_assignment_10"] = "NA"
        adata.obs["X_scar_assignment_10_n_BC"] = "NA"
        adata.obs.loc[raw_counts.columns,"X_scar_assignment_10"] = sgRNAs.feature_assignment.sgRNAs.astype(str)
        adata.obs.loc[raw_counts.columns,"X_scar_assignment_10_n_BC"] = sgRNAs.feature_assignment.n_sgRNAs.astype(str)
    if batch:
        adata.obs["Batch"] = batch
    return adata.copy()

sample_paths = ["/media/chang/HDD-8/chang/cloneseq/fastqs/GW15_rep1/GW15_rep1_tube1/GW15_rep1_tube1_Solo.out/",
                "/media/chang/HDD-8/chang/cloneseq/fastqs/GW15_rep1/GW15_rep1_tube2/GW15_rep1_tube2_Solo.out/",
                "/media/chang/HDD-8/chang/cloneseq/fastqs/GW15_rep1/GW15_rep1_tube3/GW15_rep1_tube3_Solo.out/",
                "/media/chang/HDD-8/chang/cloneseq/fastqs/GW15_rep2/GW15_rep2_tube1/GW15_rep2_tube1_Solo.out/",
                "/media/chang/HDD-8/chang/cloneseq/fastqs/GW15_rep2/GW15_rep2_tube2/GW15_rep2_tube2_Solo.out/",
                "/media/chang/HDD-8/chang/cloneseq/fastqs/GW15_rep2/GW15_rep2_tube3/GW15_rep2_tube3_Solo.out/",
                "/media/chang/HDD-8/chang/cloneseq/fastqs/GW15_rep3/GW15_rep3_Solo.out/",
                "/media/chang/HDD-8/chang/cloneseq/fastqs/GW15_rep4/GW15_rep4_Solo.out/",
                "/media/chang/HDD-8/chang/cloneseq/fastqs/GW15_rep5/GW15_rep5_Solo.out/",
                "/media/chang/HDD-8/chang/cloneseq/fastqs/GW18_MGE/GW18_MGE_Solo.out/",
                "/media/chang/HDD-8/chang/cloneseq/fastqs/GW18_PFC/GW18_PFC_tube1/GW18_PFC_tube1_Solo.out/",
                "/media/chang/HDD-8/chang/cloneseq/fastqs/GW18_PFC/GW18_PFC_tube2/GW18_PFC_tube2_Solo.out/",
                "/media/chang/HDD-8/chang/cloneseq/fastqs/GW18_V1/GW18_V1_tube1/GW18_V1_tube1_Solo.out/",
                "/media/chang/HDD-8/chang/cloneseq/fastqs/GW18_V1/GW18_V1_tube2/GW18_V1_tube2_Solo.out/",
                "/media/chang/HDD-8/chang/cloneseq/fastqs/GW18_V1/GW18_V1_tube3/GW18_V1_tube3_Solo.out/"]
barcode_paths = ["/media/chang/HDD-8/chang/cloneseq/final/GW15_Rep1_Tube1_Final_Barcodes.tsv",
                 "/media/chang/HDD-8/chang/cloneseq/final/GW15_Rep1_Tube2_Final_Barcodes.tsv",
                 "/media/chang/HDD-8/chang/cloneseq/final/GW15_Rep1_Tube3_Final_Barcodes.tsv",
                 "/media/chang/HDD-8/chang/cloneseq/final/GW15_Rep2_Tube1_Final_Barcodes.tsv",
                 "/media/chang/HDD-8/chang/cloneseq/RND_STICR_GEO/tsv/GW15_6weeks_Rep2_Tube_2.tsv",
                 "/media/chang/HDD-8/chang/cloneseq/RND_STICR_GEO/tsv/GW15_6weeks_Rep2_Tube_3.tsv",
                 "/media/chang/HDD-8/chang/cloneseq/final/GW15_Rep3_Final_Barcodes.tsv",
                 "/media/chang/HDD-8/chang/cloneseq/final/GW15_Rep4_Final_Barcodes.tsv",
                 "/media/chang/HDD-8/chang/cloneseq/final/GW15_Rep5_Final_Barcodes.tsv",
                 "/media/chang/HDD-8/chang/cloneseq/final/GW18_MGE_Final_Barcodes.tsv",
                 "/media/chang/HDD-8/chang/cloneseq/final/GW18_PFC_Tube1_Final_Barcodes.tsv",
                 "/media/chang/HDD-8/chang/cloneseq/final/GW18_PFC_Tube2_Final_Barcodes.tsv",
                 "/media/chang/HDD-8/chang/cloneseq/final/GW18_V1_Tube1_Final_Barcodes.tsv",
                 "/media/chang/HDD-8/chang/cloneseq/final/GW18_V1_Tube2_Final_Barcodes.tsv",
                 "/media/chang/HDD-8/chang/cloneseq/final/GW18_V1_Tube3_Final_Barcodes.tsv"]


ad = []
for sample, bc in zip(sample_paths,barcode_paths):
    batch = sample.split("/")
    ad.append(buildAnndataFromStar(sample,bc,batch[8]))
merged = sc.concat(ad)
merged.write("/media/chang/HDD-8/chang/cloneseq/new_merged.h5ad")

