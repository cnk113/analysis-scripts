def build(filename, pca = False, umap = False):
    import anndata
    import pandas as pd
    from itertools import chain
    ad = anndata.read_mtx('combined.mtx')
    #ad.layers['sct'] = anndata.read_mtx('sct.mtx').X
    ad.obs = pd.read_csv('meta.csv',header=0)
    if umap:
        ad.obsm['X_umap'] = pd.read_csv('umap.csv',header=0).values
    if pca:
        ad.obsm['X_pca'] = pd.read_csv('pca.csv',header=0).values
    ad.var_names = list(chain.from_iterable(pd.read_csv('genes.csv', header=0).values.tolist()))
    ad.obs_names = list(chain.from_iterable(pd.read_csv('cells.csv', header=0).values.tolist()))
    ad.write(filename)
