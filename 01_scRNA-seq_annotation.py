import muon as mu
import os
import pandas as pd
import scanpy as sc
import celltypist
from celltypist import models

def create_dictionary_from_file(filename):
        data_dict = {}
        with open(filename, 'r') as file:
                for line in file:
                        key, value = line.strip().split("\t")
                        data_dict[key] = value
        return data_dict

markers = {'T' : ['CD3D', 'CD3E', 'CD3G'],
       'CD4+T' : ['CD4', 'CD28'],
       'Treg' : 'FOXP3',
       'CD8+T' : ['CD8A', 'CD8B'],
       'MAIT' : ['KLRB1', 'SLC4A10'],
       'gdT' : ['TRGV9', 'TRDV2'],
       'CD56loNK' : ['GNLY', 'GZMB'],
       'CD56hiNK' : ['NCAM1', 'NCR1'],
       'Bcell' : ['MS4A1', 'CD79A', 'IGHD', 'IGHM'],
       'Plasma' : ['IGHA1', 'IGHG1'],
       'CD14+mono' : ['CD14', 'VCAN', 'S100A12', 'CSF3R'],
       'CD16+mono' : ['FCGR3A', 'MS4A7', 'LILRB1', 'CSF1R', 'CDKN1C'],
       'cDC' : ['HLA-DQA1', 'HLA-DQB1', 'CD1C'],
       'pDC' : ['LILRA4', 'IL3RA', 'IRF7'],
       'HSPC' : ['CD34', 'SPINK2'],
       'RBC' : ['HBA2', 'HBB'],
       'Platelet' : ['PF4', 'PPBP'],
       'Plymph' : ['MKI67', 'TOP2A']}


#annotate RNA-seq modality

#initial processing and filtering

fol = "./in/"
fol_list = os.listdir(fol)

adatas = []
samp_list = []

for f in fol_list:
        samp_list = samp_list + [f.replace("_raw_feature_bc_matrix.h5", "").replace("_", "-")]
        tmp_a = sc.read_10x_h5(fol + f, gex_only=True)
        tmp_a.var_names_make_unique()
        adatas = adatas + [tmp_a]

adata = adatas[0].concatenate(adatas[1:], batch_categories = samp_list)
adata.var_names_make_unique()

sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata = adata[adata.obs.pct_counts_mt < 10, :]
sc.pp.filter_genes(adata, min_cells=3)

adata.raw = adata

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, n_top_genes = 3000)
adata = adata[:, adata.var.highly_variable]

cell_cycle_genes = [x.strip() for x in open('./regev_lab_cell_cycle_genes.txt')]
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
sc.pp.regress_out(adata, ['S_score', 'G2M_score'])

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack', n_comps=50)

sc.external.pp.bbknn(adata)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
sc.tl.umap(adata)

#Clustering
sc.tl.leiden(adata, resolution = 50, key_added = 'leiden_50')  

#Annotating
models.download_models(force_update = True)
model = models.Model.load(model = 'Immune_All_Low.pkl')
predictions = celltypist.annotate(adata, model = 'Immune_All_Low.pkl', majority_voting = True)
adata1 = predictions.to_adata()

mv = adata1.obs['majority_voting'].to_dict()
adata.obs['cells'] = adata.obs.index
adata.obs['majority_voting'] = adata.obs['cells'].map(mv)
mv = adata1.obs['predicted_labels'].to_dict()
adata.obs['predicted_labels'] = adata.obs['cells'].map(mv)

sc.pl.umap(adata, color='majority_voting', title=' ', legend_fontsize = "xx-small", frameon=False)
sc.tl.dendrogram(adata, groupby = 'majority_voting')
sc.pl.dotplot(adata, markers, groupby = 'majority_voting', dendrogram=False)
sc.pl.umap(adata, color='predicted_labels', title=' ', legend_fontsize = "xx-small", frameon=False)
sc.tl.dendrogram(adata, groupby = 'predicted_labels')
sc.pl.dotplot(adata, markers, groupby = 'predicted_labels', dendrogram=False)
sc.pl.umap(adata, color='leiden_50', title=' ', legend_fontsize = "xx-small", frameon=False)
sc.tl.dendrogram(adata, groupby = 'leiden_50')
sc.pl.dotplot(adata, markers, groupby = 'leiden_50', dendrogram=False)

#Manualy curate cell cluster annotations based on marker genes and 'predicted_labels' and 'majority_voting' output from celltypist

ct_d = create_dictionary_from_file("tmp.txt")
adata.obs['cell_type'] = adata.obs['leiden_50'].map(ct_d)
adata = adata[~adata.obs['cell_type'].isin(['Ambient/Doublet'])]

adata.write("./rna.h5ad")


#process ATAC-seq modality and annotate according to the RNA-seq

fol = "/home/omics/DATA4/lisa/projects/AS_multiomics/GEX_ATAC_scanpy_asdas_2_1/fragments/"
fol_list = os.listdir(fol)

files = [(x.replace("_atac_fragments.tsv", "").replace("_", "-"), fol + x ) for x in fol_list]

adatas = snap.pp.import_data([fl for _, fl in files], 
                             file=[name + '.h5ad' for name, _ in files],tempdir="./tmp/", 
                             chrom_sizes=snap.genome.hg38, 
                             sorted_by_barcode=False, 
                             min_num_fragments=500)

snap.metrics.tsse(adatas, snap.genome.hg38)
snap.pp.filter_cells(adatas, min_tsse=7)
snap.pp.add_tile_matrix(adatas, bin_size=5000)
snap.pp.select_features(adatas, n_features=50000)

data = snap.AnnDataSet(adatas=[(name, adata) for (name, _), adata in zip(files, adatas)],filename="atac.h5ads")

data.obsm['fragment_paired'] = data.adatas.obsm['fragment_paired']
anndata_atac = data.to_adata()

atac.obs['cell_raw'] = atac.obs.index
atac.obs['cells'] = atac.obs.apply(lambda row: row['cell_raw'] + '-' + row['sample'], axis=1)
atac.obs_names = list(atac.obs['cells'])
atac = atac[atac.obs['cells'].isin(list(set(atac.obs['cells']).intersection(set(rna.obs['cells']))))]

ct = dict(zip(rna.obs['cells'], rna.obs['cell_type']))
atac.obs['ASDAS_response_CRP'] = atac.obs['cells'].map(ASDAS_response_CRP)

atac.raw = atac

snap.pp.select_features(atac, n_features=50000)
snap.tl.spectral(atac)
snap.pp.mnc_correct(atac, batch="sample")
snap.pp.harmony(atac, batch="sample", max_iter_harmony=20)
snap.tl.umap(atac, use_rep="X_spectral_mnn")
snap.pl.umap(atac, color="cell_type", interactive=False)

atac.write("./atac.h5ad")

#combine into mudata

atac = sc.read_h5ad("./atac.h5ad")
rna = sc.read_h5ad("./rna.h5ad")

atac = atac[atac.obs['cells'].isin(list(set(atac.obs['cells']).intersection(set(rna.obs['cells']))))]
rna = rna[rna.obs['cells'].isin(list(set(atac.obs['cells']).intersection(set(rna.obs['cells']))))]

mdata = mu.MuData({"rna": rna, "atac": atac})
mdata.write("./mudata.h5mu")

