import pycisTopic
import os, sys
import pandas as pd
import scanpy as sc
import matplotlib
import scvi
import muon as mu
from muon import atac as ac
import snapatac2 as snap
import anndata as ad
import numpy as np
from scipy.sparse import issparse

cts = ['B', 'NK', 'CD14+_Monocytes', 'CD4+_T', 'CD8+_T']
ct = cts[0] #repeat for all cell types

#Initial data load

global_dir = "./Scenic/"
out_dir = global_dir + ct + "_asdas_0vs12_tp1/"
contrast = "atac:group_012"

if not os.path.exists(os.path.join(out_dir)):
    os.makedirs(os.path.join(out_dir))

mdata = mu.read("./mudata.h5mu")
mdata.obs['cell_type'] = mdata.obs.apply(lambda row: row['rna:cell_type_global'], axis=1)
mdata = mdata[mdata.obs['cell_type'].isin([ct])]
mdata = mdata[mdata.obs['atac:time_point'].isin(['1'])]
mdata['atac'].obs = mdata['atac'].obs.rename(columns = {'cell_raw' : 'barcode'})
mdata['atac'].obs['group_012'] = mdata['atac'].obs.apply(lambda row: row['group'].replace('asdas_1', 'asdas_12').replace('asdas_2', 'asdas_12'), axis=1)

group_012 = mdata['atac'].obs.set_index('cells')['group_012'].to_dict()
mdata['rna'].obs['group_012'] = mdata['rna'].obs['cells'].map(group_012)


chromsizes = pd.read_table("http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes", header = None, names = ["Chromosome", "End"])
chromsizes.insert(1, "Start", 0)

fr_d = "./fragments/"

fragments_dict = {'S077-1' : fr_d + "S077-1_atac_fragments.tsv.gz",
 'm-SpA-105-1' : fr_d + "m-SpA-105-1_atac_fragments.tsv.gz",
 'S094-2' : fr_d + "S094-2_atac_fragments.tsv.gz",
 'm-SpA-112-1' : fr_d + "m-SpA-112-1_atac_fragments.tsv.gz",
 'SpA-113-2' : fr_d + "SpA_113_2_atac_fragments.tsv.gz",
 'm-SpA-010-2' : fr_d + "m-SpA-010-2_0124_atac_fragments.tsv.gz",
 'SpA-115-2' : fr_d + "SpA_115_2_atac_fragments.tsv.gz",
 'SpA-028-2' : fr_d + "SpA_028_2_atac_fragments.tsv.gz",
 'm-SpA-010-1' : fr_d + "m-SpA-010-1_atac_fragments.tsv.gz",
 'm-SpA-098-2' : fr_d + "m-SpA-098-2_atac_fragments.tsv.gz",
 'S087-2' : fr_d + "S087-2_atac_fragments.tsv.gz",
 'SpA-113-1' : fr_d + "SpA_113_1_atac_fragments.tsv.gz",
 'm-SpA-098-1' : fr_d + "m-SpA-098-1_atac_fragments.tsv.gz",
 'o-SpA-060-2' : fr_d + "o-SpA-060-2_atac_fragments.tsv.gz",
 'm-SpA-099-1' : fr_d + "m-SpA-099-1_atac_fragments.tsv.gz",
 'm-SpA-105-2' : fr_d + "m-SpA-105-2_atac_fragments.tsv.gz",
 'm-SpA-109-1' : fr_d + "m-SpA-109-1_atac_fragments.tsv.gz",
 'o-SpA-060-1' : fr_d + "o-SpA-060-1_atac_fragments.tsv.gz",
 'S094-1' : fr_d + "S094-1_atac_fragments.tsv.gz",
 'S077-2' : fr_d + "S077-2_atac_fragments.tsv.gz",
 'SpA-028-1' : fr_d + "SpA_028_1_atac_fragments.tsv.gz",
 'm-SpA-099-2' : fr_d + "m-SpA-099-2_atac_fragments.tsv.gz",
 'm-SpA-112-2' : fr_d + "m-SpA-112-2_atac_fragments.tsv.gz",
 'm-SpA-013-1' : fr_d + "m-SpA-013-1_atac_fragments.tsv.gz",
 'm-SpA-109-2' : fr_d + "m-SpA-109-2_atac_fragments.tsv.gz",
 'SpA-115-1' : fr_d + "SpA_115_1_atac_fragments.tsv.gz",
 'S087-1' : fr_d + "S087-1_atac_fragments.tsv.gz",
 'm-SpA-013-2' : fr_d + "m-SpA-013-2_0124_atac_fragments.tsv.gz"
}

fragments_dict = {key: fragments_dict[key] for key in list(mdata['atac'].obs['sample'].drop_duplicates()) if key in fragments_dict}

#run pycisTopic

from pycisTopic.pseudobulk_peak_calling import export_pseudobulk

os.makedirs(os.path.join(out_dir, "./consensus_peak_calling"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "./consensus_peak_calling/pseudobulk_bed_files"), exist_ok = True)
os.makedirs(os.path.join(out_dir, "./consensus_peak_calling/pseudobulk_bw_files"), exist_ok = True)

bw_paths, bed_paths = export_pseudobulk(
    input_data = mdata['atac'].obs,
    variable = "group_012",
    sample_id_col = "sample",
    chromsizes = chromsizes,
    bed_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"),
    bigwig_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"),
    path_to_fragments = fragments_dict,
    n_cpu = 1,
    normalize_bigwig = True,
    temp_dir = "/home/tmp" #,
#    split_pattern = "___"
)

with open(os.path.join(out_dir, "consensus_peak_calling/bw_paths.tsv"), "wt") as f:
    for v in bw_paths:
        _ = f.write(f"{v}\t{bw_paths[v]}\n")

with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv"), "wt") as f:
    for v in bed_paths:
        _ = f.write(f"{v}\t{bed_paths[v]}\n")


bw_paths = {}
with open(os.path.join(out_dir, "consensus_peak_calling/bw_paths.tsv")) as f:
    for line in f:
        v, p = line.strip().split("\t")
        bw_paths.update({v: p})
bed_paths = {}
with open(os.path.join(out_dir, "consensus_peak_calling/bed_paths.tsv")) as f:
    for line in f:
        v, p = line.strip().split("\t")
        bed_paths.update({v: p})


from pycisTopic.pseudobulk_peak_calling import peak_calling
macs_path = "macs2"
import ray
ray.shutdown()

os.makedirs(os.path.join(out_dir, "consensus_peak_calling/MACS"), exist_ok = True)

narrow_peak_dict = peak_calling(
    macs_path = macs_path,
    bed_paths = bed_paths,
    outdir = os.path.join(os.path.join(out_dir, "consensus_peak_calling/MACS")),
    genome_size = 'hs',
    n_cpu = 10,
    input_format = 'BEDPE',
    shift = 73,
    ext_size = 146,
    keep_dup = 'all',
    q_value = 0.05,
    _temp_dir = '/home/tmp012/'
)


from pycisTopic.iterative_peak_calling import get_consensus_peaks

peak_half_width=250

consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict = narrow_peak_dict,
    peak_half_width = peak_half_width,
    chromsizes = chromsizes)

consensus_peaks.to_bed(
    path = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed"),
    keep =True,
    compression = 'infer',
    chain = False)

#run in terminal
#mkdir -p ./B_asdas_0vs12_tp1/outs/qc
#pycistopic tss get_tss \
#    --output ./B_asdas_0vs12_tp1/outs/qc/tss.bed \
#    --name "hsapiens_gene_ensembl" \
#    --to-chrom-source ucsc \
#    --ucsc hg38


regions_bed_filename = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed")
tss_bed_filename = os.path.join(out_dir, "outs/qc", "tss.bed")

pycistopic_qc_commands_filename = "pycistopic_qc_commands.txt"

# Create text file with all pycistopic qc command lines.
with open(pycistopic_qc_commands_filename, "w") as fh:
    for sample, fragment_filename in fragments_dict.items():
        print(
            "pycistopic qc",
            f"--fragments {fragment_filename}",
            f"--regions {regions_bed_filename}",
            f"--tss {tss_bed_filename}",
            f"""--output {os.path.join(out_dir, "outs/qc")}/{sample}""",
            sep=" ",
            file=fh,
        )

#run in terminal
#cat pycistopic_qc_commands.txt | parallel -j 4 {}

from pycisTopic.qc import get_barcodes_passing_qc_for_sample
sample_id_to_barcodes_passing_filters = {}
sample_id_to_thresholds = {}
for sample_id in fragments_dict:
    (
        sample_id_to_barcodes_passing_filters[sample_id],
        sample_id_to_thresholds[sample_id]
    ) = get_barcodes_passing_qc_for_sample(
            sample_id = sample_id,
            pycistopic_qc_output_dir = out_dir + "outs/qc",
            unique_fragments_threshold = None, # use automatic thresholding
            tss_enrichment_threshold = None, # use automatic thresholding
            frip_threshold = 0,
            use_automatic_thresholds = True,
    )

mdata = mu.read("./mudata.h5mu")
mdata.obs['cell_type'] = mdata.obs.apply(lambda row: row['rna:cell_type_global'], axis=1)
mdata = mdata[mdata.obs['cell_type'].isin([ct])]
mdata = mdata[mdata.obs['atac:time_point'].isin(['1'])]

my_sample_id_to_barcodes_passing_filters = {}

for s in list(mdata['atac'].obs['sample'].drop_duplicates()):
    my_sample_id_to_barcodes_passing_filters[s] = np.array(mdata['atac'][mdata['atac'].obs['sample'].isin([s])].obs['cell_raw'])

path_to_regions = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed")
#path_to_blacklist = "pycisTopic/blacklist/hg38-blacklist.v2.bed"
pycistopic_qc_output_dir = out_dir + "outs/qc"

from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
import polars as pl

cistopic_obj_list = []
for sample_id in my_sample_id_to_barcodes_passing_filters.keys(): #fragments_dict
    sample_metrics = pl.read_parquet(
        os.path.join(pycistopic_qc_output_dir, f'{sample_id}.fragments_stats_per_cb.parquet')
    ).to_pandas().set_index("CB").loc[ my_sample_id_to_barcodes_passing_filters[sample_id] ] #sample_id_to_barcodes_passing_filters
    cistopic_obj = create_cistopic_object_from_fragments(
        path_to_fragments = fragments_dict[sample_id],
        path_to_regions = path_to_regions,
        #path_to_blacklist = path_to_blacklist,
        metrics = sample_metrics,
        valid_bc = my_sample_id_to_barcodes_passing_filters[sample_id], #sample_id_to_barcodes_passing_filters
        n_cpu = 1,
        project = sample_id,
        split_pattern = '___'
    )
    cistopic_obj_list.append(cistopic_obj)
    
from pycisTopic.cistopic_class import merge
cistopic_obj = merge(cistopic_obj_list)

mdata = mu.read("./mudata.h5mu")

cell_data = mdata.obs
cistopic_obj.cell_data['cells'] = cistopic_obj.cell_data.index
cell_data['cells_cistopic'] = cell_data.apply(lambda row: row['atac:cells'].split('-')[0] + '-1___' + "-".join(row['atac:cells'].split('-')[2:]) , axis = 1)
cell_data = cell_data.set_index('cells_cistopic')

cistopic_obj.add_cell_data(cell_data, split_pattern='___')

import pickle
pickle.dump(
    cistopic_obj,
    open(os.path.join(out_dir, "cistopic_obj.pkl"), "wb")
)


cistopic_obj.cell_data['atac:group_012'] = cistopic_obj.cell_data.apply(lambda row: row['atac:group'].replace('asdas_1', 'asdas_12').replace('asdas_2', 'asdas_12'), axis=1)


os.environ['MALLET_MEMORY'] = '200G'
from pycisTopic.lda_models import run_cgs_models_mallet
# Configure path Mallet
mallet_path="Mallet-202108/bin/mallet"
# Run models
models=run_cgs_models_mallet(
    cistopic_obj,
    n_topics=[10, 20, 30, 40, 50, 60, 65, 70],
    n_cpu=5,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    tmp_path="/tmp012/mallet",
    save_path="/tmp012/mallet",
    mallet_path=mallet_path,
)

from pycisTopic.lda_models import evaluate_models
model = evaluate_models(
    models,
    select_model = 50, ### choose the model!
    return_model = True
)


#Prepare scenic+ input

mdata = mu.read("./mudata.h5mu")

atac = mdata['atac']
rna = mdata['rna']
atac.obs['cell_cistopic'] = atac.obs.apply(lambda row: row['cell_raw'] + "___" + row['sample'], axis = 1)
atac.obs_names = list(atac.obs['cell_cistopic'])
rna.obs['cell_cistopic'] = rna.obs.apply(lambda row: row['cells'].split('-')[0] + "-1___" + row['batch'], axis = 1)
rna.obs_names = list(rna.obs['cell_cistopic'])

rna = rna[rna.obs['cell_cistopic'].isin(list(cistopic_obj.cell_data.index))]
atac = atac[atac.obs['cell_cistopic'].isin(list(cistopic_obj.cell_data.index))]

atac.obs['group_012'] = atac.obs.apply(lambda row: row['group'].replace('asdas_1', 'asdas_12').replace('asdas_2', 'asdas_12'), axis=1)
group_012 = atac.obs.set_index('cell_cistopic')['group_012'].to_dict()
rna.obs['group_012'] = rna.obs['cell_cistopic'].map(group_012)
rna
mdata = mu.MuData({'rna' : rna, 'atac' : atac})


mdata.write(out_dir + "mudata_0vs12.h5mu")
rna.write(out_dir + "rna_0vs12.h5ad")
atac.write(out_dir + "atac_0vs12.h5ad")


cistopic_obj.cell_data['atac:group_012'] = cistopic_obj.cell_data.apply(lambda row: row['atac:group'].replace('asdas_1', 'asdas_12').replace('asdas_2', 'asdas_12'), axis=1)



from pycisTopic.topic_binarization import *
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)

from pycisTopic.diff_features import *

if not os.path.exists(os.path.join(out_dir + "regions_0vs12")):
    os.makedirs(os.path.join(out_dir + "regions_0vs12"))

imputed_acc_obj = impute_accessibility(cistopic_obj, selected_cells=None, selected_regions=None, scale_factor=10**6)
normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)
variable_regions = find_highly_variable_features(normalized_imputed_acc_obj, plot = False)
markers_dict = find_diff_features(cistopic_obj, imputed_acc_obj, variable=contrast, var_features=variable_regions, split_pattern = '-')

for k in list(markers_dict.keys()):
    markers_dict[k].to_csv(out_dir + 'cistopic_' + str(k)  + '_markers_0vs12.tsv', sep = '\t')

if not os.path.exists(os.path.join(out_dir, 'candidate_enhancers_0vs12')):
    os.makedirs(os.path.join(out_dir, 'candidate_enhancers_0vs12'))
import pickle
pickle.dump(region_bin_topics_otsu, open(os.path.join(out_dir, 'candidate_enhancers_0vs12/region_bin_topics_otsu.pkl'), 'wb'))
pickle.dump(region_bin_topics_top3k, open(os.path.join(out_dir, 'candidate_enhancers_0vs12/region_bin_topics_top3k.pkl'), 'wb'))
pickle.dump(markers_dict, open(os.path.join(out_dir, 'candidate_enhancers_0vs12/markers_dict.pkl'), 'wb'))

import pickle, os
region_bin_topics_otsu = pickle.load(open(os.path.join(out_dir, 'candidate_enhancers_0vs12/region_bin_topics_otsu.pkl'), 'rb'))
region_bin_topics_top3k = pickle.load(open(os.path.join(out_dir, 'candidate_enhancers_0vs12/region_bin_topics_top3k.pkl'), 'rb'))
markers_dict = pickle.load(open(os.path.join(out_dir, 'candidate_enhancers_0vs12/markers_dict.pkl'), 'rb'))

import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['topics_top_3'] = {}
region_sets['DARs'] = {}
for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for topic in region_bin_topics_top3k.keys():
    regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for DAR in markers_dict.keys():
    regions = markers_dict[DAR].index[markers_dict[DAR].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['DARs'][DAR] = pr.PyRanges(region_names_to_coordinates(regions))

if not os.path.exists(os.path.join(out_dir + "regions_0vs12/Topics_top_3k")):
    os.makedirs(os.path.join(out_dir + "regions_0vs12/Topics_top_3k"))
if not os.path.exists(os.path.join(out_dir + "regions_0vs12/Topics")):
    os.makedirs(os.path.join(out_dir + "regions_0vs12/Topics"))
if not os.path.exists(os.path.join(out_dir + "regions_0vs12/DARs")):
    os.makedirs(os.path.join(out_dir + "regions_0vs12/DARs"))

dir_dic = {
    'topics_otsu' : out_dir + "regions_0vs12/Topics/",
    'topics_top_3' : out_dir + "regions_0vs12/Topics_top_3k/",
    'DARs' : out_dir + "regions_0vs12/DARs/"
}
for key, value in region_sets.items():
    for key1, value1 in region_sets[key].items():
        bed_df = region_sets[key][key1].as_df()
        bed_df.to_csv(dir_dic[key] + key1 + '.bed', sep='\t', header=False, index=False)


if not os.path.exists(os.path.join(out_dir, 'scplus_pipeline_0vs12')):
    os.makedirs(os.path.join(out_dir, 'scplus_pipeline_0vs12'))

os.chdir(out_dir)

#in terminal
#scenicplus init_snakemake --out_dir scplus_pipeline_0vs12
#snakemake --cores 20


