# Single-Cell Multiomic Analysis of TNF-α Inhibitor Response in Axial Spondyloarthritis
This repository contains the main analysis code for the manuscript:
"Single-Cell Multiomic Profiling of Response to TNF-α Inhibitor in Axial Spondyloarthritis"
Elizaveta Ignatova, Jae Soon Park, Jung Kyoon Choi, Hoon-Suk Cha, Seulkee Lee

01_scRNA-seq_annotation.py:
Preprocessing, filtering, and annotation of scRNA-seq and scATAC-seq data. Includes cell type assignment using CellTypist, integration with SnapATAC2, and creation of a MuData object for multiomic downstream analysis.

02_scRNA-seq_ATAC-seq_Scenic+.py:
Pipeline for eRegulon inference using pycisTopic and SCENIC+. 

03_target_gene-driven_GRN_construction_and_metric.R
Constructs target gene-driven gene regulatory networks (GRNs) using scRank on scRNA-seq data and calculates network metrics.
