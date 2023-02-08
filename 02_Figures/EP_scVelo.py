# Script for scVelo

import sys
import scanpy as sc
import scvelo as scv
import numpy as np
import matplotlib.pyplot as plt

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualization

adata = scv.read("/home/mazhuo/Rproject/fetal_panc/EP.h5ad")
adata
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata, mode="stochastic")
scv.tl.velocity_graph(adata, n_jobs=16)
scv.pl.velocity_embedding_stream(adata, basis='umap',color="cell_type",
                                 groups=["EP early","EP mid","EP late","EP alpha","Alpha/PP","Epsilon","EP beta","Beta","Delta"],
                                 dpi=300,palette=["#c6dbef","#6baed6","#3182bd","#fdae6b","#e6550d","#636363","#a1d99b","#31a354","#756bb1"])
