import joblib
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
#!pip install pysal
import pysal
from pysal.explore import esda
import pysal.lib as lps
from esda.moran import Moran, Moran_Local, Moran_BV, Moran_Local_BV
import splot
from splot.esda import moran_scatterplot, plot_moran, lisa_cluster, plot_moran_bv_simulation, plot_moran_bv, plot_local_autocorrelation
from libpysal.weights.contiguity import Queen
from libpysal import examples
import numpy as np
import os
import scanpy as sc
import scipy
import stlearn as st
gene_matrix=np.transpose(pd.read_csv('gene_uTAR_mat.txt', sep="\t", header=0, index_col=0))
spatial=pd.read_csv('spatial_anndata', sep="\t", header=0)
data = st.create_stlearn(count=gene_matrix,spatial=spatial,library_id="MelD", scale=1,background_color="white")
spatial.iloc[1:5,:]
spatial=spatial.set_index('BC')
data.obs=spatial.iloc[:,[0,1,2,3,4]]
data.obsm["gpd"] = gpd.GeoDataFrame(spatial.iloc[:,[1,2,3,4]],geometry=gpd.points_from_xy(y=spatial.imagerow,x=spatial.imagecol))
#spatial
## CALC MORAN'S INDEX
epr_df = data.to_df()
#cols=['genes_sig', 'uTARs']
#genes_sig=['genes_sig']
#uTARs=['uTARs']
flist=pd.read_csv('gene_list', sep=",", header=0)
flist2=pd.read_csv('spatDE_KidneyA_sig_uTARs.csv', sep=",", header=0)
gene1_name = flist['g'].tolist()
gene2_name = flist2['g'].tolist()
gene1_name=tuple(gene1_name)
gene2_name=tuple(gene2_name)
appended_data = []
for i, j in [(i,j) for i in gene1_name for j in gene2_name]:
    data.obsm["gpd"][i] = epr_df[i]
    data.obsm["gpd"][j] = epr_df[j]
    w = Queen.from_dataframe(data.obsm["gpd"])
    x=np.array(epr_df[i]).astype(int) #need to set type to int or float
    y=np.array(epr_df[j]).astype(int)

    #calculate moran index
    moran = Moran(y,w)
    moran_loc_bv = Moran_Local_BV(y, x, w)
    
    ## ADD RESULTS TO DF
    # Assign pseudo P-values to `db`
    data.obs['p-sim'] = moran_loc_bv.p_sim
    # `1` if significant (at 5% confidence level), `0` otherwise
    sig = 1 * (moran_loc_bv.p_sim < 0.05)
    # Assign significance flag to `db`
    data.obs['sig'] = sig
    # Print top of the table to inspect
    #data.obs[['sig','p-sim']].head()
    #data.obs[['sig','p-sim']].tail()
      ## ADD LABELS TO DF
    # Pick as part of a quadrant only significant polygons, 
    # assign `0` otherwise (Non-significant polygons)
    spots = moran_loc_bv.q * sig
    # Mapping from value to name (as a dict)
    spots_labels = {
        0: 'Non-Significant', 1:'HH', 2: 'LH', 3:'LL', 4: 'HL'
    }
    # Create column in `db` with labels for each polygon
    data.obs['labels'] = pd.Series(
        # First initialise a Series using values and `db` index
        spots, index=data.obs.index
    # Then map each value to corresponding label based 
    # on the `spots_labels` mapping
    ).map(spots_labels)
    # Print top for inspection
    #data.obs['labels'].head()

    
    t=data.obs['labels'].value_counts()
    t=t.sort_index(ascending=True)
    t=t.to_frame()
    #t=t.rename(columns={'labels': i +':'+ j})
    t['sample']= i +':'+ j
    appended_data.append(t)
# see pd.concat documentation for more info
appended_data = pd.concat(appended_data)

    #t.to_csv('/Users/ritanya/Downloads/PhD/Analysis/ST/QUAN/PolyA_FF/Vis9_Skin/9B/spatial/moran'+i+':'+j+'.txt', sep="\t")
appended_data.to_csv('KA_cancergenes_topuTARs_moran.txt', sep="\t")
   
