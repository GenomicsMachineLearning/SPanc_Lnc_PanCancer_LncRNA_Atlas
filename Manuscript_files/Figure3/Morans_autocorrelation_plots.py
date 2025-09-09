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
import stlearn as st

gene_matrix=np.transpose(pd.read_csv('/Users/ritanya/Downloads/PhD/Analysis/ST/QUAN/PolyA_FF/Vis15_HN/C1/genes_uTAR_mat.txt', sep="\t", header=0, index_col=0))
spatial=pd.read_csv('/Users/ritanya/Downloads/PhD/Analysis/ST/QUAN/PolyA_FF/Vis15_HN/C1/spatial_anndata', sep="\t", header=0)
data = st.create_stlearn(count=gene_matrix,spatial=spatial,library_id="HNC", scale=1,background_color="white")

#data.obs["in_tissue"]=1
#data.obs["array_row"]=
#data.obs["array_col"]=
#data.obs.append(spatial, ignore_index=True)
spatial.iloc[1:5,:]
spatial=spatial.set_index('BC')
data.obs=spatial.iloc[:,[0,1,2,3,4]]


from typing import Literal

library_id='HNC'
data.obs.iloc[1:5,:]
_QUALITY = Literal["fulres", "hires", "lowres"]
quality: _QUALITY = "hires"
scale = data.uns["spatial"][library_id]["scalefactors"][
            "tissue_" + quality + "_scalef"]

image_coor = data.obsm["spatial"] * scale

spatial["image_col"] = image_coor[:, 5]*0.115
spatial["image_row"] = image_coor[:, 4]*0.115


data.obsm["gpd"] = gpd.GeoDataFrame(spatial.iloc[:,[1,2,3,4]],geometry=gpd.points_from_xy(y=spatial.image_row,x=spatial.image_col))

## CALC MORAN'S INDEX
epr_df = data.to_df()


gene1_name="S100A10"#"NDRG1"#"SAT1"#"PRDX5"#"S100A10" #"ELF3" #"CKS1B"# "ID3" #"HNRNPD"#"FGFR3" #"SERPINB5"
#"VHL"#"PTEN"#"TMSB10"#"ITM2B" #"PKM" #"ANXA2"  #CCND1:chr19_9687799_9691349_-_3392_0
gene2_name="chr4_182814299_182820899_+_4324_0"#"chr18_54606899_54607299_-_2526_0"#"chr21_43784599_43787199_+_1362_0"#"chr18_54606899_54607299_-_2526_0"#"chr13_29654899_29655899_+_2535_0" #"chr21_43777049_43779999_-_705_0"#"chr20_51170099_51187099_+_11139_0"#"chr6_6676549_6678399_-_765_0"#"chr16_35802399_35802799_+_1642_0"#"chr4_182814299_182820899_+_4324_0" 

data.obsm["gpd"][gene1_name] = epr_df[gene1_name]
#y = pd.DataFrame(y[:,:1])
data.obsm["gpd"][gene2_name] = epr_df[gene2_name]
w = Queen.from_dataframe(data.obsm["gpd"])
x=np.array(epr_df[gene1_name]).astype(int) #need to set type to int or float
y=np.array(epr_df[gene2_name]).astype(int)

#calculate moran index
moran = Moran(y,w)
moran_loc_bv = Moran_Local_BV(y, x, w)

fig, ax = plt.subplots(figsize=(5,5))
moran_scatterplot(moran_loc_bv, p=0.05, ax=ax)
ax.set_xlabel('Gene {}'.format(gene1_name))
ax.set_ylabel('Gene {}'.format(gene2_name))
plt.tight_layout()
plt.show()

data.obsm["gpd"].to_csv('moran_spots.txt', sep="\t")

import matplotlib.image as mpimg

# plot categorical values for 2 genes on a tissue 
def plot_choropleth(gdf, 
                    attribute_1,
                    attribute_2,
                    bg_img,
                    alpha=0.5,
                    scheme='Quantiles', 
                    cmap='YlGnBu', 
                    legend=True):
    
    fig, axs = plt.subplots(2,1, figsize=(5, 8),
                            subplot_kw={'adjustable':'datalim'})
    
    # Choropleth for attribute_1
    gdf.plot(column=attribute_1, scheme=scheme, cmap=cmap,
             legend=legend, legend_kwds={'loc': 'upper left',
                                         'bbox_to_anchor': (0.92, 0.8)},
             ax=axs[0], alpha=alpha, markersize=2)
    
    axs[0].imshow(bg_img)
    axs[0].set_title('choropleth plot for {}'.format(attribute_1), y=0.8)
    axs[0].set_axis_off()
    
    # Choropleth for attribute_2
    gdf.plot(column=attribute_2, scheme=scheme, cmap=cmap,
             legend=legend, legend_kwds={'loc': 'upper left',
                                         'bbox_to_anchor': (0.92, 0.8)},
             ax=axs[1], alpha=alpha, markersize=2)
    
    axs[1].imshow(bg_img)
    axs[1].set_title('choropleth plot for {}'.format(attribute_2), y=0.8)
    axs[1].set_axis_off()
    plt.tight_layout()
    return fig, ax

img = mpimg.imread("/Users/ritanya/Downloads/PhD/Analysis/ST/QUAN/PolyA_FF/Vis15_HN/C1/spatial/tissue_hires_image.png")
plot_choropleth(data.obsm["gpd"], 
                gene1_name,
                gene2_name,
                img)
plt.show()

fig, axs = plt.subplots(1,1, figsize=(5, 8),
                            subplot_kw={'adjustable':'datalim'})

# Main spatial correlation plot 
lisa_cluster(moran_loc_bv, 
             data.obsm["gpd"], 
             p=0.05,
             figsize = (9,9),
             markersize=3,
             **{"alpha":1},
             ax = axs,
            legend_kwds={'loc': 'upper left',
                         'bbox_to_anchor': (0.92, 0.8)}
            )
axs.imshow(img)
axs.set_title(gene1_name + ' vs ' + gene2_name, y=0.75)
axs.set_axis_off()
plt.show()
fig.savefig('/Users/ritanya/Downloads/PhD/Analysis/ST/QUAN/PolyA_FF/Vis15_HN/C1/Spatial_correlation_' + gene1_name + 'vs' + gene2_name  +'.pdf', dpi=300)


