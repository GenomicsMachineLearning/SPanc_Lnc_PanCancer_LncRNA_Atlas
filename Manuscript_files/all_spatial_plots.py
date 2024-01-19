import scipy
import scanpy as sc
import pandas as pd
import numpy as np
import scvelo as scv
import diopy
import stlearn as st
import matplotlib.pyplot as plt
sc.pl.spatial(anndata_object,  alpha_img=0.6, color =['cuTARID1','cuTARID2'],size=0.8, cmap='inferno', vmax=2,save="spatial_plot.png")
