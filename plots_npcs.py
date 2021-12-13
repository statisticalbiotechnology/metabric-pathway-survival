import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pickle


survival_cross_genes = pickle.load(open('results/metabric_gene_cross.p', 'rb'))
survival_cross_pathways = pickle.load(open('results/metabric_path_cross.p', 'rb'))
survival_cross_pathways_2 = pickle.load(open('results/metabric_path_cross_2pcs.p', 'rb'))
survival_cross_pathways_5 = pickle.load(open('results/metabric_path_cross_5pcs.p', 'rb'))


survival_cross_genes['mean'] = np.mean(survival_cross_genes.values, axis=1)
survival_cross_pathways['mean'] = np.mean(survival_cross_pathways.values, axis=1)
survival_cross_pathways_2['mean'] = np.mean(survival_cross_pathways_2.values, axis=1)
survival_cross_pathways_5['mean'] = np.mean(survival_cross_pathways_5.values, axis=1)

bins=np.histogram(np.hstack((survival_cross_genes['mean'], survival_cross_pathways_2['mean'],survival_cross_pathways['mean'], survival_cross_pathways_5['mean'])), bins=50)[1]

# fig, (ax1,ax2) = plt.subplots(2,1, sharey=True, sharex=True, figsize=(7,10))


plt.hist(survival_cross_genes['mean'], bins, density=True, alpha = 0.5, histtype='step', fill=True)
plt.hist(survival_cross_pathways['mean'], bins, density = True, alpha = 0.5, histtype='step', fill=True)
plt.hist(survival_cross_pathways_2['mean'], bins, density=True, alpha = 0.5, histtype='step', fill=True)
plt.hist(survival_cross_pathways_5['mean'], bins, density=True, alpha = 0.5, histtype='step', fill=True)
plt.legend(['Transcripts','1pc','2pc', '5pc'])
plt.xlabel('Concordance Index')
plt.ylabel('Density')

