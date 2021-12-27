import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import pickle

from porch_npcs import qvalue

survival_cross_genes = pickle.load(open('results/bolt/metabric_gene_cross.p', 'rb'))
survival_cross_pathways = pickle.load(open('results/bolt/metabric_path_cross.p', 'rb'))
survival_cross_pathways_2 = pickle.load(open('results/bolt/metabric_path_cross_2pcs.p', 'rb'))
survival_cross_pathways_3 = pickle.load(open('results/bolt/metabric_path_cross_3pcs.p', 'rb'))
survival_cross_pathways_4 = pickle.load(open('results/bolt/metabric_path_cross_4pcs.p', 'rb'))
survival_cross_pathways_5 = pickle.load(open('results/bolt/metabric_path_cross_5pcs.p', 'rb'))
survival_cross_pathways_10 = pickle.load(open('results/metabric_path_cross_10pcs.p', 'rb'))
survival_cross_pathways_50 = pickle.load(open('results/metabric_path_cross_50pcs.p', 'rb'))
survival_cross_pathways_100 = pickle.load(open('results/metabric_path_cross_100pcs.p', 'rb'))


survival_cross_genes['mean'] = np.mean(survival_cross_genes.values, axis=1)
survival_cross_pathways['mean'] = np.mean(survival_cross_pathways.values, axis=1)
survival_cross_pathways_2['mean'] = np.mean(survival_cross_pathways_2.values, axis=1)
survival_cross_pathways_3['mean'] = np.mean(survival_cross_pathways_3.values, axis=1)
survival_cross_pathways_4['mean'] = np.mean(survival_cross_pathways_4.values, axis=1)
survival_cross_pathways_5['mean'] = np.mean(survival_cross_pathways_5.values, axis=1)
survival_cross_pathways_10['mean'] = np.mean(survival_cross_pathways_10.values, axis=1)
survival_cross_pathways_50['mean'] = np.mean(survival_cross_pathways_50.values, axis=1)
survival_cross_pathways_100['mean'] = np.mean(survival_cross_pathways_100.values, axis=1)

bins=np.histogram(np.hstack((
    survival_cross_genes['mean'], 
    survival_cross_pathways['mean'], 
    survival_cross_pathways_2['mean'],
    survival_cross_pathways_3['mean'],
    survival_cross_pathways_4['mean'],
    survival_cross_pathways_5['mean'], 
    survival_cross_pathways_10['mean'], 
    survival_cross_pathways_50['mean'],
    survival_cross_pathways_100['mean']
    )), bins=50)[1]

# fig, (ax1,ax2) = plt.subplots(2,1, sharey=True, sharex=True, figsize=(7,10))


plt.hist(survival_cross_genes['mean'], bins, density=True, alpha = 0.5, histtype='step', fill=False)
plt.hist(survival_cross_pathways['mean'], bins, density = True, alpha = 0.5, histtype='step', fill=False)
plt.hist(survival_cross_pathways_2['mean'], bins, density=True, alpha = 0.5, histtype='step', fill=False)
plt.hist(survival_cross_pathways_3['mean'], bins, density=True, alpha = 0.5, histtype='step', fill=False)
plt.hist(survival_cross_pathways_4['mean'], bins, density=True, alpha = 0.5, histtype='step', fill=False)
plt.hist(survival_cross_pathways_5['mean'], bins, density=True, alpha = 0.5, histtype='step', fill=False)
plt.hist(survival_cross_pathways_10['mean'], bins, density=True, alpha = 0.5, histtype='step', fill=False)
plt.hist(survival_cross_pathways_50['mean'], bins, density=True, alpha = 0.5, histtype='step', fill=False)
plt.hist(survival_cross_pathways_100['mean'], bins, density=True, alpha = 0.5, histtype='step', fill=False)
plt.legend(['Transcripts','1sv','2sv','4sv', '50sv'], loc = 'upper left')
plt.xlabel('Concordance Index')
plt.ylabel('Density')

plt.show()

sns.set_context('talk')
plt.figure(figsize=(10, 8))

sns.histplot(data = survival_cross_pathways['mean'], bins = bins, stat = 'density', element='step', fill=False)
sns.histplot(data = survival_cross_pathways_2['mean'], bins = bins, stat = 'density', element='step', fill=False)
sns.histplot(data = survival_cross_pathways_3['mean'], bins = bins, stat = 'density', element='step', fill=False)
sns.histplot(data = survival_cross_pathways_4['mean'], bins = bins, stat = 'density', element='step', fill=False)
sns.histplot(data = survival_cross_pathways_5['mean'], bins = bins, stat = 'density', element='step', fill=False)

plt.legend(['1sv','2sv','5sv'], loc = 'upper left')
plt.xlabel('Concordance Index')
plt.ylabel('Density')

plt.savefig('cross_npc.png')
plt.show()

#######

survival_genes = pickle.load(open('results/metabric_gene_survival.p', 'rb'))
survival_ssgsea = pickle.load(open('results/ssgsea_survival.p', 'rb'))
survival = pickle.load(open('results/metabric_path_survival.p', 'rb'))
survival_2 = pickle.load(open('results/metabric_path_survival_2pcs.p', 'rb'))
survival_3 = pickle.load(open('results/metabric_path_survival_3pcs.p', 'rb'))
survival_4 = pickle.load(open('results/metabric_path_survival_4pcs.p', 'rb'))
survival_5 = pickle.load(open('results/metabric_path_survival_5pcs.p', 'rb'))


bins=np.histogram(np.hstack(( 
    survival_genes['p'],
    survival['p'], 
    survival_2['p'], 
    survival_3['p'], 
    survival_4['p'], 
    survival_5['p']
    )), bins=50)[1]

sns.set_context('talk')

plt.hist(survival_genes['p'], bins, density=True, alpha = 0.5, histtype='step', fill=False)
plt.hist(survival['p'], bins, density=True, alpha = 0.5, histtype='step', fill=False)
plt.hist(survival_2['p'], bins, density=True, alpha = 0.5, histtype='step', fill=False)
plt.hist(survival_3['p'], bins, density=True, alpha = 0.5, histtype='step', fill=False)
plt.hist(survival_4['p'], bins, density=True, alpha = 0.5, histtype='step', fill=False)
plt.hist(survival_5['p'], bins, density=True, alpha = 0.5, histtype='step', fill=False)

plt.legend(['T','1','2','3', '4', '5'])
plt.xlabel('p value')
plt.ylabel('Density')

plt.show()

###########

survival = pickle.load(open('results/metabric_path_survival.p', 'rb'))
survival_2 = pickle.load(open('results/metabric_path_survival_2pcs_sep.p', 'rb'))
survival_3 = pickle.load(open('results/metabric_path_survival_3pcs_sep.p', 'rb'))
survival_4 = pickle.load(open('results/metabric_path_survival_4pcs_sep.p', 'rb'))
survival_5 = pickle.load(open('results/metabric_path_survival_5pcs_sep.p', 'rb'))


bins=np.histogram(np.hstack(( 
    survival['p'], 
    survival_2['p'], 
    survival_3['p'], 
    survival_4['p'], 
    survival_5['p']
    )), bins=50)[1]


plt.hist(survival['p'], bins, density=True, alpha = 0.5, histtype='step', fill=False)
plt.hist(survival_2['p'], bins, density=True, alpha = 0.5, histtype='step', fill=False)
plt.hist(survival_3['p'], bins, density=True, alpha = 0.5, histtype='step', fill=False)
plt.hist(survival_4['p'], bins, density=True, alpha = 0.5, histtype='step', fill=False)
plt.hist(survival_5['p'], bins, density=True, alpha = 0.5, histtype='step', fill=False)

plt.legend(['1sv','2sv_sep','3sv_sep', '4sv_sep', '5sv_sep'])
plt.xlabel('pval')
plt.ylabel('Density')
plt.title('Each sv regressed separatelly')

plt.show()



########


survival_genes['q'] = qvalue.qvalues(pd.DataFrame(survival_genes['p']))['q']
survival_ssgsea['q'] = qvalue.qvalues(pd.DataFrame(survival_ssgsea['p']))['q']
survival['q'] = qvalue.qvalues(pd.DataFrame(survival['p']))['q']
survival_2['q'] = qvalue.qvalues(pd.DataFrame(survival_2['p']))['q']
survival_3['q'] = qvalue.qvalues(pd.DataFrame(survival_3['p']))['q']
survival_4['q'] = qvalue.qvalues(pd.DataFrame(survival_4['p']))['q']
survival_5['q'] = qvalue.qvalues(pd.DataFrame(survival_5['p']))['q']


survival_genes = survival_genes.sort_values('p')
survival_ssgsea = survival_ssgsea.sort_values('p')
survival = survival.sort_values('p')
survival_2 = survival_2.sort_values('p')
survival_3 = survival_3.sort_values('p')
survival_4 = survival_4.sort_values('p')
survival_5 = survival_5.sort_values('p')

survival_genes['frac'] = np.arange(1,len(survival_genes['p'])+1)/len(survival_genes['p'])
survival_ssgsea['frac'] = np.arange(1,len(survival_ssgsea['p'])+1)/len(survival_ssgsea['p'])
survival['frac'] = np.arange(1,len(survival['p'])+1)/len(survival['p'])
survival_2['frac'] = np.arange(1,len(survival_2['p'])+1)/len(survival_2['p'])
survival_3['frac'] = np.arange(1,len(survival_3['p'])+1)/len(survival_3['p'])
survival_4['frac'] = np.arange(1,len(survival_4['p'])+1)/len(survival_4['p'])
survival_5['frac'] = np.arange(1,len(survival_5['p'])+1)/len(survival_5['p'])


survival = survival.loc[survival['q'] <= 0.1]
survival_2 = survival_2.loc[survival_2['q'] <= 0.1]
survival_3 = survival_3.loc[survival_3['q'] <= 0.1]
survival_4 = survival_4.loc[survival_4['q'] <= 0.1]
survival_5 = survival_5.loc[survival_5['q'] <= 0.1]


sns.set_context('talk')
plt.figure(figsize=(10, 8))
sns.lineplot(x='q', y='frac', data = survival)
sns.lineplot(x='q', y='frac', data = survival_ssgsea)
sns.lineplot(x='q', y='frac', data = survival_2)
sns.lineplot(x='q', y='frac', data = survival_3)
sns.lineplot(x='q', y='frac', data = survival_4)
sns.lineplot(x='q', y='frac', data = survival_5)
sns.lineplot(x='q', y='frac', data = survival_genes)


plt.legend(['1','2','3', '4', '5', 'Transcripts'])
plt.legend(['1','2','3', '4', '5'])
plt.legend(['Eigengene','ssGESEA'])
plt.xlabel('q value')
plt.ylabel('Fraction below threshold')
plt.title('A')

plt.xlim([None,0.1])
plt.savefig('p_q_multiple_sv.png')
plt.show()

##########

survival = pickle.load(open('results/metabric_path_survival_DiseasesMitoticCellCycle.p', 'rb'))
survival_2 = pickle.load(open('results/metabric_path_survival_DiseasesMitoticCellCycle_2pcs.p', 'rb'))
survival_3 = pickle.load(open('results/metabric_path_survival_DiseasesMitoticCellCycle_3pcs.p', 'rb'))
survival_4 = pickle.load(open('results/metabric_path_survival_DiseasesMitoticCellCycle_4pcs.p', 'rb'))
survival_5 = pickle.load(open('results/metabric_path_survival_DiseasesMitoticCellCycle_5pcs.p', 'rb'))

survival['q'] = qvalue.qvalues(pd.DataFrame(survival['p']))['q']
survival_2['q'] = qvalue.qvalues(pd.DataFrame(survival_2['p']))['q']
survival_3['q'] = qvalue.qvalues(pd.DataFrame(survival_3['p']))['q']
survival_4['q'] = qvalue.qvalues(pd.DataFrame(survival_4['p']))['q']
survival_5['q'] = qvalue.qvalues(pd.DataFrame(survival_5['p']))['q']


survival = survival.sort_values('p')
survival_2 = survival_2.sort_values('p')
survival_3 = survival_3.sort_values('p')
survival_4 = survival_4.sort_values('p')
survival_5 = survival_5.sort_values('p')

survival['frac'] = np.arange(1,len(survival['p'])+1)/len(survival['p'])
survival_2['frac'] = np.arange(1,len(survival_2['p'])+1)/len(survival_2['p'])
survival_3['frac'] = np.arange(1,len(survival_3['p'])+1)/len(survival_3['p'])
survival_4['frac'] = np.arange(1,len(survival_4['p'])+1)/len(survival_4['p'])
survival_5['frac'] = np.arange(1,len(survival_5['p'])+1)/len(survival_5['p'])

survival = survival.loc[survival['q'] <= 0.1]
survival_2 = survival_2.loc[survival_2['q'] <= 0.1]
survival_3 = survival_3.loc[survival_3['q'] <= 0.1]
survival_4 = survival_4.loc[survival_4['q'] <= 0.1]
survival_5 = survival_5.loc[survival_5['q'] <= 0.1]


plt.figure(figsize=(10, 8))
sns.set_context('talk')
sns.lineplot(x='q', y='frac', data = survival)
sns.lineplot(x='q', y='frac', data = survival_2)
sns.lineplot(x='q', y='frac', data = survival_3)
sns.lineplot(x='q', y='frac', data = survival_4)
sns.lineplot(x='q', y='frac', data = survival_5)


plt.legend(['1','2','3', '4', '5'])
plt.xlabel('q value')
plt.ylabel('Fraction below threshold')
plt.title('Controled for Diseases of Mitotic Cell Cycle')
plt.title('B')

plt.xlim([None,0.1])
plt.savefig('p_q_multiple_sv_controlled.png')
plt.show()

###############

survival_cross_sep = pickle.load(open('results/metabric_path_cross_5pcs_sep.p', 'rb'))


survival_cross_sep['mean'] = np.mean(survival_cross_sep.values, axis=1)

survival_cross_pathways = survival_cross_sep.loc[[x for x in survival_cross_sep.index if x.endswith(('1'))]]
survival_cross_pathways_2 = survival_cross_sep.loc[[x for x in survival_cross_sep.index if x.endswith(('1', '2'))]]
survival_cross_pathways_3 = survival_cross_sep.loc[[x for x in survival_cross_sep.index if x.endswith(('1', '2', '3'))]]
survival_cross_pathways_4 = survival_cross_sep.loc[[x for x in survival_cross_sep.index if x.endswith(('1', '2', '3', '4'))]]
survival_cross_pathways_5 = survival_cross_sep.loc[[x for x in survival_cross_sep.index if x.endswith(('1', '2', '3', '4', '5'))]]



bins=np.histogram(np.hstack((
    survival_cross_pathways['mean'], 
    survival_cross_pathways_2['mean'],
    survival_cross_pathways_3['mean'],
    survival_cross_pathways_4['mean'],
    survival_cross_pathways_5['mean'], 
    )), bins=50)[1]


sns.set_context('talk')
plt.figure(figsize=(10, 8))

sns.histplot(data = survival_cross_pathways['mean'], bins = bins, stat = 'density', element='step', fill=False)
sns.histplot(data = survival_cross_pathways_2['mean'], bins = bins, stat = 'density', element='step', fill=False)
sns.histplot(data = survival_cross_pathways_3['mean'], bins = bins, stat = 'density', element='step', fill=False)
sns.histplot(data = survival_cross_pathways_4['mean'], bins = bins, stat = 'density', element='step', fill=False)
sns.histplot(data = survival_cross_pathways_5['mean'], bins = bins, stat = 'density', element='step', fill=False)

plt.legend(['1sv','2sv','5sv'], loc = 'upper left')
plt.xlabel('Concordance Index')
plt.ylabel('Density')

plt.savefig('cross_npc_sep.png')
plt.show()


