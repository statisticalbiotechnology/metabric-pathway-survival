import pandas as pd
import numpy as np
import sklearn as skl
import porch
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn  import preprocessing
from sklearn.metrics import pairwise_distances
import pickle
import lifelines
import qvalue

from run_analysis import load_metabric, illumina2ensembl_dictionary, get_reactome_illumina, return_pathway

def return_PC1(pathway_data):
    pca = PCA(n_components = 1)
    pca.fit(pathway_data.dropna())
    return pca.components_[0], pca.explained_variance_ratio_[0]

def check_stability(pathway_data, nsamples = 10, fraction = 0.5):
    pcs = []
    for i in range(nsamples): 
        sample = pathway_data.sample(frac=fraction)
        PC1, exp_var = return_PC1(sample)
        pcs.append(PC1)
    distances = pairwise_distances(pcs, metric = 'cosine')
    correct_direction = np.vectorize(lambda x: np.min([x, 2-x]))
    distances = correct_direction(distances)
    return [np.mean(distances), np.std(distances)]

if __name__ == "__main__":
    ## Load data

    metabric_path = '../../data/metabric'
    illumina2ensembl_path = 'data/illumina2ensembl.txt'

    data = load_metabric(metabric_path)
    expression_df = data.iloc[8:,:]
    metadata_df = data.iloc[:8,:]

    illumina_reactome_df = get_reactome_illumina(illumina2ensembl_path)

    survival = pickle.load(open('results/metabric_path_survival.p', 'rb'))
    survival.index = [x.replace('_','-') for x in survival.index]

    genes_cox = pickle.load(open('results/metabric_gene_survival.p', 'rb'))

    ## Calculate the mean pairwise distances

    stability_check_df = pd.DataFrame(columns = ['mean', 'std', 'exp_var', 'ngenes'])
    for pathway in np.unique(illumina_reactome_df['reactome_id']):
        print('Stability test' + pathway)
        pathway_data = return_pathway(data, illumina_reactome_df, pathway)
        ngenes = pathway_data.shape[1]
        pc1, exp_var = return_PC1(pathway_data)
        if pathway_data.shape[1] > 0:
            try:
                stability_check = check_stability(pathway_data, nsamples = 10, fraction = 0.2)
                stability_check_df.loc[pathway] = [stability_check[0], stability_check[1], exp_var, ngenes]
            except:
                continue


    ## Save the results of the stability check
    pickle.dump(stability_check_df, open("results/stability_check_df.p", "wb"))

    ## Load the previously saved results
    stability_check_df = pickle.load(open('results/stability_check_df.p', 'rb'))
    stability_check_df['p'] = survival['p']


    ## Plot the results of the stability check

    fig, ax = plt.subplots()

    ax.plot(stability_check_df['mean'], stability_check_df['p'], '.', alpha = 0.5)
    ax.set_yscale('log')
    ax.set_ylabel('Cox p-value')
    ax.set_xlabel('Mean pairwise cosine distance')
    plt.savefig('plots/stability_p.png')

    ## In here we try a new kind of plot to show survival and gene correlation

    metadata = metadata_df
    reactome_df = illumina_reactome_df
    reactome_id = 'R-HSA-196757' 

    pathway_data = return_pathway(data, reactome_df, reactome_id)
    vec, exp = return_PC1(pathway_data)
    sort_df = pd.DataFrame(vec, index=pathway_data.columns)
    sort_seq = sort_df.sort_values(0).index.tolist()

    x = preprocessing.StandardScaler().fit_transform(pathway_data)
    pathway_data.loc[:,:] = x
    genes = pathway_data.columns
    pathway_data['T'] = metadata.T['T']
    pathway_data['E'] = (metadata.T['last_follow_up_status'] != 'a')*1

    pathway_data = pathway_data.where(pathway_data['E'] == 1).dropna()

    pathway_data_quantiles = pathway_data.groupby(pd.cut(pathway_data['T'], pathway_data['T'].quantile(np.arange(0.1,1,0.1)))).mean()

    pathway_data_quantiles['id'] = ['(' + str(np.round(x.left/365, 1)) + ', ' + str(np.round(x.right/365, 1)) + ']' for x in pathway_data_quantiles.index]

    pathway_data_long = pd.wide_to_long(pathway_data_quantiles, stubnames='ILMN', sep = '_', i='id', j='probe')
    pathway_data_long = pathway_data_long.reset_index(level=1)
    pathway_data_long['T'] = [np.round(x.mid) for x in pathway_data_long.index]
    pathway_data_long['T'] = pathway_data_long.index


    order = [int(x[5:]) for x in sort_seq]
    name = reactome_df.where(reactome_df['reactome_id'] == reactome_id).dropna()['reactome_name'].iloc[0]

    chart = sns.pointplot(x = 'probe', y = 'ILMN', hue = 'T', data=pathway_data_long, order = order, palette='vlag', scale = 0.5, errwidth=0)
    plt.legend(bbox_to_anchor=(1, 1), title = 'Survival time (years)')
    plt.ylabel('Mean normalized expression')
    plt.xlabel('Gene')
    chart.set_xticklabels(chart.get_xticklabels(), rotation=90)
    plt.plot(sort_df.sort_values(0).values, color = 'y', linestyle = '--', linewidth = 2)
    plt.title(name)

    ax = plt.gca()
    labels = ['ILMN_' + x.get_text() for x in ax.get_xticklabels()]
    dictionary_pd = illumina2ensembl_dictionary('data/illumina2hugo.txt').set_index('probe')
    labels_gene = [np.unique(dictionary_pd.loc[x,'gene'])[0] for x in labels]

    chart.set_xticklabels(labels_gene)

    plt.tight_layout()
    plt.savefig('plots/survial_eigenpatient.png') 

    #### We compare patways against their most predictive gene

    genes_cox['q'] = qvalue.qvalues(genes_cox)

    for path in survival.index:
        print('Compare' + path)
        genes = return_pathway(genes_cox,illumina_reactome_df,path).T
        max_p = genes.sort_values('p').iloc[0]['p']
        max_q = genes.sort_values('q').iloc[0]['q']
        ngenes = genes.shape[0]
        survival.loc[path,'max_probe_p'] = max_p
        survival.loc[path,'max_probe_q'] = max_q
        survival.loc[path,'ngenes'] = ngenes

    survival['q'] = qvalue.qvalues(survival, pcolname='p')

    plot = sns.scatterplot(x = survival['p'],
            y = survival['max_probe_p'],
            alpha = 0.5)
    plot.set(xscale = 'log', yscale = 'log', xlabel = 'Pathway p value', ylabel = 'Minimum probe p value')
    plt.plot([1e-24,1], [1e-24,1], linewidth=1, color = 'r')
    plt.savefig('plots/pathway_gene_p.png')


    plot = sns.scatterplot(x = survival['q'],
            y = survival['max_probe_q'],
            alpha = 0.5)
    plot.set(xscale = 'log', yscale = 'log', xlabel = 'Pathway $q$ value', ylabel = 'Minimum probe $q$ value')
    plt.plot([1e-21,1], [1e-21,1], linewidth=1, color = 'r')
    plt.savefig('plots/pathway_gene_q.png')


    ### compare concordance indexes

    survival_cross_pathways = pickle.load(open('results/metabric_path_cross.p', 'rb'))
    survival_cross_genes = pickle.load(open('results/metabric_gene_cross.p', 'rb'))


    survival_cross_pathways['mean'] = np.mean(survival_cross_pathways.values, axis=1)
    survival_cross_genes['mean'] = np.mean(survival_cross_genes.values, axis=1)

    bins=np.histogram(np.hstack((survival_cross_genes['mean'],survival_cross_pathways['mean'])), bins=50)[1]

    plt.hist(survival_cross_pathways['mean'], bins, density = True, alpha = 0.5)
    plt.hist(survival_cross_genes['mean'], bins, density=True, alpha = 0.5)
    plt.legend(['Pathways','Transcripts'])
    plt.xlabel('Concordance Index')
    plt.ylabel('Density')
    plt.savefig('plots/concodance.png')

    #### Permutation test

    permutation_results = pickle.load(open('results/permutation_results.p', 'rb'))
    cox_results = pickle.load(open('results/metabric_path_survival.p', 'rb'))

    for pathway, row in permutation_results.iterrows():
        z = np.abs(row['base'])
        perms = np.abs(row['perms'])
        num_higher = sum(x >= z for x in perms)
        p = num_higher/len(perms)
        cox_results.loc[pathway,'p_perms'] = p

    cox_results = cox_results.dropna()

    plt.scatter(cox_results['p'],cox_results['p_perms'], alpha = 0.5)
    plt.xlabel('Cox derived p value')
    plt.ylabel('Permutation derived p value')
    plt.yscale('log')
    plt.xscale('log')
    lim = np.min([np.min(cox_results['p']), np.min(cox_results['p_perms'])])
    plt.plot([1,lim],[1,lim], color='r', alpha = 0.5)
    plt.savefig('plots/permutation.png')
