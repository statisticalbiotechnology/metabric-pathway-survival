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

    duplicates = ["MB-0025", "MB-0196", "MB-0326", "MB-0329", "MB-0330", "MB-0335", "MB-0355", "MB-0407", "MB-0433", "MB-0547", "MB-2720", "MB-6206"]

    data = data.drop(duplicates, axis = 1)

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


    stability_check_random_df = pd.DataFrame(columns = ['mean', 'std', 'ngenes'])
    for pathway in np.unique(illumina_reactome_df['reactome_id']):
        pathway_data = return_pathway(data, illumina_reactome_df, pathway)
        ngenes = pathway_data.shape[1]
        print("Random size " + str(ngenes))
        random_data = pd.DataFrame(np.random.rand(pathway_data.shape[0], pathway_data.shape[1]))
        if pathway_data.shape[1] > 1:
            stability_check_random = check_stability(random_data, nsamples = 10, fraction = 0.2)
            print(stability_check_random)
            stability_check_random_df.loc[pathway] = [stability_check_random[0], stability_check_random[1], ngenes]


    bins=np.histogram(np.hstack((stability_check_df['mean'],stability_check_random_df['mean'])), bins=50)[1]

    plt.hist(stability_check_df['mean'], bins, density = True, alpha = 0.5)
    plt.hist(stability_check_random_df['mean'], bins, density=True, alpha = 0.5)
    plt.legend(['Sampled eigenpatients','Random'])
    plt.xlabel('Mean pairwise cosine distance')
    plt.ylabel('Density')
    plt.savefig("plots/stability_vs_random.png", bbox_inches='tight')


    stability_check_df['p'] = survival['p']


    ## Plot the results of the stability check

    fig, ax = plt.subplots()

    ax.plot(stability_check_df['mean'], stability_check_df['p'], '.', alpha = 0.5)
    ax.set_yscale('log')
    ax.set_ylabel('Cox p-value')
    ax.set_xlabel('Mean pairwise cosine distance')
    plt.savefig('plots/stability_p.png', bbox_inches='tight')

    ## In here we try a new kind of plot to show survival and gene correlation

    # metadata = metadata_df
    # reactome_df = illumina_reactome_df
    # reactome_id = 'R-HSA-196757' 

    # pathway_data = return_pathway(data, reactome_df, reactome_id)
    # vec, exp = return_PC1(pathway_data)
    # sort_df = pd.DataFrame(vec, index=pathway_data.columns)
    # sort_seq = sort_df.sort_values(0).index.tolist()

    # x = preprocessing.StandardScaler().fit_transform(pathway_data)
    # pathway_data.loc[:,:] = x
    # genes = pathway_data.columns
    # pathway_data['T'] = metadata.T['T']
    # pathway_data['E'] = (metadata.T['last_follow_up_status'] != 'a')*1

    # pathway_data = pathway_data.where(pathway_data['E'] == 1).dropna()

    # pathway_data_quantiles = pathway_data.groupby(pd.cut(pathway_data['T'], pathway_data['T'].quantile(np.arange(0.1,1,0.1)))).mean()

    # pathway_data_quantiles['id'] = ['(' + str(np.round(x.left/365, 1)) + ', ' + str(np.round(x.right/365, 1)) + ']' for x in pathway_data_quantiles.index]

    # pathway_data_long = pd.wide_to_long(pathway_data_quantiles, stubnames='ILMN', sep = '_', i='id', j='probe')
    # pathway_data_long = pathway_data_long.reset_index(level=1)
    # pathway_data_long['T'] = [np.round(x.mid) for x in pathway_data_long.index]
    # pathway_data_long['T'] = pathway_data_long.index


    # order = [int(x[5:]) for x in sort_seq]
    # name = reactome_df.where(reactome_df['reactome_id'] == reactome_id).dropna()['reactome_name'].iloc[0]

    # chart = sns.pointplot(x = 'probe', y = 'ILMN', hue = 'T', data=pathway_data_long, order = order, palette='vlag', scale = 0.5, errwidth=0)
    # plt.legend(bbox_to_anchor=(1, 1), title = 'Survival time (years)')
    # plt.ylabel('Mean normalized expression')
    # plt.xlabel('Gene')
    # chart.set_xticklabels(chart.get_xticklabels(), rotation=90)
    # plt.plot(sort_df.sort_values(0).values, color = 'y', linestyle = '--', linewidth = 2)
    # plt.title(name)

    # ax = plt.gca()
    # labels = ['ILMN_' + x.get_text() for x in ax.get_xticklabels()]
    # dictionary_pd = illumina2ensembl_dictionary('data/illumina2hugo.txt').set_index('probe')
    # labels_gene = [np.unique(dictionary_pd.loc[x,'gene'])[0] for x in labels]

    # chart.set_xticklabels(labels_gene)

    # plt.tight_layout()
    # plt.savefig('plots/survial_eigenpatient.png', bbox_inches='tight') 

    #### We compare patways against their most predictive gene

    # genes_cox['q'] = qvalue.qvalues(genes_cox)

    # for path in survival.index:
    #     print('Compare' + path)
    #     genes = return_pathway(genes_cox,illumina_reactome_df,path).T
    #     max_p = genes.sort_values('p').iloc[0]['p']
    #     max_q = genes.sort_values('q').iloc[0]['q']
    #     ngenes = genes.shape[0]
    #     survival.loc[path,'max_probe_p'] = max_p
    #     survival.loc[path,'max_probe_q'] = max_q
    #     survival.loc[path,'ngenes'] = ngenes

    # survival['q'] = qvalue.qvalues(survival, pcolname='p')

    # plot = sns.scatterplot(x = survival['p'],
    #         y = survival['max_probe_p'],
    #         alpha = 0.5)
    # plot.set(xscale = 'log', yscale = 'log', xlabel = 'Pathway p value', ylabel = 'Minimum probe p value')
    # plt.plot([1e-24,1], [1e-24,1], linewidth=1, color = 'r')
    # plt.savefig('plots/pathway_gene_p.png')


    # plot = sns.scatterplot(x = survival['q'],
    #         y = survival['max_probe_q'],
    #         alpha = 0.5)
    # plot.set(xscale = 'log', yscale = 'log', xlabel = 'Pathway $q$ value', ylabel = 'Minimum probe $q$ value')
    # plt.plot([1e-21,1], [1e-21,1], linewidth=1, color = 'r')
    # plt.savefig('plots/pathway_gene_q.png', bbox_inches='tight')


    ### compare concordance indexes

    survival_cross_pathways = pickle.load(open('results/metabric_path_cross.p', 'rb'))
    survival_cross_genes = pickle.load(open('results/metabric_gene_cross.p', 'rb'))
    ssgsea_cross = pickle.load(open('results/ssgsea_cross.p', 'rb'))


    survival_cross_pathways['mean'] = np.mean(survival_cross_pathways.values, axis=1)
    survival_cross_genes['mean'] = np.mean(survival_cross_genes.values, axis=1)
    ssgsea_cross['mean'] = np.mean(ssgsea_cross.values, axis=1)

    bins=np.histogram(np.hstack((survival_cross_genes['mean'],survival_cross_pathways['mean'],ssgsea_cross['mean'])), bins=50)[1]

    # fig, (ax1,ax2) = plt.subplots(2,1, sharey=True, sharex=True, figsize=(7,10))


    plt.hist(survival_cross_pathways['mean'], bins, density = True, alpha = 0.5, histtype='step', fill=True)
    plt.hist(survival_cross_genes['mean'], bins, density=True, alpha = 0.5, histtype='step', fill=True)
    plt.legend(['Eigengenes','Transcripts'])
    plt.xlabel('Concordance Index')
    plt.ylabel('Density')

    plt.savefig('plots/concordance_A.png', bbox_inches='tight')

    plt.hist(survival_cross_pathways['mean'], bins, density = True, alpha = 0.5, histtype='step', fill=True)
    plt.hist(ssgsea_cross['mean'], bins, density=True, alpha = 0.5, histtype='step', fill=True)
    plt.legend(['Eigengenes','ssGSEA'])
    plt.xlabel('Concordance Index')
    plt.ylabel('Density')

    plt.savefig('plots/concordance_B.png', bbox_inches='tight')

    ### correlation plot

    activities = pickle.load(open('results/metabric_path_activities.p', 'rb'))
    activities = activities.iloc[:,:-2].T

    reactome_id = 'R-HSA-196757'

    metadata = metadata_df
    reactome_df = illumina_reactome_df

    pathway_data = return_pathway(data, reactome_df, reactome_id)
    vec, exp = return_PC1(pathway_data)
    sort_df = pd.DataFrame(vec, index=pathway_data.columns)
    sort_seq = sort_df.sort_values(0).index.tolist()

    genes = pathway_data.columns
    pathway_data['T'] = metadata.T['T']
    pathway_data['E'] = (metadata.T['last_follow_up_status'] != 'a')*1
    pathway_data['Eigengene'] = activities[reactome_id] 

    pathway_data_dead=pathway_data.where(pathway_data['E'] == 1).dropna()

    #x = preprocessing.PowerTransformer(method='box-cox').fit_transform(pathway_data)
    x = preprocessing.StandardScaler().fit_transform(pathway_data_dead)
    pathway_data_dead.loc[:,:] = x

    pathway_data_dead.sort_values(inplace=True,by=['T'])
    pathway_data_dead = pathway_data_dead[['Eigengene','T']+sort_seq]
    dictionary_pd = illumina2ensembl_dictionary('data/illumina2hugo.txt'
            ).set_index('probe')
    labels_gene = [np.unique(dictionary_pd.loc[x,'gene'])[0] for x in pathway_data_dead.columns if 'ILMN_' in x]
    pathway_data_dead.columns = [ x for x in pathway_data_dead.columns if 'ILMN_' not in x ] + labels_gene
    #pathway_data_dead
    #pathway_data_dead_time = pathway_data_dead.copy()

    pathway_data_dead.rename(columns={"T":"Survival Time"}, inplace=True)

    corr = pathway_data_dead.corr()
    mask = np.triu(np.ones_like(corr, dtype=bool))

    sns.set_context("talk")
    f, ax = plt.subplots(figsize=(11, 9))

    # Generate a custom diverging colormap
    cmap = sns.diverging_palette(230, 20, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    #sns.heatmap(corr, mask=mask, cmap=cmap, vmax=.3, center=0,
    #            square=True, linewidths=.5, cbar_kws={"shrink": .5})
    sns.heatmap(corr, cmap=cmap, center=0,
            square=True, linewidths=.5, cbar_kws={"shrink": .5})

    plt.show()

    f.savefig("plots/R-HSA-196757-Corr.png", bbox_inches='tight')



    #### Permutation test

    # permutation_results = pickle.load(open('results/permutation_results.p', 'rb'))
    # cox_results = pickle.load(open('results/metabric_path_survival.p', 'rb'))

    # for pathway, row in permutation_results.iterrows():
    #     z = np.abs(row['base'])
    #     perms = np.abs(row['perms'])
    #     num_higher = sum(x >= z for x in perms)
    #     p = num_higher/len(perms)
    #     cox_results.loc[pathway,'p_perms'] = p

    # cox_results = cox_results.dropna()

    # plt.scatter(cox_results['p'],cox_results['p_perms'], alpha = 0.5)
    # plt.xlabel('Cox derived p value')
    # plt.ylabel('Permutation derived p value')
    # plt.yscale('log')
    # plt.xscale('log')
    # lim = np.min([np.min(cox_results['p']), np.min(cox_results['p_perms'])])
    # plt.plot([1,lim],[1,lim], color='r', alpha = 0.5)
    # plt.savefig('plots/permutation.png, bbox_inches='tight'')


    #### We compare patways against their most predictive gene

    compare_df = pd.DataFrame(index = survival_cross_pathways.index, columns = ['path','probes_max', 'probes_mean', 'probes_median','ngenes'])
    for path in survival_cross_pathways.index:
        print('Compare ' + path)
        genes = return_pathway(survival_cross_genes,illumina_reactome_df,path).T
        mean_genes = np.mean(genes, axis=1)
        ngenes = genes.shape[0]
        compare_df.loc[path] = [np.mean(survival_cross_pathways.loc[path]), np.max(mean_genes),np.mean(mean_genes), np.median(mean_genes), ngenes]


    fig, (ax1,ax2) = plt.subplots(1,2, sharey=True, sharex=True, figsize=(10,4))
    sns.scatterplot(ax=ax1, y = compare_df['path'],
            x = compare_df['probes_mean'],
            alpha = 1, s=5)
    ax1.set(ylabel = 'Pathway C-index', xlabel = 'Average transcript C-index')
    ax1.plot([0.45,0.67], [0.45,0.67], linewidth=1, color = 'r')

    sns.scatterplot(ax=ax2, y = compare_df['path'],
            x = compare_df['probes_max'],
            alpha = 1, s=5)
    ax2.set(ylabel = 'Pathway C-index', xlabel = 'Best transcript C-index')
    ax2.plot([0.45,0.67], [0.45,0.67], linewidth=1, color = 'r')

    fig.savefig('plots/pathway_gene_C-index.png', bbox_inches='tight')

    #######

    path_c = pickle.load(open('results/metabric_path_cross.p', 'rb'))
    gene_c = pickle.load(open('results/metabric_gene_cross.p', 'rb'))

    gene_lists = illumina_reactome_df.groupby('reactome_id')['probe'].unique()

    path_c['mean'] = np.mean(path_c, axis=1)
    gene_c['mean'] = np.mean(gene_c, axis=1)

    n=50

    top_paths = path_c.sort_values('mean', ascending=False).head(n).index
    top_genes = gene_c.sort_values('mean', ascending=False).head(n).index

    top_path_genes = np.unique(np.concatenate(gene_lists[top_paths].values))

    len([x for x in top_genes if x in top_path_genes])
