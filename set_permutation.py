import numpy as np
import pandas as pd
import pickle
import multiprocessing
from joblib import Parallel, delayed

from run_analysis import load_metabric, illumina2ensembl_dictionary, get_reactome_illumina, return_pathway
import porch

def permute_set(set_size):
    try:
        set_genes = np.random.choice(genes, set_size)
        setname, set_annot, set_size, activity, eigen_sample_dict  = porch.porch_proc('setname', 'annot', set_genes, expression_df)
        activity_df = pd.DataFrame(data=activity, index=expression_df.columns, columns= ['set' + str(set_size)]).T
        rowresult = porch.survival(activity_df.iloc[0], metadata_df, duration_col = 'T', event_col = 'E')
        return rowresult['z']
    except:
        return 'error'


if __name__ == "__main__":
    metabric_path = '../../data/metabric'
    illumina2ensembl_path = '../../data/reactome/illumina2ensembl.txt'

    data = load_metabric(metabric_path)
    expression_df = data.iloc[8:,:]
    metadata_df = data.iloc[:8,:]

    metadata_df.loc['E'] = (metadata_df.loc['last_follow_up_status'] == 'd-d.s.')*1

    genes = expression_df.index

    activities = pickle.load(open('results/metabric_path_activities.p', 'rb'))
    set_sizes = np.unique(activities['set_size'])

    num_cores = multiprocessing.cpu_count() - 1
    num_perm = 10**4

    results = pd.DataFrame(columns=np.arange(num_perm))
    
    i = 0
    for set_size in set_sizes:
        i+=1
        print('Permutations for set size ' + str(set_size) + ', ' + str(i) + '/' + str(len(set_sizes)))
        permutaions = Parallel(n_jobs = num_cores, prefer="threads")(delayed(permute_set)(set_size) for i in range(num_perm))
        results.loc[set_size] = permutaions 
        pickle.dump(results, open('results/set_permutation_results.p', 'wb'))

###
    
    survival = pickle.load(open('results/metabric_path_survival.p', 'rb'))
    survival.index = [x.replace('_','-') for x in survival.index]
    survival['ngenes'] = activity['set_size']


    for pathway, row in survival.iterrows():
        z = np.abs(row['z'])
        perms = perm_results.loc[row['ngenes']]
        perms = np.abs(perms[perms != 'error'])
        num_higher = sum(x > z for x in perms)
        p = num_higher/len(perms)
        survival.loc[pathway, 'p_perms'] = p

    survival['p_perms'] = np.round(survival['p_perms'], decimals=4)
    survival['annotation'] = activity['annotation'] 
    survival[['annotation', 'p_perms']].sort_values('p_perms').to_csv('permutation_p_results.csv')


