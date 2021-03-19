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
from joblib import Parallel, delayed
import multiprocessing


from run_analysis import load_metabric, illumina2ensembl_dictionary, get_reactome_illumina, return_pathway


def permute_pathway(activity):
    shuffle.columns = np.random.permutation(metadata_df.columns) 
    rowresult = porch.survival(activity, shuffle, duration_col = 'T', event_col = 'E')
    return rowresult['coef']

def parallel_permutations(activity_row, metadata_df, n, num_cores = 7):
    activity = activity_row[:-2]
    results = Parallel(n_jobs = num_cores)(delayed(permute_pathway)(activity) for i in range(n))
    return results


if __name__ == "__main__":

    ## Load data

    metabric_path = '../../data/metabric'
    illumina2ensembl_path = '../../data/reactome/illumina2ensembl.txt'

    data = load_metabric(metabric_path)
    expression_df = data.iloc[8:,:]
    metadata_df = data.iloc[:8,:]
    # reactome_df = get_reactome_illumina(illumina2ensembl_path)
    activity_df = pickle.load(open('metabric_path_activities.p', 'rb'))
    activity_df.index = [x.replace('-','_') for x in activity_df.index]
    survival = pickle.load(open('metabric_path_survival.p', 'rb'))

    metadata_df.loc['E'] = (metadata_df.loc['last_follow_up_status'] == 'd-d.s.')*1
    shuffle = metadata_df.loc[['E','T']].copy()

    num_cores = multiprocessing.cpu_count() - 1

    n_tests = 100
    p_min = -5
    # p_list = np.arange(0,p_min,p_min/n_tests)
    p_list = np.random.uniform(p_min,0,n_tests)
    p_list = np.sort(10**p_list)[::-1]

    permutation_results = pd.DataFrame(columns=['base','perms'])

    for p_val in p_list:
        nearest = np.abs(survival['p'] - p_val).sort_values().index[0]
        activity_row = activity_df.loc[nearest]

        n_pemutation = np.ceil((1/survival.loc[nearest, 'p']) * 10).astype(int)

        print('Permutation test for ' + str(nearest) + ', number of permutations:' + str(n_pemutation))
        perm_results = parallel_permutations(activity_row, metadata_df, n_pemutation, num_cores=num_cores)

        permutation_results.loc[nearest, 'base'] = survival.loc[nearest,'coef']
        permutation_results.loc[nearest, 'perms'] = perm_results

        pickle.dump(permutation_results, open("results/permutation_results.p", "wb"))
