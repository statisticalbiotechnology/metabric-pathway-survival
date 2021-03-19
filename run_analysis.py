import pandas as pd
import numpy as np
import sklearn as skl
from sklearn.decomposition import PCA
from sklearn import preprocessing
from sklearn.metrics import pairwise_distances
import pickle
import lifelines
import patsy

import porch

def load_metabric(path, relevant_columns = ['METABRIC_ID', 'age_at_diagnosis', 'last_follow_up_status','Treatment','T', 'Pam50Subtype', 'ER.Expr', 'Her2.Expr', 'PR.Expr']):
     
    discovery_ex_path = path + '/discovery_ExpressionMatrix.txt'
    validation_ex_path = path + '/validation_ExpressionMatrix.txt'
    normals_ex_path = path + '/normals_ExpressionMatrix.txt'
    
    discovery_cl_path = path + '/table_S2_revised.txt'
    validation_cl_path = path + '/table_S3_revised.txt'
    
    discovery_ex = pd.read_csv(discovery_ex_path, sep = '\t')
    validation_ex = pd.read_csv(validation_ex_path, sep = ' ')
    normals_ex = pd.read_csv(normals_ex_path, sep = '\t')
    
    discovery_cl = pd.read_csv(discovery_cl_path, sep = '\t')
    validation_cl = pd.read_csv(validation_cl_path, sep = '\t')
    
    # relevant_columns = ['METABRIC_ID', 'age_at_diagnosis', 'last_follow_up_status','Treatment','T', 'Pam50Subtype', 'ER.Expr', 'Her2.Expr', 'PR.Expr']
    
    discovery_cl = discovery_cl[relevant_columns].set_index('METABRIC_ID')
    validation_cl = validation_cl[relevant_columns].set_index('METABRIC_ID')
    
    expression = discovery_ex.join(validation_ex)
    clinical = discovery_cl.append(validation_cl)
    
    data = clinical.join(expression.T).T
     
    return data

def illumina2ensembl_dictionary(path):
    """
    This function needs a file with Illumina probe IDs and Ensembl IDs downloaded from Biomart
    """
    dictionary_pd = pd.read_csv(path, sep = '\t').dropna()
    dictionary_pd.columns = ['gene','probe']
    # dictionary = dictionary_pd.groupby(['gene','probe']).size()
    return dictionary_pd

def get_reactome_illumina(illumina2ensembl_path, organism = 'HSA'):

    reactome_df = porch.get_reactome_df(organism = "HSA", gene_anot = "Ensembl")
    illumina2ensembl = illumina2ensembl_dictionary(illumina2ensembl_path)
    ill2reac = illumina2ensembl.groupby('gene')

    reactomeIllumina = []
    for index, row in reactome_df.iterrows():
        gene = row['gene']
        try:
            probes = ill2reac.get_group(gene)['probe']
        except:
            continue
        for probe in probes.tolist():
            new_row = [probe, row['reactome_id'], row['reactome_name']]
            reactomeIllumina.append(new_row)

    reactomeIllumina_df = pd.DataFrame(reactomeIllumina, columns=['probe', 'reactome_id', 'reactome_name'])

    return reactomeIllumina_df

def return_pathway(data, reactome_df, reactome_id):
    genes = np.unique(reactome_df.where(reactome_df['reactome_id'] == reactome_id)['probe'].dropna().tolist())
    return data.loc[genes].T

def porch_metabric(expression_df, reactome_df):
    exp_df = expression_df.copy()
    reactome_df = get_reactome_illumina(illumina2ensembl_path) 
    activity, eigen_samples, untested  = porch.porch_single_process(exp_df, reactome_df, gene_column = 'probe', set_column = 'reactome_id')
    return activity

def metabric_survival(activity_df, metadata_df):
    metadata_df['E'] = (metadata_df['last_follow_up_status'] == 'd-d.s.')*1
    result = pd.DataFrame()
    for index, row in activity_df.iterrows():
        print("Cox's regression: " + str(index))
        rowresult = porch.survival(row, metadata_df, duration_col = 'T', event_col = 'E')
        result = result.append(rowresult)
    return result

def metabric_survival_control(activity_df, metadata_df, control):
    metadata_df['E'] = (metadata_df['last_follow_up_status'] == 'd-d.s.')*1
    metadata_df[control] = metadata_df[control].astype(float)
    result = pd.DataFrame()
    for index, row in activity_df.iterrows():
        print("Cox's regression: " + str(index) + ", controlled for " + control)
        try:
            rowresult = porch.survival(row, metadata_df, duration_col = 'T', event_col = 'E', other_cols=[control])
        except:
            rowresult = pd.Series(np.repeat(1,10), index=rowresult.index, name=index)
        result = result.append(rowresult)
    return result

def survival_cross_validation(row, phenotype_df, duration_col = 'T', event_col = 'E', other_cols = []):
    """
    duration_col: survival time
    event_col: whether an event (death or other) has ocured or not. 0 for no, 1 for yes
    other_cols: other variables to consider in the regression
    """
    phenotype_df = phenotype_df.join(row.astype(float))
    phenotype_df[duration_col] = phenotype_df[duration_col].astype(float)
    phenotype_df[event_col] = phenotype_df[event_col].astype(int)

    # The following lines deal with char conflicts in patsy formulas
    duration_col = duration_col.replace(' ','_').replace('.','_').replace('-','_')
    event_col = event_col.replace(' ','_').replace('.','_').replace('-','_')
    other_cols = [x.replace(' ','_').replace('.','_').replace('-','_') for x in other_cols]
    row.name = row.name.replace(' ','_').replace('.','_').replace('-','_')
    phenotype_df.columns = [x.replace(' ','_').replace('.','_').replace('-','_') for x in phenotype_df.columns]

    formula = row.name + ' + ' + duration_col + ' + ' + event_col
    if not not other_cols:
        other_cols = [x.replace(' ','_').replace('.','_') for x in other_cols]
        formula = formula + ' + ' + ' + '.join(other_cols)
    X = patsy.dmatrix(formula_like = formula, data = phenotype_df, return_type = 'dataframe')
    X = X.drop(['Intercept'], axis = 1)
    cph = lifelines.CoxPHFitter()
    scores = lifelines.utils.k_fold_cross_validation(cph, X, duration_col=duration_col, event_col=event_col, k=5, scoring_method='concordance_index')
    # cph.fit(X, duration_col = duration_col, event_col = event_col)
    # result = cph.summary.loc[row.name]
    return scores

def metabric_cross_validation(activity_df, metadata_df):
    metadata_df['E'] = (metadata_df['last_follow_up_status'] == 'd-d.s.')*1
    result = pd.DataFrame(columns=['1','2','3','4','5'])
    for index, row in activity_df.iterrows():
        print("Cox cross validaton: " + str(index)) 
        rowresult = survival_cross_validation(row, metadata_df, duration_col = 'T', event_col = 'E')
        print(rowresult)
        result.loc[index] = rowresult
    return result


if __name__ == "__main__":
    metabric_path = '../../data/metabric'
    illumina2ensembl_path = '../../data/reactome/illumina2ensembl.txt'

    data = load_metabric(metabric_path)
    expression_df = data.iloc[8:,:]
    metadata_df = data.iloc[:8,:].T

    illumina_reactome_df = get_reactome_illumina(illumina2ensembl_path)

    activity = porch_metabric(expression_df, illumina_reactome_df)
    pickle.dump(activity, open("results/metabric_path_activities.p", "wb"))

    survival = metabric_survival(activity.iloc[:,:-2], metadata_df)
    pickle.dump(survival, open("results/metabric_path_survival.p", "wb"))

    activity = pickle.load(open("results/metabric_path_activities.p", 'rb'))
    metadata_df['DiseasesMitoticCellCycle'] = activity.iloc[:,:-2].loc['R-HSA-9675126'].T
    survival_DiseasesMitoticCellCycle = metabric_survival_control(activity.iloc[:,:-2], metadata_df, 'DiseasesMitoticCellCycle')
    pickle.dump(survival_DiseasesMitoticCellCycle, open("results/metabric_path_survival_DiseasesMitoticCellCycle.p", "wb"))

    survival_genes = metabric_survival(expression_df, metadata_df)
    pickle.dump(survival_genes, open("results/metabric_gene_survival.p", "wb"))

    survival_cross_pathways = metabric_cross_validation(activity.iloc[:,:-2], metadata_df)
    pickle.dump(survival_cross_pathways, open("results/metabric_path_cross.p", "wb"))

    survival_cross_genes = metabric_cross_validation(expression_df, metadata_df)
    pickle.dump(survival_cross_genes, open("results/metabric_genes_cross.p", "wb"))
