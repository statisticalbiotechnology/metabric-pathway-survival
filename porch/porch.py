# -*- coding: utf-8 -*-
import multiprocessing
import numpy as np
import pandas as pd
from numpy.linalg import svd
from sklearn.preprocessing import StandardScaler
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
from wpca import WPCA
import os.path
import sys
import patsy
from lifelines import CoxPHFitter
from typing import *

import porch.cache as cache

def porch_single_process(expression_df: pd.DataFrame,
                                          geneset_df: pd.DataFrame,
                                          gene_column: str = "gene",
                                          set_column: str = "pathway",
                                          annot_column: str = "reactome_name",
                        ) -> Tuple[pd.DataFrame, Dict[str,Dict[str,float]], List]:
    """
    Calculates pathway activities from the expression values of analytes,
    with a grouping given by a pathway definition.
    This call is functional equivalent to the porch function, with the difference that it is not using
    parallel processing. The function is  mostly intended for debugging purposes.

    Args:
        expression_df (pd.DataFrame): The DataFrame of the expression values we analyse. These values are logtransformed and subsequently standardized before analysis
        geneset_df (pd.DataFrame): The DataFrame of the pathway definitions.
        gene_column (str): The name of the column within geneset_df containing names of analytes.
        set_column (str): The name of the column within geneset_df containing names of pathways.
        annot_column (str): The name of the column within geneset_df containing annotation of pathways.

    Returns:
        tuple(pd.DataFrame, list): tuple containing:
            - **activity_df** (*pd.DataFrame*): A pandas DataFrames activity_df, containing the pathway activity values for each sample and pathway.
            - **untested** (*list*): a list of the pathway that were not possible to decompose, due to shortage of data in expression_df.
    """
    # expression_df = expression_df[phenotype_df.columns]
    results_df = pd.DataFrame()
    set_df = geneset_df[[gene_column, set_column,annot_column]]
    set_of_all_genes = set(expression_df.index)
    setnames, set_annots, set_sizes, untested, activities, eigen_samples = [], [], [], [], [], {}
    for setname, geneset in set_df.groupby([set_column]):
        genes = list(set(geneset[gene_column].tolist()) & set_of_all_genes)
        annot = (geneset[annot_column].iloc[0] if len(geneset.index)>0 else "Default")
        setname, set_annot, set_size, activity, eigen_sample_dict  = porch_proc(setname, annot, genes, expression_df)
        if activity is None:
            untested += [setname]
        else:
            setnames += [setname]
            set_annots += [set_annot]
            set_sizes += [set_size]
            activities += [activity]
            eigen_samples[setname] = eigen_sample_dict

    activity_df = pd.DataFrame(data=activities, columns=expression_df.columns, index=setnames)
    activity_df["annotation"] = set_annots
    activity_df["set_size"] = set_sizes
    return activity_df, eigen_samples, untested


def porch(expression_df: pd.DataFrame,
                geneset_df: pd.DataFrame,
                gene_column: str = "gene",
                set_column: str = "pathway",
                annot_column: str = "reactome_name",
                ) -> Tuple[pd.DataFrame, Dict[str,Dict[str,float]], List]:
    """
    Calculates pathway activities from the expression values of analytes,
    with a grouping given by a pathway definition.

    Args:
        expression_df (pd.DataFrame): The DataFrame of the expression values we analyse. These values are logtransformed and subsequently standardized befor analysis
        geneset_df (pd.DataFrame): The DataFrame of the pathway definitions.
        gene_column (str): The name of the column within geneset_df containing names of analytes.
        set_column (str): The name of the column within geneset_df containing id:s of pathways.
        annot_column (str): The name of the column within geneset_df containing annotation of pathways.

    Returns:
        Tuple(pd.DataFrame, Dict[str,Dict[str,float]], list): tuple containing:
            - **activity_df** (*pd.DataFrame*): A pandas DataFrames activity_df, containing the pathway activity values for each sample and pathway.
            - **eigen_samples** (*Dict[str,Dict[str,float]]*): A dictonary of the pathways eigen samples, i.e. represenative pattern of analyte expression values
            - **untested** (*list*): a list of the pathway that were not possible to decompose, due to shortage of data in expression_df.
    """
    set_df = geneset_df[[gene_column, set_column,annot_column]]
    set_of_all_genes = set(expression_df.index)
    call_args = []
    for setname, geneset in set_df.groupby([set_column]):
        genes = list(set(geneset[gene_column].tolist()) & set_of_all_genes)
        annot = (geneset[annot_column].iloc[0] if len(geneset.index)>0 else "Default")
        call_args += [(setname, annot, genes, expression_df)]
    print("Processing with {} parallel processes".format(os.cpu_count()), file=sys.stderr)
    setnames, set_annots, set_sizes, untested, activities, eigen_samples = [], [], [], [], [], {}
    with multiprocessing.Pool() as executor:
        for setname, set_annot, set_size, activity, eigen_sample_dict in executor.starmap(porch_proc,  call_args):
            if activity is None:
                untested += [setname]
            else:
                setnames += [setname]
                set_annots += [set_annot]
                set_sizes += [set_size]
                activities += [activity]
                eigen_samples[setname] = eigen_sample_dict
    activity_df = pd.DataFrame(data=activities, columns=expression_df.columns, index=setnames)
    activity_df["annotation"] = set_annots
    activity_df["set_size"] = set_sizes
    return activity_df, eigen_samples, untested


def porch_proc(setname, set_annotation, genes, expression_df,keep_feature_stdv=True):
    """ Core processing node of porch. Takes the analysis from expression values to significance testing. """
    # print("Decomposing " + setname, file=sys.stderr)
    expr = expression_df.loc[genes]
    expr.dropna(axis=0, how='any', inplace=True)
    expr = expr.loc[~(expr<=0.0).any(axis=1)]
    if expr.shape[0]>2:
        try:
            standardizer = StandardScaler(with_std=keep_feature_stdv)
            log_data = np.log(expr.values.T.astype(float))
            standard_log_data = standardizer.fit_transform(log_data).T
            eigen_genes, eigen_samples = decomposition_method(standard_log_data)
            eigen_sample_dict = dict(zip(expr.index,eigen_samples))
            return setname, set_annotation, len(genes), eigen_genes, eigen_sample_dict
        except ValueError:
            pass
#   print("Not enough data to evaluate " + setname, file=sys.stderr)
    return setname, set_annotation, len(genes), None, None

def porch_reactome(expression_df: pd.DataFrame,
                                 organism: str = "HSA",
                                 gene_anot: str = "Ensembl") -> Tuple[pd.DataFrame, List]:
    """
    Download the Reactome database and subsequently call porch.

    Args:
        expression_df (pd.DataFrame): The DataFrame of the expression values we analyse. These values are logtransformed and subsequently standardized befor analysis
        organism (str): The three letter reactome abriviation of organism, e.g. HSA or MMU
        gene_anot (str): Reactome name of row annotation, e.g. Ensembl or ChEBI

    Returns:
        tuple(pd.DataFrame, list): tuple containing:
            - **activity_df** (*pd.DataFrame*): A pandas DataFrames activity_df, containing the pathway activity values for each sample and pathway.
            - **untested** (*list*): a list of the pathway that were not possible to decompose, due to shortage of data in expression_df.
    """
    reactome_df = get_reactome_df(organism, gene_anot)
#    return porch_single_process(expression_df, reactome_df, "gene", "reactome_id")
    return porch(expression_df, reactome_df,
        "gene", "reactome_id", "reactome_name")

def porch_multi_reactome(expression_df,list_of_expression_annotations):
    """ Download the Reactome database and subsequently call porch. """
    reactome_df = None
    for organism, gene_anot in list_of_expression_annotations:
        r_df = get_reactome_df(organism, gene_anot)
        if reactome_df is None:
            reactome_df = r_df
        else:
            reactome_df.append(r_df)
    return porch_single_process(expression_df, reactome_df, "gene", "reactome_id", "reactome_name")
#    return porch(expression_df, reactome_df,
#        "gene", "reactome_id")

def wpca_decomposition(data):
    weights = 0. + np.isfinite(data)
    kwds = {'weights': weights}
    pca = WPCA(n_components=1).fit(data, **kwds)
    eigen_samples = pca.transform(data)[:,0]
    eigen_genes = pca.components_[0,:]
    return eigen_genes, eigen_samples

def svd_decomposition(data):
    U, S, Vt = svd(data, full_matrices=False)
    eigen_genes = (Vt.T)[:,0]
    eigen_samples = U[:,0]
    return eigen_genes, eigen_samples

decomposition_method = svd_decomposition
#decomposition_method = wpca_decomposition

def linear_model(test,activity_df,phenotype_df):
    """
    Applies a linear model, test, that is a function of variables in phenotype_df,
    to each row in the activity_df.

    Args:
        activity_df (pd.DataFrame): The DataFrame of the pathway activity values we analyse.
        phenotype_df (pd.DataFrame): The DataFrame containing any sample oriented variables that are included in the model.
        test (str): linear model that should be tested. The model should contain the variable Pathway, that will be replaced with each pathway's activity.

    Returns:
        results_df
        A pandas DataFrames results_df, containing the output of the significance tests
    """
    expression_df = activity_df.copy()
    phenotype_df =  phenotype_df[[ col for col in phenotype_df.columns  if col in expression_df.columns]]
    expression_df = expression_df[phenotype_df.columns]
    significance_df = expression_df.apply(applicable_linear_model,
        axis=1, result_type='reduce',
        args=(test,phenotype_df))
    significance_df["annotation"] = activity_df["annotation"]
    significance_df["set_size"] = activity_df["set_size"]
    return significance_df

def applicable_linear_model(row,test,phenotype_df):
    phenotype_df.loc["Pathway"] = row.values
    lm = ols(test, phenotype_df.T).fit()
    try:
        pvals = anova_lm(lm)["PR(>F)"].T.iloc[:-1]
        #pvals.rename(row.name)
    except ValueError:
        pvals = None
    return pvals

reactome_fn = "2Reactome_All_Levels.txt"
#reactome_fn = "UniProt2Reactome_All_Levels.txt"
reactome_url = "https://reactome.org/download/current/"

def get_reactome_df(organism = "HSA", gene_anot = "Ensembl"):
    fn = gene_anot + reactome_fn
    path = os.path.join(cache.get_cache_path(),fn)
    url = reactome_url + fn
    reactome_df = pd.read_csv(cache.download_file(path, url),
                        sep='\t',
                        header=None,
                        usecols=[0,1,3],
                        names=["gene","reactome_id","reactome_name"])
    organism = "R-" + organism
    reactome_df = reactome_df[reactome_df["reactome_id"].str.startswith(organism) ]
    return reactome_df

def read_triqler(file_name):
    """Code for reading a protein.tsv file from triqler"""
    pid_col, first_dat_col = 2, 7
    proteins, data = [], []
    with open(file_name) as infile:
        header = infile.readline().split('\t')
        last_dat_col = len(header) - 1
        col_names = [w.split(':')[2] for w in header[first_dat_col:last_dat_col]]
        phen_values = [[int(w.split(':')[0]) for w in header[first_dat_col:last_dat_col]]]
        for line in infile.readlines():
            words = line.split('\t')
            proteins += [words[pid_col]]
            data += [[np.exp2(float (w)) for w in words[first_dat_col:last_dat_col]]]
    values_df = pd.DataFrame(index=proteins, columns=col_names, data=data)
    phenotype_df = pd.DataFrame(index=["SampleGroup"], columns=col_names, data=phen_values)
    return values_df, phenotype_df

def survival(row, phenotype_df, duration_col = 'T', event_col = 'E', other_cols = []):
    """
    duration_col: survival time
    event_col: whether an event (death or other) has ocured or not. 0 for no, 1 for yes
    other_cols: other variables to consider in the regression
    """
    phenotype_df = phenotype_df.T
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
    cph = CoxPHFitter()
    cph.fit(X, duration_col = duration_col, event_col = event_col)
    result = cph.summary.loc[row.name]
    return result
