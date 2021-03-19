import pandas as pd
import numpy as np
import networkx as nx
from networkx.readwrite import json_graph
import json
import pickle
import qvalue


data_path = "../../data/"
sub_dir = "reactome/74/"

relation_file = data_path + sub_dir + "ReactomePathwaysRelation.txt"

def generate_tree(relation_file = relation_file):
    
    rel_df = pd.read_csv(relation_file, sep = "\t", header = None)

    rel_df.columns = ['parentId', 'id']
    cut = rel_df['parentId'].str.contains('HSA') & rel_df['id'].str.contains('HSA')
    rel_df = rel_df.loc[cut]

    G = nx.DiGraph()
    G.add_edges_from(rel_df.values)
    roots = [n for n,d in G.in_degree() if d==0]

    roots_df = pd.DataFrame(np.transpose([np.repeat('HomoSapiens',len(roots)),roots]), columns = [['parentId', 'id']])
    # roots_df['id'] = roots
    # roots_df['parentId'] = 'HomoSapiens'

    roots_df = pd.DataFrame(roots_df.values, columns = ['parentId', 'id'])
    rel_df = pd.DataFrame(rel_df.values, columns = ['parentId', 'id'])

    tree = roots_df.append(rel_df)
    return tree

rel_df = generate_tree()



def default(o):
     if isinstance(o, np.integer): return int(o)
     raise TypeError


def sunburst(in_df, outname = 'sun_tree.json'):
    # in_df has as indexes the pathway reactome ID, 'value' as the colors to be plot, 'ngenes' as the width, and 'descr' as the description of the pathway

    in_df = in_df.loc[[x for x in in_df.index if 'HSA' in x]]

    topPaths = rel_df.loc[(rel_df['parentId'] == 'HomoSapiens'), 'id']
    homoNgenes = np.sum(in_df.loc[[x in topPaths.tolist() for x in in_df.index],'ngenes'])
    homoNode = pd.DataFrame([[1,homoNgenes,"Homo Sapiens"]], columns = ["q", "ngenes", "description"]).xs(0)
    homoNode.name = 'HomoSapiens'

    in_df = in_df.append(homoNode)

    topDict = in_df.to_dict()

    pathways = in_df.index

    n_path = len(pathways)

    subset_vec = [x in pathways for x in rel_df.iloc[:,0]] and [x in pathways for x in rel_df.iloc[:,1]] 
    sub_rel_df = rel_df[subset_vec] 

    G = nx.DiGraph()

    G.add_nodes_from(pathways)
    G.add_edges_from(sub_rel_df.values)

    tree = nx.algorithms.dag.dag_to_branching(G)

    secondDict = nx.get_node_attributes(tree,'source')
    
    thirdDict = {'value':{}, 'ngenes':{}, 'description': {}, 'timelineHigh': {},'timelineLow': {}, 'high': {}, 'low': {}}
    for key, value in secondDict.items():
        thirdDict['value'].update({key : topDict['q'][value]})
        thirdDict['ngenes'].update({key : topDict['ngenes'][value]})
        thirdDict['description'].update({key : topDict['description'][value]})

    nx.set_node_attributes(tree, thirdDict['value'], name = 'value')
    nx.set_node_attributes(tree, thirdDict['ngenes'], name = 'ngenes')
    nx.set_node_attributes(tree, thirdDict['description'], name = 'description')

    root = [v for v, d in tree.in_degree() if d == 0][0]
    out_json = json_graph.tree_data(tree, root)

    with open(outname, 'w') as outfile:
        json.dump(out_json, outfile, default=default)

if __name__ == "__main__":

    survival = pickle.load(open('results/metabric_path_survival.p', 'rb'))
    survival.index = [x.replace('_','-') for x in survival.index]
    activities = pickle.load(open('results/metabric_path_activities.p', 'rb'))
    survival['q'] = qvalue.qvalues(survival)

    in_df = pd.DataFrame(survival['q'])
    in_df['ngenes'] = activities['set_size'].astype(int)
    in_df['description'] = activities['annotation']
    in_df = in_df.dropna()

    sunburst(in_df, outname='docs/results.json')

    survival = pickle.load(open('results/metabric_path_survival_DiseasesMitoticCellCycle.p', 'rb'))
    survival.index = [x.replace('_','-') for x in survival.index]
    activities = pickle.load(open('results/metabric_path_activities.p', 'rb'))
    survival['q'] = qvalue.qvalues(survival)

    in_df = pd.DataFrame(survival['q'])
    in_df['ngenes'] = activities['set_size'].astype(int)
    in_df['description'] = activities['annotation']
    in_df = in_df.dropna()

    sunburst(in_df, outname='docs/results_DiseasesMitoticCellCycle.json')
