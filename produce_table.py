import pickle
import pandas as pd
import qvalue

results = pickle.load(open('results/metabric_path_survival.p','rb'))
results_proliferation = pickle.load(open('results/metabric_path_survival_DiseasesMitoticCellCycle.p','rb'))
results_perm = pickle.load(open('set_perm_statistics.p','rb'))


results_proliferation.index = [x.replace('_','-') for x in results_proliferation.index]
results.index = [x.replace('_','-') for x in results.index]

activity = pickle.load(open('results/metabric_path_activities.p', 'rb'))

results['q'] = qvalue.qvalues(results)
out1 = results[['p','q']].copy()
out1['Name'] = activity['annotation']
out1.to_csv('results/ST1.tsv', sep = '\t')


out2 = results_proliferation[['p','q']].copy()
out2['Name'] = activity['annotation']
out2.to_csv('results/ST2.tsv', sep = '\t')

results_perm['q'] = qvalue.qvalues(results_perm, pcolname='p_perms')
out3 = results_perm[['p_perms','q']].copy()
out3.columns = ['p','q']
out3['Name'] = activity['annotation']
out3.to_csv('results/ST3.tsv', sep = '\t')
