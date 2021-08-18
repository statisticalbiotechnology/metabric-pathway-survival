# metabric-pathway-survival

* run_analysis.py will load the metabric dataset, perform the preliminary analysis and save the results. Edit file for the correct location of the dataset.
* permutation_test.py and set_permutation.py will do the same, but for the phenotype permutation and gene set permutation respectively. Since they are computationally intensive, they were left on separate scripts.
* plots.py will load the results of the analysis and produce the plots. Edit to comment out the permutation plot if these results are not present.
* generate_sunburst.py will also load the results and produce a two json files that are then used on the interactive sunburst plots.
