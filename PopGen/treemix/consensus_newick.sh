

source
conda activate /ibex/scratch/projects/c2016/SANTOSJ_DIR/biopython_env

newick_replicates_dir=run1/
format=newick
output=majority_tree.nex


python -u treemix_consenus.py --rundir $newick_replicates_dir \
-o $output \
-m $format
