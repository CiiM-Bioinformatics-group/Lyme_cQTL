#Annotate tables and add loci numbers to the manhattan file

Rscript resources/Table_genes.R
echo 'Table_genes finished'

Rscript resources/DE_genes.R
echo 'DE_genes finished'

Rscript resources/extract_pval_gw.R
echo 'extract pval done'

Rscript resources/finemapping.R
echo 'finemapping done'

Rscript resources/Group_output_manhattan.R
echo 'manhattan file done'

sh code_figures/Figures.sh
