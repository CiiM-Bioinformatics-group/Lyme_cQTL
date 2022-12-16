#Running all the figures

cd ../code_figures
Rscript Figures.R
echo 'Figures.R done'

Rscript loci_dotplot.R
echo 'loci_dotplot.R done'

Rscript Locuszoom.R
echo 'Locuszoom.R done'

Rscript Consistency.R
echo 'Consistency.R'

Rscript Disease_relevance.R
echo 'Disease_relevance.R done'
