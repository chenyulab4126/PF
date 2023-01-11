# Description of DEG analysis

Differentially expressed genes were called using DESeq2 package (v1.26.0). 
Considering that COVID-19/SLE and IPF samples have different sequencing depths, 
we used different empirical parameters for these two data sets. 

For COVID-19/SLE data, we used FC (Fold Change) > 3, padj (adjusted P-value) < 0.05 for identifying differentially expressed genes. 
And for IPF data, the following thresholds are used: FC = 2, padj = 0.05. 

