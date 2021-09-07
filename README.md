# pfenningsummerproject

Hello! This repo contains the scripts I've created analyzing single cell data from the Rada Lab and figures created from these scripts. 

**Data**

The data for the NA-9 Sample is in the folder labeled "NA9". It is the output from cellranger count and specifically the "filtered_feature_bc_matrix” folder. The data for the I-A10 Sample is in the folder named "IA10" \
The "differentialexpressedgenes.csv" contains the 106 genes that were differtially expressed between the two samples along with their log fold change and p-values. 

**R Scripts**

fastqc.R has the function used to generate fastqc figures. \
NA9.R has the code for analysis done to just the NA9 sample. \
IA10.R has the code for analysis done to just the IA10 sample.\
integrated.R has code after I integrated the two datasets as well as the differential expression analysis between the two datasets :D

**Figures**

The "NA9 figures" folder has figures created from just analyzing the NA9 data. \
The "IA10 figures" folder has figures created from just analyzing the IA10 data. \
The "integrated data figures" folder has figures created from after I integrated the NA9 and IA10 togehter. It has differential expression analysis figures between the two samples here. \
The "qc htmlfiles for each sample" folder has the report from the fastqc function for each of the samples. \
The "UMAPs patchworked together" folder just has some UMAP figures that I conjoined for aesthetics on the presentation. \

