# SingleViper
Accurate annotation of cell types from single-cell RNA sequencing (scRNA-seq) data remains
challenging under conditions of low transcript depth and high technical noise. Traditional anno-
tation pipelines rely on transcript abundance, which poorly reflects cellular state when genes are
sparsely expressed. To address this, we present SingleViper, a workflow that integrates transcrip-
tional regulatory network inference and protein activity estimation to classify cell types based on
regulatory dynamics instead of simply gene expression alone. Using ARACNe3, we infer a regulon
from reference expression matrices, which is then used with the VIPER algorithm to transform
query and reference gene expression data into protein activity scores. The resulting activity matri-
ces serves as inputs for supervised classification using either SingleR or a custom machine learning
classifier. To test the pipeline’s robustness, we downsample the query dataset and benchmark
purity, accuracy, and per-class of SingleViper against standard gene-expression–based annotation
methods.
