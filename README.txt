Software requirement: R and MATLAB

INPUT files: all the inputs are required in csv format.
1. miRNA Expression: features in rows and samples in col (mirna.csv)
 
5.926	1.9351	4.58	5.2786
5.5466	3.6695	6.0834	5.9058
2.8668	1.484	2.585	2.4994
2.3181	0.7524	0.5025	2.6848
4.687	4.4346	2.5512	3.8461

2. mRNA Expression: features in rows and samples in col (mrna.csv)
 
3.5329	2.4308	1.7269	0
5.5684	5.5339	3.9319	5.8608
3.2185	3.3134	1.8526	3.7722
1.3281	0.8353	0.4872	1.1967
5.0155	4.9153	4.3776	5.4042
10.753	10.5206	7.9547	10.8958

3. miRNA_names.csv

x
hsa-mir-675
hsa-mir-660
hsa-mir-671
hsa-mir-615
hsa-mir-7-1

4. mRNA_names.csv 

x
ARHGEF10L
HIF3A
RNF17
RNF10
RNF11
RNF13

5. samples.csv

x
TCGA-4J-AA1J-01
TCGA-DS-A1OA-01
TCGA-ZJ-AAXF-01
TCGA-IR-A3LF-01

Execution: Please keep all the input files along with main.R and robust_weight.m in one folder. 
Follow the instructions in main.R for cancer subtypes identification.

OUTPUT: weight_miRNA.csv, weight_mRNA.csv, cancer_subtypes.csv
