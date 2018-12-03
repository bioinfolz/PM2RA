# PM2CA
Profile Monitoring for Microbiome Community Alteration

Version 1.0

Updated date: 2018.12.03

Cite us:

＃Author and License

Author: Zhi Liu

Email: liuzhi@njum.edu.cn

Licensed under the GNU Affero General Public License version 3 or later

##DEPENDENCIES

R 3.5.1(Other version not tested)
R packages:
gtools,MASS,sfsmisc,pracma,DMwR,qcc,stringr,ggplot2,igraph,RColorBrewer,pheatmap,reshape2,

##USAGE
This PM2CA script is developed with two functional modules.
1. 2D scanning(2DScan):  scanning for pairwise association change among the microbiome community between two conditions, and the two-component community alterations form a community alteration(CA) network, where the nodes represent microbiota and the edges represent the measure of the community alteration (or association changes) of the two connected microbiota by PM score. Hub microbiotas in the CA network are the ones with extensively altered associations between two compared conditions.

2. Multiple-dimensional search(MultiSearch):  where the PM score of any defined sub-community with two or more microbiotas could be calculated.

#2DScan usage
Six parameters should be set for 2DScan function
Rscript PM2CA.R abundance.table.csv ConditionA-Label ConditionB-Label prevalence.filter project.name 2DScan
1. abundance table: the abundance table of microbiota in csv format, where the first colunm indicates the labels for samples from two conditions.
label	Bifidobacterium_animalis	Bifidobacterium_breve	Bacteroides_ovatus	Parabacteroides_distasonis
H	1.64E-05	0.000212566	0.318188189	0.005378391
H	0.000106283	0.001448253	0.142495246	0.011803262
H	0.000205558	0.000363231	0.022839183	0.023493233
H	5.72E-05	0.001193641	0.053880851	0.082181349
H	0.000629523	0.001410879	0.027530822	0.030052417
H	0.000157673	0.002483053	0.005934333	0.009705631
H	0.000196215	0.000908662	0.016441175	0.018664944
D	0.00030133	0.001388688	0.019018832	0.010586262
D	0.000416957	0.000828074	0.045962177	0.043511827
D	0.000153001	0.001358321	0.042005176	0.030979767
D	0.00013665	0.003584426	0.032091651	0.039965943
D	0.000140154	0.000798875	0.041534494	0.01880393

2.ConditionA-Label: label for samples from conditionA(e.g. H,disease,A)
3.ConditionB-Label: label for samples from conditionB(e.g. D,control,B)
4.prevalence filter: remove microbiotas with prevalence less than the prevalence filter(e.g. 0.1)
5.project.name: prefix attached to results files
6.function type: 2DScan

Run Example
Rscript PM2CA.R Example.csv H D 0.1 2DscanExample 2DScan


＃MultiSearch Usase
Seven  parameters should be set for MultiSearch function. Except for those set for 2DScan, a microbiota list should be provided for MultiSearch.

Rscript PM2CA.R abundance.table.csv ConditionA-Label ConditionB-Label prevalence.filter project.name MultiSearch micriobita.list
1. abundance table: the abundance table of microbiota in csv format, where the first colunm indicates the labels for samples from two conditions.
label	Bifidobacterium_animalis	Bifidobacterium_breve	Bacteroides_ovatus	Parabacteroides_distasonis
H	1.64E-05	0.000212566	0.318188189	0.005378391
H	0.000106283	0.001448253	0.142495246	0.011803262
H	0.000205558	0.000363231	0.022839183	0.023493233
H	5.72E-05	0.001193641	0.053880851	0.082181349
H	0.000629523	0.001410879	0.027530822	0.030052417
H	0.000157673	0.002483053	0.005934333	0.009705631
H	0.000196215	0.000908662	0.016441175	0.018664944
D	0.00030133	0.001388688	0.019018832	0.010586262
D	0.000416957	0.000828074	0.045962177	0.043511827
D	0.000153001	0.001358321	0.042005176	0.030979767
D	0.00013665	0.003584426	0.032091651	0.039965943
D	0.000140154	0.000798875	0.041534494	0.01880393

2.ConditionA-Label: label for samples from conditionA(e.g. H,disease,A)
3.ConditionB-Label: label for samples from conditionB(e.g. D,control,B)
4.prevalence filter: remove microbiotas with prevalence less than the prevalence filter(e.g. 0.1)
5.project.name: prefix attached to results files
6.function type: MultiSearch
7.microbita.list: a plain text file provides a list of micriobita，where each line contains one microbiota.







