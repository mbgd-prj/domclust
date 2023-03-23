# domclust

DomClust is a hierarchical clustering tool for orthologous grouping
in multiple genomes. The method takes as input all-against-all
similarity data and classifies genes based on the traditional
hierarchical clustering algorithm UPGMA. In the course of clustering,
the method detects domain fusion or fission events, and splits
clusters into domains if required. The subsequent procedure splits
the resulting trees such that intra-species paralogous genes are
divided into different groups so as to create plausible orthologous
groups. As a result, the procedure can split genes into the domains
minimally required for ortholog grouping. 

DomClust outputs a set of hierarchical clustering trees,
but these trees may overlap with each other. The overlapping
trees actually result from the domain fusion/fission event,
and are the salient feature of the DomClust program.

DomClust has been used for classifying more than 200 microbial
genomes in MBGD (Microbial genome database for comparative analysis;
http://mbgd.genome.ad.jp/).
