# snakemake_for_tassel
usage: 
  - set the correct parameters in the file config.yaml
  - On a sge environnement, you can use the file cluster.yaml
  - command: snakemake --cluster-config cluster.yaml --cluster 'qsub -V -b y -q {cluster.queue} -l mem_free={cluster.memory} -N {cluster.name}' -j

Caution: some of tassel parameters are not configure yet in the snakefile, but it's fonctionnal for most of the cases.
