# CAMAMED
```diff
-A pipeline for composition-aware mapping-based analysis of metagenomic data
```

This pipeline can analyze metagenomic samples at both taxonomic ‎and functional profiling levels. Using this pipeline, metagenome sequences can be mapped to non-‎redundant gene catalogs and the gene frequency in the samples are obtained. Cumulative sum-scaling ‎algorithm for compositional data analysis is used in our pipeline. Additionally, by mapping the genes to the ‎KEGG database, annotations related to each gene can be extracted at the level of KEGG ortholog (KO) ‎groups, Enzyme Commission (EC) numbers, and reactions. Furthermore, the pipeline enables the user to ‎identify potential biomarkers in case-control metagenomic samples. Using the CAMAMED pipeline, not ‎only one can easily analyze metagenome data at taxonomical (taxon) and functional (gene) level, but also ‎can go further by analyzing the potential functional differences at other functional levels, that is, KO, EC ‎number and reaction.

```diff
-To run CAMAMED, first run the following command
-     python camamed_init
```

****For more information on how to run this software, see the CAMAMED-manual.pdf file.***

Also, this software is available through two Docker images called the camamed_pipeline (without MetaPhlAn2 databases) and the camamed_pipeline_db (with MetaPhlAn2 databases) at www.hub.docker.com.
(Please refer to https://hub.docker.com/r/camamed/camamed_pipeline_db and https://hub.docker.com/r/camamed/camamed_pipeline).
Note: There is no need to install any software dependencies if you use the Docker image.


