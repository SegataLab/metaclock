# MetaClock

MetaClock is an integrated framework for reconstructing time-resolved evolutionary history for microbial genomes using large-scale (meta)genomics data from ancient and contemporary samples.<br />

![MetaClock](https://github.com/SegataLab/metaclock/blob/master/images/MetaClock_Logo.png "MetaClock")<br />
## Main features

* Automated reconstruction of whole-genome alignment from large-scale metagenomics data using multiple reference genomes
* Automated authentication for ancient DNA origin, damaged sites removal, SNV analysis etc.
* Rich utilities for manual curation to enhance the quality of phylogenetic and molecular dating analysis

 


# Installation

## Bioconda

You can install MetaClock using conda as follows:

~~~Bash
The package will be publically available once the manuscript is submitted
~~~

## Private conda channel

Temporarily, you can install MetaClock thorough a private conda channel. But this lead will be closed once bioconda channel is ready.

~~~Bash
conda create -n metaclock -c khuang16 -c compbiocore -c kantorlab -c bioconda -c conda-forge metaclock
~~~

Due to the dependency bowtie2 is encountering tbb-related issue, after metaclock installation is completed please manually install tbb library as described below:

~~~
conda activate metaclock
~~~

~~~
conda install tbb=2020.2
~~~

## Repository

You can clone the MetaClock repository from GitHub:

~~~Bash
git clone https://github.com/SegataLab/metaclock.git
~~~
Dependencies:



# Tutorials

* [MetaClock User manual](https://github.com/SegataLab/metaclock/wiki/User-manual)
* [MetaClock Tutorials](https://github.com/SegataLab/metaclock/wiki)


# Support


# Citation
