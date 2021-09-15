# MetaClock

MetaClock is an integrated framework for reconstructing time-resolved evolutionary history for microbial genomes using large-scale (meta)genomics data from ancient and contemporary samples.<br />

![MetaClock](https://github.com/SegataLab/metaclock/blob/master/images/MetaClock_Logo.png "MetaClock")<br />
## Main features

* Automated reconstruction of whole-genome alignment from large-scale metagenomics data using multiple reference genomes
* Automated authentication for ancient DNA origin, damaged sites removal, SNV analysis etc.
* Rich utilities for manual curation to enhance the quality of phylogenetic and molecular dating analysis

 


# Installation

## Bioconda

You can install MetaClock in an isolated environment using conda as follows:

~~~Bash
conda create -n mc -c bioconda metaclock
~~~

~~~Bash
conda activate mc
~~~

## Private conda channel

Temporarily, you can install MetaClock through a private conda channel. But this option will be replaced once bioconda channel is ready.

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
Please post your issues in our [metaclock issues section](https://github.com/SegataLab/metaclock/issues)

# Citation
If you think our tool is helpful for your research please cite https://github.com/SegataLab/metaclock. Note: citation source will be updated soon when manuscript is published. Thank you for your patience!  
