# sbg2004-cutseq
Preprocessing and cna calling pipelines used for the sbg2004 manuscript

## Preprocessing dependencies
Tools:
* alfred (v0.2.1)

Python:
* pysam (v0.16.0.1)
* pandas (v1.2.4)
* argparse (v1.4.0)

## Preprocessing
Since after basespace demultiplexing the `fastq` files still contain multiple CUTseq samples we further demultiplex using a `python` script. This script needs to be build first and this can be done in the following way: Build the python/cython library using `$ python setup.py build_ext --inplace`, dependencies for demultiplexing are: pandas, argparse and pysam.

After building the cython library you can edit the `config.yaml` that is present alongside the `snakefile`, making sure that all the paths to the directories, executables and the reference genome are correct. 

## GATK copynumber calling
First of, make sure that gatk 4.1.4.0 is installed and present in your $PATH to be able to run the snakefile. If you have multiple/other versions of gatk installed you can change the `gatk` to your prefered gatk version in the snakefile. Following this, change the paths in `config.yaml` to the correct directories and human reference (hg19) and its dictionary file and specify the required bin size (100kb in our case). After succesfully running the copynumber pipeline you will end up with the segmented copy number files ready for plotting and further analyses. 

## Data
In the `data` directory you can find all processed data files (`fastq` files are not available due to patient confidentiality) to perform and visualize all analysis used in the manuscript. There are two subdirectories, one for CNA analyses and one for TIL analyses.

## Analysis
For further analysis of copynumber calling output and TIL analysis scripts are available in the `Analysis` directory. These scripts are used to plot copy number profiles, SCNA frequencies and perform the analyses of TILs.  The output of this can be used to generate all the plots in the manuscript.


If you have any questions or if you are missing files feel free to email me at: luuk.harbers@scilifelab.se