# Ovarian cancer Visium Spatial transcriptomics
All of the code needed to recreate the analysis for the paper titled [Small cell carcinoma of the ovary hypercalcemic type (SCCOHT): A review and novel case with dual germline SMARCA4 and BRCA2 mutations](https://doi.org/10.1016/j.gore.2022.101077)

This includes a `snakemake` pipeline to run the first steps of the analysis on an lsf cluster and R scripts (numbered based on the analysis order) that can be run locally. R scripts were run in `R` version 4.0.3 and all packages needed are in the `renv.lock` file and can be loaded directly using the `renv` package.

All scripts in this repository were writen by Kristen Wells. Any use of these scripts should cite

```
1.	Sanders B.E., Wolsky R., Doughty E.S., Wells K.L., Ghosh D., Ku L., Pressey J. G., Bitler B.B., Brubaker L.W. Small cell carcinoma of the ovary hypercalcemic type (SCCOHT): A review and novel case with dual germline SMARCA4 and BRCA2 mutations. Gynecologic Oncology Reports. 2022; 44:101077
```

The raw data can be downloaded from GEO [GSE213699](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213699)

To use the `snakemake` pipeline:

1. Download all scripts from this repository:

```bash
git clone https://github.com/kwells4/visium_ovarian_cancer.git
```

2. Download and install miniconda3: For Linux
```bash
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh bash Miniconda3-latest-Linux-x86_64.sh
```

3. Install packages (I like using conda):

`Mamba` is a very fast package manager that relies on the `conda` framework

```bash
conda install -c conda-forge mamba
```

`Snakemake` install via `mamba`
```bash
mamba install snakemake -c bioconda -c conda-forge
```

Install [`spaceranger`](https://support.10xgenomics.com/spatial-gene-expression/software/downloads/latest) and add it to your path

4. Update the config file (config.yaml) 
>* RAW_DATA: Location of the raw data (downloaded from GEO [GSE213699](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE213699))
>* IMAGES: Path to the images. The images are included in this repository, so you shouldn't need to change this path unless you move the images.
>* SAMPLES: the list of samples you want to test. This is the name that will be in the output files. The samples needs to be a layered dictionary that contains the slide information: First position is the slide name, the second is the position. These have already been updated for the samples analyzed here, but you can change this if you add additional samples.
>* RESULTS: Path to the output directory
>* TRANSCRIPTOME: Path to the cellranger genome reference (download from [10x genomics](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest))
>* PROBE_SET: The path to the probe set downloaded from [10x geomics](https://support.10xgenomics.com/spatial-gene-expression-ffpe/probe-sets/overview) Human Version 1.0 was used for this project.
>* SLIDE_PATH: Path to the slide files. These are included in the repository and this should not be changed unless you moved the slide files.
>* LSF_TEMPLATE: Path to an LSF template. One is included in this git repo.

4. Update snakecharmer.sh to your specific cluster specs. 
>* change the -q argument to the queue you want to use 

5. submit the job using `bsub < snakecharmer.sh`

6. I highly recommend looking at the csv files that are generated and passed to cell ranger to ensure that the correct fastq files have been detected for each sample.
