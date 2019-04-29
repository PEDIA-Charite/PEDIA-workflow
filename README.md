# Project PEDIA
Prioritization of Exome Data by Image Analysis (PEDIA) investigates the value of computer-assisted analysis
of medical images and clinical features in the diagnostic workup of patients with rare genetic disorders.
We integrate the facial gestalt analysis from DeepGestalt approach
(https://www.nature.com/articles/s41591-018-0279-0) in Face2Gene with the other phenotypic analysis
tools and deleteriousness scores from the molecular level to prioritize potential disease genes.

This tool is already integrated in Face2Gen LAB and https://pedia-study.org.
For the user who would like to utilize it to analyze your patients,
please contact Prof. Peter Krawitz (pkrawitz@uni-bonn.de)
in Institute for Genomic Statistics and Bioinformatics for more details.

## Contents
* [General information](#general-information)
* [Download dataset](#download-dataset)
* [Environment setup](#environment-setup)
  * [Configuration](#configuration)
  * [Required external files](required-external-files)
  * [Description of the subfolders](description-of-subfolders)
* [Running PEDIA](#running-pedia)
  * [Activate environment](#activate-environment)
  * [Usage](#usage)
  * [Example](#example)
  * [Preprocessing (Phenotypic level)](#preprocessing)
  * [CADD annotation (Genomic level)](#cadd-annotation)
    * [VCF annotation](#vcf-annotation)
	* [Variant filtering](#variant-filtering)
  * [Classification](#classification)
  * [Results](#results)

## General Information

It requires the **Key and Secret** to fetch the gestalt score and patient data from
face2gene platform. Therefore, you are not able to run the PEDIA-workflow without the Key and Secret.
However, we provide the PEDIA datasets which is described in PEDIA paper. You could run the pipeline
by downloading the datasets from [pedia-study.org](https://pedia-study.org/pedia_services/download)
on the test case below.

The whole workflow of the PEDIA project uses [snakemake](https://snakemake.readthedocs.io/) to run a pipeline together with [conda/bioconda](https://bioconda.github.io/) to install the necessary programs. So pelase get familiar with both BEFORE starting the workflow. A good start is the [snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html). 

### Requirment
   * python version >= 3.6
   * miniconda
   * snakemake
   * storage >= 30 GB (reference genomoe and external database such as ExAC and CADD)
   * **Key and secrect of your Face2Gene LAB (You are not able to use PEDIA without key and secret)**


### Install miniconda

Go have a look at the [miniconda website](https://conda.io/miniconda.html). Be sure that you choose the right version depending on your python version. 

### Install necessary software via conda
Now we will generate a shell environment with all necessary programs for downloading and process files. The software needed for the downloading is in the `environment.yaml` file. Conda can read it:

```
conda env create -f environment.yaml
```

Now we created an enviroment called `pedia`. We can activate it and we should have snakemake installed.

```
source activate pedia
snakemake -h
```
We can deactivate the environment using `source deactivate`. The command `conda env list` will list you all environments.

Now lets run the pedia download workflow. We can make a "dry run" using `snakemake -n` to see what snakemake would probably do.
```
cd data
snakemake -p -n all
```

### Setup classifier submodule
classifier is the submodule, so we have to clone it by the following command.
```
git submodule update --recursive
```

## Download dataset
Please download the PEDIA datasets we used in PEDIA paper in the following link
[https://pedia-study.org/pedia_services/download](https://pedia-study.org/pedia_services/download).
The download link requires registration to pedia-study.org.
The PEDIA datasets contain the following two data sets.
Please find more details about PEDIA dataset in wiki page ([PEDIA dataset](https://github.com/PEDIA-Charite/PEDIA-workflow/wiki/PEDIA-datasets)).
 * PEDIA cohort
 * Deep-Gestalt publication test set

## Environment setup
### Configuration

Most configuration options are in a `config.ini` file, with options commented.
A `config.ini.SAMPLE` in the project directory can be used as reference for
creating an own configuration.

#### Setup Key and Secret of your Face2Gene LAB
To fetch the patient data from your Face2Gene LAB, please put the secret and key in the following setting in config.ini
```
[your_lab_name]
; not necessary to be the same name in Face2Gene LAB
lab_id = your lab id in Face2Gene LAB
key = your key in Face2Gene LAB 
secret = your secrect in Face2Gene LAB
download_path = process/lab/bonn (the folder you would like to save the downloaded JSON files)
```

#### Setup Phenomizer account (Optional)
It requires username and password for using phenomizer.
You could still use Feature-Match from Face2Gene if you don't fill up this part.
```
[phenomizer]
; addition of phenomization scores based on HPO terms
url = 
user = 
password = 
```

### Required external files
   * Go to data folder, and run 'snakemake all' to download all necessary files such as reference genome, population data.
   ```
   cd data
   source activate pedia
   snakemake all
   ```
   * Copy corrected JSON files to process/correct/ (optional)
   * Copy hgvs_errors.json to project folder (optional)

### Description of the subfolders

* **3_simulation** - 
	Scripts and pipelie about generation background populations and simulationg the data.
	In addition the spike in is made and molecular pathogenicity scores are added to the jsons.
* **data** - 
	It contains the files from external database (dbSNP, reference genome, CADD, ExAC). Besides, data/PEDIA/json/phenomized contains the JSON files after QC and phenomization. 
* **classifier** - 
 It is a submodule of PEDIA-workflow. The repository is here. (https://github.com/PEDIA-Charite/classifier). 
 * rest of the files and folders belong to preprocessing such as lib, helper, test
	All About the quality check and phenomization

### HGVS Error dict (optional)

HGVS variant overrides are specified in `hgvs_errors.json`. Which is per-default
searched for in the project root.
The hgvs version is specified in `lib/constants.py` and will cause an error if
an hgvs errors file of not at least the specified version is found.
The number can be lowered manually to accept older hgvs error files.
A version of 0 will accept no hgvs_errors file.

### Activate environment
```
source activate pedia
```
Please check if you already have the following files before you run the pipeline.
   * Go to data folder, and run 'snakemake all' to download all necessary files such as reference genome, population data.
   * You could download the training data we used in PEDIA paper in the following link (https://pedia-study.org/pedia_services/download)
   * Copy jsons folder to 3_simulation. 3_simulation/jsons/1KG/CV/* .json will be used for training data.


## Running PEDIA

There are the following three steps in PEDIA pipeline.
The whole workflow is connected by snakemake.
You are able to get the PEDIA results by one command. Please check the [example](#example).
1. Preprocessing (perform phenomizer and map syndrome to gene)
   * Without the authentication from Face2Gene, you are not able to run this step. Please go to [Example](#example) directly.
1. CADD annotation (annotate CADD and merge with the phenotype scores)
1. Classification

### Usage
You could run PEDIA analysis on the patient in your Face2Gene LAB by the following command.
By specifying the Lab ID and Lab cases ID, you can fetch the case with phenotypic information such as
gestalt and feature-match score and HPO terms from your Face2Gene LAB.
With argument -v, you could specify the VCF file of the patient.

```
# run for a single file on whole PEDIA workflow with -v and specify the VCF file
python3 preprocess.py -l lab_name_in_config.ini --lab-case-id the_lab_case_id_of_your_case -v your_vcf_file
```

If you already downloaded the case from Face2Gene, you could use argument -s to specify the case.

```
python3 preprocess.py -s PATH_TO_FILE -v your_vcf_file
```

### Example
You could use the example in tests/data/cases/123.json and tests/data/vcfs/123.vcf.gz.
By excuting the command below, you will find the PEDIA results in classifier/output/test/1KG/123.

```
python3 preprocess.py -s tests/data/cases/123.json -v tests/data/vcfs/123.vcf.gz
```

### Preprocessing
This step is the processing on **phenotype level**.
It will donwload the data from Face2Gene LAB and perform the following:
* Convert JSON format to PEDIA format
* Perform phenomizer
* Map syndrome to gene

Since some steps depend on the existence of API keys, running the preprocess.py
script without a configuration file will **not work**.

The **preprocess.py** script contains most information necessary for running a
conversion of json files from your Face2Gene LAB to PEDIA format.

If you add ```-v your_vcf_file```, it will automatically trigger the whole workflow.

```
# get a list of usable options
./preprocess.py -h

# run complete process on all the cases in your lab
./preprocess.py -l lab_name_in_config.ini

# run complete process on a single case in your lab
./preprocess.py -l lab_name_in_config.ini --lab-case-id the_lab_case_id_of_your_case

# run for a single file (specifying output folder is beneficial)
./preprocess.py -s PATH_TO_FILE -o OUTPUT_FOLDER
```

**Output of preprocessing**
   * config.yml contains the cases passed quality check. SIMPLE_SAMPLES is the case with disease-causing mutation but without real VCF file. VCF_SAMPLES is the case with real VCF file. TEST_SAMPLE is the case with real VCF but without disease-causing mutation.
   * process/lab/lab_name is the folder of cases downloaded from LAB (in config.ini).
   * data/PEDIA/jsons/phenomized is the folder which contains the JSON files passed QC.
   * data/PEDIA/mutations/case_id.vcf.gz  is the VCF file which contains disease-causing mutations of all cases.
   * data/PEDIA/vcfs/original is the folder which contains the VCF files. In mapping.py, we rename the filename of VCF files to case_id.vcf.gz and store to ../data/PEDIA/vcfs/original/. The new filename is added in vcf field of the JSON file. For example,
   ```
   "vcf": [
           "28827.vcf.gz"
       ],
   ```

### CADD annotation
This step is the processing on **genomic level**.

To obtain the CADD scores of variants, we have to annotate the VCF files first and
further retrieve the CADD score and append it to the geneList in JSON file.

Therefore, we could further separate this steps into two parts:
* Annotating VCF file (CADD and allele frequency)
* Retrieve CADD score (Perform [variant filtering](https://github.com/PEDIA-Charite/PEDIA-workflow/wiki/Variants-filtering) and retrieve CADD score)

#### VCF annotation
The working directory for VCF annotation is in data/PEDIA/vcfs.
Please find data/PEDIA/vcfs/Snakefile for more detail.

To run the annotation, please run the following command.
```
snakemake annotated_vcfs/{case_id}_annotated.vcf.gz
```

#### Variant filtering
We filter out the variants with high allele frequency or without phenotype mapping in OMIM.
Please find [variant filtering](https://github.com/PEDIA-Charite/PEDIA-workflow/wiki/Variants-filtering) for more details.
To get the JSON file with scores from all 5 methods, please run the following command.
```
snakemake jsons/real/unknown_test/{case_id}.json
```

**Output**
* The final JSON files are in 3_simulation/json_simulation folder.
    * 3_simulation/jsons/real/unknown_test is the folder for the cases with real VCF file.
 
### Classification
1. Go to classifier folder to classify the patient. 
Train with all cases and test on patient with **unknown diagnosis**. You will find the results in output/test/1KG/case_id/.
Please find more detail in (https://github.com/PEDIA-Charite/classifier).
   ```
   snakemake output/test/1KG/123/run.out
   ```

1. Cross-validation evaluation
   * Perform 10 fold cross-validation on 1KG simulation data, please execute this command

   ```
   snakemake output/cv/CV_1KG/run.log
   ```


## Results
You will find the results in the output dir you specified in the command.

```
ls output_dir/cv_0/
# *.csv contain all five scores and pedia score for each gene in csv format
# *.json contain the PEDIA score in JSON format
# count_*.csv list the number of cases in each rank
# rank_*.csv list the rank of each case
```

**45254.csv**

Here, we listed the top ten genes in 45254.csv. You will find the five scores and PEDIA score.
The label indicates whether this gene is disease-causing gene or not.
In this case, ARID1B has the highest PEDIA score and it is the disease-causing gene of this case.

```
gene_name gene_id pedia_score feature_score cadd_score gestalt_score boqa_score pheno_score label
ARID1B    57492   4.509       0.836         25.0       0.721         0.0        0.9982      1
ARID1A    8289    1.238       0.836         0.001      0.721         0.0        0.9982      0
SMARCB1   6598    1.238       0.836         0.001      0.721         0.0        0.9982      0
SOX11     6664    1.238       0.836         0.001      0.721         0.0        0.9982      0
SMARCE1   6605    1.238       0.836         0.001      0.721         0.0        0.9982      0
SMARCA4   6597    1.238       0.836         0.001      0.721         0.0        0.9982      0
FIG4      9896    0.942       0.738         38.0       0.0           0.0        0.0         0
CYP26C1   340665  0.074       0.0           24.0       0.273         0.0        0.0         0
RFT1      91869   0.0207      0.0           35.0       0.0           0.0        0.0         0
VEGFC     7424    -0.110      0.0           34.0       0.0           0.0        0.0         0
```

