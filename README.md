# Project PEDIA

## Description of the subfolders

* **3_simulation** - 
	Scripts and pipelie about generation background populations and simulationg the data.
	In addition the spike in is made and molecular pathogenicity scores are added to the jsons.
* **data** - 
	It contains the files from external database (dbSNP, reference genome, CADD, ExAC). Besides, data/PEDIA/json/phenomized contains the JSON files after QC and phenomization. 
* **classifier** - 
 It is a submodule of PEDIA-workflow. The repository is here. (https://github.com/PEDIA-Charite/classifier). 
 * rest of the files and folders belong to preprocessing such as lib, helper, test
	All About the quality check and phenomization
	
## General Information

The whole workflow of the PEDIA project uses [snakemake](https://snakemake.readthedocs.io/) to run a pipeline together with [conda/bioconda](https://bioconda.github.io/) to install the necessary programs. So pelase get familiar with both BEFORE starting the workflow. A good start is the [snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html). 

### Install miniconda

Go have a look at the [miniconda website](https://conda.io/miniconda.html). Be sure that you choose the right version depending on your python version. To find out what python version you have please type in

```
 python --version
```
### Install necessary software via conda

Let's go into the data folder to download external files.

```
cd data
```

Now we will generate a shell environment with all necessary programs for downloading and process files. The software needed for the downloading is in the `environment.yaml` file. Conda can read it:

```
conda env create -f environment.yaml
```

Now we created an enviroment called `pedia_download`. We can activate it and we should have snakemake installed.

```
source activate pedia_download
snakemake -h
```
We can deactivate the environment using `source deactivate`. The command `conda env list` will list you all environments.

Now lets run the pedia download workflow. We can make a "dry run" using `snakemake -n` to see what snakemake would probably do.

## Setup for preprocessing

* Create a python virtualenv preferably using the venv module
```
python3 -m venv <dir of choice, eg env>
./<env>/bin/activate (command can vary depending on used shell)
```
* Install dependencies
  * Python 3.6 is currently the only tested version

```
pip install -r requirements.txt
```
* Optional: Run tests
    * Populate the tests/data directory. The case with the specified case ID and the specified genomic entry is needed.
        ```
        tests/data/cases/51702.json
        tests/data/genomics_entries/2669.json
        tests/data/config.ini -- necessary for API keys
        ```
    * Run tests
        ```
        python3 -m unittest discover
        ```

* Obtain additional data files - put them into the specified locations
  * `data/mim2gene.txt` - publicly available on omim.org
    * [https://omim.org/downloads/][OMIM download page]
  * `data/morbidmap.txt` - requires API-key access to download
  * note that if new versions of the files are obtained MD5 checksums might have
    to be regenerated

## Usage instructions

Files inside `scripts/` can be used for inspiration for own usecases. All
scripts should be run from the project base directory to automatically include
the lib package containing the actual program code.

### Configuration

Most configuration options are in a `config.ini` file, with options commented.
A `config.ini.SAMPLE` in the project directory can be used as reference for
creating an own configuration.

### HGVS Error dict

HGVS variant overrides are specified in `hgvs_errors.json`. Which is per-default
searched for in the project root.

The hgvs version is specified in `lib/constants.py` and will cause an error if
an hgvs errors file of not at least the specified version is found.

The number can be lowered manually to accept older hgvs error files.

A version of 0 will accept no hgvs_errors file.

### Required external files
   * Go to data folder, and run 'snakemake all' to download all necessary files such as reference genome, population data.
   * Copy dbsnp files to data/dbSNP/b147
   * Copy IRAN_trio files to 3_simulation/background/data/IRAN_trio/
   * Copy mim_to_ps.json, mim2gene.txt, mimTitles.txt, morbidmap.txt, omim_deprecated_replacement.json, phenoptypicSeries.txt to data/
   * Copy corrected JSON files to process/correct/
   * Copy hgvs_errors.json to project folder
   * Copy config.ini to project folder
   * Copy genemap2.txt to 3_simulation/OMIM/ folder

### Running preprocessing

Since some steps depend on the existence of API keys, running the preprocess.py
script without a configuration file will **not work**.

Keep in mind that the virtual environment needs to be enabled for script
execution.

The **preprocess.py** script contains most information necessary for running a
conversion of new json files into the old format necessary for conversion.

```
# do not forget to activate the previously created virtual environment

# get a list of usable options
./preprocess.py -h

# run complete process with AWS synchronization
./preprocess.py

# run for a single file (specifying output folder is beneficial)
./preprocess.py -s PATH_TO_FILE -o OUTPUT_FOLDER
```

### Run PEDIA pipeline
There are three steps to run pipeline.
1. Environment setup
   * Go to data folder, and run 'snakemake all' to download all necessary files such as reference genome, population data.
   ```
   source activate pedia_download
   snakemake all
   ```
   * Copy dbsnp files to data/dbSNP/b147
   * Copy IRAN_trio files to 3_simulation/background/data/IRAN_trio/
   * Copy mim_to_ps.json, mim2gene.txt, mimTitles.txt, morbidmap.txt, omim_deprecated_replacement.json, phenoptypicSeries.txt to data/
   * Copy corrected JSON files to process/correct/
   * Copy hgvs_errors.json to project folder
   * Copy config.ini to project folder
   * Copy genemap2.txt to 3_simulation/OMIM/ folder

1. Download cases and perform preprocessing

   ```
   python3 preprocess.py
   ```
   * config.yml contains the cases passed quality check. SIMPLE_SAMPLES is the case with disease-causing mutation but without real VCF file. VCF_SAMPLES is the case with real VCF file. TEST_SAMPLE is the case with real VCF but without disease-causing mutation.
   * process/aws_dir is the folder of cases downloaded via aws.
   * data/PEDIA/jsons/phenomized is the folder which contains the JSON files passed QC.
   * data/PEDIA/mutations/variants.vcf.gz  is the VCF file which contains disease-causing mutations of all cases.
   * data/PEDIA/vcfs/original is the folder which contains the VCF files. In mapping.py, we rename the filename of VCF files to case_id.vcf.gz and store to ../data/PEDIA/vcfs/original/. The new filename is added in vcf field of the JSON file. For example,
   ```
   "vcf": [
           "28827.vcf.gz"
       ],
   ```

1. Get JSON files of simulated cases and real cases

    To obtain the CADD scores of variants, we need to annotate the VCF files and retrieve the CADD score and append it to the geneList in JSON file. Now, we go to 3_simulation folder and activate simulation environment. 
   
    Note: you could skip this step by running the experiemnt in classifier. The classifier will trigger this subworkflow to generate JSON files.
   
    Before we start, we would like to explain the two experiments we want to conduct in this study. First one is that we want to perform cross-validation on all cases to evaluate the performance among three simulation samples (1KG, ExAC and IRAN). The second one is that we want to train the model with simulated cases and test on the real cases. To achieve these two goals, we have the following command to perform simulation and generate the final JSON files.


    * To perform the CV experiment, we run the following command to obtain the JSON files simulated from 1KG, ExAC and IRAN data. You could replace 1KG with ExAC and IRAN
    ```
    snakemake performanceEvaluation/data/CV/1KG.csv
    ```
   
    * To peform the second experiemnt, we run the following command to obtain the training and testing data sets. Generate the JSON files of **real cases** the output will be in 3_simulation/json_simulation/real/test
    ```
    snakemake createCsvRealTest
    ```
   
    * Generate the JSON files of **simulated cases** the output will be in 3_simulation/json_simulation/real/train/1KG. You could replace 1KG with ExAC and IRAN
    ```
    snakemake performanceEvaluation/data/Real/train_1KG.csv
    ```
    
    * The final JSON files are in 3_simulation/json_simulation folder.
        * 3_simulation/json_simulation/1KG is the folder for all cases simulated by 1KG.
        * 3_simulation/json_simulation/ExAC is the folder for all cases simulated by ExAC.
        * 3_simulation/json_simulation/IRAN is the folder for all cases simulated by IRAN.
        * 3_simulation/json_simulation/real/train is the folder for the cases without simulated by 1KG, ExAC or IRAN. We also have three folder under this folder.
        * 3_simulation/json_simulation/real/test is the folder for the cases with real VCF file.
      
1. Cross-validation evaluation
   * Go to classifier folder.  Run 'source activate classifier' to enable the environment. If you haven't created the environment, please execute 'conda env create -f environment.ymal'.
   * Perform 10 times 10 fold cross-validation on all 3 simulation population by running 'snakemake CV_all'. You can add --cores 3 to run it on parallel. The output will be in output/cv/CV_1KG, output/cv/CV_ExAC and output/cv/CV_IRAN. The classifier will trigger the simulation and phenomization if the files haven't been generated. It takes a long time for running the first time due to the process of simulation population data.

   ```
   snakemake -p --cores 3 CV_all
   ```
   * If you only want to run on 1KG simulation data, please execute this command 'snakemake ../output/cv/CV_1KG/run.log'. Please remind the working directory of classifier is in scripts folder, so to run on 1KG simulation you need to specify the output file '../output/cv/CV_1KG/run.log' instead of 'output/cv/CV_1KG/run.log '.

   ```
   snakemake ../output/cv/CV_1KG/run.log
   ```
   
1. Train and test evaluation
   * Training set is in 3_simulation/json_simulation/real/train/1KG (1KG, ExAC and IRAN)
   * Testing set is in 3_simulation/json_simulation/real/test
   ```
   snakemake ../output/real_test/1KG/run.log
   ```

1. Train with all cases and test on patient with unknown diagnosis
   ```
   snakemake ../output/test/1KG/21147/21147.csv
   ```

1. How to read the PEDIA results?
   * case_id.csv contains all genes with corresponding pedia scores in this case
   * manhattan_case_id.png is the manhattan plot of the case
   * manhattan_all.png is the manhattan plot of all cases
   * rank_gene_1KG.csv contains the case_id and the predicted rank of disease-causing gene of each case.
   * rank_1KG is the performance overview file. You will find the number of cases in each rank.
