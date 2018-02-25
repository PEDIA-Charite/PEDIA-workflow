# Project PEDIA

This folder contains the project data of the PEDIA study

## Description of the subfolders

* 1_qualityCheck
	All About the quality check. Scripts, pipeline, data etc.
* 2_phenomization
	Scripts about running the phenomizer/boqa/... priorization using HPO terms
* 3_simulation
	Scripts and pipelie about generation background populations and simulationg the data.
	In addition the spike in is made and molecular pathogenicity scores are added to the jsons.
* data
	On the one side databases (ans their download scripts) needed for several tasks in the whole project. On the other side the read PEDIA data with samples etc.


## General Information

The whole workflow of the PEDIA project uses [snakemake](https://snakemake.readthedocs.io/) to run a pipeline together with [conda/bioconda](https://bioconda.github.io/) to install the necessary programs. So pelase get familiar with both BEFORE starting the workflow. A good start is the [snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html).

Anyway here are the major thinsg you have to do to run the workflow.

### Install miniconda

Go have a look at the [miniconda website](https://conda.io/miniconda.html). Be sure that you choose the right version depending on your python version. To find out what python version you have please type in

```
 python --version
```

### Install necessary software via conda

For example we want to run the quality check pipeline. This will require python version 2.7. But snakemake needs at least 3.5. So in general this is problematic. BUT: snakemake can handel this. So we need first an enviroment with snakemake. Let's go into the quality folder.

```
cd 1_qualityCheck
```

Now we will generate a shell environment with all necessary programs for the quality check. The software needed for the quality check is in the `environment.yml` file. Conda can read it:

```
conda env create -f environment.yml
```

Now we created an enviroment called `pedia_quality`. We can activate it and we should have snakemake installed.

```
source activate pedia_quality
snakemake -h
```
We can deactivate the environment using `source deactivate`. The command `conda env list` will list you all environments.

### Use snakemake

Now lets run the quality workflow. We can make a "dry run" using `snakemake -n` to see what snakemake would probably do.

Because of this special case (python 2.7 and 3.5) we have to tell snakemake to use conda.

```
 snakemake --use-conda
```

Using this command conda will first create an enviroment with python 2.7 (defined in the file `envs/python2.7.yml`)
