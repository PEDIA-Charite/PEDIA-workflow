# Project PEDIA

## UPDATE: Preprocessing update

The preprocessing steps for quality control of json files have been rewritten.

The new preprocessing includes multiple entrypoints for file overrides and fixing of faulty json files.

## UPDATE: Current todo

* Transcript search for completion of HGVS codes
* Creation of VCF Files from hgvs
* Running the simulation

## UPDATE: Goals

* Single Files should always be runnable
* Externalized configuration in central config.ini file
* Future adaptations to Face2Gene format changes should be easier
* More complicated mapping/per disease/per mutation is completely supported and more transparent in the preprocessing procedure

* Complete packaging for distribution

## Installation instructions

Python 3.6 and pip. Using virtualenvs is recommended. Python includes a builtin venv module for this. Install requirements.txt to virtualenv and all steps up to phenomization should work.

Config.ini should be adapted from config.ini.SAMPLE, some API Keys might be needed to access OMIM, Phenomizer and Face2Gene data.
