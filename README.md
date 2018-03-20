# Project PEDIA

##Preprocessing update

The preprocessing steps for quality control of json files have been rewritten.

## Current todo

* Transcript search for completion of HGVS codes
* Creation of VCF Files from hgvs
* Running the simulation
* Complete packaging for distribution

## Installation instructions

Python 3.6 and pip. Using virtualenvs is recommended. Python includes a builtin venv module for this. Install requirements.txt to virtualenv and all steps up to phenomization should work.

* Create a python virtualenv preferably using the venv module
```
python3 -m venv <dir of choice, eg env>
./<env>/bin/activate (command can vary depending on used shell)
```
* Install dependencies

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

## Usage instructions

Files inside `scripts/` can be used for inspiration for own usecases. All
scripts should be run from the project base directory to automatically include
the lib package containing the actual program code.

Keep in mind that the virtual environment needs to be enabled for script
execution.
