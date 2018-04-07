# Project PEDIA

##Preprocessing update

## Current todo

* Transcript search for completion of HGVS codes
* Creation of VCF Files from hgvs
* Running the simulation
* Complete packaging for distribution

## Setup

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

### Running

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

