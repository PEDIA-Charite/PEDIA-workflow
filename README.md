# Project PEDIA

##Preprocessing update

## Current todo

* Transcript search for completion of HGVS codes
* Creation of VCF Files from hgvs
* Running the simulation
* Complete packaging for distribution

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

### Configuration

Most configuration options are in a `config.ini` file, with options commented.
A `config.ini.SAMPLE` in the project directory can be used as reference for
creating an own configuration.

Since some steps depend on the existence of API keys, running the preprocess.py
script without a configuration file will **not work**.
