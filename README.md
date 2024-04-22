# RNA2way_pred

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

## Install

To install RNA2way_pred 

```shell
git clone https://github.com/Sanduni522/RNA2way_pred
cd RNA2way_pred
python -m pip install .
```


## How to run 

```shell
RNA2way_pred --help
Usage: RNA2way_pred [OPTIONS]

  A program that predicts DMS reactivity using structural parameters

Options
  --help            Show this message and exit.
```

```shell
RNA2way_pred
```

## Important information

1) The content in Important_files should be placed where the PDB files are to run this. If you have separate DMS reactivity values stored in a CSV file, replace 'Reactivity_values.csv' with it. 

2) The package must be executed in a virtual environment. You can read below to learn how to make a virtual environment.

  Create a virtual environment.
```shell
$ brew install virtualenv
```
  Create a directory (./pythonenv) to hold it.
```shell
$ virtualenv --system-site-packages -p python3 ./pythonenv
$ cd ./pythonenv
```
  Activate the virtual environment.
```shell
$ source bin/activate
```

3) You must have the Rosetta rna_denovo package installed on your computer.

4) The package must be executed where the PDB files are located.


