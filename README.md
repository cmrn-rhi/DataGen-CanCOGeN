# DataGen-CanCOGeN v1.7

## Table of contents
* [General Info](#General-Info)
* [Technologies](#Technologies)
* [How-To](#How-To)
* [Files](#Files)

## General Info
For generation of CanCOGeN metadata to test applications.

## Technologies
Project is created with:
* Python version: 3.8.3

## How-To
Run the following in cmd when files are in current directory:

`$ python main.py <file-name> <int for number of rows to generate> <'comma' or 'tab' delimiter>`

<'tab'> delimiter will save file as <.csv>

<'comma'> delimiter will save file as <.csv>

Defaults to 'comma' delimiter if no delimiter is specified.

*In-progress: How-to specify error requests for the error-grid output.*

## Files
* **main.py** - for generating delimited csv files.
* **library.py** - library of functions for use in main.py.
* **vocab-lists.csv** - lists CanCOGeN COVID19 picklist vocabularies and additional lists for imitation data generation.
* **variables.csv** - lists CanCOGeN COVID19 variable names.
