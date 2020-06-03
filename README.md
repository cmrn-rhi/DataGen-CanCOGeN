# DataGen-CanCOGeN v1.1

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

`$ python main.py <file-name.csv> <int for number of rows to generate> <'comma' or 'tab' deliminter>`

## Files
* main.py - for generating delimited csv files.
* library.py - library of functions for use in main.py.
* vocab-lists.csv - lists CanCOGeN COVID19 picklist vocabularies and additional lists for imitation data generation (yellow header).
* variables.csv - lists CanCOGeN COVID19 variable names.
