# spectra-cluster-cli [![Build Status](https://travis-ci.org/spectra-cluster/spectra-cluster-cli.svg?branch=master)](https://travis-ci.org/spectra-cluster/spectra-cluster-cli)

This is a stand-alone implementation of the new updated 
[PRIDE Cluster](https://www.ebi.ac.uk/pride/cluster) algorithm. 
It is based on the [spectra-cluster](https://github.com/spectra-cluster/spectra-cluster) 
API and uses a highly similar logic as the Hadoop implementation 
[spectra-cluster-hadoop](https://github.com/spectra-cluster/spectra-cluster-hadoop) 
used to build the [PRIDE Cluster](https://www.ebi.ac.uk/pride/cluster) resource.

You can find an overview over all our clustering related tools at
https://spectra-cluster.github.io.

__WARNING__: This software is still in beta phase. We do expect it to still have several bugs. Bug reports are highly welcome! Just submit in [issue](https://github.com/spectra-cluster/spectra-cluster-cli/issues) to let us know.

## Changelog

### Version 1.1.2

* Fixed a bug that caused a crash when re-binning very small datasets
  with scarce m/z regions.

### Version 1.1.1

* Fixed a bug when running on systems with non-english locale
* Fixed a bug in the merging of adjacent bins: Spectra were not re-binned correctly
* Adapted binning procedure to create fewer files

### Version 1.1.0

* Updated to spectra-cluster API version 1.1 Therefore, the following new features are available
  * Automatically detect the number of comparison to use. This no longer has to be
    set by the user.
  * A new command line option was added (`-x_min_adapt_comparison`) that sets a minimum number of
    comparisons that will always be used (as a correction factor) in case the number learned from the
    data is smaller. This feature is only recommended when clustering very small datasets.
* New command line option `-add_scores` was added. This then outputs additional scores to the generated
  .clustering file. These additional properties are stored as JSON encoded strings.
* Added command line options to only cluster identified or only cluster unidentified spectra.
* For a complete list of command line options, please refer to the output of the `-help` option

### Version 1.0.4

* Updated to spectra-cluster API version 1.0.11. Changes include:
  * New version of the .clustering file format which stores how often a consensus spectrum peak 
    was observed.
  * Option to add debug information to spectra such as the score when the spectrum was merged to a cluster.
* Added support for **incremental clustering**. The spectra-cluster-cli tool can
  now process MGF files and .clustering files. This only works with .clustering files
  created with the spectra-cluster-cli tool version >= 1.0.4.

### Version 1.0.3

* Changed to spectra-cluster API version 1.0.10 which includes contaminant peak filters
* Added CLI options to remove contaminant peaks from spectra.

### Version 1.0.2

* Changed to spectra-cluster API version 1.0.9. This fixes issues caused
  by unknown charge states (as well as 
  [#2](https://github.com/spectra-cluster/spectra-cluster-cli/issues/2)).

### Version 1.0.1

* updated to new .clustering format version
  * .clustering files now includes complete reference to original spectrum using the same indexing system
    as the PSI standard file formats
* added feature to learn the cumulative distribution function (CDF) from a given dataset and then use
  this newly learned CDF
* added a list of experimental / advanced parameters (all starting with "-x_"). For
  more information simply launch the application with the "-help" parameter
* __Note:__ For small datasets (<100 MS runs), the option `-x_min_comparisons` should be set to 10,000.
* fixed [#3](https://github.com/spectra-cluster/spectra-cluster-cli/issues/3) see 1db57ae

## Installation
The spectra-cluster-cli application is written in Java and therefore runs on Windows, Linux, and Mac OS X. 

__Java__: needs to be installed on your system for the spectra-cluster-cli to work. You can download the latest Java version for your system [here](https://www.java.com).

To install the spectra-cluster-cli simply download the latest [release](https://github.com/spectra-cluster/spectra-cluster-cli/releases) and extract the zip file.

## Usage
The spectra-cluster-cli tool is a command-line only tool. If you prefer to use a graphical use interface, please use the [spectra-cluster-gui](https://github.com/spectra-cluster/spectra-cluster-gui) tool. A detailed tutorial on how to prepare your files for clustering can be found at https://spectra-cluster.github.io.

To use the spectra-cluster-cli tool, follow these steps: 

Open a command line and navigate to the folder where you extracted the spectra-cluster-cli to. The following example assumes that you are already in this folder.

This command launches the clustering job using the default values (as used for the [PRIDE Cluster](https://www.ebi.ac.uk/pride/cluster) resource) all available CPU cores and writes the results to my_clustering_result.clustering.

__Note__: You need to replace the spectra-cluster-cli-1.1.0.jar with the name of the downloaded version.

```bash
$ java -jar spectra-cluster-cli-1.1.0.jar -filter immonium_ions -output_path my_clustering_result.clustering C:\my_first_file.mgf C:\my_second_file.mgf
```

To improve the clustering accuracy in small datasets (< 100 MS runs) the default value for `-x_min_comparisons` was changed to 10,000 (changed in version 1.0.3). When clustering a repository scale dataset, a value of 5,000 is used (default value in the Hadoop version). 

Additionally, we recommend to always either use the `immonium_ions` filter, or even filter all peaks below 150 m/z ("-filter mz_150") or even below 200 m/z ("-filter mz_200").

The full list of options is printed through the -help parameter:

```bash
$ java -jar spectra-cluster-cli-1.1.0.jar -help
```


## Using the clustering results

The spectra-cluster-cli generates a .clustering file to store the clustering results.
A specification of this format can be found at the 
[clustering-file-reader page](http://github.com/spectra-cluster/clustering-file-reader)

The 
[spectra-cluster-py](https://github.com/spectra-cluster/spectra-cluster-py)
Python library contains a collection of tools to analyse the clustering results.
You can find a detailed documentation of this library at
http://spectra-cluster-py.readthedocs.io. Additionally, this library contains
many classes that should help in writing your own Python scripts to analyse your
clustering results.

Additionally, we provide a 
[Java API](http://github.com/spectra-cluster/clustering-file-reader)
that can be used to develop Java software that reads the .clustering file format.

## Getting help

### Technical issues

Should you have any problems when running the spectra-cluster-cli tool, please do
not hesitate to report this problem using the
[issue tracker](https://github.com/spectra-cluster/spectra-cluster-cli/issues).
 
### "Bioinformatic" issues

In case you have any other questions don't hesitate to post a question at 
[http://qa.proteomics-academy.org](http://qa.proteomics-academy.org).
 
These can include

* questions about what to do with the results
* questions about whether the clustering algorithm can be used for a given analysis problem

## The PRIDE Cluster resource and citation

We are using the [Hadoop version](https://github.com/spectra-cluster/spectra-cluster-hadoop) 
version of the spectra-cluster algorithm to cluster the complete 
[PRIDE Archive repository](http://www.ebi.ac.uk/pride). Thereby, we were able to
recognize millions of consistently unidentified spectra across thousands of submission.
These clustering results are presented in the
[PRIDE Cluster](http://www.ebi.ac.uk/pride/cluster) resource.

For more information see the recent paper
[Griss et al., Recognizing millions of consistently unidentified 
spectra across hundreds of shotgun proteomics datasets., 
Nat. Meth. 2016 Aug;13(8):651-656](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4968634/).

Additionally, if you are able to use our algorithm for your own project, please
cite the above reference.
