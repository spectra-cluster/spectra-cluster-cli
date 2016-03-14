# spectra-cluster-cli
This is a stand-alone implementation of the [PRIDE Cluster](https://www.ebi.ac.uk/pride/cluster) algorithm. It is based on the [spectra-cluster](https://github.com/spectra-cluster/spectra-cluster) API and uses a highly similar logic as the Hadoop implementation [spectra-cluster-hadoop](https://github.com/spectra-cluster/spectra-cluster-hadoop).

__WARNING__: This software is still in beta phase. We do expect it to still have several bugs. Bug reports are highly welcome! Just submit in [issue](https://github.com/spectra-cluster/spectra-cluster-cli/issues) to let us know.

## Installation
The spectra-cluster-cli application is written in Java and therefore runs on Windows, Linux, and Mac OS X. 

__Java__: needs to be installed on your system for the spectra-cluster-cli to work. You can download the latest Java version for your system [here](https://www.java.com).

To install the spectra-cluster-cli simply download the latest [release](https://github.com/spectra-cluster/spectra-cluster-cli/releases) and extract the zip file.

## Usage
The spectra-cluster-cli currently only provides a command line interface (CLI). 

Open a command line and navigate to the folder where you extracted the spectra-cluster-cli to. The following example assumes that you are already in this folder.

This command launches the clustering job using the default values (as used for the [PRIDE Cluster](https://www.ebi.ac.uk/pride/cluster) resource) all available CPU cores and writes the results to my_clustering_result.clustering.

__Note__: You need to replace the spectra-cluster-cli-1.0-SNAPSHOT.jar with the name of the downloaded version.

```bash
$ java -jar spectra-cluster-cli-1.0-SNAPSHOT.jar -output_path my_clustering_result.clustering C:\my_first_file.mgf C:\my_second_file.mgf
```

The full list of options is printed through the -help parameter:

```bash
$ java -jar spectra-cluster-cli-1.0-SNAPSHOT.jar -help
```

