# CNCA Standalone
 
MANUAL


CNCA aligns multiple small, closely related annotated genomes (<50kb, typically virus) using information from protein alignments.


This version represents a quick adaptation of the CNCA webtool to a standalone version. 
For instance, it only offers the run with default parameters. To explore all the possibilities of CNCA, please visit https://cnca.ijm.fr/

Currently, this standalone version is complicated to execute and is only intended for sharing the code with individuals who are interested in examining the algorithm's details. 
The actual tool corresponds to the web version.

## Preliminary steps

1. Install Python and the required libraries.

2. Ensure access to the MAFFT software in the environment.

3. Install R and the necessary packages.

4. Set up the required dependencies and configurations.

Furthermore, running this pipeline requires the installation of the following Python libraries:

In addition, access to MAFFT within the environment is required.

And finally, R along with its various packages is also needed.

INPUT
The submitted input must consist of genbank files (https://www.ncbi.nlm.nih.gov/genbank/). These genbank files must present closely related species (ANI >= 90%, to test your data: https://jspecies.ribohost.com/jspeciesws/). This proximity constraint is particularly important, as the organization of the genomes and in particular the order of the coding sequences must be comparable for the use of this tool to make sense (see Protein alignment part).

PROCESS
Loading and parsing the user data
From these genbank files, the CDS of each genome are retrieved and organized in temporary files allowing direct access to the full genome for the nucleotide alignment and the protein sequence for the protein alignment. 

Nucleotide alignment
The nucleotide alignment process is performed with Mafft. 

Protein alignment
For protein alignment, CDS are recovered in the order of appearance of the genomes. For each genome a temporary file with the 


# tmp dir creation
the tmp directory must be owned by "www-data".
To do it, create the directory as root with mkdir and change the owner with:
sudo chown www-data tmp