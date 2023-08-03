# CNCA Standalone
 
MANUAL


CNCA aligns multiple small, closely related annotated genomes (<50kb, typically virus) using information from protein alignments.


This version represents a quick adaptation of the CNCA webtool to a standalone version (on linux). 
For instance, it only offers the run with default parameters. To explore all the possibilities of CNCA, please visit https://cnca.ijm.fr/

Currently, this standalone version is complicated to execute and is only intended for sharing the code with individuals who are interested in examining the algorithm's details. 
The actual tool corresponds to the web version.

## Minimal Requirement

- Python 3.9.12

Python Package:

|Package                |Version |                                                                                      
|---------------------- |-------- |                                                                                    
|argparse               |1.1|
|biopython              |1.80|  
|json                   |2.0.9|                                                                                                                                                                                 
|numpy                  |1.23.5|                                                                                        
|pandas                 |2.03|
|pathlib                |*    |                                                                                                                                                                
|python-dateutil        |2.8.2 |                                                                                        
|python-dotenv          |1.0.0|
|subprocess	       |*|
|shutil                | *|

	
                                                                                                                                                                          
- R 4.3.1

R packages should install autonomously, if there is any troubleshooting during the installation of R packages, here is a list of the required packages along with their version:

|Package                |Version |                                                                                      
|---------------------- |-------- |                                                                            
|seqinr		       |4.2-8      |                                                                                                                                                                           
|optparse               |1.7.1|
|stringr                |1.5.0|
|tibble                 |3.1.5|
|ape                    |5.6-2|
|msa                    |1.24.0|


- MAFFT within the environment is required (https://mafft.cbrc.jp/alignment/software/linux.html).

## Run command
To run the pipeline

``` shell
python3 cnca_standalone_run.py -i [path_to_gbk_files] -w [output directory]
```

The full documentation is avalaible here : https://cnca.ijm.fr/CNCA_DOC.pdf
