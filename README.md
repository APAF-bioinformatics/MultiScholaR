# ProteomeRiver

Tools for Proteomics and Phosphoproteomics analysis. 

## Getting Started

git clone repository 
```
git clone git@bitbucket.org:cmri-bioinformatics/proteomeriver.git
```

### Prerequisites

R version 3.6.3 (>=). 

## Installation

Edit `config.mk` file. Change `REXEC` and `PREFIX` as necessary. 

Install the R package and R scripts using make:
```
make install 
```
Add the scripts `PATH` as indicated in the command line. See all make options:
``` 
make help
```

## Usage
ProteomeRiver consists of several command line scripts. You can create a config.ini script with all the parameters and run the script sequentially:

```
clean_proteins_cmd.R -c config.ini
de_analysis_cmd.R -c config.ini
annot_proteins_cmd.R -c config.ini
```

Typical directory structure of your project should look like:

``` 

# General results folder 
cache
Reactome
UniProt

# Proteomics results folder 
clean_proteins
de_proteins
annot_proteins
prot_publication_graphs
proteins_pathways_enricher

# Phosphoproteomics results folder 
clean_phos
de_phos
annot_phos
phos_publication_graphs
norm_phos_by_prot_abundance
phos_pathways_enricher
phos_kinswingr_ST

```

## Annotations Download Procedures (keep here for now)

## Download UniProt data table
wget -O data.xml.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&download=true&format=xml&query=%28%28proteome%3AUP000005640%29%29"
# stream?compressed=true&download=true&format=xml&query=((proteome:UP000005640))
gunzip ata.xml.gz
mv data data.xml
tail -n +2 data.xml > data_updated.xml

## Edit the parameters in python file
python /home/ignatius/PostDoc/2021/proteomeriver/Source/Python/parse_go_terms.py "/home/ignatius/PostDoc/2022/pml_apex_tony_cesare_pml_20220822/Data/UniProt/data_updated.xml" \
 '/home/ignatius/PostDoc/2022/pml_apex_tony_cesare_pml_20220822/Results/UniProt/go_terms_table_python_all.tab'
 
wget -O data.tab.gz  "https://rest.uniprot.org/uniprotkb/stream?compressed=true&download=true&fields=accession%2Cid%2Corganism_id%2Creviewed%2Cannotation_score%2Cprotein_name%2Cgene_names%2Clength%2Ckeyword&format=tsv&query=%28%28proteome%3AUP000005640%29%29"
gunzip data.tab.gz
mv data data.tab

## Download Reactome
/home/ignatius/PostDoc/2021/proteomeriver/Source/R/download_reactome_table.R

## Download KEGG 
 /home/ignatius/PostDoc/2021/proteomeriver/Source/Rmd/retrive_kegg_annotation.Rmd


## Installing Anaconda
Installation instructions based on: 
https://phoenixnap.com/kb/how-to-install-anaconda-ubuntu-18-04-or-20-04

wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh
sha256sum Anaconda3-2022.10-Linux-x86_64.sh 
sudo apt install libgl1-mesa-glx libegl1-mesa libxrandr2 libxrandr2 libxss1 libxcursor1 libxcomposite1 libasound2 libxi6 libxtst6

bash Anaconda3-2022.10-Linux-x86_64.sh 
 
source ~/.bashrc 
conda info
conda update conda
conda update anaconda
conda create --name proteomeriver python=3
conda activate proteomeriver
conda deactivate

## History

* First release 01/10/2021.

## Credits

Authors: 

Ignatius Pang

Pablo Galaviz 

Contact:  cmri-bioinformatics@cmri.org.au


**Children’s Medical Research Institute, finding cures for childhood genetic diseases**  

## License

Edit license statement or refer to file. 

ProteomeRiver is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
any later version.

ProteomeRiver is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with ProteomeRiver.  If not, see <http://www.gnu.org/licenses/>.
