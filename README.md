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
ProteomeRiver consists of several command line scripts. 
### Clean phosphorylation data  
```
clean_phos_cmd.R \
 --output-dir=$PHOS_DIR \
 --fasta=$DATA_DIR/ALPK1-set1MusMusculus20201226CanIso.fasta \
 --fasta-save=$PHOS_DIR/aa_seq_tbl.RDS  \
 --raw-counts=$DATA_DIR/ALPK1-set1evidence.txt \
 --site-prob=0.75 \
 --recover-prob=0.5 \
 --column-pattern="Reporter intensity corrected" \
 --add-columns="experiment"
```

## Annotations Download Procedures (keep here for now)

## Download UniProt data table
wget -O data.xml.gz "https://rest.uniprot.org/uniprotkb/stream?compressed=true&download=true&format=xml&query=%28%28proteome%3AUP000005640%29%29"
# stream?compressed=true&download=true&format=xml&query=((proteome:UP000005640))
gunzip ata.xml.gz
mv data data.xml
tail -n +2 data.xml > data_updated.xml
# cd /home/ignatius/PostDoc/2022/Embryology_BMP_14/Source/UniProt
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
