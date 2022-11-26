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
