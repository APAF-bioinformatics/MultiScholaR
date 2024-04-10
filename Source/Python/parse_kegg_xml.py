#!/usr/bin/env python3

library(reticulate)


```{python}
from lxml import etree, objectify
import pandas as pd
import sys

xmlFile =  "/home/ignatius/PostDoc/2021/Podocyte_2021/Data/KEGG/hsa00600.xml"  #  sys.argv[1]  # ".xml"
output_file = "/home/ignatius/Desktop/hsa00600_reactions.txt" # sys.argv[2] # '.tab'


# print(f"Path of input XML file      : {sys.argv[1]}")
# print(f"Path of output file     : {sys.argv[2]}")


with open(xmlFile) as fobj:
	xml = fobj.read()

## need to remove UTF declaration string from the file first
root = etree.fromstring(xml, parser=etree.XMLParser(huge_tree=True)) 
```

```{python}
namespace_x = {'x':"https://www.kegg.jp/kegg/xml/docs/KGML_v0.7.2_.dtd.html" }

entries = root.xpath(".//x:entry", namespaces=namespace_x )


def parseOneGoTerm(one_go_term):
	namespace_x = {'x':"http://uniprot.org/uniprot" }

	has_go_id = one_go_term.xpath(".//@id", namespaces=namespace_x)

	if len(has_go_id ) == 0:
		return( "NA" + ";" + "NA" + ";" + "NA" + 
				";" + "NA" +  ";" + "NA" )

	go_id = has_go_id[0]

	go_term_combined = one_go_term.xpath(".//x:property[@type='term']", 
	namespaces=namespace_x )[0].xpath(".//@value", 
	namespaces=namespace_x )[0]

	go_type = go_term_combined.split(":")[0]

	go_term = go_term_combined.split(":")[1]

	go_evidence = one_go_term.xpath(".//x:property[@type='evidence']", 
	namespaces=namespace_x )[0].xpath(".//@value", 
	namespaces=namespace_x )[0]

	go_project = one_go_term.xpath(".//x:property[@type='project']", 
	namespaces=namespace_x  )[0].xpath(".//@value", 
	namespaces=namespace_x  )[0]

	return( go_id + ";" + go_type + ";" + go_term + 
	";" + go_evidence +  ";" + go_project)



def getGOTermFromOneUniprotEntry(one_entry):
	namespace_x = {'x':"http://uniprot.org/uniprot" }

	uniprot_acc = one_entry.xpath(".//x:accession", 
		namespaces=namespace_x )[0].text
	
	list_of_go_terms = one_entry.xpath(".//x:dbReference[@type='GO']", 
		namespaces=namespace_x )
	
	if len(list_of_go_terms) == 0: 
		return( [ uniprot_acc + ";NA;NA;NA;NA;NA"] )
		
	parsed_go_terms = list( map(parseOneGoTerm, list_of_go_terms) )
	
	accession_and_go_terms = list( map(lambda x: uniprot_acc + ";" + x, parsed_go_terms) )
	
	return( accession_and_go_terms)





## Go through one example for parsing 
#one_go_term = list_of_go_terms[0]

#parseOneGoTerm(one_go_term)

## Go through one entry for parsing
#one_entry = entries[0]

# getGOTermFromOneUniprotEntry( one_entry)

## Run all entries


all_go_terms = list( map( getGOTermFromOneUniprotEntry, entries))

flat_list = [item for sublist in all_go_terms for item in sublist]


df = pd.DataFrame(flat_list, columns=['go_term_string'])

df[['uniprot_acc',
'go_id', 
'go_type', 
'go_term', 
'go_evidence', 
'go_project']] = df.go_term_string.str.split(";",expand=True)

df.drop('go_term_string', axis=1, inplace=True)

df.to_csv(output_file, 
	sep="\t",
	index=False)
```

