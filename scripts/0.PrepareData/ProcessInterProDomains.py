"""
Process InterPro domains
"""
docker run -it --gpus all -v $PWD:$PWD -w $PWD gnn_python:1.0

import json
import sys
import os
from urllib.error import HTTPError
from urllib.request import urlopen
import pandas as pd
import re
import gzip
from tqdm import tqdm
import xml.etree.ElementTree as ET

# Load proteins in mouse
## Define patterns
id_pattern = re.compile(r'^AC\s+(\w+)(?:; (\w+);)?\s+')
os_pattern = re.compile(r'^OS\s+(.+)')

## Initialize lists to store ID and OS data
uniprot_ids = []
organism_list = []
uniprot_id = ""
## Open the UniProt DAT file for reading
with gzip.open('data/Uniprot/UP000000589_10090.dat.gz', 'rt') as file:
    # Iterate over each line in the file
    for line in file:
        # Match the ID pattern
        id_match = id_pattern.match(line)
        if id_match:
            uniprot_id = id_match.group(1)
        # Match the OS pattern
        os_match = os_pattern.match(line)
        if os_match:
            organism_species_text = os_match.group(1)
            # OS field might be spread across multiple lines
            # In such cases, continue reading lines until the full OS field is obtained
            while line.strip().endswith(';'):
                next_line = next(file).strip()
                organism_species_text += " " + next_line[:-1]
            uniprot_ids.append(uniprot_id.split("; "))
            organism_list.append(organism_species_text)

## Check all proteins are from mouse
set(organism_list)

## Get unique list of uniprot ids
flattened_uniprotids = [item for sublist in uniprot_ids for item in sublist]
unique_uniprotids = set(flattened_uniprotids)

## Load mapping between MGI and uniprot ids
gene_map = pd.read_csv("results/preprocess/MGI_uniprot_map.tsv", sep='\t')

gene_map_filt = gene_map[gene_map['UniProtIDs'].isin(unique_uniprotids)]

def process_xml_by_blocks(filename):
    protein_data = {}
    with gzip.open(filename, 'rb') as f:
        xml_iter = ET.iterparse(f, events=('start', 'end'))
        ## Initialize protein_id
        protein_id = ""
        for event, element in xml_iter:
            if event == 'start' and element.tag == 'protein':
                if element.get('id') in unique_uniprotids:
                    protein_id = element.get('id')
                    print(protein_id)
                    prot_len = int(element.get('length'))
                    domain_names = []
                    interpro_ids = []
                    domain_props = []
                    domain_Ns = []
                    domain_dbs = []
            elif protein_id != "":
                if event == 'start' and element.tag == 'match':
                    domain_id = element.get("id")
                    domain_db = element.get("dbname")
                    domain_size = 0
                    domain_N = 0
                    domain_name = ""
                elif event == 'start' and element.tag == 'ipr':
                    if element.get("type") == "Domain" and domain_id != "":
                        domain_name = element.get('id')
                elif event == 'start' and element.tag == 'lcn' and domain_name != "":
                    domain_size += (int(element.get('end')) - int(element.get('start')))
                    domain_N += 1
                elif event == 'end' and element.tag == 'match':
                    if domain_name != "":
                        domain_names.append(domain_id)
                        domain_dbs.append(domain_db)
                        interpro_ids.append(domain_name)
                        domain_props.append(domain_size/prot_len)
                        domain_Ns.append(domain_N)
                    domain_id = ""
                    domain_name = "" 
            if event == 'end' and element.tag == 'protein':
                if protein_id != "":
                    out = pd.DataFrame({"UniprotID" : protein_id, 
                        "Domain_ID" : interpro_ids,
                        "Domain_DB" : domain_dbs,
                        "Domain_Name" : domain_names,
                        "Prop" : domain_props,
                        "N": domain_Ns
                    })
                    protein_id = ""
                    yield out
                element.clear()


prot_domains = [df for df in process_xml_by_blocks("data/Interpro/match_complete.xml.gz")]
prot_domains_df = pd.concat(prot_domains)

prot_domains_comb = pd.merge(gene_map_filt, prot_domains_df, 
                             left_on = "UniProtIDs", right_on = "UniprotID",
                             how = "inner")
prot_domains_comb = prot_domains_comb.drop("UniProtIDs", axis = 1)
prot_domains_comb.to_csv("results/preprocess/UniprotAllDomains.tsv",
                         sep = "\t", index=False
                         )


def getUniqueMarkerIDs(df, db):
    df_filt = df[df['Domain_DB'] == db]
    n_markers = len(df_filt.MGIMarkerAccessionID.unique())
    n_domains = len(df_filt.Domain_ID.unique())
    return n_markers, n_domains

db_names = ["PFAM", "PROFILE", "SMART", "CDD", "PROSITE"]
for db in db_names:
    print(getUniqueMarkerIDs(prot_domains_comb, db))

## Select only PFAM domains
pfam_domains_comb = prot_domains_comb[prot_domains_comb['Domain_DB'] == "PFAM"]
pfam_domains_comb.to_csv("results/preprocess/UniprotPFAMDomains.tsv",
                         sep = "\t", index=False
                         )

dom_counts = pfam_domains_comb.Domain_ID.value_counts()



def getUniqueMarkerIDs2(df, n):
    dom_counts = df.Domain_ID.value_counts()
    df_filt = df[df['Domain_ID'].isin(dom_counts[dom_counts > n].index)]
    n_markers = len(df_filt.MGIMarkerAccessionID.unique())
    n_domains = len(df_filt.Domain_ID.unique())
    return n_markers, n_domains

dom_freq = [0, 1, 5, 10, 20, 40, 60]
for n in dom_freq:
    print(n, getUniqueMarkerIDs2(prot_domains_comb, n))

## Convert format form long to wide (domains in columns)
##df.pivot(index = "Protein", columns = "Domain_ID", values = "Prop")