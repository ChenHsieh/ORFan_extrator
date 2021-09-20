# orthoFinder parser
# generate the orphan gene from different level

# set up file path
species_file_name_table = 'file_list.csv'
orthogroup_path = 'Orthogroups.tsv'
peptide_file_dir = 'fasta_files'
# set up focal groups in list
focal_list = [
"PdeltoidesWV94_445_v2.1.protein_primaryTranscriptOnly.fa",
"PtremulaxPopulusalbaaltv2.1.primaryTrs.pep.fa",
"PtremulaxPopulusalbav2.1p.primaryTrs.pep.fa",
"Ptrichocarpa_444_v3.1.protein.fa",
"PtrichocarpaStettler14_532_v1.1.protein.fa",
"Ptrichocarpav4.1g.primaryTrs.pep.fa"
]

# import packages and read in files
from Bio import SearchIO
from Bio import SeqIO
from pathlib import Path
import argparse
import json
import pandas as pd 
import sys
import seaborn as sns
import matplotlib.pyplot as plt
from tqdm import tqdm
import itertools
from Bio.SeqFeature import SeqFeature, FeatureLocation
from multiprocessing import Pool
import time
# helper functions
def s2l(species_id):
    return species_csv.set_index('species_id').loc[species_id]['file_name.fa'][:-3]

# read files
orthogroup = pd.read_csv(orthogroup_path, delimiter='\t')
species_csv = pd.read_csv(species_file_name_table)
# read fasta files

peptide_file_dir = Path.cwd() / Path(peptide_file_dir)

# read fasta files
peptide_file_dir = Path.cwd() / Path(peptide_file_dir)
species_pep_dict = {}

print("reading fasta files")
start = time.time()
for index, row in tqdm(species_csv.iterrows()):
    species = row['species_id']
    pep_seq_file = peptide_file_dir / row['file_name.fa']
    seq_dict = SeqIO.to_dict(SeqIO.parse(pep_seq_file, "fasta"))
    species_pep_dict[species] = seq_dict
print(f"fasta reading: {time.time()-start} used")

orthogroup = orthogroup.set_index("Orthogroup")

## start parsing
orphan_dict = {}
subset_dict = {}
subset_sn = 0
count = 0
print("start parsing")
for L in tqdm(range(0, len(focal_list)+1)):
    for subset in itertools.combinations(focal_list, L):
        subset_sn += 1
        if len(subset) == 0: continue
        # print(list(subset))
        
        # trim OG containing non focal species
        # prepare a total_filter that is defo not in the focal list
        total_filter = orthogroup[s2l('Acomosus_321_v3.protein.fa')].isnull()
        for species_id in species_csv.set_index('species_id').index:
            # OG have focal sp not null
            if species_id in subset:
                total_filter = total_filter & orthogroup[s2l(species_id)].notnull()
                continue
            # every not focal species should be mull
            species_filter = orthogroup[s2l(species_id)].isnull()
            total_filter = total_filter & species_filter
            
        current_table = orthogroup[total_filter]
        orphan_dict[subset] = current_table.dropna(axis=1)
        count += current_table.shape[0]
        print(f"{list(subset)}: {current_table.shape[0]}")

        gene_list = []
        reference_list = []
        seq_record_list = []
        subset_str = ""
        for sp in subset:
            subset_str = subset_str + sp + "_"
            for orphan_list in [ og.split(', ') for og in current_table[s2l(sp)] ]:
                for orphan in orphan_list:
                    gene_list.append(orphan)
                    reference_list.append(sp)
                    seq_dict = species_pep_dict[sp]
                    seq_record_list.append(seq_dict[orphan])
        file_name = f"orphan_{subset_sn}_{len(subset)}_{len(gene_list)}.fasta"
        SeqIO.write(seq_record_list, file_name, "fasta")
        subset_dict[subset_sn] = (subset, current_table.shape[0])
import json
with open('subset_info.json', 'w') as fp:
    json.dump(subset_dict, fp)
print(f"A total of {count} orphan genes were identified")




