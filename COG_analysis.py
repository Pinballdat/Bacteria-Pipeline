import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
import pathlib 
import argparse

parser = argparse.ArgumentParser(
    prog = "GO_barchart",
    usage = "%(prog)s [options]",
    add_help = True,
    description="Gene ontology clasification",
    )

parser.add_argument("-i", "--input", type=str, help="Input path of the %(prog)s program", required=True)
parser.add_argument("-o", "--output", type=str, help="Output path of the %(prog)s program", required=True)
args = parser.parse_args()

COG_dict = {
    'A': ['RNA processing and modification', '#4E79A7'],
    'B': ['Chromatin structure and dynamics', '#F28E2B'],
    'C': ['Energy production and conversion', '#E15759'],
    'D': ['Cell cycle control, cell division, chromosome partitioning', '#76B7B2'],
    'E': ['Amino acid transport and metabolism', '#59A14F'],
    'F': ['Nucleotide transport and metabolism', '#EDC948'],
    'G': ['Carbohydrate transport and metabolism', '#B07AA1'],
    'H': ['Coenzyme transport and metabolism', '#FF9DA7'],
    'I': ['Lipid transport and metabolism', '#9C755F'],
    'J': ['Translation, ribosomal structure and biogenesis', '#BAB0AC'],
    'K': ['Transcription', '#1F77B4'],
    'L': ['Replication, recombination and repair', '#2CA02C'],
    'M': ['Cell wall/membrane/envelope biogenesis', '#FF7F0E'],
    'N': ['Cell motility', '#D62728'],
    'O': ['Posttranslational modification, protein turnover, chaperones', '#9467BD'],
    'P': ['Inorganic ion transport and metabolism', '#8C564B'],
    'Q': ['Secondary metabolites biosynthesis, transport and catabolism', '#E377C2'],
    'R': ['General function prediction only', '#7F7F7F'],
    'S': ['Function unknown', '#BCBD22'],
    'T': ['Signal transduction mechanisms', '#17BECF'],
    'U': ['Intracellular trafficking, secretion, and vesicular transport', '#AEC7E8'],
    'V': ['Defense mechanisms', '#98DF8A'],
    'W': ['Extracellular structures', '#FFBB78'],
    'Y': ['Nuclear structure', '#C5B0D5'],
    'Z': ['Cytoskeleton', '#C49C94']
}

COG_analysis_input = args.input
COG_df = pd.read_csv(COG_analysis_input, 
                    sep="\t", 
                    comment="#",
                    names=["query", "seed_ortholog", "evalue", "score", "eggNOG_OGs", "max_annot_lvl",	"COG_category",	"Description", "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "KEGG_TC", "CAZy", "BiGG_Reaction", "PFAMs"],
                    usecols = ["query", "COG_category"],
                    index_col = 0)

number_chars = defaultdict(int)
for i in COG_df.index:
    if COG_df.at[i, "COG_category"] != "-":
        for char in COG_df.at[i, "COG_category"]:
            number_chars[char] += 1

for key in COG_dict.keys():
    if key in number_chars:
        COG_dict[key].append(number_chars[key])
    else:
        COG_dict[key].append(0)

COG_table = pd.DataFrame.from_dict(COG_dict, orient="index", columns=["Description", "Color", "Counts"])

COG_output = args.output
COG_table.to_csv(COG_output, sep="\t", index=True, index_label="Letter")