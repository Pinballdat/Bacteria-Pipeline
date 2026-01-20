import pandas as pd
from collections import defaultdict
from itertools import islice
import argparse

parser = argparse.ArgumentParser(prog="GO_analysis.py", description="Program for analysing Gene Ontology function in Eggnog_mapper")
parser.add_argument("-i", "--input", type=str, help="Input table of %(prog)s program", required=True)
parser.add_argument("-o", "--output", type=str, help="Output table of %(prog)s program", required=True)
parser.add_argument("-db", "--database", type=str, help="Basic database of %(prog)s program", required=True)

args = parser.parse_args()
GO_analysis_database = args.database

# GO_analysis_database = "/data18tb/datnguyen/LuanChu/test/go-basic.obo"
with open(GO_analysis_database, "r") as hd:
    content = hd.read().split("\n\n[Term]\n")[1:]
# print(content)
bio_pro = {}
mol_fun = {}
cel_com = {}
Name_sp_db = {}
for c in content:
    if len(c.split("\n")) >= 3:
        Name_sp_db[c.split("\n")[0].split(": ")[1]] = c.split("\n")[2].split(": ")[1]
        if c.split("\n")[2].split(": ")[1] == "biological_process":
            bio_pro[c.split("\n")[0].split(": ")[1]] = c.split("\n")[1].split(": ")[1]
        elif c.split("\n")[2].split(": ")[1] == "molecular_function":
            mol_fun[c.split("\n")[0].split(": ")[1]] = c.split("\n")[1].split(": ")[1]
        elif c.split("\n")[2].split(": ")[1] == "cellular_component":
            cel_com[c.split("\n")[0].split(": ")[1]] = c.split("\n")[1].split(": ")[1]

GO_analysis_input = args.input
# GO_analysis_input = "/data18tb/datnguyen/LuanChu/03.Results/12.eggnog_mapper/BS1.emapper.emapper.annotations"
# with open(GO_analysis_input, "r") as hd:
#     for line in hd:
#         if "##" not in line:
#             print(line)
GO_df = pd.read_csv(GO_analysis_input, 
                    sep="\t", 
                    comment="#",
                    names=["Query", "Seed_ortholog", "Evalue", "Score", "EggNOG_OGs", "Max_annot_lvl", "COG_category", "Description", "Preferred_name", "GOs", "EC", "KEGG_ko", "KEGG_Pathway", "KEGG_Module", "KEGG_Reaction", "KEGG_rclass", "BRITE", "BRITE_TC", "CAZy", "BiGG_Reaction", "PFAMs"],
                    usecols = ["Query", "GOs"],
                    )
# print(GO_df)
GO_dict = {}
for i in GO_df.index:
    GO_dict[GO_df.loc[i,"Query"]] = GO_df.loc[i,"GOs"].split(",")
# print(GO_dict)
Filtered_list = []
for value in GO_dict.values():
    for v in value:
        if v != "-":
            Filtered_list.append(v)
# print(Filtered_list)
BioPro_dict = defaultdict(int)
MolFun_dict = defaultdict(int)
CelCom_dict = defaultdict(int)
for value in Filtered_list:
    if value in bio_pro.keys():
        BioPro_dict[bio_pro[value]] += 1
    if value in mol_fun.keys():
        MolFun_dict[mol_fun[value]] += 1
    if value in cel_com.keys():
        CelCom_dict[cel_com[value]] += 1
sorted_BioPro_dict = dict(sorted(BioPro_dict.items(), key=lambda item: item[1]))
sorted_MolFun_dict = dict(sorted(MolFun_dict.items(), key=lambda item: item[1]))
sorted_CelCom_dict = dict(sorted(CelCom_dict.items(), key=lambda item: item[1]))
# print(sorted_BioPro_dict)
merged_dict = {**sorted_BioPro_dict, **sorted_MolFun_dict, **sorted_CelCom_dict}
sorted_merged_dict = dict(sorted(merged_dict.items(), key=lambda item: item[1]))

final_df = pd.DataFrame.from_dict(sorted_merged_dict, orient="index", columns=["Number of genes"])
# print(final_df)
for i in final_df.index:
    if i in sorted_BioPro_dict.keys():
        final_df.loc[i, "Category"] = "Biological Process"
    if i in sorted_MolFun_dict.keys():
        final_df.loc[i, "Category"] = "Molecular Function"
    if i in sorted_CelCom_dict.keys():
        final_df.loc[i, "Category"] = "Cellular Component"

sorted_final_df = final_df.sort_values("Number of genes", ascending=True)

GO_analysis_output = args.output
sorted_final_df.to_csv(GO_analysis_output, sep="\t", index_label="Namespace")

