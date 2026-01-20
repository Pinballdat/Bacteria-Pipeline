import pandas as pd
# from GO_analysis import sorted_final_df
import matplotlib.pyplot as plt 
import seaborn as sns 
import argparse
# final_BioPro_table["Category"] = "Biological Process"
# final_MolFun_table["Category"] = "Molecular Function"
# final_CelCom_table["Category"] = "Cellular Component"

# merged_df = pd.concat([final_BioPro_table, final_MolFun_table, final_CelCom_table])
plt.figure(figsize=(40,20))

parser = argparse.ArgumentParser(prog="GO_chart.lv2.py", description="Program for drawing Gene Ontology level 1 graph")
parser.add_argument("-i", "--input", type=str, help="Input table of the %(prog)s program", required=True)
parser.add_argument("-o", "--output", type=str, help="Output graph of the %(prog)s program", required=True)
args = parser.parse_args()

GO_analysis_output = args.input

content = pd.read_csv(GO_analysis_output)
content = content.set_index(content.columns[0])
# print(content)
colors = ["#3274a1", "#e1812c", "#3a923a"]
fig = sns.barplot(data=content, x=content["Number of genes"], y=content.index, orient = 'h', hue="Category", palette=colors, dodge=False)

# fig = sns.barplot(data=final_MolFun_table, x="Molecular Function", y=final_MolFun_table.index, orient = 'h', hue="Category", palette="Oranges")

# fig = sns.barplot(data=final_CelCom_table, x="Cellular Component", y=final_CelCom_table.index, orient = 'h', hue="Category", palette='Greens')
fig.bar_label(fig.containers[0], size=25)
fig.bar_label(fig.containers[1], size=25)
fig.bar_label(fig.containers[2], size=25)


# colors = [mcolors.to_hex(bar.get_facecolor()) for bar in fig.patches]
# print(colors)
fig.set_ylabel("Function", fontsize=30, fontweight="bold")
fig.set_xlabel("Number of genes", fontsize=30, fontweight="bold")
fig.tick_params(labelsize=30)

# plt.legend(title="Function", fontsize=30, title_fontsize=30, loc='lower left')
plt.legend(loc="upper right", fontsize=30)
plt.title("GO enrichment analysis", fontsize=30, fontweight="bold")
plt.tight_layout()

GO_lv2_graph = args.output
plt.savefig(GO_lv2_graph, dpi=400)  