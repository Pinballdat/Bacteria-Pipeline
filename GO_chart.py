import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse

plt.figure(figsize=(30,10))

parser = argparse.ArgumentParser(prog="GO_chart.py", description="Program for drawing Gene Ontology level 1 graph")
parser.add_argument("-i", "--input", type=str, help="Input table of the %(prog)s program", required=True)
parser.add_argument("-o", "--output", type=str, help="Output graph of the %(prog)s program", required=True)
args = parser.parse_args()

# GO_table_input = "/data18tb/datnguyen/LuanChu/test/BS1_table.go.tsv"
GO_table_input = args.input

content = pd.read_csv(GO_table_input, sep="\t")
# print(content)
df = pd.DataFrame(content.groupby(["Category"])["Number of genes"].sum())
# print(df)
ax = sns.barplot(data=df, x="Number of genes", y=df.index, orient = 'h', palette=["#3274a1", "#e1812c", "#3a923a"])
ax.bar_label(ax.containers[0], fontsize=25)

plt.xlabel("Number of genes", fontsize=25, fontweight="bold")
plt.ylabel("Function", fontsize=25, fontweight="bold")
plt.title("GO Enrichment Analysis", fontsize=30, fontweight="bold")
plt.xticks(fontsize=20)  
plt.yticks(fontsize=20)  
# fig.tick_params(labelsize=20)

plt.tight_layout()
# plt.savefig("test.png")

GO_figure_output = args.output
plt.savefig(GO_figure_output, dpi=400) 
