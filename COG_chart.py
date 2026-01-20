import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import argparse
import pathlib


parser = argparse.ArgumentParser(
    prog = "GO_analysis",
    usage = "%(prog)s [options]",
    add_help = True,
    description="Gene ontology clasification",
    )

parser.add_argument("-i", "--input", type=str, help="Input table of the %(prog)s program", required=True)
parser.add_argument("-o", "--output", type=str, help="Output graph of the %(prog)s program", required=True)
args = parser.parse_args()

COG_table = args.input
content = pd.read_csv(COG_table, sep="\t")


# Draw barchart
plt.figure(figsize=(30, 15))
ax = sns.barplot(data=content, x="Letter", y="Counts", palette=content["Color"])
ax.bar_label(ax.containers[0], fontsize=25)

# Title and Labels
plt.xlabel("Functional Category", fontsize=30, fontweight="bold")
plt.ylabel("Number of genes", fontsize=30, fontweight="bold")
plt.title("COG Functional Classification", fontsize=40, fontweight="bold")
plt.xticks(fontsize=25)  
plt.yticks(fontsize=25)       
plt.tight_layout()

# Legend
# plt.legend(
#     labels=content["Description"],
#     loc="center left",
#     bbox_to_anchor=(1, 0.5),
#     fontsize=16,
#     title="Functional Categories",   
#     title_fontsize=20,
#     frameon=False
# )
# Save figure
COG_barchart = args.output
plt.savefig(COG_barchart, 
            dpi=400)
