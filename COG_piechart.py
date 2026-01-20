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

# Filter and sort value 
df_filtered = content[content["Counts"] > 0].sort_values(by="Counts", ascending=False).copy()

# Calculate sum of counts
total = df_filtered["Counts"].sum()
letters = df_filtered["Letter"].tolist()

# Print % if percent > 3%
def make_autopct(letters, threshold):
    def my_autopct(pct):
        index = my_autopct.index
        letter = letters[index]
        my_autopct.index += 1
        if pct >= threshold:
            return f"{pct:.1f}%"
        else:
            return ""
    my_autopct.index = 0
    return my_autopct

# Draw piechart
plt.figure(figsize=(23, 10))
plt.pie(
    df_filtered["Counts"],
    labels=None,
    colors=df_filtered["Color"],
    autopct=make_autopct(letters, threshold=3),
    startangle=90,
    pctdistance=0.85,
    counterclock=False,
    textprops={'fontsize': 20, 'weight': 'bold'}
)

# Add legend
plt.legend(
    labels=df_filtered["Description"],
    loc="center left",
    bbox_to_anchor=(1, 0.5),
    fontsize=16,
    title="Functional Categories",   
    title_fontsize=20,
    frameon=False
)

# Title
plt.title("COG Functional Classification", fontsize=30, fontweight='bold')
plt.tight_layout()


# Save figure
COG_piechart = args.output
plt.savefig(args.output, dpi=400)