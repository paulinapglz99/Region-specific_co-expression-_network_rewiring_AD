"""Rule graph of the pipeline, with the rules colored by block.

Not part of the Snakemake workflow: run it by hand, from the repository root,
with the driver environment (workflow/envs/ad-networks-smk.yaml) and Graphviz
(`dot`) installed:
    python workflow/scripts/pipeline_rulegraph.py

Output: docs/figures/pipeline_rulegraph.pdf and .png
Rules not listed in BLOCKS are drawn in grey and reported, so add new rules here.
"""

import os
import re
import subprocess
import sys

OUT = "docs/figures/pipeline_rulegraph"

# Color of each block (same colors as docs/figures/experimental_design) and its rules
BLOCKS = {
    "data and subsampling": (
        "#8C8C8C",
        [
            "download_synapse",
            "select_samples",
            "split_expression",
            "plot_demographics",
            "subsample_matched",
            "subsample_equal_n",
            "subsample_curve",
        ],
    ),
    "ARACNe-AP networks": (
        "#3C5488",
        [
            "install_aracne",
            "aracne_matrix",
            "aracne_threshold",
            "aracne_bootstrap",
            "aracne_consolidate",
        ],
    ),
    "A: rewiring AD vs CT (1-3, 6)": (
        "#4DBBD5",
        ["network_metrics", "rewiring_index", "plot_network_metrics"],
    ),
    "B and C: networks vs DE (4, 5, 7)": (
        "#00A087",
        [
            "gene_rewiring",
            "rewiring_vs_de",
            "go_rewired",
            "hubs",
            "risk_genes",
            "plot_networks_de",
        ],
    ),
}
DEFAULT_COLOR = "#000000"  # rule all and any rule not listed above

colour_of = {rule: colour for colour, rules in BLOCKS.values() for rule in rules}

# Rule graph of the default target (all). The Snakefile asks for the Synapse
# token, which is not needed to draw the graph
env = {**os.environ, "SYNAPSE_AUTH_TOKEN": os.environ.get("SYNAPSE_AUTH_TOKEN", "unused")}
dot = subprocess.run(
    ["snakemake", "--rulegraph"], capture_output=True, text=True, check=True, env=env
).stdout

# Color of each rule
unlisted = []


def recolor(match):
    rule = match.group(1)
    if rule not in colour_of and rule != "all":
        unlisted.append(rule)
    return f'label = "{rule}", color = "{colour_of.get(rule, DEFAULT_COLOR)}"'


dot = re.sub(r'label = "([A-Za-z0-9_]+)", color = "[^"]+"', recolor, dot)
dot = dot.replace(
    "graph[bgcolor=white, margin=0];",
    "graph[bgcolor=white, margin=0.2, rankdir=TB, nodesep=0.3, ranksep=0.4];",
)

# Legend: one box per block, in a frame next to the graph
legend_nodes = "\n".join(
    f'        legend{i}[label = "{name}", color = "{colour}"];'
    for i, (name, (colour, _)) in enumerate(BLOCKS.items())
)
legend_order = " -> ".join(f"legend{i}" for i in range(len(BLOCKS)))
legend = f"""
    subgraph cluster_legend {{
        label = "Rules"; fontname = sans; fontsize = 10; color = grey80;
{legend_nodes}
        {legend_order} [style = invis];
    }}
}}"""
dot = dot.rstrip().removesuffix("}") + legend + "\n"

for fmt in ("pdf", "png"):
    subprocess.run(
        ["dot", f"-T{fmt}", "-Gdpi=200", "-o", f"{OUT}.{fmt}"],
        input=dot,
        text=True,
        check=True,
    )
print(f"Saved {OUT}.pdf and .png")

if unlisted:
    print("Rules drawn in black (add them to BLOCKS):", ", ".join(unlisted), file=sys.stderr)
