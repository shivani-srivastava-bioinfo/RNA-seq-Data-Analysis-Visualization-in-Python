# RNA-seq Data Analysis & Visualization in Python

Python scripts for RNA-seq visualization and DEG exploration- heatmap, PCA, volcano, and violin plots- using pandas, numpy, seaborn, scikit-learn, and scipy on *E. coli* tetracycline vs. control data.

This is my Python-based downstream analysis for RNA-seq data- I take the gene count table (the output from my RNA-seq pipeline built earlier in Bash/R) and use Python to explore it further: heatmaps, PCA plots, finding differentially expressed genes (DEGs), a volcano plot, and a violin plot for individual genes.

**Dataset used here:** 4 samples total- 2 Control (`SRR22578513`, `SRR22578515`) and 2 Tetracycline-treated (`SRR22578517`, `SRR22578535`), same samples as used in the main R/DESeq2 pipeline- columns appear in that order in `count_matrix.csv`.

## What this covers

1. Loading and exploring the count data
2. Basic stats- mean and median expression per sample
3. Heatmap of gene expression
4. PCA plot- to see how similar/different my samples are
5. Finding DEGs
6. Volcano plot of the results
7. Violin plot for a single gene of interest

## Python Libraries I Used

pandas, numpy, matplotlib, seaborn, scikit-learn, and scipy are all **Python libraries** (also called **packages**)- not built into plain Python, they're separate collections of pre-written code that I install (`pip install ...`) and then `import` to add specific functionality like data tables, math, plotting, or statistics.

- **pandas**– a Python library for working with tables of data (rows and columns), similar to Excel but inside code. I use it to load my `count_matrix.csv` file, look at rows/columns, calculate averages, and filter genes based on conditions.

- **numpy** – a library for doing math on large sets of numbers quickly. I use it mainly for the log2 transform (`np.log2()`), which compresses very large and very small expression numbers so they're easier to compare and plot.

- **matplotlib** – the base Python library for making plots and charts (lines, scatter plots, titles, axis labels, etc.). Most other plotting libraries, including seaborn, are actually built on top of it. I use it directly for the volcano plot and PCA plot, and to control titles/labels/saving on all my figures.

- **seaborn** – built on top of matplotlib, but makes certain plots (like heatmaps and violin plots) much easier to create with nicer default styling. I use it for the heatmap and the violin plot.

- **scikit-learn (sklearn)** – a machine learning library. I don't use it for machine learning here- just two of its tools: `StandardScaler` (which rescales all genes so they're on a comparable scale before PCA) and `PCA` (Principal Component Analysis, which reduces complex gene expression data down to 2 dimensions so I can plot and see how similar or different my samples are).

- **scipy** – a scientific computing library. I use one specific function from it, `ttest_ind`, which runs a t-test- a statistical test that checks whether the average expression of a gene is really different between two groups (Control vs Treated), or if the difference could just be random chance.

## Folder structure

```
rnaseq_python_analysis/
├── count_matrix.csv        # gene count table (input, from the RNA-seq pipeline)
├── analysis.py               # the main script
└── DEG_results.csv           # output: significant genes found by the t-test
```

> **Note:** as written, the script displays each plot on screen with `plt.show()`- it doesn't save the plots as image files. If you want to keep copies of the plots (e.g. to show on GitHub), see [Saving the plots as image files](#saving-the-plots-as-image-files-optional) below.

## Requirements

```bash
pip install pandas numpy matplotlib seaborn scikit-learn scipy
```

---

## Step 1: Load the count table

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.stats import ttest_ind

# read the gene count table, use the gene ID column as the row index
counts = pd.read_csv("count_matrix.csv", index_col=0)

# check the first few rows to make sure it loaded correctly
print(counts.head())
```

* This just opens my count table (genes as rows, samples as columns) into a table I can work with in Python.*

## Step 2: Basic stats- mean and median per sample

```python
# all sample columns (since index_col=0 already removed the gene ID column,
# I use "counts" directly here- not a slice of it)
sample_counts = counts

print(sample_counts.mean())
print(sample_counts.median())
```

*Simple explanation: mean/median tells me the average gene expression level in each sample- useful for a quick sanity check (samples with very different averages might need extra normalization).*

> **Bug I fixed:** in my first draft I wrote `sample_counts = counts.iloc[:, 1:]`, which accidentally dropped my first sample column (since `index_col=0` already removes the gene ID column- there's nothing extra left to skip). Fixed by just using `counts` directly, so all 4 samples are included in every step below.

## Step 3: Log-transform the data

```python
# log2 transform makes expression values easier to compare and visualize
# (+1 avoids errors from taking log of zero)
log_counts = np.log2(sample_counts + 1)
print(log_counts.head())
```

*Simple explanation: raw counts can range from 0 to millions, which makes plots hard to read. Log2 transforming compresses that range so patterns are easier to see.*

## Step 4: Heatmap of gene expression

```python
plt.figure(figsize=(10, 8))
sns.heatmap(log_counts.iloc[:50], cmap="viridis")
plt.title("E. coli Heatmap of 50 Genes")
plt.show()
```

*Simple explanation: a heatmap shows expression of many genes across all samples at once, colored by intensity- darker/lighter patches make it easy to spot genes that behave differently between samples.*

## Step 5: PCA plot- how similar are my samples?

```python
# transpose so each row is a sample (not a gene)- PCA needs samples as rows
X = log_counts.T

# scale the data so all genes are treated equally regardless of their range
scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# reduce everything down to 2 dimensions for plotting
pca = PCA(n_components=2)
pc = pca.fit_transform(X_scaled)

plt.figure(figsize=(8, 6))
plt.scatter(pc[:, 0], pc[:, 1], s=100)

# label each point with its sample name
for i, sample in enumerate(X.index):
    plt.text(pc[i, 0], pc[i, 1], sample)

plt.xlabel("PC1")
plt.ylabel("PC2")
plt.title("E. coli PCA Plot")
plt.show()
```

*Simple explanation: PCA plots each sample as a dot based on its overall gene expression pattern. Samples that are biologically similar (e.g. all the Control samples) should cluster close together, and Treated samples should cluster separately- a good visual check that my experiment worked as expected.*

## Step 6: Find differentially expressed genes (DEGs) with a t-test

```python
# split samples into Control and Treated groups
# (first 2 columns = Control, next 2 = Treated, based on my sample order:
# SRR22578513, SRR22578515 = Control | SRR22578517, SRR22578535 = Treated)
control = log_counts.iloc[:, 0:2]
treated = log_counts.iloc[:, 2:4]

# average fold-change: how much higher/lower is Treated vs Control
log2FC = treated.mean(axis=1) - control.mean(axis=1)

# run a t-test for every gene to check if the difference is statistically significant
pvalues = []
for gene in log_counts.index:
    stat, p = ttest_ind(control.loc[gene], treated.loc[gene])
    pvalues.append(p)

results = pd.DataFrame({
    "Gene": log_counts.index,
    "log2FC": log2FC,
    "Pvalue": pvalues
})

# keep only genes that are both statistically significant AND changed a meaningful amount
deg = results[
    (results["Pvalue"] < 0.05) &
    (abs(results["log2FC"]) > 1)
]

print(deg.head())
deg.to_csv("DEG_results.csv", index=False)
```

*Simple explanation: for every gene, I compare its expression in Control vs Treated samples using a t-test, which tells me if the difference is likely real or just random noise.*

> **Worth knowing:** with only 2 samples per group, a plain t-test has very little statistical power- it's easy to miss real changes or flag noisy ones. This is exactly why the main pipeline uses **DESeq2** (Step 11 in the main pipeline README) instead of a plain t-test for the "official" DEG list- DESeq2 is built to handle small sample sizes properly. This Python t-test version is a good learning/exploration exercise and quick visual sanity check, but I treat `Significant_DEGs.csv` from DESeq2 as the primary result, not `DEG_results.csv` from this script.

**What each condition in the filter means:**
- `Pvalue < 0.05` → there's less than a 5% chance this gene's difference happened by random chance (i.e., it's statistically significant).
- `abs(log2FC) > 1` → the gene's expression changed by at least 2-fold (doubled or halved) between groups- a meaningful biological change, not just a tiny wobble.
- Both conditions together (`&`) → only genes that pass **both** tests are called significant DEGs- significant *and* biologically meaningful.

## Step 7: Volcano plot

```python
results["minuslog10P"] = -np.log10(results["Pvalue"])

plt.figure(figsize=(8, 6))
plt.scatter(results["log2FC"], results["minuslog10P"], s=15)
plt.xlabel("log2 FoldChange")
plt.ylabel("-Log10(P-value)")
plt.title("Volcano Plot")
plt.show()
```

*Simple explanation: a volcano plot shows every gene at once- how much it changed (x-axis) vs. how significant that change is (y-axis). Genes in the top-left and top-right corners are the most interesting (big change + highly significant).*

## Step 8: Violin plot for one specific gene

```python
gene = "b0002"
expression = counts.loc[gene]

violin_df = pd.DataFrame({
    "Expression": expression.values,
    "Condition": ["Control", "Control", "Treated", "Treated"]
})

plt.figure(figsize=(6, 5))
sns.violinplot(x="Condition", y="Expression", data=violin_df)
sns.stripplot(x="Condition", y="Expression", data=violin_df, color="black", size=5)

plt.title(f"Expression of {gene}")
plt.xlabel("Condition")
plt.ylabel("Expression")
plt.show()
```

*Simple explanation: once I have a gene I care about (e.g. from the DEG list), a violin plot shows its exact expression value in every sample, grouped by Control vs Treated- useful for double-checking a specific result visually.*

> **Note:** with only 2 samples per group, the "violin" shape itself won't show much of a distribution- the black dots (from `stripplot`) are really doing the work here, showing the 2 actual data points per condition.

> **Note:** change `gene = "b0002"` to any gene ID from `DEG_results.csv` that you want to look at closely.

---

## Output

As written, this script only produces one saved file:

| File | What it is |
|---|---|
| `DEG_results.csv` | genes that passed both the p-value and fold-change filters |

The heatmap, PCA plot, volcano plot, and violin plot are only **displayed on screen** (`plt.show()`)- they aren't saved as image files, so they close once you close the plot window and don't appear anywhere in the project folder.

## Saving the plots as image files (optional)

If I want to keep the plots as files (e.g. to show them on GitHub), I add one line before each `plt.show()`:

```python
plt.savefig("heatmap.png")
plt.show()
```

The same pattern (`plt.savefig("<name>.png")` right before `plt.show()`) works for the PCA plot, volcano plot, and violin plot too- just give each one a different filename. This isn't in my original script, so it's an optional addition if I decide I want image files later.

## Things to remember if I re-run this

- This script assumes the first 2 sample columns are Control and the next 2 are Treated- if the column order changes, `control`/`treated` slicing in Step 6 needs to be updated too.
- Make sure `count_matrix.csv` is in the same folder as the script, or update the path.

## License

MIT- see [LICENSE](LICENSE) for details.
