########################
# Importing Dependencies
########################

import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

##########################################################
# Relevant genes for each of the four diseases of interest
##########################################################

disease_genes = {
    "AD": ["GRIN2B", "KCNMA1", "PLCB1", "PPP2R2B", "NTRK2", "APBB2", "GSK3B", "CLOCK", "FYN", "APP", "SORL1"],
    "ALS": ["NRG1", "STMN2", "TIA1", "GSK3B", "FUS", "ADARB1", "DPP6", "ERBB4"],
    "PD": ["DNAJC6", "PRKN", "HTT", "SYNJ1", "ATXN2", "VPS13C"],
    "SCZ": ["GRIN2B", "NRG1", "ANK3", "CTNND2", "CNTNAP2", "FOXP1", "GSK3B", "FYN", "GRIA4", "SYN2", "KIF5C", "GRIN2A"]
}

#################################################################
# Function to plot the frequency of genes in AD, ALS, PD, and SCZ
#################################################################

def generate_gene_frequency_plot():
    genes = [gene for disease_genes in disease_genes.values() for gene in disease_genes]
    gene_count = {}

    for gene in genes:
        if gene in gene_count:
            gene_count[gene] += 1
        else:
            gene_count[gene] = 1
    
    gene_count = {k: v for k, v in sorted(gene_count.items(), key=lambda item: item[1], reverse=True)}
    gene_diseases = {}

    for gene in gene_count.keys():
        diseases = []
        diseases.append(1 if gene in disease_genes["AD"] else 0)
        diseases.append(1 if gene in disease_genes["ALS"] else 0)
        diseases.append(1 if gene in disease_genes["PD"] else 0)
        diseases.append(1 if gene in disease_genes["SCZ"] else 0)
        gene_diseases[gene] = diseases
    
    plt.style.use('seaborn')
    
    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(2, 1, height_ratios=[4, 1], hspace=0.05)
    ax1 = fig.add_subplot(gs[0])
    colors = sns.color_palette("viridis", max(gene_count.values()))
    x = np.arange(len(gene_count))
    bars = ax1.bar(x, gene_count.values(), 
                   edgecolor='white', linewidth=1.5, width=0.8)
    
    for i, bar in enumerate(bars):
        bar.set_facecolor(colors[gene_count[list(gene_count.keys())[i]]-1])
    
    ax1.set_xticks([])
    ax1.set_ylabel("Number of Diseases", fontsize=12, weight="bold", labelpad=10)
    ax1.set_title("Gene Frequency Across Neurological Disorders", 
                  fontsize=14, weight="bold", pad=20)
    
    for i, (gene, count) in enumerate(gene_count.items()):
        ax1.text(i, count, str(count), ha='center', va='bottom', fontsize=10, 
                weight='bold', color='black')
    
    ax1.grid(axis='y', linestyle='--', alpha=0.3)
    ax2 = fig.add_subplot(gs[1])
    ax2.set_facecolor('#f8f9fa')
    marker_size = 120

    for i, (gene, diseases) in enumerate(gene_diseases.items()):
        for j, present in enumerate(diseases):
            if present:
                ax2.scatter(i, j, color='#2c3e50', marker='o', s=marker_size, 
                          alpha=0.9, edgecolors='white', linewidth=1)
            else:
                ax2.scatter(i, j, color='none', marker='o', s=marker_size, 
                          edgecolors='#95a5a6', linewidth=1)
                
    for y in [-0.5, 0.5, 1.5, 2.5, 3.5]:
        ax2.axhline(y=y, color='#e0e0e0', linestyle='-', linewidth=0.5, zorder=0)
    
    ax2.set_xticks(x)
    ax2.set_xticklabels(gene_count.keys(), rotation=45, ha='right', fontsize=10)
    ax2.set_yticks(range(4))
    ax2.set_yticklabels(['AD', 'ALS', 'PD', 'SCZ'], fontsize=10, weight='bold')
    ax2.set_axisbelow(True)

    for x_pos in np.arange(-0.5, len(gene_count), 1):
        ax2.axvline(x=x_pos, color='#e0e0e0', linestyle='-', linewidth=0.5, zorder=0)
    
    ax1.set_xlim(-0.6, len(gene_count) - 0.4)
    ax2.set_xlim(-0.6, len(gene_count) - 0.4)
    
    sns.despine(ax=ax1)
    sns.despine(ax=ax2)
    plt.subplots_adjust(bottom=0.2, left=0.1, right=0.98)
    
    return plt

######
# Main
######

if __name__ == "__main__":
    generate_gene_frequency_plot()
    plt.show()