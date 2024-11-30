########################
# Importing Dependencies
########################

import requests
import pandas as pd

#########################################################################
# Function to convert mouse genes to human homologs using Ensembl BioMart
#########################################################################

def get_human_homologs(mouse_genes):
    # BioMart URL and query
    server = "http://www.ensembl.org/biomart/martservice?"
    query = f"""<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="0" count="" datasetConfigVersion="0.6">
        <Dataset name="mmusculus_gene_ensembl" interface="default">
            <Attribute name="external_gene_name"/>
            <Attribute name="hsapiens_homolog_associated_gene_name"/>
            <Attribute name="hsapiens_homolog_orthology_confidence"/>
        </Dataset>
    </Query>
    """
    
    try:
        # Get all homology data
        r = requests.get(server + "query=" + query)
        r.raise_for_status()
        
        # Create DataFrame with all homology data
        all_homologs = pd.DataFrame([line.split('\t') for line in r.text.strip().split('\n')], columns=['mouse_gene', 'human_gene', 'confidence'])
        
        # Clean up the data
        all_homologs['human_gene'] = all_homologs['human_gene'].replace('', 'Not found')
        
        # Filter for requested genes
        mouse_genes_upper = [g.upper() for g in mouse_genes]
        result_df = all_homologs[all_homologs['mouse_gene'].str.upper().isin(mouse_genes_upper)].copy()
        
        # Add missing genes
        found_genes = set(result_df['mouse_gene'].str.upper())
        missing_genes = [(g, 'Not found', 'N/A') for g in mouse_genes if g.upper() not in found_genes]
        
        if missing_genes:
            missing_df = pd.DataFrame(missing_genes, columns=['mouse_gene', 'human_gene', 'confidence'])
            result_df = pd.concat([result_df, missing_df], ignore_index=True)
        
        # Ensure original order
        result_df['order'] = result_df['mouse_gene'].str.upper().map({g.upper(): i for i, g in enumerate(mouse_genes)})
        result_df = result_df.sort_values('order').drop('order', axis=1)
        
        return result_df
        
    except Exception as e:
        print(f"Error: {str(e)}")
        return pd.DataFrame({
            'mouse_gene': mouse_genes,
            'human_gene': ['Error'] * len(mouse_genes),
            'confidence': ['Error'] * len(mouse_genes)
        })
    
######
# Main 
######

if __name__ == "__main__":

    # Loading mouse DEGs from mRNA and single cell data
    mice_mrna_degs_full = pd.read_csv('data/mice_mrna_degs_full.csv')
    mice_single_cell_degs_full = pd.read_csv('data/mice_single_cell_degs_full.csv')

    # Processing data
    mice_mrna_degs_full = mice_mrna_degs_full[mice_mrna_degs_full['pvalue'] <= 0.05]
    mice_genes = list(mice_mrna_degs_full['gene_name'].append(mice_single_cell_degs_full['gene_name']))
    mice_genes_common = list(set(mice_mrna_degs_full['gene_name']).intersection(set(mice_single_cell_degs_full['gene_name'])))
    mice_genes_without_duplicates = list(set(mice_genes))

    # Printing some statistics before querying
    print("\n#########################\n" + 
          "Statistics for mouse DEGs\n" +
          "#########################\n")
    print(f"\t> Number of genes in mRNA DEGs: {len(mice_mrna_degs_full)}")
    print(f"\t> Number of genes in single cell DEGs: {len(mice_single_cell_degs_full)}")
    print(f"\t> Number of genes before removing duplicates: {len(mice_genes)}")
    print(f"\t> Number of genes after removing duplicates: {len(mice_genes_without_duplicates)}")
    print(f"\t> Number of common genes: {len(mice_genes_common)}")
    print(f"\t> Names of common genes: {mice_genes_common}")

    # Querying mouse genes
    print(f"\n---> Querying {len(mice_genes_without_duplicates)} mouse genes...")
    result = get_human_homologs(mice_genes_without_duplicates)

    # Printing the result
    print("\n#############################\n" + 
          "Human homologs of mouse genes\n" +
          "#############################\n")
    print(result)

    # Printing some statistics after querying
    print("\n########################\n" + 
          "Statistics of the result\n" +
          "########################\n")
    print(f"\t> Number of genes not found: {len(result[result['human_gene'] == 'Not found'])}")
    print(f"\t> Number of genes found with a confidence score of 0: {len(result[(result['human_gene'] != 'Not found') & (result['confidence'] == '0')])}")
    print(f"\t> Number of genes found with a confidence score of 1: {len(result[(result['human_gene'] != 'Not found') & (result['confidence'] == '1')])}")

    # Saving the result
    result.to_csv('data/mouse_to_human_homologs.csv', index=False)
