########################
# Importing Dependencies
########################

import pandas as pd
import matplotlib.pyplot as plt

######################
# Function Definitions
######################

# Function to get data and labels for the piechart
def get_piechart_requirements(homology_data):
    not_found = len(homology_data[homology_data['human_gene'] == 'Not found'])
    confidence_0 = len(homology_data[(homology_data['human_gene'] != 'Not found') & (homology_data['confidence'] == 0)])
    confidence_1 = len(homology_data[(homology_data['human_gene'] != 'Not found') & (homology_data['confidence'] == 1)])
    data = [not_found, confidence_0, confidence_1]
    labels = ['Not found', 'Confidence 0', 'Confidence 1']
    return data, labels

# Custom function to format labels with count and percentage
def autopct_format(pct, all_vals):
    absolute = int(round(pct/100.0 * sum(all_vals)))
    return f'{absolute}\n({pct:.1f}%)'

# Function to generate the piechart
def generate_piechart(data, labels):
    plt.figure(figsize=(10, 10))
    wedges, texts, autotexts = plt.pie(
        data, 
        labels=labels, 
        autopct=lambda pct: autopct_format(pct, data),
        textprops=dict(color="black"),
        startangle=140,
        colors=['#ff9999','#66b3ff','#99ff99']
    )

    for text in texts:
        text.set_fontsize(12)
    for autotext in autotexts:
        autotext.set_fontsize(10)

    plt.title('Mouse Genes to Human Homologs', fontsize=16)
    plt.tight_layout()

###### 
# Main 
######

if __name__ == "__main__":

    # Loading homology data
    homology_data = pd.read_csv('data/mouse_to_human_homologs.csv')

    # Getting data and labels for the piechart
    data, labels = get_piechart_requirements(homology_data)

    # Generating the piechart
    generate_piechart(data, labels)
    plt.show()
