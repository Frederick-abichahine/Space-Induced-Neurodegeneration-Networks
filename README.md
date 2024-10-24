# Research Proposal: Causal Network Analysis of Space-Induced Neurodegeneration in Mice with Human Gene Mapping

## Background
Spaceflight exposes organisms to unique environmental stressors, including microgravity and cosmic radiation, which can impact neurological health. While general effects have been studied, the specific causal mechanisms leading to neurodegenerative processes in space remain unclear. This study aims to delve deeper into these mechanisms using advanced causal inference techniques and explore their potential relevance to human health across multiple neurodegenerative diseases.

## Our Goal
Investigate the causal mechanisms of space-induced neurodegeneration using advanced causal inference techniques, with applications to human health.

## Research Objectives
1. Construct a network of gene interactions related to neurodegeneration in mice exposed to spaceflight conditions.
2. Focus on a set of neurodegenerative diseases within the broader analysis.
3. Use the PC algorithm to infer causal relationships between gene expression changes and neurodegenerative processes.
4. Map findings from mice to potential human impacts through comparative genomics.

## Data Source
NASA's Open Science Data Repository: OSD-352 dataset 
- Reference: https://osdr.nasa.gov/bio/repo/data/studies/OSD-352
- Includes gene expression data from mice exposed to spaceflight conditions
- Provides ground control data for comparison

## Methodology

### 1. Data Preprocessing
- Clean and normalize the OSD-352 dataset

### 2. Analysis of Differentially Expressed Genes (DEGs)
- Perform differential gene expression analysis between spaceflight and control mice
- Identify significantly altered genes

### 3. Functional Enrichment Analysis
- Perform functional enrichment analysis between spaceflight and control mice
- Identify significantly altered pathways

### 4. Extract DEGs Relevant to the Neurodegenerative Diseases
- Utilization of tools to extract the DEGs of the neurological disease of interest

### 5. PC Algorithm Utilization & Network Analysis
- Utilize the PC algorithm on the extracted DEGs to construct and analyze the causal network
- Identify key hub genes and central pathways in the causal network
- Visualize the network to highlight critical nodes and edges

### 6. Human Gene Mapping and Comparative Analysis
- Map mouse genes in the causal network to their human orthologs
- Conduct comparative analysis to identify conserved pathways between mice and humans
- Assess the potential relevance of findings to human astronauts

## Novelty and Significance
- Utilizing the PC algorithm to space neurodegeneration data provides a new perspective on causal mechanisms
- Focus on multiple neurodegenerative diseases offers insights into various neurological processes in space
- Mapping to human relevance enhances the impact for astronaut health and long-duration space missions

## Expected Outcomes
1. Identification of key genes and pathways related to various neurodegenerative processes in space
2. A comprehensive causal network of space-induced neurodegeneration in mice
3. Hypotheses about conserved neurodegenerative mechanisms between mice and humans in space environments

## Potential Challenges (Limitations & Future Studies) and Solutions
- Limited Data: Limited statistical power due to small sample size (n=6 female mice; 3 ground control vs. 3 spaceflight-exposed) and single time point data collection, restricting temporal analysis of spaceflight effects
- Confounding Variables: Potential confounding variables (environmental stressors, dietary variations) coupled with limited sample size challenge the ability to isolate direct spaceflight effects from other factors
