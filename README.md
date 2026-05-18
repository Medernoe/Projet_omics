# Projet_omics

DEGO est une application interactive open-source développée en R (Shiny), dédiée à l'analyse et à la visualisation de données transcriptomiques. Conçue pour les chercheurs et biologistes, elle permet d'explorer facilement les gènes différentiellement exprimés (DEG) et de réaliser des analyses d'enrichissement fonctionnel (GO, KEGG, Reactome) via des approches ORA et GSEA.

## Sommaire

- À propos du projet  
- Arborescence du projet  
- Prérequis et format des données  
- Fonctionnalités principales  
- Architecture technique et dépendances  
- Auteurs et remerciements  
- Références  

## À propos du projet

Ce projet a été développé dans le cadre du Master 2 Bioinformatique et Modélisation (BIMS) à l’Université de Rouen Normandie.

L’outil vise à simplifier l’analyse des données transcriptomiques en permettant :
- l’identification des termes GO significativement enrichis  
- la mise en évidence des voies biologiques pertinentes  
- une exploration intuitive des résultats via une interface Shiny  

## Arborescence du projet

Structure modulaire de l’application :

```
Projet_omics/
├── global.R # Initialisation de l'environnement, chargement des bibliothèques et scripts
├── ui.R # Interface utilisateur (Shiny Dashboard)
├── server.R # Logique réactive du serveur
├── fonctions.R # Logique métier (ORA/GSEA, nettoyage, calculs, graphiques)
├── www/ # Fichiers statiques
│ ├── styles.css
│ ├── logo-modified.png
│ ├── logouni.png
│ ├── logomaster.png
│ ├── img_deg.png
│ ├── img_go.png
│ ├── img_pathway.png
│ └── musique_attente.mp3
└── README.md
```


## Prérequis et format des données

### Organismes supportés

- Homo sapiens (humain)  
- Mus musculus (souris)  
- Drosophila melanogaster (drosophile)  

### Préparation du fichier

Les données issues de DESeq2 ou edgeR doivent être exportées en CSV avec séparateur point-virgule (;).

Le fichier doit contenir :

| Colonne  | Description |
|----------|-------------|
| GeneName | Symbole du gène (ex: TP53, BRCA1) |
| log2FC   | Log2 Fold Change |
| pval     | p-value ou p-value ajustée |

## Fonctionnalités principales

### Inspection des données (DEG)

- Volcano plot interactif avec seuils ajustables  
- Tableau interactif synchronisé avec le plot  

### Enrichissement ORA

Bases :
- Gene Ontology (BP, CC, MF)
- KEGG
- Reactome

Visualisations :
- Dotplot
- Barplot
- Cnetplot
- Emapplot
- Upsetplot
- Heatplot

### Enrichissement GSEA

Analyse globale basée sur le classement des gènes.

Visualisations :
- Ridgeplot
- GSEAplot2
- GSEArank

## Architecture technique et dépendances

### Packages CRAN

- shiny, shinydashboard, waiter, shinyalert
- ggplot2, plotly, ggarchery, qqman
- DT, dplyr  

### Packages Bioconductor

- clusterProfiler, ReactomePA, DOSE, enrichplot
- org.Hs.eg.db, org.Mm.eg.db, org.Dm.eg.db  

## Auteurs et remerciements

- Noé Méderlet (noe.mederlet@univ-rouen.fr)  
- Mehdi Tachekort (mehdi.tachekort@univ-rouen.fr)  
- Mathieu Cartier (mathieu.cartier@univ-rouen.fr)  
- Valentin Fourdraine (valentin.fourdraine@univ-rouen.fr)  

Encadrement :
- Hélène Dauchel  
- Solène Pety  

## Références

- Nguyen T. M. et al. (2019). Over-representation analysis: a comprehensive review. BMC Bioinformatics. DOI: 10.1186/s12859-019-2710-2  
- Subramanian A. et al. (2005). Gene set enrichment analysis. PNAS. DOI: 10.1073/pnas.0506580102  
- Wu T. et al. (2021). clusterProfiler 4.0. The Innovation. DOI: 10.1016/j.xinn.2021.100141  
- Yu G & He QY (2016). ReactomePA. Molecular BioSystems. DOI: 10.1039/C5MB00663E  
