# HEATraN

## About

HEATraN, litteraly ***H**yper-**E**xpression **A**nalysis **T**ool **ra**mpantly developed in **N**ormandy* is a bioinformatics analysis tool dedicated to transcriptomic analysis. 
It was developed as part of a student project in the Bioinformatics Master of Rouen Normandy University.

<img src='./www/logo.png' width='256' alt='HEATraN logo' style='display:block;margin-left: auto;margin-right: auto;'>


You can find its last version on its [GitHub](https://github.com/TheLokj/HEATraN).

## Requirements

| **Library** | **Version** |
|-------------|-------------|
| R           | 4.3.3       |
| shiny       | 1.9.1       |
| shinyalert  | 3.1.0       |
| shinyjs     | 2.1.0       |
| ggplot2     | 3.5.1       |
| readxl      | 1.4.3       |
| fresh       | 0.2.1       |

## Repository structure

```
HEATraN/
├── app.R                # Main script                
├── bin/                 # Directory containing secondary scripts
│   ├── server.R         # Shiny Server script
│   └── ui.R             # Shiny UI script
├── data/                # Data directory
│   └── example.tsv      # Differentially expressed results example file
└── www/                 # Additional graphics resources directory
    ├── logo.png         # HEATraN logo 
    └── style.css        # HEATraN custom CSS stylesheet
```

The provided example is from a [Galaxy training by Maria Doyle](https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-viz-with-volcanoplot/tutorial.html). 