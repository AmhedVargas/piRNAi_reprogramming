
<h1 align="center">piRNA-mediated interference (piRNAi)</h1>
<p align="center">
Main repository containing pipelines, data, and software described in: <br><b>Reprogramming an endogenous small RNA pathway for multiplexed and transgenerational gene silencing in <i>C. elegans</i></b>
</p>

## Abstract
Short guide RNAs allow easy CRISPR/Cas targeting to specific DNA or RNA sequences and has enabled a wealth of genetic tools. In this project, we repurpose an endogenous small RNA pathway for gene silencing in the germline of *C. elegans* by expressing short guide piRNAs (21 nucleotides) from extra-chromosomal DNA. piRNA-mediated interference ("piRNAi") allows multiplexed gene silencing, is more efficient than RNAi, and can be reversibly controlled by auxin. We show in our manuscript that piRNAi decreases mRNA levels, increases secondary small interfering RNAs, and results in repressive chromatin modifications at the target locus. Gene silencing is highly specific and can be induced with short (300 bp) fragments amplified from array-based oligo pools, making piRNAi compatible with large-scale libraries. We use piRNAi to study transgenerational silencing of (*oma-1*) and two other genes (*him-5* and *him-8*). Epigenetic silencing is inherited for a limited time (4-6 generations) when gene-specific piRNAs are lost whereas inherited silencing becomes permanent (RNAe) when the entire piRNA pathway is impaired. We propose piRNAi as a novel method to study germline processes and epigenetic inheritance.

## Repository description
This github repository contains the code used to map and analize short RNA libraries as well as DNA libraries. For information regarding our shiny app and the data generated for its development please see below.

### Folder structure
This repository is structured into two main folders:
1. DNA libraries
Folder structured into three subfolders where DNA sequences were mapped to quantify the content of extra chromosomal arrays. See shell script (.sh) for more info.
2. RNA libraries
Folder containing script used to map short RNA sequences to quantify the expression of endogenous and synthetic piRNA sequences. See shell script (.sh) for more info.

### Data download
We have deposited DNA and RNA libraries under the [GEO](https://www.ncbi.nlm.nih.gov/geo/) accesion number [GSE165210](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165210)

## Additional information
We have build a [shiny web app](https://wormbuilder.dev/piRNAi/) to ease the construction of piRNAi fragments. 
It's code can be found in this repository: https://github.com/AmhedVargas/piRNAi_dev

Originally, the data to produce the app was processed using BLAST search as a proxy to identify thew number of mis-matches between a piRNAi guide and any other part of the genome. How the data was processed can be found here: https://github.com/AmhedVargas/piRNAi-DB

However, in the latest version of the piRNAi app we opted to use an exhaustive algorithm based on the calculation of the Hamming distance. Our c++ implementation can be found here: https://github.com/AmhedVargas/CelegansHammingAlignments

## Contact us
Please feel free to send any question or comment to [me](mailto:avargas0lcg@gmail.com) or [Christian](mailto:cfjensen@kaust.edu.sa)
