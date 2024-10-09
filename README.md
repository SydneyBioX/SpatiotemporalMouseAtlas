---
title: "README"
author: "Shila Ghazanfar"
date: "09/07/2024"
output:
  html_document: default
  word_document: default
---

<center>

#### A Spatiotemporal Atlas of Mouse Gastrulation and Early Organogenesis to Explore Axial Patterning and Project In Vitro Models onto In Vivo Space:

<a href="TBA" target="_blank">bioRxiv Link (to come)</a>

</center>

<br>

#### Abstract

At the onset of murine gastrulation, pluripotent epiblast cells migrate through the primitive streak, generating mesodermal and endodermal precursors, while the ectoderm arises from the remaining epiblast. Together, these germ layers establish the body plan, defining major body axes and initiating organogenesis. Although comprehensive single cell transcriptional atlases of dissociated mouse embryos across embryonic stages have provided valuable insights during gastrulation, the spatial context for cell differentiation and tissue patterning remain underexplored. In this study, we employed spatial transcriptomics to measure gene expression in mouse embryos at E6.5 and E7.5 and integrated these datasets with previously published E8.5 spatial transcriptomics1 and a scRNA-seq2 atlas spanning E6.5 to E9.5. This approach resulted in a comprehensive spatiotemporal atlas, comprising over 150,000 cells with 88 refined cell type annotations as well as genome-wide transcriptional imputation during mouse gastrulation and early organogenesis. The atlas facilitates exploration of gene expression dynamics along anterior-posterior and dorsal-ventral axes at cell type, tissue, and organismal scales, revealing insights into mesodermal fate decisions within the primitive streak. Moreover, we developed a bioinformatics pipeline to project additional scRNA-seq datasets into a spatiotemporal framework and demonstrate its utility by analysing cardiovascular models of gastrulation3. To maximise impact, the atlas is publicly accessible via a user-friendly web portal empowering the wider developmental and stem cell biology communities to explore mechanisms of early mouse development in a spatiotemporal context.

<br>

#### Experimental Design

<center><img src="Figure1.png" width="75%"/></center>

Figure 1: Exploring gene expression patterns during mouse gastrulation using <a href="https://dx.doi.org/10.1038%2Fnmeth.2892" target="_blank">seqFISH</a> 
(a, c) 3D illustrations of Theiler stage (TS)10 (a) and TS11 (c) mouse embryos, adapted from eMouseAtlas. Dotted red lines mark the estimated positions of sagittal optical tissue sections shown in (b) and (d). Orientation abbreviations: D, distal; V, ventral; R, right; L, left; A, anterior; P, posterior; PR, proximal.
(b, d) Tile scans of 4-μm sagittal sections from two independently sampled E6.5 embryos (b) and one E7.5 embryo (d), imaged using seqFISH with DAPI nuclear staining (white).
(e) Schematic outlining the seqFISH pipeline. 
(f, g) UMAP projections generated from integrated seqFISH expression data. In (f), cells are coloured by their embryo of origin, and in (g), by cell type.
(h) Spatial maps of E6.5 and E7.5 embryos, with cells coloured according to their cell types. The black dotted line indicates the extraembryonic (ExE) embryonic (Em) boundary. BI, blood islands; EPC, ectoplacental cone.
(i) Representative visualisation of normalised log expression counts of selected genes, measured by seqFISH, to validate performance in E6.5 (top) and E7.5 (bottom) embryos.
(j) Stacked bar chart showing the proportion of cell types per seqFISH embryo.
(k) Dot plot displaying the average gene expression for marker genes across different cell types identified in the E6.5 and E7.5 seqFISH embryos.

<br>

#### Download links

-   <a href="https://maths.usyd.edu.au/u/sheilag/SpatioTemporalMouseAtlas/seqFISH.h5ad" target="_blank">seqFISH observed data (450MB)</a>.
-   <a href="https://maths.usyd.edu.au/u/sheilag/SpatioTemporalMouseAtlas/imputed.h5ad" target="_blank">imputed logcounts for seqFISH cells (2.7GB)</a>.
-   <a href="https://maths.usyd.edu.au/u/sheilag/SpatioTemporalMouseAtlas/joint_annotation.csv" target="_blank">joint annotation for seqFISH and scRNA-seq data (64MB)</a>.

<br>

#### Contact

Please contact Shila Ghazanfar (shila.ghazanfar-at-sydney.edu.au) with any web app-related queries or suggestions.
