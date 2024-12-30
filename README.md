# Genetic Map Interpolator (GMI)

Bioconductor version          
Version 1.5 February 10, 2016

Originally written by Xinyu Tang and Nandita Mukhopadhyay under the supervision of Daniel E. Weeks. Updated by Robert V. Baron.

Copyright © 2010-24 by Xinyu Tang, Nandita Mukhopadhyay, Robert V. Baron, Daniel E. Weeks, and the University of Pittsburgh.

Supported by NIH grant R01GM076667 "Mega2: Manipulation Environment for Genetic Analyses" (PI Daniel E. Weeks) and the University of Pittsburgh.

# Overview

Many statistical algorithms for analyzing genetic data require genetic maps of the markers, which specify the recombination rates between adjacent markers. However, while it is relatively easy to extract physical map positions from the online databases, it is more difficult to extract genetic map positions. Furthermore, some of the genetic map tools require that one input not only the marker ID, but also the marker's physical position in a specified older build of the genome.

The Genetic Map Interpolator (GMI) package is designed to create interpolated genetic maps of single nucleotide polymorphism (SNP) markers. Starting from a list of single nucleotide polymorphism (SNP) rs numbers, GMI uses the R packages biomaRt and RMySQL to fetch the most up-to-date SNP and microsatellite physical positions from Ensembl, and then combines these with the Rutgers combined genetic and physical map (Matise et al., 2007; Kong et al., 2004; Kong and Matise, 2004) to estimate the corresponding Kosambi genetic positions by linear interpolation for these SNPs. The resulting information is then output in map files formatted to be read in by our data-reformatting program, Mega2 (Mukhopadhyay et al., 2005; http://watson.hgen.pitt.edu/register/).

SNPs were assigned rs IDs at the time of discovery, and, later if another SNP was found to be identical to a previously discovered SNP, their rs IDs were merged to one representative rs ID (usually the one with the lower number). If a SNP is not found in Ensembl, GMI then queries NCBI’s Entrez SNP database using eutils to see if this SNP has been renamed; if so, GMI then searches Ensembl a second time with the new SNP ID. The user has the option to use either the old or the new ID in the output map.

For each of these database queries, the list of SNPs being queried is broken up into smaller chunks so that each individual query can be completed within a few seconds. This has been done to prevent premature timing out of web connections, and in order to obey the limit guidelines required for Entrez queries.
 
GMI is a Unix-based program written in Perl and R, and is portable over most Unix platforms. We have used and tested it extensively on both Intel and PPC-based Macs. Setup is done via a text-based user-interface within a Unix shell. Subsequently, queries can be executed from the Unix command line.
GMI has a number of useful features (see Figure 1 in the manual 'GMI_manual_v1.5.pdf'):

* GMI automatically looks up and uses the physical positions for each marker from the most recent Ensembl build.
* For each SNP that is not initially found in Ensembl, GMI automatically checks to see if that SNP has been assigned a new rs number.
* GMI automatically figures out which chromosome each marker is on, so knowledge
of a marker's chromosome is not required to run GMI.

Please note that GMI relies crucially upon the Rutgers Combined Linkage-Physical Map, created by Tara Matise and colleagues (Matise et al., 2007; Kong et al., 2004; Kong and Matise, 2004), and we would like to thank them for making this important well-validated framework map available to the scientific community. GMI would not be possible without it.

# Documentation

Please see the PDF file [GMI_manual_v1.5.pdf](https://github.com/DanielEWeeks/GMI/blob/main/GMI_bioc_v1.5/GMI_manual_v1.5.pdf) for complete documentation.

# Contact

Daniel E. Weeks       
weeks@pitt.edu        
