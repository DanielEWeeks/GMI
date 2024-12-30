##   GMI: Genetic Map Interpolator
##   Copyright (C) 2010 Nandita Mukhopadhyay, Xinyu Tang, Daniel E. Weeks
## 
##   This file is part of the GMI program, which is free software; you
##   can redistribute it and/or modify it under the terms of the GNU
##   General Public License as published by the Free Software Foundation;
##   either version 2 of the License, or (at your option) any later
##   version.
## 
##   GMI is distributed in the hope that it will be useful, but WITHOUT
##   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
##   FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
##   for more details.
## 
##   You should have received a copy of the GNU General Public License
##   along with this program; if not, write to the Free Software
##   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
## 
##   For further information contact:
##       Daniel E. Weeks
##       e-mail: weeks@pitt.edu
##
## ===========================================================================
library(biomaRt)
snpmart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")

#Obtain the current Ensembl version

lmrts <- listMarts()
ever <- as.character(lmrts[lmrts[,1] == "ENSEMBL_MART_SNP", 2])

ver <- strsplit(ever, "[ ]+")[[1]][3]
write(file="version.txt", as.integer(ver))
q()

