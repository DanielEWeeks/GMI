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

## Functions for querying Ensembl and Entrez
## Chromosome codes:
##            X, Y, MT = X, Y and MT (Y and MT have physical bit not genetics maps)
##            M = Missing from Ensembl (No maps)
##            D = Duplicates in Ensembl (Will be assigned Unknown chromsome)
##            UM = Unknown (no known physical position in Ensembl)
##            U = internal code (will be assigned to M, D, and UM SNPs)
##

gmi.bioc.scale.negative = TRUE
gmi.bioc.debug = FALSE

snp.bioc.mapping <- function(snpfile, output,
                             num, mpath, create.plots, use.aliases,
                             request.from)
{
  ## # Define the names of output files
  output.ann <- paste(output, "_annot.txt", sep="")
  output.non.ann <- paste(output, "_nonannot.txt", sep="")
  output.pdf <- paste(output, ".pdf", sep="")
  output.warnings <- paste(output, "_log.txt", sep="")
  
  ## Remove files with the same name as output files before those files are created
  if (file.exists(output.ann)) {
    file.remove(output.ann)
  }
  if (file.exists(output.non.ann)) {
    file.remove(output.non.ann)
  }
  if (file.exists(output.pdf)) {
    file.remove(output.pdf)
  }

  options(warn=-1)
  snplists <- snp.bioc.querysnps(snpfile, request.from,
                                 output.warnings, use.aliases)
  if (is.null(snplists)) {
    return(FALSE)
  }

  log.query.problems(snplists, output.warnings, use.aliases)

  retval <-
    do.bioc.mapping(snplists, output, num, mpath, create.plots, use.aliases)
  
  return(retval)
  
}

snp.bioc.querysnps <- function(snpfile, request.from, output.warnings, use.aliases)
{
  library(biomaRt)
  library(XML)

  snpmart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
  
  ## allocate the maximum needed, then we will shorten the array
  attr <- c("Name", "Ensembl.chr", "Ensembl.p")
  
  num.returned <- 0
  chunk.size <- 5000
  snplists <- NULL
  
  ## #--------------------#
  ## #                    #
  ## #    Read in Data    #
  ## #                    #
  ## #--------------------#
  
  ## # Read in the original snplist
  
  if (file.exists(snpfile) == TRUE) {
    
    snplist <- read.table(snpfile, header=F, stringsAsFactors=F)

  } else {
    cat(paste("ERROR: File", snpfile, "does not exist.", sep=" "))
    return(NULL)
    
  }
  
  if(dim(snplist)[2] == 3) {
    names(snplist) <- c("Name", "User.chr", "User.p")
  }
  else if (dim(snplist)[2] == 1) {
    names(snplist) <- "Name"
  }

  ##  snplist$Name <- tolower(snplist$Name)

  a <- data.frame(matrix(nrow=dim(snplist)[1], ncol=length(attr)))
  names(a) <- attr
  dup <- NULL
  
  ##--------------------##
  ##                    ##
  ## Query from Ensembl ##
  ##    	        ##
  ##--------------------##
  
  ## Do the querying in chunks
  
  num.chunks <- ceiling(length(snplist$Name)/chunk.size)
  
  for (i in 1:num.chunks) {
    start.index <- (i-1)*chunk.size + 1
    end.index <- min(length(snplist$Name), i*chunk.size)
    cat(paste("Fetching from Ensembl ",
              as.character(round(end.index/length(snplist$Name) * 100),0),
              "%\r", sep=""))
    snplist.subset.names <- snplist$Name[start.index:end.index]
##rr  a1 <- getBM(attributes=c("refsnp_id","chr_name","chrom_start"),
##rr              filters=c("snp_filter"),values=snplist.subset.names,
##rr              mart=snpmart)
    a1 <- getBM(attributes=c("refsnp_id","chr_name","chrom_start"),
                filters=c("snp_filter", "chr_name"),
                values=list(snplist.subset.names, c(1:22, "X")),
                mart=snpmart)

    cat("Fetching from Ensembl ")
    if (gmi.bioc.debug)
      browser()
    if (! is.null(a1) && nrow(a1) > 0) {
      names(a1) <- attr
      a1$Ensembl.chr[which(a1$Ensembl.chr == "" | is.na(a1$Ensembl.chr))] <- "U"
      if (length(which(duplicated(a1$Name))) > 0) {
        ## safe for merge as both dr1 and dr2 are non empty
        if (length(which(duplicated(a1$Name))) > 0) {
          dr2 <- a1[which(duplicated(a1$Name)), ]          
          a1 <- a1[-which(duplicated(a1$Name)), ]
          dr1 <- a1[which(a1$Name %in% dr2$Name), ]
          a1[which(a1$Name %in% dr2$Name), 2] <- "D"
          names(dr2) <- c("Name", "Ensembl.chr.dupl", "Ensembl.p.dupl")
          dup <- rbind(dup, merge(dr1, dr2, by="Name"))
        }
      }
      ## replace the names in a1 with names in snplist if they only differ by case
      ## This should maintain rsIDs with lower-case rs since snplist already had
      ## lowercase rs IDs
      ## b[which(tolower(b) %in% tolower(a))] <- a[which(tolower(a) %in% tolower(b))]
      a1 <- a1[order(tolower(a1$Name)),]
      snplist.subset.names <- snplist[order(tolower(snplist.subset.names)),1]
      
      a1$Name[which(tolower(a1$Name) %in% tolower(snplist.subset.names))]  <-
        snplist.subset.names[which(tolower(snplist.subset.names) %in% tolower(a1$Name))]
      a[(num.returned + 1):(num.returned + dim(a1)[1]), ] <- a1          
      num.returned <- num.returned + dim(a1)[1]
    }
    cat(paste(as.character(round(end.index/length(snplist$Name) * 100),0),
              "%\r", sep=""))
  }
  
  cat(paste("Total number of SNPs queried: ",
            as.character(length(snplist$Name)), ".\n", sep=""))
  cat(paste("Total number of SNPs found: ",
            as.character(num.returned), ".\n", sep=""))

  ## Search SNPs not placed the first time for aliases inside NCBI

  ## names that can be searched in Entrez
  if (length(which(snplist$Name == "rs")) > 0 ||
      length(which(sub("rs", "", snplist$Name) <= 0)) > 0) {
      
    mis.a <- snplist[snplist$Name == "rs" ||
                     sub("rs", "", snplist$Name) <= 0, ]
  }
  else {
    mis.a <- NULL
  }
  
  if (length(which(! snplist$Name %in% a$Name &
                 substr(snplist$Name,1,2) == "rs")) > 0) {
    mis <- snplist[which(! snplist$Name %in% a$Name &
                 substr(snplist$Name,1,2) == "rs"), ]
    if (length(snplist) == 3) {
      names(mis) <- names(snplist)
    }
    else {
      mis <-  data.frame(matrix(unlist(mis), ncol=1, byrow=TRUE),
                         row.names = NULL, stringsAsFactors = F)
      
      names(mis) <- "Name"
    }
    searchID <- sub("rs", "", mis$Name)
  }
  else {
    mis <- NULL
    searchID <- NULL
  }

  ## names that cannot be searched in Entrez,
  ## will be missing throughout

  if (length(which(! snplist$Name %in% a$Name &
                   ! substr(snplist$Name,1,2) == "rs")) > 0) {
    mis.a <- rbind(mis.a,
                   snplist[! snplist$Name %in% a$Name &
                           ! substr(snplist$Name,1,2) == "rs",])

    if (length(snplist) == 1) {
      mis.a <-  data.frame(matrix(unlist(mis.a), ncol=1, byrow=TRUE),
                           row.names = NULL, stringsAsFactors = F)
      names(mis.a) <- "Name"
    }
    else {
      names(mis.a) <- names(snplist)
    }
  }
  else {
    mis.a <- NULL
  }

  num.mis <- length(searchID)

  if (num.mis == 0) {
    if (num.returned > 0) {
      a <- a[1:num.returned, ]
    }
    else {
      a <- NULL
    }
    if (!is.null(dup)) {
      if (length(which(duplicated(dup$Name))) > 0) {
        dup <- dup[-which(duplicated(dup$Name)), ]
      }
    }

    snplists$snplist <- snplist
    snplists$a <-  a
    snplists$dup <- dup
    ## Every SNP not found could not be searched in NCBI either
    snplists$mis <- mis.a
    snplists$alias.list <- NULL
    snplists$more.aliases <- NULL
    
    return(snplists)
  }

  ## Weed out SNPs named rs, rs0, rs00, rs000 etc.
  b.snplist <-  NULL
  searchID <- as.character(searchID[!is.na(searchID) & searchID > 0])

  search.rsID <- paste("rs", searchID, sep="")

  ##   write.table(file="saved.notfound1", mis, quote=F)
  more.aliases <- NULL
  
  ncbi.chunk.size <- 100
  ## alias.notfound is a list of searched IDs which either did not
  ## return a new name, or did not return anything at all
  
  alias.notfound <- NULL

  ## As per NCBI guidelines, chunk these out into 100 SNP lists
  if (num.mis > 0) {
    ncbi.num.chunks <- ceiling(num.mis/ncbi.chunk.size)
  
    reverse.aliases <- NULL
    ## See if NCBI ever returns a new name but no reverse alias
    new.no.reverse <- NULL
    for (j in 1:ncbi.num.chunks) {
      ncbi.start.index <- (j-1)*ncbi.chunk.size + 1
      ncbi.end.index <- min(num.mis, j*ncbi.chunk.size)
      ## print(c(ncbi.start.index, ncbi.end.index))
      if (length(snplist) == 3) {
        oldID <- as.vector(searchID[ncbi.start.index : ncbi.end.index])
      }
      else  {
        oldID <- searchID[ncbi.start.index : ncbi.end.index]
      }
      old.rsID <- paste("rs", oldID, sep="")
      idlist <- paste(oldID,collapse=",")
      cat(paste("Resolving aliases ",
                as.character(round(ncbi.end.index/num.mis * 100),0),
                "%\r", sep=""))
      URL <-
        paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?",
              "tool=GMI&email=", request.from, "&db=snp&report=XML&id=",
              idlist,"&retmode=xml",sep="")
      
      z <- xmlTreeParse(URL,useInternalNode=T)
      ns <- getDefaultNamespace(xmlRoot(z))
      namespaces <- c(ns=ns)

      ## Have to make sure that b.snplist is the same length as reverse.aliases
      ## otherwise we get weird results
      ## So, if a reverse alias is not returned, then fill in the old name
      ## There may be two possibilities:
      ## 1) merge history is empty
      ## 2) merge history is not empty, but does not contain the searched name(!)
        
      if (length(namespaces) > 0) {
        els <- getNodeSet(z,"//ns:Rs",namespaces)      
        rsID <- sapply(els,xmlGetAttr,"rsId")

        ## Found some merged SNP names
        if (length(which(oldID %in% rsID)) > 0) {

          ## NCBI returned aliases that was the same as the oldID
          ## no point in searching these SNPs again
          alias.notfound <- append(alias.notfound,
                                   paste("rs", oldID[which(oldID %in% rsID)], sep=""))
          ## cat(paste("Returned same alias as query name ",
          ##               as.character(length(which(rsID %in% oldID))), "\n",
          ##               sep=""))
          
          old.subset.rsID <- old.rsID[which(!oldID %in% rsID)]
        }
        else {
          old.subset.rsID <- old.rsID
        }
        if (length(which(!rsID %in% oldID)) > 0) {
          rsID <- rsID[which(!rsID %in% oldID)]
        }
        else {
          rsID <- NULL
        }
        
        num.reverse <- 0
        for (rs in rsID) {
          ## Create a combined list of all merge histories
          
          ## Look for reverse aliases
          qu <- paste("//ns:ExchangeSet/ns:Rs[@rsId=", rs, "]/ns:MergeHistory",
                      sep="")
          els = getNodeSet(z, qu, namespaces)
          
          if (length(els) > 0) {
            
            ## Each merged SNP name may have more than one reverse alias.
            ## Now find which old snp name is in the reverse alias list
            
            rv <- sapply(els, xmlGetAttr, "rsId")
            if (length(which(rv %in% oldID)) > 0) {
              reverse.aliases <-
                append(reverse.aliases, paste("rs", rv[which(rv %in% oldID)], sep=""))
              num.reverse <- num.reverse + 1
            }
            else {
              ## NCBI returned the a new SNP name, but no reverse alias(does this happen?)
              new.no.reverse <-
                append(new.no.reverse, paste("rs", rs, sep=""))              
            }
            ##  print(c(rs, rv[which(rv %in% oldID)]))
          }
          else {
            ## Some SNPs do not have merge history in Entrez
            reverse.aliases <- append(new.no.reverse, paste("rs", rs, sep=""))
          }
        }

        ##           cat(paste("Merged SNPs ", as.character(length(rsID)),
        ##                     ", reverse aliases ", as.character(num.reverse), "\n",
        ##                     sep=""))

        ## All new merged SNPs appear to return reverse aliases
        ## Therefore it is reasonable to consider those old ids not in reverse aliases
        ## to be missing
        
        b.snplist <- append(b.snplist, paste("rs", rsID, sep=""))
        
        if (length(which(!old.subset.rsID %in% reverse.aliases)) > 0) {
          alias.notfound <- append(alias.notfound,
                                   old.subset.rsID[which(!old.subset.rsID %in% reverse.aliases)])

          ## cat(paste("Queried ", as.character(length(old.rsID)),
##                     "Entrez returned ", as.character(length(rsID)),
##                     " added to missing ",
##                     as.character(length(which(!old.subset.rsID %in% reverse.aliases))),
##                     "\n",
##                     sep=""))
          
        }
      }
      else {
        ## Not a single searched name was found 
        alias.notfound <- append(alias.notfound, old.subset.rsID)
      }
    }
    cat(paste("2nd round query on ", as.character(length(b.snplist)), " SNPs\n",
              sep=""))
    
    if (!is.null(reverse.aliases)) {
      alias.list <- data.frame(t(rbind(reverse.aliases, b.snplist)),
                               stringsAsFactors=F)
      ## Now remove those whose new and old IDs are the same
      names(alias.list) <- c("Old", "New")
      
      alias.list <- alias.list[alias.list$Old != alias.list$New, ]
      ## alias.no.reverses <- alias.list[alias.list$Old == alias.list$New]

      if (length(which(duplicated(alias.list))) > 0) {        
        alias.list <- alias.list[-which(duplicated(alias.list)), ]
      }
      
##       if (length(which(duplicated(alias.no.reverses))) > 0) {        
##         alias.no.reverses <- alias.no.reverses[-which(duplicated(alias.no.reverses)), ]
##       }
    }
    else {
      alias.list <- NULL
##      alias.no.reverses <-  NULL
    }

    ## Now also create the mis.aliases matrix
    ## from those for which we did not find a new merged name
    ## in Entrez. or for some reason the record returned by ENtrez does not
    ## have our searched SNP

    if (!is.null(alias.notfound)) {
      mis.aliases <- mis[which(mis$Name %in% alias.notfound),]
      if (length(mis) == 3) {
        names(mis.aliases) <- names(mis)
      }
      else {
        mis.aliases <-  data.frame(matrix(mis.aliases,
                                          nrow=length(mis.aliases), ncol=1))
        names(mis.aliases) <- "Name"
      }
    }
    else {
      mis.aliases <- NULL
    }

    ## Now query Ensembl with new SNP names as previously
    ## SNPs not found will be stored in mis.b
    cat("\n")

    if (!is.null(b.snplist)) {
      b <- data.frame(matrix(nrow=length(b.snplist), ncol=length(attr)))
      names(b) <- attr
      b.dup <- NULL
      num.returned.b <- 0
      num.chunks <- ceiling(length(b.snplist)/chunk.size)
      for (i in 1:num.chunks) {
        start.index <- (i-1)*chunk.size + 1
        end.index <- min(length(b.snplist), i*chunk.size)
        cat("Fetching from Ensembl \r")
        b1 <- getBM(attributes=c("refsnp_id","chr_name","chrom_start"),
                    filters=c("snp_filter"),values=b.snplist[start.index:end.index],
                    mart=snpmart)
        
        if (! is.null(b1) && nrow(b1) > 0) {
          names(b1) <- attr
          b1$Ensembl.chr[which(b1$Ensembl.chr == "" | is.na(b1$Ensembl.chr))] <- "U"
          if (length(which(duplicated(b1$Name))) > 0) {
            ## OK for merge
            dr2 <- b1[which(duplicated(b1$Name)), ]
            b1 <- b1[-which(duplicated(b1$Name)), ]
            dr1 <- b1[which(b1$Name %in% dr2$Name), ]
            b1[which(b1$Name %in% dr2$Name), 2] <- "D"
            names(dr2) <- c("Name", "Ensembl.chr.dupl", "Ensembl.p.dupl")
            b.dup <- rbind(b.dup, merge(dr1, dr2, by="Name"))
          }
          b[(num.returned.b + 1):(num.returned.b + dim(b1)[1]), ] <- b1 
          num.returned.b <- num.returned.b + dim(b1)[1]
          
        }
        cat(paste(as.character(round(end.index/length(b.snplist) * 100),0), "%\r", sep=""))
      }
      
      ## Look for merged SNPs not found in Ensembl
      ## Create a list mis.b with old or new names as necessary.
      
      if (num.returned.b > 0) {
        b <- b[1:num.returned.b, ]
        ## SNPs missing from the 2nd search
        if (length(which(!b.snplist %in% b$Name)) > 0) {
          b.mis.names <- b.snplist[which(!b.snplist %in% b$Name)]
        }
        else {
          b.mis.names <- NULL
        }        
      }
      else {
        ## None of the merged SNPs were found, so all queried SNPs are missing
        b <- NULL
        b.mis.names <- b.snplist
      }
      
      if (!is.null(b.mis.names)) {
        ## Create a mis.b with new names and old user info if they exist
        ## if old names are required, create one or more rows per alias

        mis.new <-  b.mis.names[which(b.mis.names %in% alias.list$New)]
        mis.b <- mis[which(mis$Name %in% alias.list$Old), ]
        ## If old names are output. no other modification is necessary
        
        if (use.aliases) {
          for (new in mis.new) {
            mis.indexes <- which(alias.list$Old[alias.list$New == new])
            if (length(mis.indexes) > 1) {
              ## We have more than one set of user info mapping to one SNP
              ## Mark these ambiguous
              mis.b$Name[mis.indexes[1]] <- new
              if (length(mis) == 3) {
                ## set the chr and pos to unknowns
                mis.b[old.indexes, 2:3] <- c("U", "NA")
              }
              ## Now delete the other rows 
              mis.b <- mis.b[-mis.indexes[2:length(mis.indexes)], ]
            }
          }
        }
        
      }
      else {
        mis.b <- NULL
      }
      
      cat(paste("After 2nd query, total number of SNPs found: ",
                as.character(num.returned + num.returned.b), ".\n", sep=""))

      ## Replace new names inside b with old if necessary

      if (length(alias.list) == 0) {
        alias.list <- NULL
      }
      else {        
        ## Replace merged names with their input names in b if newsnps option not set
        ## This may end up adding rows, if two or more SNPs merged
        
        add.rows <- NULL

        for (new in unique(alias.list$New)) {
          ## Get the list of old names for each new name
          old <- alias.list$Old[alias.list$New == new]            
          if (!is.null(dim(b)) && dim(b)[1] > 0) {          
            b.indexes <- which(b$Name == new)
          }
          else {
            b.indexes <- integer(0)
          }           
          if (length(b.indexes) > 0) {
            if (use.aliases) {
              ## Do not replace the merged SNP names in b and with old names
              ## If a already has the merged SNP name, remove this from B,
              ## otherwise, do nothing, b only has a single incidence of a
                            
              if (length(which(new %in% a$Name)) > 0) {
                ## remove all occurrences of b$Name
                more.aliases <- rbind(more.aliases,
                                      alias.list$Old[alias.list$New %in% b$Name[b.indexes]])
                b <- b[-b.indexes, ]
              }
            }
            else {
              copy.b.record <- b[b$Name == new, 2:3]
              b$Name[b$Name == new] <- old[1]
              if (length(old) > 1) {
                for (i in 2:length(old)) {                  
                  ## Add input SNP names
                  add.rows <- rbind(add.rows, (cbind(old[i], copy.b.record)))
                }
              }
            }
          }
        }
          
        if (!is.null(add.rows)) {
          names(add.rows) <- names(b)
          b <- rbind(b, add.rows)
        }
        
        if (!is.null(dim(b)) && dim(b)[1] > 0) {
          num.returned.b <- dim(b)[1]
          a[(num.returned + 1) : (num.returned + num.returned.b), ] <-  b
          num.returned <- num.returned + num.returned.b
        }
        
        ## End of replacing aliases in b

        ## Now replace SNPs inside snplist if necessary

        if (use.aliases) {
          for (new in unique(alias.list$New)) {
            ## Get the list of old names for each new name
            old <- alias.list$Old[alias.list$New == new]
            snplist.indexes <-  which(snplist$Name %in% c(new, old))

            if (length(snplist.indexes) > 1) {
              ## replace, and discard
              if (length(snplist) == 3) {
                snplist <- snplist[-snplist.indexes[2:length(snplist.indexes)], ]
              }
              else {
                snplist <-
                  data.frame(matrix(snplist[-snplist.indexes[2:length(snplist.indexes)],],
                                    ncol=1, byrow=T),
                             row.names = NULL, stringsAsFactors = F)
                names(snplist) <- "Name"
              }
            }
            else if (length(snplist.indexes) == 1) {
              ## Simply replace
              snplist$Name[snplist.indexes] <- new
            }
          }
        }
      }  
    }

    dup <- rbind(dup, b.dup)
    if (num.returned > 0) {
      a <- a[1:num.returned, ]
    }
    else {
      a <- integer(0)
    }
  }
  else {
    cat("Found 0 SNps after searching Entrez SNPs!\n");
    a <- integer(0)
    alias.list <- NULL
  }


  if (length(which(duplicated(dup$Name))) > 0) {
    dup <- dup[-which(duplicated(dup$Name)), ]
  }

  mis <-  data.frame(rbind(mis.a,  mis.aliases, mis.b),
                     row.names=NULL, stringsAsFactors = F)
  if (length(mis) > 0) {
    names(mis) <-  names(snplist)
    if (length(which(duplicated(mis$Name))) > 0) {
      mis <- mis[-which(duplicated(mis$Name)), ]
    }
  }
  else {
    mis <- NULL
  }

  snplists$snplist <- snplist
  snplists$a <- a
  snplists$dup <- dup
  snplists$alias.list <- alias.list
  snplists$mis <- mis
  snplists$more.aliases <- more.aliases
  
  return(snplists)
  
}

log.query.problems <- function(snplists, output.warnings, use.aliases)

{
  aliases <- snplists$alias.list
  mis <- snplists$mis
  dup <- snplists$dup
  more.aliases <- snplists$more.aliases
  
  write("  WARNINGS FROM INTERPOLATION", file=output.warnings, append=T)
  line.dash(file=output.warnings)
  
  if (!is.null(aliases)) {
    write("   Names of the following SNPs have merged:",
          file=output.warnings, append=T)
    line.dash(file=output.warnings)
    write.table(aliases, file=output.warnings, append=T, quote=F, row.names=F, sep="\t")
    line.dash(file=output.warnings)
    line.blank(file=output.warnings)
  }

  ##    print(more.aliases)
  if (!is.null(more.aliases)) {
    line.dash(file=output.warnings)
    write("   The following SNPs excluded from output because", 
          file=output.warnings, append=T)
    write("   they have merged with other SNPs also present in input:",
          file=output.warnings, append=T)
    line.dash(file=output.warnings)
    write.table(more.aliases, file=output.warnings, append=T, quote=F,
                col.names=F, row.names=F, sep="\t")
  }

##   if (!is.null(really.bad.snps)) {
##     line.dash(file=output.warnings)
##     write("   The following SNPs excluded from output because", 
##           file=output.warnings, append=T)
##     write("   their names have merged and have duplicate Ensembl entries:",
##           file=output.warnings, append=T)
##     line.dash(file=output.warnings)
##     write.table(really.bad.snps, file=output.warnings, append=T, quote=F,
##                 col.names=F, row.names=F, sep="\t")
##   }


  ## SNPs not found in Ensembl

  if(!is.null(dim(mis)) && (dim(mis)[1] > 0)){
    line.dash(file=output.warnings)
    write("   SNPs not found in Ensembl:", file=output.warnings, append=T)
    line.dash(file=output.warnings)
    if (ncol(mis) == 3) {
      write.table(mis, file=output.warnings, append=T, quote=F,
                  row.names=F, sep="\t")
    }
    else {
      write.table(mis, file=output.warnings, append=T, quote=F,
                  col.names = F, row.names=F, sep="\t")
    }
      
    line.dash(file=output.warnings)
    line.blank(file=output.warnings)
  }

  ## Bad SNPs with multiple records in Ensembl
    
  if(length(dup$Name) > 0){
    ## Assign NAs to physical positions
    chrstr <- paste(dup$Ensembl.chr, dup$Ensembl.chr.dupl, sep="")
    dupxy <- dup[chrstr == "XY" | chrstr == "YX", ]
    dupauto <- dup[chrstr != "XY" & chrstr != "YX", ]
    
    if (length(dupauto$Name) > 0) {
      line.dash(file=output.warnings)
      write("    Duplicate SNPs Found in Ensembl:", file=output.warnings, append=T)
      line.dash(file=output.warnings)
      line.blank(file=output.warnings)
      write.table(dupauto, file=output.warnings, append=T, quote=F, row.names=F, sep="\t")
    }
    
    if (length(dupxy$Name) > 0) {
      line.dash(file=output.warnings)
      write("   SNPS mapped to both X and Y:", file=output.warnings, append=T)
      line.dash(file=output.warnings)
      write.table(dupxy, file=output.warnings, append=T, quote=F, row.names=F, sep="\t")
      line.dash(file=output.warnings)
      line.blank(file=output.warnings)
    }

  }
  return
}

gmi.bioc.read.rutgers.map = function(mpath, file) {
  map <- read.table(paste(mpath, file, sep=""),
                    header = T, stringsAsFactors=F)
  map <- map[!is.na(map$Ensembl_map_physical_position),]
  map <- map[duplicated(map$Ensembl_map_physical_position) == FALSE,]
  return (map);

}

gmi.bioc.interpolate = function(map, bppos, sel.gen) {

  ## Mapping region
  int <- which(bppos >= min(map$Ensembl_map_physical_position) & bppos <= max(map$Ensembl_map_physical_position))

  ## Extrapolation region
  min.ext <- which(bppos < min(map$Ensembl_map_physical_position))
  max.ext <- which(bppos > max(map$Ensembl_map_physical_position))

  map.n <- length(map$Markers_name)
  ans = NULL

  ## Linear interpolation for internal points
  for (sel in sel.gen) {
    map.gen  <- rep(NA, length(bppos))

    if (!is.na(sel)) {
      gen = map[, sel]
      a.fun <- approxfun(map$Ensembl_map_physical_position, gen)

      ## Ensembl records mapping
      map.gen[int] <- round(a.fun(bppos[int]), 4)

      ## Use first slope (2 snps)
      lm.min <- lm(gen[1:2]~map$Ensembl_map_physical_position[1:2])
      map.gen[min.ext] <- round(lm.min$coefficient[1] + lm.min$coefficient[2] * bppos[min.ext], 4)

      ## Use last slope (2 snps)
      lm.max <- lm(gen[(map.n-1):map.n]~map$Ensembl_map_physical_position[(map.n-1):map.n])
      map.gen[max.ext] <- round(lm.max$coefficient[1] + lm.max$coefficient[2] * bppos[max.ext], 4)

      ## map.gen[snp.e[[i]]$Name %in% dup$Name] <- NA
      ## snp.e[[i]]$Ensembl.chr[snp.e[[i]]$Name %in% dup$Name] <- "U"
    }

    ans = cbind(ans, map.gen)

  }

  if (gmi.bioc.scale.negative) {
    del = min(ans, na.rm=T)
    if (del < 0)
      ans = ans - del
  }
  
  return (cbind(ans, as.numeric(1:length(map.gen) %in% c(min.ext, max.ext) )) )

}

do.bioc.mapping <- function(snplists, output, num, mpath, create.plots, use.aliases)
{

  ## #--------------------#
  ## #                    #
  ## #    Initialize some #
  ## #       variables    #
  ## #                    #
  ## #--------------------#
  
  cat("Finished querying, now interpolating genetic positions ...\n")
  if (gmi.bioc.debug)
    browser()
  snplist <- snplists$snplist
  a <- snplists$a

  aliases <- snplists$alias.list
  more.aliases <- snplists$more.aliases
  mis <- snplists$mis

  if (!is.null(mis) && dim(mis)[2] > 1)
    mis = mis[!(mis$User.chr %in% c("Y", "MT", "U", "XY")),]
  
  ## # Define the names of output files
  output.ann <- paste(output, "_annot.txt", sep="")
  output.non.ann <- paste(output, "_nonannot.txt", sep="")
  output.pdf <- paste(output, ".pdf", sep="")
  output.warnings <- paste(output, "_log.txt", sep="")

  ## SNPs on Chromosome Y, MT or Unmapped Chrmosomes
  if (! is.null(dim(a))) {
    snp.bad <- a[a$Ensembl.chr %in% c("Y", "MT", "U", "XY"), ]
    if(length(snp.bad$Name) != 0){
      snp.bad$Ensembl.p[snp.bad$Ensembl.chr == "U"] <- NA   
      line.dash(file=output.warnings)
      write("  SNPs on XY, Y, MT or unknown chromosomes:", file=output.warnings, append=T)
      line.dash(file=output.warnings)
      write.table(snp.bad, file=output.warnings, append=T, quote=F, row.names=F, sep="\t")
      line.dash(file=output.warnings)
      line.blank(file=output.warnings)
    }
    else {
      snp.bad <- NULL
    }

    dup <- a[a$Ensembl.chr == "D", ]
    if (length(dup$Name) > 0) {
      dup$Ensembl.chr <- "U"
      dup$Ensembl.p <- NA
    }
    else {
      dup <- NULL
    }

  }
  else {
    dup <- NULL
    snp.bad <- NULL
  }

  ## Define build number
  build.e <- paste("Ensembl", num, ".p", sep="")
  build.woe <- paste("Ensembl", num, sep="")
  index <- length(snplist)
  chr.reg <- c("1","2","3","4","5","6","7","8","9","10","11","12",
               "13","14","15","16","17","18","19","20","21","22","X")

  attr <- c("Name", "Ensembl.chr", "Ensembl.p")

  plot.created <- 0
  if (create.plots > 0) {
    pdf(file=output.pdf)
  }

  if(index == 3) {
        
    ##--------------------##
    ##  Diagnostic file   ##
    ##   ("bad" SNPs)     ##
    ##      PartI	  ##
    ##--------------------##
    
    ## Clean up the ensembl SNPs
    ## SNPs on Chromosome Y or Irregular Chrmosomes (Ensemble File)

    user.bad <- snplist[snplist$User.chr %in% c("Y", "MT", "U", "XY"), ]

    if (length(user.bad$Name) > 0) {
      if (use.aliases && !is.null(aliases)) {
        aliased <- which(user.bad$Name %in% aliases$Old)
        for (ai in aliased) {
          user.bad$Name[ai] <- aliases$New[aliases$Old == user.bad$Name[ai]]
        }
      }      
      line.dash(file=output.warnings)
      write("  SNPs on XY, Y, MT or unknown chromosomes (User records):", file=output.warnings, append=T)
      line.dash(file=output.warnings)
      write.table(user.bad, file=output.warnings, append=T, quote=F, row.names=F, sep="\t")
      line.dash(file=output.warnings)
      line.blank(file=output.warnings)
    }
    else {
      user.bad <- NULL
    }

    work.u <- snplist[snplist$User.chr %in% chr.reg,]

    ## Split and order the ensembl working data

    if (! is.null(dim(a))) {
      work.e <- a[a$Ensembl.chr %in% chr.reg, ]   
      work.ord.e <- work.e[order(work.e$Ensembl.p),]
      snp.e <- split(work.ord.e, work.ord.e$Ensembl.chr)
      chr.e <- names(snp.e)
      chr.e <- chr.e[chr.e %in% chr.reg]
            
      ##-----------------------##
      ##                       ##
      ## Interpolated map file ##
      ##                       ##
      ##-----------------------##
      
      
      ## Interpolation using Ensembl records
      snp.mapped.e <- NULL
      if (gmi.bioc.debug) {
        print("snp.e"); browser()
      }
      for (i in chr.e) {
        if (i == "X") {
          ## Read in the map file for interpolation
          map = gmi.bioc.read.rutgers.map(mpath, "/chr23.comb")
          gen = gmi.bioc.interpolate(map, snp.e[[i]]$Ensembl.p, c(NA, "Smoothed_Female_map_position", NA))
          colnames(gen) = c("map.ave.e", "map.female.e", "map.male.e", "ind.e")
          snp.e[[i]] <- cbind(snp.e[[i]], gen)

          ## Plot the Ensembl mapping results
          p.min.e <- min(snp.e[[i]]$Ensembl.p, na.rm=T)
          p.max.e <- max(snp.e[[i]]$Ensembl.p, na.rm=T)
          g.min.e <- min(snp.e[[i]]$map.female.e, na.rm=T)
          g.max.e <- max(snp.e[[i]]$map.female.e, na.rm=T)
          
          if (create.plots > 0 &&
              length(which(is.finite(c(p.min.e, p.max.e, g.min.e, g.max.e)))) == 4) {
            if (p.min.e == p.max.e) {
              p.min.e <- p.min.e - 5000000;
              p.max.e <- p.max.e + 5000000;
            }
            if (g.min.e == g.max.e) {
              g.min.e <- g.min.e - 5;
              g.max.e <- g.max.e + 5;
            }
            
            ptitle.e = paste("Chromosome", i, "Genetic vs. Physical Positions", "(Ensembl Records)")
            plot(snp.e[[i]]$Ensembl.p, snp.e[[i]]$map.female.e,
                 main=ptitle.e, xlab="physical position",
                 ylab="female genetic position", xlim=c(p.min.e, p.max.e),
                 ylim=c(g.min.e, g.max.e), type="p", pch=21, col="blue")
            
            points(map$Ensembl_map_physical_position,
                   map$Smoothed_Female_map_position, pch=20, col="blue")
            rug(jitter(map$Ensembl_map_physical_position, amount=0.01), side = 1, col = "black")
            plot.created <- 1
          }
          
        } else {

          ## Read in the map file for interpolation
          map = gmi.bioc.read.rutgers.map(mpath, paste("/chr", i, ".comb", sep=""))

          gen = gmi.bioc.interpolate(map, snp.e[[i]]$Ensembl.p,
                                     sel.gen=c("Smoothed_sex.averaged_map_position",
                                               "Smoothed_Female_map_position",
                                               "Smoothed_male_map_position"))
          colnames(gen) = c("map.ave.e", "map.female.e", "map.male.e", "ind.e")
          snp.e[[i]] <- cbind(snp.e[[i]], gen)
          
          ## Plot the Ensembl mapping results
          p.min.e <- min(snp.e[[i]]$Ensembl.p, na.rm=T)
          p.max.e <- max(snp.e[[i]]$Ensembl.p, na.rm=T)
          g.min.e <- min(snp.e[[i]]$map.ave.e, snp.e[[i]]$map.female.e, snp.e[[i]]$map.male.e, na.rm=T)
          g.max.e <- max(snp.e[[i]]$map.ave.e, snp.e[[i]]$map.female.e, snp.e[[i]]$map.male.e, na.rm=T)
          
          if (create.plots > 0 &&
              length(which(is.finite(c(p.min.e, p.max.e, g.min.e, g.max.e)))) == 4) {
            if (p.min.e == p.max.e) {
              p.min.e <- p.min.e - 5000000;
              p.max.e <- p.max.e + 5000000;
            }
            if (g.min.e == g.max.e) {
              g.min.e <- g.min.e - 5;
              g.max.e <- g.max.e + 5;
            }
            
            ptitle.e = paste("Chromosome", i, "Genetic vs. Physical Positions", "(Ensembl Records)")
            
            plot(snp.e[[i]]$Ensembl.p, snp.e[[i]]$map.ave.e, main=ptitle.e,
                 xlab="physical position", ylab="genetic position",
                 xlim=c(p.min.e, p.max.e), ylim=c(g.min.e, g.max.e),
                 type="p", pch=21, col="red")
            points(snp.e[[i]]$Ensembl.p, snp.e[[i]]$map.female.e,
                   pch=21, col="blue")
            points(snp.e[[i]]$Ensembl.p, snp.e[[i]]$map.male.e,
                   pch=21, col="black")
            
            points(map$Ensembl_map_physical_position,
                   map$Smoothed_sex.averaged_map_position, pch=".", col="red")
            points(map$Ensembl_map_physical_position,
                   map$Smoothed_Female_map_position, pch=".", col="blue")
            points(map$Ensembl_map_physical_position,
                   map$Smoothed_male_map_position, pch=".", col="black")
            legend(p.min.e, g.max.e,
                   c("average(queried)", "female", "male",
                     "average(Rutgers)", "female", "male"),
                   lwd=0, col = c("red","blue","black","red","blue","black"),
                   text.col = "black",
                   cex=0.75, pch = c(21, 21, 21, 20, 20, 20), merge = TRUE)
            
            rug(jitter(map$Ensembl_map_physical_position, amount=0.01), side = 1, col = "black")
            
            
            ## Plot the ladder
            xtitle = paste("Chromosome", i, "Ladder Plot (Ensembl Records)")
            ##           print(length(snp.e[[i]]$map.ave.e))
            ##           print(length(snp.e[[i]]$map.female.e))
            ##           print(length(snp.e[[i]]$map.male.e))
            x <- data.frame(snp.e[[i]]$map.ave.e, snp.e[[i]]$map.female.e, snp.e[[i]]$map.male.e)
            names(x) <- c("Map.k.a", "Map.k.f", "Map.k.m")
            xx<-stack(x)
            with(xx, stripchart(values~ind, xlim=c(1,3), pch=19, main=xtitle, ylab="values", vertical=TRUE, col=c("red", "green", "blue")))
            if (length(x[,1]) == 1) lines(xx$value) else apply(x,1,lines)
            plot.created <- 1
          }        
        }
        
        snp.mapped.e <- rbind(snp.mapped.e, snp.e[[i]])
      }
    }
    else {
      snp.mapped.e <- NULL
      chr.e <- NULL
    }
  
    
    ## Interpolation using user records

    ## Split and order the user working data
    work.ord.u <- work.u[order(work.u$User.p),]
    snp.u <- split(work.ord.u, work.ord.u$User.chr)
    chr.u <- names(snp.u)
    chr.u <- chr.u[chr.u %in% chr.reg]
    snp.mapped.u <- NULL
    
    if (gmi.bioc.debug) {
      print("snp.u"); browser()
    }
    for (i in chr.u) {
      if (i == "X") {
        
        ## Read in the map file for interpolation
        map = gmi.bioc.read.rutgers.map(mpath, "/chr23.comb")
        gen = gmi.bioc.interpolate(map, snp.u[[i]]$User.p, c(NA, "Smoothed_Female_map_position", NA))
        colnames(gen) = c("map.ave.u", "map.female.u", "map.male.u", "ind.u")
        snp.u[[i]] <- cbind(snp.u[[i]], gen)

        ## Plot the user mapping results
        if (create.plots > 0 && dim(snp.u[[i]])[1] > 0) {
          p.min.u <- min(snp.u[[i]]$User.p, na.rm=T)
          p.max.u <- max(snp.u[[i]]$User.p, na.rm=T)
          g.min.u <- min(snp.u[[i]]$map.female.u, na.rm=T)
          g.max.u <- max(snp.u[[i]]$map.female.u, na.rm=T)
          if (length(which(is.finite(c(p.min.u, p.max.u, g.min.u, g.max.u)))) == 4) {
            if (p.min.u == p.max.u) {
              p.min.u <- p.min.u - 5000000;
              p.max.u <- p.max.u + 5000000;
            }
            if (g.min.u == g.max.u) {
              g.min.u <- g.min.u - 5;
              g.max.u <- g.max.u + 5;
            }
            
            ptitle.u = paste("Chromosome", i, "Genetic vs. Physical Positions", "(User Records)")
            plot(snp.u[[i]]$User.p, snp.u[[i]]$map.female.u,
                 main=ptitle.u, xlab="physical position",
                 ylab="female genetic position", xlim=c(p.min.u, p.max.u),
                 ylim=c(g.min.u, g.max.u), type="p", pch=23, col="blue")
            points(map$Ensembl_map_physical_position, map$Smoothed_Female_map_position, pch=20, col="blue")
            rug(jitter(map$Ensembl_map_physical_position, amount=0.01), side = 1, col = "black")
            plot.created <- 1
          }
        }              
      }
      else {
        
        ## Read in the map file for interpolation
        map = gmi.bioc.read.rutgers.map(mpath, paste("/chr", i, ".comb", sep=""))

        gen = gmi.bioc.interpolate(map, snp.u[[i]]$User.p,
                                   sel.gen=c("Smoothed_sex.averaged_map_position",
                                             "Smoothed_Female_map_position",
                                             "Smoothed_male_map_position"))
        colnames(gen) = c("map.ave.u", "map.female.u", "map.male.u", "ind.u")
        snp.u[[i]] <- cbind(snp.u[[i]], gen)

        ## Plot the user mapping results
      
        if (create.plots > 0 && dim(snp.u[[i]])[1] > 0) {

          p.min.u <- min(snp.u[[i]]$User.p, na.rm=T)
          p.max.u <- max(snp.u[[i]]$User.p, na.rm=T)
          g.min.u <- min(snp.u[[i]]$map.ave.u, snp.u[[i]]$map.female.u, snp.u[[i]]$map.male.u, na.rm=T)
          g.max.u <- max(snp.u[[i]]$map.ave.u, snp.u[[i]]$map.female.u, snp.u[[i]]$map.male.u, na.rm=T)

          if (length(which(is.finite(c(p.min.u, p.max.u, g.min.u, g.max.u)))) == 4) {
            if (p.min.u == p.max.u) {
              p.min.u <- p.min.u - 5000000;
              p.max.u <- p.max.u + 5000000;
            }
            if (g.min.u == g.max.u) {
              g.min.u <- g.min.u - 5;
              g.max.u <- g.max.u + 5;
            }
            
            ptitle.u = paste("Chromosome", i, "Genetic vs. Physical Positions", "(User Records)")
            plot(snp.u[[i]]$User.p, snp.u[[i]]$map.ave.u, main=ptitle.u,
                 xlab="physical position", ylab="genetic position",
                 xlim=c(p.min.u, p.max.u), ylim=c(g.min.u, g.max.u),
                 type="p", pch=21, col="red")
            points(snp.u[[i]]$User.p, snp.u[[i]]$map.female.u,
                   pch=21, col="blue")
            points(snp.u[[i]]$User.p, snp.u[[i]]$map.male.u,
                   pch=21, col="black")
            
            points(map$Ensembl_map_physical_position,
                   map$Smoothed_sex.averaged_map_position, pch=".", col="red")
            points(map$Ensembl_map_physical_position,
                   map$Smoothed_Female_map_position, pch=".", col="blue")
            points(map$Ensembl_map_physical_position,
                   map$Smoothed_male_map_position, pch=".", col="black")
            legend(p.min.u, g.max.u,
                   c("average(queried)", "female", "male","average(Rutgers)", "female", "male"),
                   col = c("red","blue","black","red","blue","black"),
                   cex=0.75, text.col = "black",
                   lwd=0, pch = c(21, 21, 21, 20, 20, 20), merge = TRUE)
            rug(jitter(map$Ensembl_map_physical_position, amount=0.01), side = 1,
                col = "black")
            ## Plot the ladder
            xtitle <- paste("Chromosome", i, "Ladder Plot (User Records)")
            x <- data.frame(snp.u[[i]]$map.ave.u,
                            snp.u[[i]]$map.female.u,
                            snp.u[[i]]$map.male.u)
            names(x) <- c("Map.k.a", "Map.k.f", "Map.k.m")
            xx <- stack(x)
            with(xx, stripchart(values~ind, xlim=c(1,3), pch=19, main=xtitle,
                                ylab="values", vertical=TRUE,
                                col=c("red", "green", "blue")))
            if (length(x[,1]) == 1)
              lines(xx$value)
            else
              apply(x,1,lines)

            plot.created <- 1

          }
        }  
      }
      snp.mapped.u <- rbind(snp.mapped.u, snp.u[[i]])
    }

    ## Merge ensembl interpolated file and user interpolated file
    ## 
    if (gmi.bioc.debug) {
      print("merge u&e"); browser()
    }
    snp.mapped.names = c("Name", "User.chr", "User.p",
                         "User.Map.k.a", "User.Map.k.f", "User.Map.k.m",
                         "User.Extrapolation", "Ensembl.chr", "Ensembl.p",
                         "Ensembl.Map.k.a", "Ensembl.Map.k.f", "Ensembl.Map.k.m",
                         "Extrapolation")
    if (! is.null(snp.mapped.e)) {
      ## snp.mapped.u is not NULL, as we took care of that when
      ## creating the snplist file for gmi
      snp.mapped <- merge(snp.mapped.u, snp.mapped.e, by="Name")
      snp.mapped$User.chr[is.na(snp.mapped$User.chr)] <- "U"
      snp.mapped$Ensembl.chr[is.na(snp.mapped$Ensembl.chr)] <- "U"
      snp.mapped <- snp.mapped[order(snp.mapped$User.chr, snp.mapped$User.p),]

      ## snp.mapped$Ensembl.p[snp.mapped$Name %in% dup$Name] <- NA
      ## snp.mapped[snp.mapped$Name %in% dup$Name, 10:12] <- NA

      names(snp.mapped) <- snp.mapped.names      
    }

    else {
      snp.mapped <- NULL
      chr <- NULL
    }

    
    ## Extract SNP rows where user provided known position and chromosome info
    ## If unknown chromsome, these will end up in user.bad, and get added on at
    ## the end.

    ## Need to add in map.a, map.m, map.f, extrapolation, Ensembl.p, Ensembl.chr

    if (gmi.bioc.debug) {
      print("missing"); browser()
    }
    if (!is.null(dim(mis)) && (dim(mis)[1] > 0)) {
      if (any( mis[,1] %in% snp.mapped.u$Name ))
      {

        mis <- merge(mis, snp.mapped.u, by.x=1, by.y="Name")

        ## mis looks like this if user provided known chromosomes
        ##      Name User.chr.x User.p.x User.chr.y User.p.y map.ave.u map.female.u
        ##     1 hCV16170970          1   111111          1   111111   -1.1638      -0.9975
        ##     2  rs99010163          2   222222          2   222222    0.1479       0.1858
        ##       map.male.u ind.u
        ##     1    -1.5795     1
        ##     2     0.2258     0

        snp.uf <- data.frame(mis[,1:3],
                             mis[,6:9],
                             rep("U", length(mis[,1])), rep(NA, length(mis[,1])),
                             rep(NA, length(mis[,1])), rep(NA, length(mis[,1])),
                             rep(NA, length(mis[,1])), rep(0, length(mis[,1])),
                             stringsAsFactors=F)
      }
      else {
        ## mis only has 3 columns
        snp.uf <- data.frame(mis,
                             rep(NA, length(mis[,1])), rep(NA, length(mis[,1])),
                             rep(NA, length(mis[,1])), rep(0, length(mis[,1])),
                             rep("U", length(mis[,1])), rep(NA, length(mis[,1])),
                             rep(NA, length(mis[,1])), rep(NA, length(mis[,1])),
                             rep(NA, length(mis[,1])), rep(0, length(mis[,1])),
                             stringsAsFactors=F)
      }                             
      names(snp.uf) <- snp.mapped.names
      snp.uf$User.chr[is.na(snp.uf$User.chr)] <- "U"

    }
    else {
      snp.uf <- NULL
    }

    ## snp.mapped.u looks like this (same as user.bad)
    ##         Name User.chr   User.p map.a, map.f, map.m, interp
    ## dup looks like thus
    ##         Name Ensembl.chr Ensembl.p
    ## user.dup has columns of snp.mapped.u followed by columns of dup
    ##

    if (gmi.bioc.debug) {
      print("dup"); browser()
    }
    if (!is.null(dup) && dim(dup)[1] > 0) {
      user.dup <- merge(dup, snp.mapped.u, by="Name")
      snp.dup <- data.frame(user.dup[,1], user.dup[,4:9],
                            rep("U",length(user.dup[,1])),
                            rep(NA,length(user.dup[,1])),
                            rep(NA,length(user.dup[,1])),
                            rep(NA,length(user.dup[,1])),
                            rep(NA,length(user.dup[,1])),
                            rep(0,length(user.dup[,1])),
                            stringsAsFactors=F);
      names(snp.dup) <- snp.mapped.names
    }
    else {
      snp.dup <- NULL
    }

    if (gmi.bioc.debug) {
      print("user.bad"); browser()
    }
    if (!is.null(user.bad)) {
      ## if Ensembl returned different chromosomes, we need to take that
      ## into account
      if (!is.null(snp.mapped.e)) {
        user.bad <- merge(user.bad, snp.mapped.e, all.x=T, by="Name")
        ## This adds the columns Ensembl.chr, Ensembl.p, and 3 map columns
        ## and the extrapolation flag from snp.mapped.e.
      }
      else  {
        l = length(user.bad[,1])
        user.bad=cbind(user.bad,
                       array(rep(NA,l), dim=c(l,5)),
                       array(rep(0,l), dim=c(l,1)),
                       stringsAsFactors=F)
        names(user.bad) <- c(names(user.bad)[1:3],
                             c("Ensembl.chr",  "Ensembl.p",
                               "map.ave.e", "map.female.e",
                               "map.male.e", "ind.e"))
      }

      if (!is.null(snp.bad)) {
        ## Union of user.bad and snp.bad
        ## This adds two more columns Ensembl.p.y and Ensembl.chr.y
        ## and renames the previous columns to Ensembl.p.x and Ensembl.chr.x
        user.bad <-
          merge(user.bad, snp.bad, by="Name", all.x =T, all.y=T)
        user.bad$Ensembl.chr.x[is.na(user.bad$Ensembl.chr.x)] <-
          user.bad$Ensembl.chr.y[is.na(user.bad$Ensembl.chr.x)]
        user.bad$Ensembl.p.x[is.na(user.bad$Ensembl.p.x)] <-
          user.bad$Ensembl.p.y[is.na(user.bad$Ensembl.p.x)]
        user.bad$Ensembl.chr.x[is.na(user.bad$Ensembl.chr.x)] <- "U"
        user.bad <- user.bad[,1:9]
      }

      ## user.bad looks like this when Ensembl also found a SNP on MT
      ##     Name User.chr User.p Ensembl.chr.x Ensembl.p.x map.ave.e map.female.e
      ##     rs28358286       MT  10000             U          NA        NA           NA
      ##     map.male.e ind.e Ensembl.chr.y Ensembl.p.y
      ##       NA    NA            MT       11674
      
      ## user.bad looks like this when Ensembl did not find the MT SNP
      ## Name User.chr User.p Ensembl.chr Ensembl.p map.ave.e map.female.e
      ##     rs74586330       MT  20000        <NA>        NA        NA           NA
      ##     map.male.e ind.e
      ##             NA    NA
      
      if (gmi.bioc.debug) {
        print("user.bad ymt"); browser()
      }
      user.bad$User.chr[is.na(user.bad$User.chr)] <- "U"
      snp.ymt <- data.frame(user.bad[,1:3],
                            rep(NA,length(user.bad[,1])),
                            rep(NA,length(user.bad[,1])),
                            rep(NA,length(user.bad[,1])),
                            rep(0,length(user.bad[,1])),
                            user.bad[,4:9],
                            stringsAsFactors=F);
      
      names(snp.ymt) <- snp.mapped.names
      snp.ymt <- snp.ymt[order(snp.ymt$Ensembl.p),]
    }
    
    else {
      snp.ymt <- NULL
    }
##     print("snp mapped")
##     print(snp.mapped)
##     print("snp YMT")
##     print(snp.ymt)
##     print("snp UF")
##     print(snp.uf)
##     print("snp DUP")
##     print(snp.dup)

    if (gmi.bioc.debug) {
      print("Output"); browser()
    }
    snp.mapped.out <- rbind(snp.mapped, snp.ymt, snp.uf, snp.dup)
    snp.mapped.out <-
      snp.mapped.out[order(snp.mapped.out$User.chr, snp.mapped.out$User.p),]        
    ## Reorder and rename columns of snp.mapped.out
    snp.mapped.out <-
      snp.mapped.out[,c(2,4,1,6,5,3,7,8,10,12,11,9,13)]
    names(snp.mapped.out) <- c("Chromosome", "Map.k.a", "Name",
                               "Map.k.m", "Map.k.f", "User.p", "X.Extrapolation",
                               "X.Chromosome.e", "X.Map.k.a.e", "X.Map.k.m.e",
                               "X.Map.k.f.e", paste("X.", build.e, sep=""),
                               "X.Extrapolation.e")

    write.table(snp.mapped.out, file=output.ann, quote=F, col.names=T, row.names=F, sep="\t")

    ## Non-annotated format
    names(snp.mapped.out) <- c("Chromosome", "Kosambi", "Name", "Male",
                               "Female", "User.p", "Extrapolation",
                               "Chromosome.E", "Kosambi.e", "Male.e",
                               "Female.e", build.woe, "Extrapolation.e")
    
    snp.mapped.out$Chromosome[which(snp.mapped.out$Chromosome == "X")] <- "23"
    snp.mapped.out$Chromosome[which(snp.mapped.out$Chromosome == "U")] <- "999"
#oldsnp.mapped.out$Chromosome[which(snp.mapped.out$Chromosome == "Y")] <- "25"
#oldsnp.mapped.out$Chromosome[which(snp.mapped.out$Chromosome == "XY")] <- "24"
    snp.mapped.out$Chromosome[which(snp.mapped.out$Chromosome == "Y")] <- "24"
    snp.mapped.out$Chromosome[which(snp.mapped.out$Chromosome == "XY")] <- "25"
    snp.mapped.out$Chromosome[which(snp.mapped.out$Chromosome == "MT")] <- "26"
    write.table(snp.mapped.out, file=output.non.ann, quote=F, col.names=T,
                row.names=F, sep="\t", na="-99.00")

    
    ##--------------------##
    ##                    ##
    ##  Diagnostic file   ##
    ##      PartII	  ##
    ##--------------------##
    
    chr.m <- c(chr.u, chr.e)
    chr.m <- chr.m[!duplicated(chr.m)]
    
    if (! is.null(snp.mapped)) {
      snp.mapped.p <- snp.mapped[which(snp.mapped$User.chr!="U" & snp.mapped$Ensembl.chr!="U"),]
      if (length(which(! is.na(snp.mapped.p$Ensembl.p))) > 0) {
        diff.p <- snp.mapped.p$Ensembl.p - snp.mapped.p$User.p
        if (create.plots > 0) {
          for (i in chr.m) {        
            chr.ind <- which(snp.mapped$User.chr==i & snp.mapped$Ensembl.chr==i)
            pdiff <- diff.p[chr.ind]
            pdif <- pdiff[!is.na(pdiff)]
            if (length(pdiff) != 0) {
              # print(pdiff)                          
              plot(pdiff, main=paste("Chromosome ", i, " (Ensembl_", num, ")", sep=""), xlab="",
                   ylab="Physical Difference", col="red")
              abline(h=0)
              plot.created <- 1
            }
          }
        }
        ##  inconsistent user and ensembl physical positions
        inc.p <- snp.mapped[which(diff.p != 0),]
        if(length(inc.p$Name) != 0){
          line.dash(file=output.warnings)
          write("   SNPs with inconsistent Ensembl and User positions:",
                file=output.warnings, append=T)
          line.dash(file=output.warnings)
          write.table(inc.p, file=output.warnings, append=T, quote=F, col.names=T, row.names=F, sep="\t")
          line.dash(file=output.warnings)
          line.blank(file=output.warnings)
          
        }
      }
      else {
        inc.p <- NULL
      }
    }
    else {
      inc.p <- NULL
    }
  
    
    ##--------------------##
    ##                    ##
    ##  Diagnostic file   ##
    ##     PartIII        ##
    ##                    ##
    ##--------------------##
    
    ## Inconsistent map order
    inc.u.a <- NULL
    inc.e.a <- NULL
    ## SNPs with extrapolated map positions
    extra.u.a <- NULL
    extra.e.a <- NULL
    
    ## Ensembl records
    for(i in chr.e){
      if(i == "X"){
        snp.e[[i]]$diff.a.e <- rep(NA, length(snp.e[[i]]$map.female.e))
        snp.e[[i]]$diff.f.e <- c(abs(snp.e[[i]][1,5]), diff(snp.e[[i]]$map.female.e, 1))
        snp.e[[i]]$diff.m.e <- rep(NA, length(snp.e[[i]]$map.female.e))
        snp.e[[i]] <- snp.e[[i]][!is.na(snp.e[[i]]$diff.f.e),]
        inc.e <- snp.e[[i]][which(snp.e[[i]]$diff.f.e<0), ]
        extra.e <- snp.e[[i]][which(snp.e[[i]]$ind.e==1), c(1:7)]
      }
      else {
        
        snp.e[[i]]$diff.a.e <- c(abs(snp.e[[i]][1,4]), diff(snp.e[[i]]$map.ave.e, 1))
        snp.e[[i]]$diff.f.e <- c(abs(snp.e[[i]][1,5]), diff(snp.e[[i]]$map.female.e, 1))
        snp.e[[i]]$diff.m.e <- c(abs(snp.e[[i]][1,6]), diff(snp.e[[i]]$map.male.e, 1))
        snp.e[[i]] <- snp.e[[i]][!is.na(snp.e[[i]]$diff.f.e),]
        inc.e <- snp.e[[i]][which(snp.e[[i]]$diff.a.e<0 | snp.e[[i]]$diff.f.e<0 | snp.e[[i]]$diff.m.e<0), ]
        extra.e <- snp.e[[i]][which(snp.e[[i]]$ind.e==1), c(1:7)]
      }

      if (dim(inc.e)[1] > 0) {
        inc.e.a <- rbind(inc.e.a, inc.e)
      }
      if (dim(extra.e)[1] > 0) {
        extra.e.a <- rbind(extra.e.a, extra.e)
      }
    }
    
    ## User records
    for(i in chr.u){
      if(i == "X"){
        
        snp.u[[i]]$diff.a.u <- rep(NA, length(snp.u[[i]]$map.female.u))
        snp.u[[i]]$diff.f.u <- c(abs(snp.u[[i]][1,5]), diff(snp.u[[i]]$map.female.u, 1))
        snp.u[[i]]$diff.m.u <- rep(NA, length(snp.u[[i]]$map.female.u))
        snp.u[[i]] <- snp.u[[i]][!is.na(snp.u[[i]]$diff.f.u),]
        inc.u <- snp.u[[i]][which(snp.u[[i]]$diff.f.u<0), ]
        extra.u <- snp.u[[i]][which(snp.u[[i]]$ind.u==1), c(1:7)]        
      }
      else {
        
        snp.u[[i]]$diff.a.u <- c(abs(snp.u[[i]][1,4]), diff(snp.u[[i]]$map.ave.u, 1))
        snp.u[[i]]$diff.f.u <- c(abs(snp.u[[i]][1,5]), diff(snp.u[[i]]$map.female.u, 1))
        snp.u[[i]]$diff.m.u <- c(abs(snp.u[[i]][1,6]), diff(snp.u[[i]]$map.male.u, 1))
        snp.u[[i]] <- snp.u[[i]][!is.na(snp.u[[i]]$diff.f.u),]
        inc.u <- snp.u[[i]][which(snp.u[[i]]$diff.a.u<0 | snp.u[[i]]$diff.f.u<0 | snp.u[[i]]$diff.m.u<0), ]        
        extra.u <- snp.u[[i]][which(snp.u[[i]]$ind.u==1), c(1:7)]
      }

      if (dim(inc.u)[1] > 0) {
        inc.u.a <- rbind(inc.u.a, inc.u)
      }
      if (dim(extra.u)[1] > 0) {
        extra.u.a <- rbind(extra.u.a, extra.u)
      }      
    }
    
    if(!is.null(inc.u.a)){
      names(inc.u.a) <- c("Name", "Ensembl.chr", "Ensembl.p",
                        "Ensembl.Map.k.a", "Ensembl.Map.k.f", "Ensembl.Map.k.m",
                        "Extrapolation", "Ensembl.diff.a", "Ensembl.diff.f",
                        "Ensembl.diff.m")

      line.dash(file=output.warnings)
      write("   SNPs with Inconsistent Map Orders (User records):",
            file=output.warnings, append=T)
      line.dash(file=output.warnings)
      write.table(inc.u.a, file=output.warnings, append=T, quote=F, row.names=F, sep="\t")
      line.dash(file=output.warnings)
      line.blank(file=output.warnings)
    }
    
    if(!is.null(inc.e.a)){
      names(inc.e.a) <- c("Name", "Ensembl.chr", "Ensembl.p",
                          "Ensembl.Map.k.a", "Ensembl.Map.k.f", "Ensembl.Map.k.m",
                          "Extrapolaltion", "Ensembl.diff.a", "Ensembl.diff.f",
                          "Ensembl.diff.m")

      line.dash(file=output.warnings)
      write("   SNPs with Inconsistent Map Orders (Ensembl records):",
            file=output.warnings, append=T)
      line.dash(file=output.warnings)
      write.table(inc.e.a, file=output.warnings, append=T, quote=F, row.names=F, sep="\t")
      line.dash(file=output.warnings)
      line.blank(file=output.warnings)
    }	

    if(!is.null(extra.e.a)){
      names(extra.e.a) <-
        c("Name", "Chromosome", "Ensembl.p", "Map.k.a", "Map.k.f", "Map.k.m", "Extrapolation")
      line.dash(file=output.warnings)
      write("   SNPs with extrapolated map positions (Ensembl records):",
            file=output.warnings, append=T)
      line.dash(file=output.warnings)
      write.table(extra.e.a, file=output.warnings, append=T, quote=F, row.names=F, sep="\t")
      line.dash(file=output.warnings)
      line.blank(file=output.warnings)
    }

    if(!is.null(extra.u.a)) {
      names(extra.u.a) <- c("Name", "Chromosome", "User.p", "User.Map.k.a", "User.Map.k.f", "User.Map.k.m", "User.Extrapolation")

      line.dash(file=output.warnings)
      write("   SNPs with extrapolated map positions (User records):", 
            file=output.warnings, append=T)
      line.dash(file=output.warnings)
      write.table(extra.u.a, file=output.warnings, append=T, quote=F, row.names=F, sep="\t")
      line.dash(file=output.warnings)
      line.blank(file=output.warnings)
    }  
  
    ##--------------------##
    ##                    ##
    ##  Warning summary   ##
    ##                    ##
    ##--------------------##

    if(!is.null(aliases) |
       length(mis$Name) != 0 | length(user.bad$Name) != 0 |
       length(snp.bad$Name) != 0 | length(dup$Name) != 0 |
       length(inc.p$Name) != 0 | length(extra.u.a$Name) != 0 |
       length(extra.e.a$Name) != 0 | length(inc.u.a$Name) != 0 |
       length(inc.e.a$Name) != 0) {
      
      line.dash(file=output.warnings)
      write("   SUMMARY OF WARNINGS FROM INTERPOLATION", file=output.warnings, append=T)
      line.dash(file=output.warnings)
      
      ## Aliases
      
      if (!is.null(aliases)) {
        write("  -> Some SNP names have changed in Ensembl.",
              file=output.warnings, append=T)
        if (use.aliases) {
          write("     Map files contain new SNP names.",
                file=output.warnings, append=T)
        }
        else {
          write("     Map files contain old SNP names.",
                file=output.warnings, append=T)
        }
        if (!is.null(more.aliases)) {
          write("     Some SNPs omitted from output as mutliple SNPs have merged.",
                file=output.warnings, append=T)
        }
      }
      
      ## SNPs not placed
      
      if(length(mis$Name) != 0)	
        write("  -> One or more SNP(s) not found in Ensembl.",
              file=output.warnings, append=T)

      ## SNPs on Chromosome Y or Irregular Chrmosomes or Unknown Chromosomes (Users records)
      if(length(user.bad$Name) != 0)
        write("  -> Found SNPs on XY, Y, MT or unknown chromosomes in user records.",
              file=output.warnings, append=T)
      
      ## SNPs on Chromosome Y or Irregular Chrmosomes (Ensembl records)
      if(length(snp.bad$Name) != 0)
        write("  -> Found SNPs on XY, Y, MT or unknown chromosomes in Ensembl records.",
              file=output.warnings, append=T)
      
      ## Bad SNPs with multiple records in Ensembl
      
      if(length(dup$Name) > 0)
        write("  -> Found SNPs with multiple records in Ensembl.",
              file=output.warnings, append=T)
      
      ## Inconsistent user and ensembl physical positions
      if(length(inc.p$Name) != 0)
        write("  -> Some SNPs have inconsistent Ensembl and user-supplied physical positions.",
              file=output.warnings, append=T)
      
      ## SNPs with extrapolated map positions (User records)
      for(i in chr.u){
        extra.u <- snp.u[[i]][which(snp.u[[i]]$ind.u==1), ]
        if(length(extra.u$Name) != 0){			
          etitle <- paste("  ->", length(extra.u$Name), "Markers on chromosome", i, "have extrapolated map positions (User records).", sep=" ")
          write(etitle, file=output.warnings, append=T)
        }
      }
      
      ## SNPs with extrapolated map positions (Ensembl records)
      for(i in chr.e){
        extra.e <- snp.e[[i]][which(snp.e[[i]]$ind.e==1), ]
        if(length(extra.e$Name) != 0){			
          etitle <- paste("  ->", length(extra.e$Name), "Markers on chromosome", i, "have extrapolated map positions (Ensembl records).", sep=" ")
          write(etitle, file=output.warnings, append=T)
        }
      }

      ## Inconsistent map order (User records)
      if(length(inc.u.a$Name) != 0)
        write("  -> Found SNPs out of map order in user records.",
              file=output.warnings, append=T)
      
      ## Inconsistent map order (Ensembl records)
      if(length(inc.e.a$Name) != 0)
        write("  -> Found SNPs out of map order in Ensembl records.",
              file=output.warnings, append=T)

    }
  }

  if(index == 1) {
    if (gmi.bioc.debug) {
      print("index == 1"); browser()
    }
    if (! is.null(dim(a))) {
      work <- a[a$Ensembl.chr %in% chr.reg,]
    
      ##Split and order the working data
      work.ord <- work[order(work$Ensembl.p),]
      snp.e <- split(work.ord, work.ord$Ensembl.chr)
      chr <- names(snp.e)
      chr <- chr[chr %in% chr.reg]
      
      ##-----------------------##
      ##                       ##
      ## Interpolated map file ##
      ##                       ##
      ##-----------------------##
      
      snp.mapped <- NULL
      
      for (i in chr) {

        if (i == "X") {
          
          ## Read in the map file for interpolation 
          map = gmi.bioc.read.rutgers.map(mpath, "/chr23.comb")
          gen = gmi.bioc.interpolate(map, snp.e[[i]]$Ensembl.p, c(NA, "Smoothed_Female_map_position", NA))
          colnames(gen) = c("map.ave", "map.female", "map.male", "ind")
          snp.e[[i]] <- cbind(snp.e[[i]], gen)

          if (create.plots > 0 && dim(snp.e[[i]])[1] > 0) {
            ## Plot the Ensembl mapping results
            p.min <- min(snp.e[[i]]$Ensembl.p, na.rm=T)
            p.max <- max(snp.e[[i]]$Ensembl.p, na.rm=T)
            g.min <- min(snp.e[[i]]$map.female, na.rm=T)
            g.max <- max(snp.e[[i]]$map.female, na.rm=T)
            if (length(which(is.finite(c(p.min, p.max, g.min, g.max)))) == 4) {
              if (p.min == p.max) {
                p.min <- p.min - 5000000;
                p.max <- p.max + 5000000;
              }
              if (g.min == g.max) {
                g.min <- g.min - 5;
                g.max <- g.max + 5;
              }
              
              ptitle = paste("Chromosome", i, "Genetic vs. Physical Positions")
              
              plot(snp.e[[i]]$Ensembl.p, snp.e[[i]]$map.female,
                   main=ptitle, xlab="physical position", ylab="female genetic position", xlim=c(p.min, p.max),
                   ylim=c(g.min, g.max), type="p", pch=23, col="blue")
              
              points(map$Ensembl_map_physical_position, map$Smoothed_Female_map_position,
                     pch=20, col="blue")
              
              rug(jitter(map$Ensembl_map_physical_position, amount=0.01), side = 1, col = "black")
              plot.created <- 1
            }
          }
          
        } else {
          
          ## Read in the map file for interpolation 
          map = gmi.bioc.read.rutgers.map(mpath, paste("/chr", i, ".comb", sep=""))
          
          gen = gmi.bioc.interpolate(map, snp.e[[i]]$Ensembl.p,
                                     sel.gen=c("Smoothed_sex.averaged_map_position",
                                               "Smoothed_Female_map_position",
                                               "Smoothed_male_map_position"))
          colnames(gen) = c("map.ave", "map.female", "map.male", "ind")
          snp.e[[i]] <- cbind(snp.e[[i]], gen)
          
          if (create.plots > 0 && dim(snp.e[[i]])[1] > 0) {
            ## Plot the Ensembl mapping results
            p.min <- min(snp.e[[i]]$Ensembl.p, na.rm=T)
            p.max <- max(snp.e[[i]]$Ensembl.p, na.rm=T)
            g.min <- min(snp.e[[i]]$map.ave, snp.e[[i]]$map.female, snp.e[[i]]$map.male, na.rm=T)
            g.max <- max(snp.e[[i]]$map.ave, snp.e[[i]]$map.female, snp.e[[i]]$map.male, na.rm=T)
            if (length(which(is.finite(c(p.min, p.max, g.min, g.max)))) == 4) {
              
              if (p.min == p.max) {
                p.min <- p.min - 5000000;
                p.max <- p.max + 5000000;
              }
              if (g.min == g.max) {
                g.min <- g.min - 5;
                g.max <- g.max + 5;
              }
              
              ptitle = paste("Chromosome", i, "Genetic vs. Physical Positions")
              plot(snp.e[[i]]$Ensembl.p, snp.e[[i]]$map.ave, main=ptitle, xlab="physical position", ylab="genetic position",
                   xlim=c(p.min, p.max), ylim=c(g.min, g.max), type="p",
                   pch=21, col="red")
              points(snp.e[[i]]$Ensembl.p, snp.e[[i]]$map.female,
                     pch=21, col="blue")
              points(snp.e[[i]]$Ensembl.p, snp.e[[i]]$map.male,
                     pch=21, col="black")                          
              
              points(map$Ensembl_map_physical_position,
                     map$Smoothed_sex.averaged_map_position, pch=".", col="red")
              points(map$Ensembl_map_physical_position,
                     map$Smoothed_Female_map_position, pch=".", col="blue")
              points(map$Ensembl_map_physical_position,
                     map$Smoothed_male_map_position, pch=".", col="black")
              
              legend(p.min, g.max,
                     c("average(queried)", "female", "male",
                       "average(Rutgers)", "female", "male"),
                     col = c("red","blue","black","red","blue","black"),
                     cex=0.75, text.col = "black", lwd=0,
                     pch = c(21, 21, 21, 20, 20, 20),
                     merge = TRUE)
              
              
              rug(jitter(map$Ensembl_map_physical_position, amount=0.01), side = 1, col = "black")
              
              ## Plot the ladder
              xtitle = paste("Chromosome", i, "Ladder Plot")
              x <- data.frame(snp.e[[i]]$map.ave, snp.e[[i]]$map.female, snp.e[[i]]$map.male)
              names(x) <- c("Map.k.a", "Map.k.f", "Map.k.m")
              xx<-stack(x)
              with(xx, stripchart(values~ind, xlim=c(1,3), pch=19, main=xtitle, ylab="values", vertical=TRUE, col=c("red", "green", "blue")))
              if (length(x[,1]) == 1) lines(xx$value) else apply(x,1,lines)
              plot.created <- 1
            }
          }
        }
        snp.mapped <- rbind(snp.mapped, snp.e[[i]])
      }
    }
    else {
      snp.mapped <- NULL
      chr <- NULL
    }

    if (length(mis$Name) > 0) {
      snp.uf <- data.frame(mis, rep("U", length(mis)),
                           rep(NA, length(mis)), rep(NA, length(mis)),
                           rep(NA, length(mis)), rep(NA, length(mis)),
                           rep(NA, length(mis)),
                           stringsAsFactors=F)
      names(snp.uf) <- c("Name", "Ensembl.chr", "Ensembl.p",
                         "map.ave", "map.female", "map.male", "ind")
    }
    else {
      snp.uf <- NULL
    }

    if (!is.null(dup) && dim(dup)[1] > 0) {
      snp.dup <- data.frame(dup, rep(NA, length(dup[,1])),
                            rep(NA, length(dup[,1])), rep(NA, length(dup[,1])),
                            rep(NA, length(dup[,1])),
                            stringsAsFactors=F)
      names(snp.dup) <- c("Name", "Ensembl.chr", "Ensembl.p",
                         "map.ave", "map.female", "map.male", "ind")
    }
    else {
      snp.dup <- NULL
    }
    ## create the rows for the Y and MT snps
    ## already have chr.p and physical-pos.p
    
    if (!is.null(snp.bad) && dim(snp.bad)[1] > 0) {
      snp.ymt <- data.frame(snp.bad, rep(NA, length(snp.bad[,1])),
                            rep(NA, length(snp.bad[,1])),
                            rep(NA, length(snp.bad[,1])),
                            rep(NA, length(snp.bad[,1])),
                            stringsAsFactors=F)
      snp.ymt <- snp.ymt[order(snp.ymt$Ensembl.p),]
      names(snp.ymt) <- c("Name", "Ensembl.chr", "Ensembl.p", "map.ave",
                          "map.female", "map.male", "ind")
    }
    else {
      snp.ymt <- NULL
    }

    ## snp.mapped is only SNPs that had proper chromsomes and distances
    ## So, now paste all the other tables with unmapped SNPs

    snp.mapped.out <- rbind(snp.mapped, snp.ymt, snp.uf, snp.dup)
    snp.mapped.out <- snp.mapped.out[,c(2,4,1,6,5,3,7)]
    
    ## Order the augmented table wrt chromosome, then position within chromosome
    snp.mapped.out <-
      snp.mapped.out[order(snp.mapped.out$Ensembl.chr, snp.mapped.out$Ensembl.p),]        

    names(snp.mapped.out) <-
      c("Chromosome", "Map.k.a", "Name", "Map.k.m", "Map.k.f", build.e, "X.Extrapolation")
    write.table(snp.mapped.out, file=output.ann, quote=F, col.names=T, row.names=F, sep="\t")

    ## Now for non-annotated output
    ## Change the header names
    names(snp.mapped.out) <-
      c("Chromosome", "Kosambi", "Name", "Male", "Female", build.woe, "Extrapolation")    
    snp.mapped.out$Chromosome[which(snp.mapped.out$Chromosome == "X")] <- "23"
    snp.mapped.out$Chromosome[which(snp.mapped.out$Chromosome == "U")] <- "999"
#oldsnp.mapped.out$Chromosome[which(snp.mapped.out$Chromosome == "Y")] <- "25"
#oldsnp.mapped.out$Chromosome[which(snp.mapped.out$Chromosome == "XY")] <- "24"
    snp.mapped.out$Chromosome[which(snp.mapped.out$Chromosome == "Y")] <- "24"
    snp.mapped.out$Chromosome[which(snp.mapped.out$Chromosome == "XY")] <- "25"
    snp.mapped.out$Chromosome[which(snp.mapped.out$Chromosome == "MT")] <- "26"
    
    write.table(snp.mapped.out, file=output.non.ann, quote=F, col.names=T,
                row.names=F, sep="\t", na="-99.00")
    
    ##--------------------##
    ##                    ##
    ##  Diagnostic file   ##
    ##     (Cont.)        ##
    ##                    ##
    ##--------------------##
    
    ## Inconsistent map order
    inc.a <- NULL
    ## SNPs with extrapolated map positions
    extra.a <- NULL
    for(i in chr){
      if(i == "X"){
        
        snp.e[[i]]$diff.a <- rep(NA, length(snp.e[[i]]$map.female))
        snp.e[[i]]$diff.f <- c(abs(snp.e[[i]][1,5]), diff(snp.e[[i]]$map.female, 1))
        snp.e[[i]]$diff.m <- rep(NA, length(snp.e[[i]]$map.female))
        snp.e[[i]] <- snp.e[[i]][!is.na(snp.e[[i]]$diff.f),]
        inc <- snp.e[[i]][which(snp.e[[i]]$diff.f<0),]
        extra <- snp.e[[i]][which(snp.e[[i]]$ind==1), c(1,2,3,4,5,6,7)]
      }
      else {
        
        snp.e[[i]]$diff.a <- c(abs(snp.e[[i]][1,4]), round(diff(snp.e[[i]]$map.ave, 1), 4))
        snp.e[[i]]$diff.f <- c(abs(snp.e[[i]][1,5]), round(diff(snp.e[[i]]$map.female, 1), 4))
        snp.e[[i]]$diff.m <- c(abs(snp.e[[i]][1,6]), round(diff(snp.e[[i]]$map.male, 1), 4))
        snp.e[[i]] <- snp.e[[i]][!is.na(snp.e[[i]]$diff.f),]
        inc <- snp.e[[i]][which(snp.e[[i]]$diff.a<0 | snp.e[[i]]$diff.f<0 | snp.e[[i]]$diff.m<0),]
        extra <- snp.e[[i]][which(snp.e[[i]]$ind==1), c(1,2,3,4,5,6,7)]
      }

      if (dim(inc)[1] > 0) {
        inc.a <- rbind(inc.a, inc)
      }
      if (dim(extra)[1] > 0) {        
        extra.a <- rbind(extra.a, extra)
      }
    }

    if (!is.null(inc.a)) {
      names(inc.a) <- c("Name", "Chromosome", "Ensembl.p", "Map.k.a", "Map.k.f", "Map.k.m", "Extrapolation", "diff.a", "diff.f", "diff.m")
      line.dash(file=output.warnings)
      write("   SNPS with inconsistent map orders:", file=output.warnings, append=T)
      line.dash(file=output.warnings)
      write.table(inc.a, file=output.warnings, append=T, quote=F, row.names=F, sep="\t")
      line.dash(file=output.warnings)
      line.blank(file=output.warnings)
    }

    if (!is.null(extra.a)) {
      names(extra.a) <- c("Name", "Chromosome", "Ensembl.p", "Map.k.a", "Map.k.f", "Map.k.m", "Extrapolation")      
      line.dash(file=output.warnings)
      write("   SNPs with extrapolated map positions:",
            file=output.warnings, append=T)
      line.dash(file=output.warnings)
      write.table(extra.a, file=output.warnings, append=T, quote=F, row.names=F, sep="\t")
      line.dash(file=output.warnings)
      line.blank(file=output.warnings)
    }

    ##--------------------##
    ##                    ##
    ##  Warning summary   ##
    ##                    ##
    ##--------------------##

    if(!is.null(aliases) |
       length(mis$Name) != 0 |
       length(snp.bad$Name) != 0 |
       length(dup$Name) != 0 |
       length(extra.a$Name) != 0 |
       length(inc.a$Name) != 0) {

      line.dash(file=output.warnings)
      write("   SUMMARY OF WARNINGS FROM INTERPOLATION", file=output.warnings, append=T)
      line.dash(file=output.warnings)

      ## Aliases

      if (!is.null(aliases)) {
        write("  -> Some SNP names have changed in Ensembl.",
              file=output.warnings, append=T)
        if (use.aliases) {
          write("     Map files contain new SNP names.",
                file=output.warnings, append=T)
        }
        else {
          write("     Map files contain old SNP names.",
                file=output.warnings, append=T)
        }
        if (!is.null(more.aliases)) {
          write("     Some SNPs omitted from output as mutliple SNPs have merged.",
                file=output.warnings, append=T)
        }
      }

      ## SNPs not placed
      if(length(mis$Name) > 0)
        write("  -> One or more SNP(s) not found in Ensembl.", file=output.warnings, append=T)
      
      ## SNPs on Chromosome Y or Irregular Chrmosomes
      if(length(snp.bad$Name) != 0)
        write("  -> Found SNPs on XY, Y, MT or unknown chromosomes.",
              file=output.warnings, append=T)
      
      ## Bad SNPs with multiple records in Ensembl
      if(length(dup$Name) > 0)
        write("  -> Found SNPs with multiple records in Ensembl.",
              file=output.warnings, append=T)
      
      ## SNPs with extrapolated map positions
      for(i in chr){
        extra <- snp.e[[i]][which(snp.e[[i]]$ind==1), ]
        if(length(extra$Name) != 0){			
          etitle <- paste("  ->", length(extra$Name), "Markers on chromosome", i, "have extrapolated map positions.", sep=" ")
          write(etitle, file=output.warnings, append=T)
        }
      }
      
      ## Inconsistent map order
      if(length(inc.a$Name) != 0)
        write("  -> Found SNPs out of map order in Ensembl records.",
              file=output.warnings, append=T)
      
    }
  }
  
  
  if (create.plots > 0) {
    if (plot.created == 0) {
      write("  -> No plots were created, as no mapped SNPs were found.",
            file=output.warnings, append=T)
      unlink(output.pdf)
    }
    else {
      dev.off()
      
    }
  }

  line.dash(file=output.warnings)
  
  return(TRUE)
  
}

line.dash <- function(file="")

{
  write("---------------------------------------------------",
        file=file, append=T)
  return
}

line.blank <- function(file="")

{
  write(" ", file=file, append=T)
  return
}

