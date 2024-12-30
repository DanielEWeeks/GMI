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
gmi.bioc.setdata <- function(mpath, rpath) {

  vfile <- paste(mpath, "/version.txt", sep="")
  
  unlink(vfile)
  
  library(biomaRt)
  snpmart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")
  
  ## Obtain the current Ensembl version
  
  lmrts <- listMarts()
  ever <- as.character(lmrts[lmrts[,1] == "ENSEMBL_MART_SNP", 2])
  
  ver <- strsplit(ever, "[ ]+")[[1]][3]
  
  write(ver, vfile)
  
  library(RMySQL)
  
  bcon <- dbConnect(MySQL(), user="anonymous", host="ensembldb.ensembl.org", port=5306)
  
  qw <- paste("SHOW DATABASES LIKE 'homo_sapiens_core_", ver, "_%'", sep="")
  
  q <- dbSendQuery(bcon, qw)
  
  qdb <- fetch(q)[1,1]
  
  con <- dbConnect(MySQL(), user="anonymous", host="ensembldb.ensembl.org", dbname=qdb, port=5306)
  
  thres <- 200
  
  chr.reg <- c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
  ## Create a directory for local database
  
  merpath <- paste(mpath, "/rutgers_merged_", ver, sep="")
  
  dir.create(merpath, showWarnings = FALSE)
  
  ## Change to directory for local database
  
  setwd(merpath)
  
  ## Process with regular chromosomes (1-22)
  
  cat("Creating merged maps for SNPs on all chromosomes ...")
  
  
  for(i in chr.reg){
    
    output.db <- paste("chr", i, ".comb", sep="")
    output.war <- paste("chr", i, ".warnings", sep="")
    
    unlink(output.db)
    unlink(output.war)
    
    ##----------------------------##
    ##                            ##
    ##    Read in Rugters Data    ##
    ##                            ##
    ##----------------------------##
    
    f <- paste(rpath, "/chr", i, ".sm.map2", sep="")
    
    if (file.exists(f) == TRUE) {
      
      r <- read.table(f, header = T)
      
    } else {
      
      warning(paste("File", f, "does not exist!", sep=" "))
      return(-1)
      
    }
    
    ##-----------------------------##
    ##                             ##
    ##    Query using Rutgers id   ##
    ##                             ##
    ##-----------------------------##
    
    e <- getBM(attributes=c("refsnp_id","chr_name","chrom_start"),filters=c("snp_filter"),values=r$Markers_name,mart=snpmart)
    
    names(e) <- c("name", "chr", "Ensembl_map_physical_position")
    
    dup.e <- subset(e, e$name %in% e$name[duplicated(e$name)])
    
    work.e <- e[!e$name %in% dup.e$name,]
    
    ##-----------------------------##
    ##                             ##
    ##  Query for microsatellites  ##
    ##                             ##
    ##-----------------------------##
    
    micro <- r$Markers_name[!r$Markers_name %in% work.e$name]
    
    n1 <- length(micro)
    
    q1 <- NULL
    
    for (j in 1:n1) {
      
      if(j == 1) q1 <- paste("select marker_id, name from marker_synonym where name in (", micro[j], ",", sep="'")
      
      if(j > 1 & j < n1) q1 <- paste(q1, micro[j], ",", sep="'")
      
      if(j == n1) q1 <- paste(q1, micro[j], ")", sep="'")
      
    }
    
    m.id <- dbGetQuery(con,q1)
    
    n2 <- length(m.id$marker_id)
    
    q2 <- NULL
    
    for(j in 1:n2) {
      
      if(j == 1) q2 <- paste("select marker_id, seq_region_id, seq_region_start, seq_region_end from marker_feature where marker_id in (", m.id$marker_id[j], ",", sep="")
      
      if(j > 1 & j < n2) q2 <- paste(q2, m.id$marker_id[j], ",", sep="")
      
      if(j == n2) q2 <- paste(q2, m.id$marker_id[j], ")", sep="")
      
    }
    
    m.ft <- dbGetQuery(con,q2)
    
    n3 <- length(m.ft$seq_region_id)
    
    q3 <- NULL
    
    for(j in 1:n3) {
      
      if(j == 1) q3 <- paste("select seq_region_id, name from seq_region where seq_region_id in (", m.ft$seq_region_id[j], ",", sep="")
      
      if(j > 1 & j < n3) q3 <- paste(q3, m.ft$seq_region_id[j], ",", sep="")
      
      if(j == n3) q3 <- paste(q3, m.ft$seq_region_id[j], ")", sep="")
      
    }
    
    m.chr <- dbGetQuery(con,q3)
    
    names(m.chr) <- c("seq_region_id", "chr")
    
    m1 <- merge(m.ft, m.chr, by.x="seq_region_id", by.y="seq_region_id", all.x=T)
    
    work.m1 <- m1[m1$chr == i,]
    
    dup.id <- subset(work.m1, work.m1$marker_id %in% work.m1$marker_id[duplicated(work.m1$marker_id)])
    dup.id.ord <- dup.id[order(dup.id$marker_id), ]
    
    work.m2 <- work.m1[!work.m1$marker_id %in% dup.id.ord$marker_id,]
    
    m2 <- merge(m.id, work.m2, by.x="marker_id", by.y="marker_id", all.x=T)
    
    work.m3 <- m2[!is.na(m2$seq_region_id),c(2,6,4)]
    
    names(work.m3) <- c("name", "chr", "Ensembl_map_physical_position")
    
    dup.name <- subset(work.m3, work.m3$name %in% work.m3$name[duplicated(work.m3$name)])
    dup.name.ord <- dup.name[order(dup.name$name), ]
    
    work.m4 <- work.m3[!work.m3$name %in% dup.name.ord$name,]
    
    ##---------------------------------------##
    ##                                       ##
    ## Dealing with Micro with multiple hits ##
    ##                                       ##
    ##---------------------------------------##
    
    uniq <- dup.name.ord[!duplicated(dup.name.ord$name),"name"]
    
    work.m5 <- NULL
    
    for(k in 1:length(uniq)) {
      
      rec <- dup.name.ord[dup.name.ord$name==uniq[k],]
      
      if(length(rec$Ensembl_map_physical_position) > 1) pdiff <- max(rec$Ensembl_map_physical_position) - min(rec$Ensembl_map_physical_position) else pdiff <- 0
      
      if(pdiff <= thres) work.m5 <- rbind(work.m5, rec[rec$Ensembl_map_physical_position==min(rec$Ensembl_map_physical_position),])
      
    }
    
    work.m6 <- work.m5[!duplicated(work.m5$name),]
    
    work.m <- rbind(work.m4, work.m6)
    
    dup <- rbind(dup.e, dup.name.ord)	
    
    work <- rbind(work.e, work.m)
    
    ##-----------------------------##
    ##                             ##
    ##    Merge two records        ##
    ##                             ##
    ##-----------------------------##
    
    w <- merge(r, work, by.x="Markers_name", by.y="name", all.x=T)
    w.ord <- w[order(w$Build36_map_physical_position), c(1:8,10)]
    
    ## Keep a copy of w.ord without Ensembl_map_physical_position being NA
    
    w.ord.c2 <- w.ord.c <- w.ord[!is.na(w.ord$Ensembl_map_physical_position), ]
    
    ## Keep a copy of w.ord with Ensembl_map_physical_position being NA
    
    w.ord.na <- w.ord[is.na(w.ord$Ensembl_map_physical_position), ]
    
    ##-----------------------------##
    ##                             ##
    ##   Check inconsistent        ##
    ##   Ensembl physical order    ##
    ##                             ##
    ##-----------------------------##
    
    ## Deal with inconsistent Ensembl physical position
    
    diff.p <- c(w.ord.c$Ensembl_map_physical_position[1], diff(w.ord.c$Ensembl_map_physical_position, 1))
    inc.p.t <- which(diff.p < 0)
    
    if(length(inc.p.t) > 0) {
      
      for(k in 1:length(inc.p.t)) {
        
        inc.p <- which(diff.p < 0)
        
        if(length(inc.p) > 0) {
          
          if(inc.p[1] == 2) ind.p.bef <- NULL else ind.p.bef <- 1:(inc.p[1]-2)
          ind.p.aft <- which(w.ord.c$Ensembl_map_physical_position >= w.ord.c$Ensembl_map_physical_position[inc.p[1]-1])
          w.ord.c <- w.ord.c[c(ind.p.bef, ind.p.aft),]
          
          diff.p <- c(w.ord.c$Ensembl_map_physical_position[1], diff(w.ord.c$Ensembl_map_physical_position, 1))
          
          k <- k+1
          
        } else { break }
        
      }
      
    }
    
    w.fin <- rbind(w.ord.c, w.ord.na)
    w.ord.fin <- w.fin[order(w.fin$Build36_map_physical_position), ]
    
    write.table(w.ord.fin, file=output.db, append=F, quote=F, col.names=T, row.names=F, sep="\t")
    
    ##-----------------------------##
    ##                             ##
    ##    Create warning file      ##
    ##                             ##
    ##-----------------------------##
    
    ## Markers with multiple records in Ensembl, these are all microsatellites with
    ## two ends of the sequence
    
    if(length(dup$name) != 0){
      
      write("Duplicate Ensembl records found (using smallest base-pair position):",
            file=output.war, append=T)
      write("Name\tChromosome\tPositions", file=output.war, append=T)

      dup.split <- split(dup, dup$name)
      for (m in unique(dup$name))  {
        if (length(unique(as.vector(as.matrix(dup.split[[m]])))) > 3) {
          write(file=output.war,
                gsub(",", "\t", toString(unique(as.vector(as.matrix(dup.split[[m]]))))),
                append=T)
        }
      }
      ##  write.table(dup, file=output.war, append=T, quote=F, col.names=F, row.names=F, sep="\t")
      
    }
    
    ## Markers with inconsistent Ensembl positions
    
    if(length(inc.p.t) > 0) {
      
      inc <- w.ord.c2[! w.ord.c2$Markers_name %in% w.ord.c$Markers_name,]
      inc.ord <- inc[order(inc$Build36_map_physical_position), ]
      
      write("Markers with inconsistent Ensembl map orders:", file=output.war, append=T)
      write.table(inc.ord, file=output.war, append=T, quote=F, col.names=F, row.names=F, sep="\t")
      
    }
    
    ## Check if there are SNPs without physical positions
    unmapped.snps <- w.ord.na[!w.ord.na$Name %in% work.m$name, ]
    
    if (dim(unmapped.snps)[1] > 0) {
      write("Found unmapped SNPs:", file=output.war, append=T)
      write.table(unmapped.snps, file=output.war, append=T, quote=F, col.names=F, row.names=F, sep="\t")
    }
    cat(i)
    cat(" ")
    
  }
  
  ## Process with X chromosome
  
  unlink("chr23.comb")
  unlink("chr23.warnings")
  
  ##----------------------------##
  ##                            ##
  ##    Read in Rugters Data    ##
  ##                            ##
  ##----------------------------##
  
  f <- paste(rpath, "/chr23.sm.map2", sep="")
  
  if (file.exists(f) == TRUE) {
    
    r <- read.table(f, header = T)
    
  } else {
    
    warning(paste("File", f, "does not exist!", sep=" "))
    return(-1)
    
  }
  
  ##-----------------------------##
  ##                             ##
  ##    Query using Rutgers id   ##
  ##                             ##
  ##-----------------------------##
  
  e <- getBM(attributes=c("refsnp_id","chr_name","chrom_start"),filters=c("snp_filter"),values=r$Markers_name,mart=snpmart)
  names(e) <- c("name", "chr", "Ensembl_map_physical_position")
  e <- e[which(e$chr == "X"),]
  
  dup.e <- subset(e, e$name %in% e$name[duplicated(e$name)])
  
  work.e <- e[!e$name %in% dup.e$name,]
  
  cat("X\n")
  
  ##-----------------------------##
  ##                             ##
  ##  Query for microsatellites  ##
  ##                             ##
  ##-----------------------------##
  
  
  micro <- r$Markers_name[!r$Markers_name %in% work.e$name]
  
  n1 <- length(micro)
  
  if (n1 > 0) {
    cat("Creating merged maps for microsatellites on all chromosomes ...")
  }
  q1 <- NULL
  
  for(j in 1:n1) {
    
    if(j == 1) q1 <- paste("select marker_id, name from marker_synonym where name in (", micro[j], ",", sep="'")
    
    if(j > 1 & j < n1) q1 <- paste(q1, micro[j], ",", sep="'")
    
    if(j == n1) q1 <- paste(q1, micro[j], ")", sep="'")
    
  }
  
  m.id <- dbGetQuery(con,q1)
  
  n2 <- length(m.id$marker_id)
  
  q2 <- NULL
  
  for(j in 1:n2) {
    
    if(j == 1) q2 <- paste("select marker_id, seq_region_id, seq_region_start, seq_region_end from marker_feature where marker_id in (", m.id$marker_id[j], ",", sep="")
    
    if(j > 1 & j < n2) q2 <- paste(q2, m.id$marker_id[j], ",", sep="")
    
    if(j == n2) q2 <- paste(q2, m.id$marker_id[j], ")", sep="")
    
  }
  
  m.ft <- dbGetQuery(con,q2)
  
  n3 <- length(m.ft$seq_region_id)
  
  q3 <- NULL
  
  for(j in 1:n3) {
    
    if(j == 1) q3 <- paste("select seq_region_id, name from seq_region where seq_region_id in (", m.ft$seq_region_id[j], ",", sep="")
    
    if(j > 1 & j < n3) q3 <- paste(q3, m.ft$seq_region_id[j], ",", sep="")
    
    if(j == n3) q3 <- paste(q3, m.ft$seq_region_id[j], ")", sep="")
    
  }
  
  m.chr <- dbGetQuery(con,q3)
  
  names(m.chr) <- c("seq_region_id", "chr")
  
  m1 <- merge(m.ft, m.chr, by.x="seq_region_id", by.y="seq_region_id", all.x=T)
  
  work.m1 <- m1[m1$chr == "X",]
  
  dup.id <- subset(work.m1, work.m1$marker_id %in% work.m1$marker_id[duplicated(work.m1$marker_id)])
  dup.id.ord <- dup.id[order(dup.id$marker_id), ]
  
  work.m2 <- work.m1[!work.m1$marker_id %in% dup.id.ord$marker_id,]
  
  m2 <- merge(m.id, work.m2, by.x="marker_id", by.y="marker_id", all.x=T)
  
  work.m3 <- m2[!is.na(m2$seq_region_id),c(2,6,4)]
  
  names(work.m3) <- c("name", "chr", "Ensembl_map_physical_position")
  
  dup.name <- subset(work.m3, work.m3$name %in% work.m3$name[duplicated(work.m3$name)])
  dup.name.ord <- dup.name[order(dup.name$name), ]
  
  work.m4 <- work.m3[!work.m3$name %in% dup.name.ord$name,]
  
  ##---------------------------------------##
  ##                                       ##
  ## Dealing with Micro with multiple hits ##
  ##                                       ##
  ##---------------------------------------##
  
  uniq <- dup.name.ord[!duplicated(dup.name.ord$name),"name"]
  
  work.m5 <- NULL
  
  for(k in 1:length(uniq)) {
    
    rec <- dup.name.ord[dup.name.ord$name==uniq[k],]
    
    if(length(rec$Ensembl_map_physical_position) > 1) pdiff <- max(rec$Ensembl_map_physical_position) - min(rec$Ensembl_map_physical_position) else pdiff <- 0
    
    if(pdiff <= thres) work.m5 <- rbind(work.m5, rec[rec$Ensembl_map_physical_position==min(rec$Ensembl_map_physical_position),])
    
  }
  
  work.m6 <- work.m5[!duplicated(work.m5$name),]
  
  work.m <- rbind(work.m4, work.m6)
  
  dup <- rbind(dup.e, dup.name.ord)	
  
  work <- rbind(work.e, work.m)
  
  ##-----------------------------##
  ##                             ##
  ##    Merge two records        ##
  ##                             ##
  ##-----------------------------##
  
  w <- merge(r, work, by.x="Markers_name", by.y="name", all.x=T)
  w.ord <- w[order(w$Build36_map_physical_position), c(1:5,7)]
  
  ## Keep a copy of w.ord without Ensembl_map_physical_position being NA
  
  w.ord.c2 <- w.ord.c <- w.ord[!is.na(w.ord$Ensembl_map_physical_position), ]
  
  ## Keep a copy of w.ord with Ensembl_map_physical_position being NA
  
  w.ord.na <- w.ord[is.na(w.ord$Ensembl_map_physical_position), ]
  
  ##-----------------------------##
  ##                             ##
  ##   Check inconsistent        ##
  ##   Ensembl physical order    ##
  ##                             ##
  ##-----------------------------##
  
  ## Deal with inconsistent Ensembl physical position
  
  diff.p <- c(w.ord.c$Ensembl_map_physical_position[1], diff(w.ord.c$Ensembl_map_physical_position, 1))
  inc.p.t <- which(diff.p < 0)
  
  if(length(inc.p.t) > 0) {
    
    for(k in 1:length(inc.p.t)) {
      
      inc.p <- which(diff.p < 0)
      
      if(length(inc.p) > 0) {
        
        if(inc.p[1] == 2) ind.p.bef <- NULL else ind.p.bef <- 1:(inc.p[1]-2)
        ind.p.aft <- which(w.ord.c$Ensembl_map_physical_position >= w.ord.c$Ensembl_map_physical_position[inc.p[1]-1])
        w.ord.c <- w.ord.c[c(ind.p.bef, ind.p.aft),]
        
        diff.p <- c(w.ord.c$Ensembl_map_physical_position[1], diff(w.ord.c$Ensembl_map_physical_position, 1))
        
        k <- k+1
        
      } else { break }
      
    }
    
  }
  
  w.fin <- rbind(w.ord.c, w.ord.na)
  w.ord.fin <- w.fin[order(w.fin$Build36_map_physical_position), ]
  
  write.table(w.ord.fin, file="chr23.comb", append=F, quote=F, col.names=T, row.names=F, sep="\t")
  
  ##-----------------------------##
  ##                             ##
  ##    Create warning file      ##
  ##                             ##
  ##-----------------------------##
  
  ## Bad Markers with multiple records in Ensembl

  output.war <- "chr23.warnings"
  if(length(dup$name) != 0){
    
    write("Duplicate Ensembl records found (using smallest base-pair position):",
          file=output.war, append=T)
    write("Name\tChromosome\tPositions", file=output.war, append=T)

    dup.split <- split(dup, dup$name)
    for (m in unique(dup$name))  {
      if (length(unique(as.vector(as.matrix(dup.split[[m]])))) > 3) {
        write(file=output.war,
              gsub(",", "\t", toString(unique(as.vector(as.matrix(dup.split[[m]]))))),
              append=T)
      }
    }
    ##  write("Duplicate Markers Found in Ensembl:", file="chr23.warnings", append=T)
    ##  write.table(dup, file="chr23.warnings", append=T, quote=F, col.names=F, row.names=F, sep="\t")
    
  }
  
  ##  Markers with inconsistent Ensembl positions
  
  if(length(inc.p.t) > 0) {
    
    inc <- w.ord.c2[! w.ord.c2$Markers_name %in% w.ord.c$Markers_name,]
    inc.ord <- inc[order(inc$Build36_map_physical_position), ]
    
    write("Markers with Inconsistent Ensembl Map Orders:", file=output.war, append=T)
    write.table(inc.ord, file=output.war, append=T, quote=F, col.names=F, row.names=F, sep="\t")
    
  }
  
  ## Check if there are SNPs without physical positions
  unmapped.snps <- w.ord.na[!w.ord.na$Name %in% work.m$name, ]
  if (dim(unmapped.snps)[1] > 0) {
    write("Found unmapped SNPs:", file=output.war, append=T)
    write.table(unmapped.snps, file=output.war, append=T, quote=F, col.names=F, row.names=F, sep="\t")
  }
  
  setwd(mpath)
  return(as.integer(ver))
  
  
}


