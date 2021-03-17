#!/usr/bin/env Rscript
library(Biostrings)
library(stringr)
library(data.table)

#args[1]: gff file
#args[2]: ref file
#args[3]: alignment (must be aligned to matching ref (and gff))
#args[4]: out directory

args = commandArgs(trailingOnly=TRUE)

#Function to split gff lines in columns
parse_gff <- function(x) {
  split <- stringr::str_split(stringr::str_split(x,';',simplify = T),'=',simplify = T)
  decoded_values <- sapply(split[,2], URLdecode)
  parsed <- data.table(t(decoded_values))
  names(parsed) <- split[,1]
  return(parsed)
}

#Read gff and only retain usefull columns
# cov_gff <- fread(args[1], skip=2)
cov_gff <- fread(args[1], skip=2)

names(cov_gff) <- c("seqid","source","type","start","end","score","strand","phase","attributes")

cov_gff <- cbind(cov_gff, cov_gff[,rbindlist(lapply(attributes,parse_gff),use.names = T, fill = T)])
cov_gff <- cov_gff[type%in%c("mature_protein_region_of_CDS","CDS")]
cov_gff[,attributes:=NULL]
cov_gff[,name:=gene]
cov_gff[is.na(name),name:=product]

#add rbd
cov_gff <- rbind(NA,(cov_gff), fill=T)
cov_gff[1,name:="S-rbd"]
cov_gff[1,start:=21563+1005]
cov_gff[1,end:=21563+1547]

#Read alignment
aln <- readDNAStringSet(args[3])
ref <- readDNAStringSet(args[2])
aln <- c(ref,aln)

#Create output dir
dir.create(args[4])

#For each gene write a fasta of the gene ORF
apply(cov_gff, 1, function(line) {
  cov_gff_region <- subseq(aln,(as.numeric(line["start"])-1),(as.numeric(line["end"])-1))
  names(cov_gff_region) <- paste(names(cov_gff_region),line["name"], sep = '|')
  writeXStringSet(cov_gff_region, paste0(args[4],"/",line["name"],".fna"))
})