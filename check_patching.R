#!/usr/bin/env Rscript

rm(list = ls())

#automatic install of packages if they are not installed already
list.of.packages <- c(
  "seqinr",
  "optparse",
  "stringr", 
  "tibble",
  "ape")

list.of.Bioconductor.packages = c(
  "msa")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

library(seqinr)
library(optparse)
library(stringr)
library(tibble)
library(ape)

#

reduce_alignment_region = function(path_alignment, start_region, end_region) {
  # load alignment
  align = read.alignment(file = path_alignment, format = "clustal")
  
  # restrict alignment on the given region
  region = start_region : end_region
  all_align_nrr = lapply(align$nam, function(virus) paste(strsplit(align$seq[[which(align$nam == virus)]][1], "")[[1]][region], collapse = ""))
  region_align = list("nam" = align$nam,
                      "seq" = all_align_nrr,
                      "nb" = align$nb,
                      "com" = align$com[align$nam %in% align$nam])
  
  return(region_align)
}

#

### PARAMETERS

option_list = list(
  make_option(c("-i", "--align_path"), type="character", help="absolute or relative path to alignment directory", metavar="character"),
  make_option(c("-n", "--version_nuc"), type="character", default="nuc", help="first version to patch [default= %default]"),
  make_option(c("-a", "--version_aa"), type="character", default="aa", help="second version to patch [default= %default]"),
  make_option(c("-m", "--method"), type="character", default="mafft", help="aligner method [default= %default]")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

align_path = "C:/Users/jean-/OneDrive/Documents/SARS-cov2/CNCA_standalone/tmp/run_2/mafft/"
vers = "nuc"
vers_aa = "aa"
method = "mafft"

###

align_path = opt$align_path
vers = opt$version_nuc
vers_aa = opt$version_aa
method = opt$method

#


path_alignment = paste(align_path, "/outfile_", method, "/", vers , "/multiple.aln", sep = "")
align = read.alignment(file = path_alignment, format = "clustal")

# retrieve patched alignment
path_patch = paste(align_path, "/patching_", vers, "-vs-", vers_aa, "/", sep = "")
files_patch = list.files(path_patch, "^.+-.+.fasta$", full.names = T)
files_patch = files_patch[order(as.numeric(sub("-.*.fasta$", "", gsub("^.*/", "", files_patch))))]
all_patch = do.call(list, lapply(files_patch, function(x) read.alignment(file = x, format = "fasta")))

names(all_patch) = sub(".fasta$", "", gsub(".*\\/", "", files_patch))

# check 1
for (region_index in 1:length(files_patch)){
  region = str_split(names(all_patch)[region_index], "-")[[1]]
  start_region = as.numeric(region[1])
  end_region = as.numeric(region[2])
  
  raw_align = reduce_alignment_region(path_alignment, start_region, end_region)
  patch_align = all_patch[[region_index]]
  
  for (seq_index in 1:length(raw_align$seq)){
    raw = gsub("-", "", raw_align$seq[[seq_index]])
    patch = gsub("-", "", patch_align$seq[[seq_index]])
    if (raw != patch){
      print(paste("region: ", region_index, ", seq: ", seq_index, sep = " "))
    }
  }
}


# patch alignment
df_align = data.frame()
for (virus in align$nam){
  seq = strsplit(align$seq[[which(align$nam == virus)]][1], "")[[1]]
  df_align = rbind(df_align, seq)
}

patched_align = align
raw_length = nchar(align$seq[[1]][1])

for (region_index in 1:length(files_patch)){
  # print(paste("region: ", region_index, " / ", length(files_patch), sep = ""))
  patch_length = ncol(df_align)
  region = str_split(names(all_patch)[region_index], "-")[[1]]
  shift = patch_length - raw_length
  start_region = as.numeric(region[1]) + shift
  end_region = as.numeric(region[2]) + shift
  seq_region = as.numeric(start_region) : as.numeric(end_region)

  # deletion on the raw region
  df_align = df_align[, -seq_region]
  
  # retrieve patched region
  raw_patch_region = all_patch[[region_index]]
  
  patch_region = data.frame()
  for (virus in raw_patch_region$nam){
    seq = strsplit(raw_patch_region$seq[[which(raw_patch_region$nam == virus)]][1], "")[[1]]
    patch_region = rbind(patch_region, seq)
  }
  
  # insert patched region
  correction = 0
  for(col_id in 1:ncol(patch_region)){
    if (shift != 0){
      correction = shift
    }
    df_align = add_column(.data = df_align, patch_region[, col_id], .after = start_region + shift + col_id - 2 - correction, .name_repair = "minimal")
  }
}

# check 2
for (seq_index in 1:length(align$seq)){
  raw = gsub("-", "", align$seq[[seq_index]])
  patch = gsub("-", "", paste(df_align[seq_index, ], collapse = ""))
  if (raw != patch){
    print(paste("seq: ", seq_index, sep = " "))
  }
}

# write patched alignment
row.names(df_align) = patched_align$nam
df_align = setNames(split(df_align, seq(nrow(df_align))), rownames(df_align))
path_out_patched_alignment = paste(align_path, "/cnca_alignment.fasta", sep = "")
write.fasta(sequences = df_align, file.out = path_out_patched_alignment, names = names(df_align))

path_out_patched_alignment = paste(align_path, "/cnca_alignment.nexus", sep = "")
write.nexus.data(df_align, path_out_patched_alignment)

# rewrite nt and prot alignment in fasta and nexus format
##nt
df_nt = setNames(split(align$seq, seq(length(align$seq))), align$nam)
write.fasta(sequences = df_nt, 
            file.out =  paste(align_path, "/nucleotide_alignment.fasta", sep = ""), 
            names = names(df_align))

write.nexus.data(df_nt, 
                 file = paste(align_path, "/nucleotide_alignment.nexus", sep = ""))

## aa
aa_align = read.alignment(file = paste(align_path, "/outfile_", method, "/aa/multiple.aln", sep = ""),
                          format = "clustal")
df_aa = setNames(split(aa_align$seq, seq(length(aa_align$seq))), aa_align$nam)
write.fasta(sequences = df_aa, 
            file.out =  paste(align_path, "/protein_alignment.fasta", sep = ""), 
            names = names(df_aa))

write.nexus.data(df_aa, 
                 file = paste(align_path, "/protein_alignment.nexus", sep = ""))

