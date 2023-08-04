#!/usr/bin/env Rscript

rm(list = ls())

#automatic install of packages if they are not installed already
list.of.packages <- c(
  "seqinr",
  "optparse",
  "rjson"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE)
}

library(seqinr)
library(optparse)
library(rjson)

#

check_phase_pos = function(coding_pos, codon, prev_phase, pos_to_remove){
  vect_cds_id = unlist(lapply(names(codon[[coding_pos]]), function(x) if (x != "coding") return(x)))
  phase = unlist(lapply(vect_cds_id, function(cds_id) codon[[coding_pos]][[cds_id]]$pos_codon))
  if ((prev_phase + 1) %in% phase){
    prev_phase = prev_phase + 1
  } else {
    pos_to_remove = c(pos_to_remove, coding_pos)
  }
  if (prev_phase == 2){
    prev_phase = -1
  }
  return(list(pos_to_remove, prev_phase))
}


is_gap = function(nuc, pos) {
  if (nuc == '-'){
    return(NA)
  } else {
    temp = which(! is.na(pos))
    if (length(temp) > 0){
      return(as.numeric(tail(pos[temp], 1)) + 1)
    }
    else {
      return(1)
    }
  }
}


insertRow <- function(existingDF, newrow, r) {
  existingDF[seq(r+1,nrow(existingDF)+1),] <- existingDF[seq(r,nrow(existingDF)),]
  existingDF[r,] <- newrow
  existingDF
}



#

### PARAMETERS

option_list = list(
  make_option(c("-i", "--align_path"), type="character", help="absolute or relative path to alignment directory", metavar="character"),
  make_option(c("-n", "--version_nuc"), type="character", default="nuc", help="first version to patch [default= %default]"),
  make_option(c("-a", "--version_aa"), type="character", default="aa", help="second version to patch [default= %default]"),
  make_option(c("-j", "--json_path"), type="character", help="absolute or relative path to json directory"),
  make_option(c("-m", "--method"), type="character", default="mafft", help="aligner method [default= %default]")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

###

align_path = opt$align_path
vers = opt$version_nuc
vers_aa = opt$version_aa
json_path = opt$json_path
method = opt$method

#

path_outfile = paste(align_path, "/", vers, "-vs-", vers_aa, "_match_align.csv", sep = "")
path_alignment_aa = paste(align_path, "/outfile_", method ,"/", vers_aa , "/multiple.aln", sep = "")
path_position_aa = paste(align_path, "/infile_", method, "/", vers_aa , "/multiple_position.fasta", sep = "")
path_alignment = paste(align_path, "/outfile_", method, "/", vers , "/multiple.aln", sep = "")

align_aa = read.alignment(file = path_alignment_aa, format = "clustal")
pos_aa = read.fasta(path_position_aa, as.string = T, strip.desc = T)
pos_aa = lapply(pos_aa, function(string) as.integer(strsplit(string, ",")[[1]]) + 1)

align = read.alignment(file = path_alignment, format = "clustal")
pos = setNames(lapply(align$nam, function(virus) 1:length(which(strsplit(align$seq[[which(align$nam == virus)]][1], "")[[1]] != '-'))), align$nam)

pos = setNames(lapply(align$nam, function(virus) 1:length(which(strsplit(align$seq[[which(align$nam == virus)]][1], "")[[1]] != '-'))), align$nam)

len_align = nchar(align$seq[[1]])

removed_non_coding = c()
removed_pos = c()
empty_align = F
species_list = align$nam

change_coding = list()
res_list = list()
cds_pos = list()
nuc_seq = list()
aa_seq = list()
frame_seq = list()
pos_to_remove = c()
for (species in species_list){
  # deal with the case where all the alignment is removed
  if (length(removed_pos) >= len_align) {
    empty_align = T
    break
  }
  
  # load json file
  path_json = paste(json_path, "/", species, "_codon.cross", sep = "")
  codon = fromJSON(file = path_json)
  
  # get aa by position
  aa_seq[[species]] = c()
  frame_seq[[species]] = c()
  p = length(aa_seq[[species]])
  p_cds = ""
  for (position in 1:length(codon)){
    if ( ! codon[[position]]$`coding`){
      aa_seq[[species]] = c(aa_seq[[species]], NA)
      frame_seq[[species]] = c(frame_seq[[species]], NA)
    } 
    else {
      all_keys = ls(codon[[position]])
      if (length(all_keys) > 2){
        for(cds in rev(all_keys)) {
          if (cds == "coding" | cds == p_cds){
            next
          }
          if ('aa' %in% ls(codon[[position]][[cds]])){
            aa_seq[[species]] = c(aa_seq[[species]], codon[[position]][[cds]]$aa)
            frame_seq[[species]] = c(frame_seq[[species]], codon[[position]][[cds]]$pos_codon)
          } else {
            aa_seq[[species]] = c(aa_seq[[species]], NA)
            frame_seq[[species]] = c(frame_seq[[species]], NA)
          }
          p_cds = ""
          break
        }
      } else {
        for( cds in all_keys){
          if (cds != "coding"){
            if ('aa' %in% ls(codon[[position]][[cds]])) {
              aa_seq[[species]] = c(aa_seq[[species]], codon[[position]][[cds]]$aa)
              frame_seq[[species]] = c(frame_seq[[species]], codon[[position]][[cds]]$pos_codon)
            } else {
              aa_seq[[species]] = c(aa_seq[[species]], NA)
              frame_seq[[species]] = c(frame_seq[[species]], NA)
            }
            p_cds = cds
            break
          }
        }
      }
      # }
    }
    if(p == length(aa_seq[[species]])){
      # no aa added for this position, pb !
      print(species)
      print(position)
    }
    p = length(aa_seq[[species]])
  }
  
  change_coding[[species]] = c()
  x = 2
  while(x < (length(aa_seq[[species]]) - 2)) {
    if (length(which(is.na((aa_seq[[species]][(x - 1) : (x + 2)])))) > 0) {
      x = x + 1
      next
    }
    if (aa_seq[[species]][x] != aa_seq[[species]][x - 1] & (aa_seq[[species]][x] != aa_seq[[species]][x + 1] | aa_seq[[species]][x + 1] != aa_seq[[species]][x + 2])){
      change_coding[[species]] = c(change_coding[[species]], c(x, x + 1))
    }
    x = x + 1
  }
  change_coding[[species]] = c(change_coding[[species]], which(aa_seq[[species]] == '*'))
  
  # find non coding position for the given species
  region = 1:len_align
  
  coding_bool = unlist(lapply(region, function(x) codon[[x]]$coding))
  coding_bool[change_coding[[species]]] = FALSE
  
  non_coding_pos = region[which(coding_bool == 0)]
  
  # remove non coding position
  if (length(non_coding_pos) > 0) {
    # check if this non-coding pos are already removed
    non_coding = unlist(lapply(non_coding_pos, function(x) which(region == x)))
  }
  
  # check if the sequence coding is phased
  prev_phase = -1
  
  for (coding_pos in region[which(coding_bool)]){
    handler = check_phase_pos(coding_pos, codon, prev_phase, pos_to_remove)
    pos_to_remove = handler[[1]]
    prev_phase = handler[[2]]
  }
  
  # edit alignment
  if (length(pos_to_remove) > 0){
    removed_pos = c(removed_pos, unique(pos_to_remove))
  }
  
  res_list[[species]] = sort(pos_to_remove)
  
  # retrieve coding coding pos in the raw sequence
  # load json file
  path_json_annot = paste(json_path, "/", species, ".txt", sep = "")
  
  annot = fromJSON(file = path_json_annot)
  
  cds_pos[[species]] = 1:nchar(annot$complete_sequence) %in% unlist(lapply(annot$cds, function(cds_id) return((cds_id$start + 1):cds_id$end)))
  nuc_seq[[species]] = unlist(strsplit(annot$complete_sequence, ''))
}

# find gap position in the aa align
aa_gap = list()
for (species in align_aa$nam){
  aa_gap[[species]] = which(unlist(strsplit(align_aa$seq[[which(align_aa$nam == species)]], '')) == '-')
}


df = data.frame()

# Merges nucleotide alignment and coding informations

for (species in align$nam){
  path_json_annot = paste(json_path, "/", species, "_codon.cross", sep = "")
  codon = fromJSON(file = path_json)
  coding_bool = unlist(lapply(names(codon), function(x) codon[[x]]$coding))
  
  gap_temp = which(unlist(strsplit(align$seq[[which(align$nam == species)]], '')) == '-')
  gap_region = split(gap_temp, cumsum(c(1, diff(gap_temp) != 1)))
  
  cds_pos_aligned = list()
  cds_pos_aligned[[species]] = cds_pos[[species]]
  
  nuc_seq_gapped = list()
  nuc_seq_gapped[[species]] = nuc_seq[[species]]
  
  nuc_seq_aligned = unlist(strsplit(align$seq[[which(align$nam == species)]], ''))
  
  gap_index_offset = -1
  for (gap in gap_region){
    cds_pos_aligned[[species]] = append(cds_pos_aligned[[species]], rep(NA, length(gap)), after = gap[1] + gap_index_offset)
    nuc_seq_gapped[[species]] = append(nuc_seq_gapped[[species]], rep(NA, length(gap)), after = gap[1] + gap_index_offset)
  }
  cds_pos_aligned[[species]][change_coding[[species]]] = FALSE
  
  pos = c()
  for (nuc in unlist(strsplit(align$seq[[which(align$nam == species)]], ''))){
    pos = c(pos, is_gap(nuc, pos))
  }
  
  temp = data.frame('pos_align_nuc' = 1:nchar(align$seq[[1]]),
                    'species_pos_nuc' = pos,
                    'species' = rep(species, length(pos)),
                    'coding' = cds_pos_aligned[[species]],
                    'nuc_aligned' = toupper(nuc_seq_aligned),
                    'nuc_species' = toupper(nuc_seq_gapped[[species]]),
                    'aa_species' = toupper(aa_seq[[species]]),
                    'frame' = frame_seq[[species]]
  )
  
  
  # align aa gestion
  
  len_aa_align = nchar(align_aa$seq[[which(align_aa$nam == species)]])
  
  handler_aa_pos = c()
  for (i in temp$species_pos){
    if (is.na(i)){
      handler_aa_pos = c(handler_aa_pos, NA)
    }
    else if (as.integer(i) %in% unlist(pos_aa[[species]])){
      handler_aa_pos = c(handler_aa_pos, i)
    } else {
      handler_aa_pos = c(handler_aa_pos, NA)
    }
  }
  
  aa_species_align = toupper(unlist(strsplit(align_aa$seq[[which(align_aa$nam == species)]], '')))
  pos_aa_aligned = c()
  gap_in_aa_align = c()
  new_row = data.frame()
  index_new_row = c()
  current_match_index = 0
  index_align_aa = 1
  codon_index = 0
  prev_gap_len = 0
  end_align_aa = FALSE
  for(index_align_nuc in 1:len_align){
    if (is.na(temp$coding[index_align_nuc])) {
      pos_aa_aligned = c(pos_aa_aligned, NA)
      next
    }
    # non coding position : no corresponding pos in the aa align
    if (! temp$coding[index_align_nuc]) {
      pos_aa_aligned = c(pos_aa_aligned, NA)
      next
    }
    # pos in aa align is a gap, must add a new line with no match with pos align nuc
    if (!end_align_aa){
      gap_len = 0
      while(aa_species_align[index_align_aa] == '-') {
        gap_len = gap_len + 1
        new_row = rbind.data.frame(new_row, c(NA, NA, species, NA, NA, NA, NA, NA, '-', index_align_aa))
        index_new_row = c(index_new_row, index_align_nuc + prev_gap_len + gap_len - 1)
        index_align_aa = index_align_aa + 1
        if (index_align_aa >= len_aa_align){
          index_align_aa = len_aa_align
          end_align_aa = TRUE
        }
      }
      prev_gap_len = prev_gap_len + gap_len
    }
    
    # matching position
    if (temp$coding[index_align_nuc]) {
      match_index = min(which(temp$aa_species[(current_match_index + 1):len_align] == aa_species_align[index_align_aa]))
      pos_aa_aligned = c(pos_aa_aligned, index_align_aa)
      codon_index = codon_index + 1
    }
    
    if (codon_index > 2) {
      codon_index = 0
      index_align_aa = index_align_aa + 1
      current_match_index = current_match_index + match_index + 2
      if (index_align_aa >= len_aa_align){
        index_align_aa = len_aa_align
        end_align_aa = TRUE
      }
    }
  }
  
  aa_aligned = aa_species_align[as.numeric(pos_aa_aligned)]
  
  temp = data.frame(temp, 
                    "aa_aligned" = aa_aligned,
                    "pos_align_aa" = pos_aa_aligned)
  i = 1
  for (index in index_new_row){
    temp = insertRow(temp, new_row[i, ], index)
    i = i + 1
  }
  
  df = rbind(df, temp)
}

# # define stop codon as coding
# df$coding[which(df$aa_species == '*')] = TRUE

write.table(df, path_outfile, col.names = T, row.names = F, quote = F, sep = '\t')
