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

insertCol <- function(existingDF, newcol, r) {
  existingDF[, seq(r+1,ncol(existingDF)+1)] <- existingDF[, seq(r,ncol(existingDF))]
  existingDF[, r] <- newcol
  existingDF
}


#

option_list = list(
  make_option(c("-i", "--align_path"), type="character", help="absolute or relative path to alignment directory", metavar="character"),
  make_option(c("-n", "--version_nuc"), type="character", default="nuc", help="first version to patch [default= %default]"),
  make_option(c("-a", "--version_aa"), type="character", default="aa", help="second version to patch [default= %default]")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# align_path = "C:/Users/jean-/OneDrive/Documents/SARS-cov2/CNCA_standalone/tmp/run_2/mafft/"
# vers = "nuc"
# vers_aa = "aa"
# method = "mafft"

###

align_path = opt$align_path
vers = opt$version_nuc
vers_aa = opt$version_aa

#

path_match_table = paste(align_path, "/", vers, "-vs-", vers_aa, "_match_align.csv", sep = "")

df = read.csv2(path_match_table, sep = '\t')


# detect position in the nucleotide alignment that need to be patch
to_patch = c() # vector storing all the position to patch
for(pos_nuc in unique(df$pos_align_nuc)) {
  # retrieve frame of each nuc if it is involved in a codon
  coding = df$frame[which(df$pos_align_nuc == pos_nuc)]
  coding = coding[! is.na(coding)]
  
  # retrieve the corresponding aa if it is involved in a codon
  aa = df$pos_align_aa[which(df$pos_align_nuc == pos_nuc)]
  aa = aa[! is.na(aa)]
  
  # if there is at least one coding nuc at this position
  if (length(coding) > 0){
    # if there is a phase mismatch, store this position
    if (length(unique(coding)) != 1){
      to_patch = c(to_patch, pos_nuc)
      # if there is aa mismatch (col position in the aa alignment), store this position
    } else if (length(unique(aa)) != 1){
      to_patch = c(to_patch, pos_nuc)
    }
  }
}

# regroup position to patch into region (successive position)
to_patch_region = split(to_patch, cumsum(c(1, diff(to_patch) != 1)))

len_align = length(unique(df$pos_align_nuc)) - 1

unresolved_region = c()

done_region = 0
# patch each region to patch
for(patch_region in 1:length(to_patch_region)){
  unresolved = F
  print(patch_region)
  
  # define region to patch on the nuc multiple alignment
  start_align_nuc = min(to_patch_region[[patch_region]])
  end_align_nuc = max(to_patch_region[[patch_region]])
  
  # find true start: nuc multiple alignment position to the left matching with the same aa for each species and not NA
  true_start_align_nuc = start_align_nuc
  while(T){
    start_frame = unique(df$frame[which(df$pos_align_nuc == as.character(true_start_align_nuc))])
    if (length(start_frame) == 1){
      if (! is.na(start_frame)){
        if (start_frame == 0){
          if(length(unique(df$pos_align_aa[which(df$pos_align_nuc == as.character(true_start_align_nuc))])) == 1){
            break
          }
        }
      }
    }
    true_start_align_nuc = true_start_align_nuc - 1
  }
  
  # find true end: nuc multiple alignment position to the rigth matching with the same aa for each species and not NA
  true_end_align_nuc = end_align_nuc
  while(T){
    end_frame = unique(df$frame[which(df$pos_align_nuc == as.character(true_end_align_nuc))])
    if (length(end_frame) == 1){
      if (! is.na(end_frame)){
        if (end_frame == 0){
          if(length(unique(df$pos_align_aa[which(df$pos_align_nuc == as.character(true_end_align_nuc))])) == 1){
            break
          }
        }
      }
    }
    true_end_align_nuc = true_end_align_nuc + 1
    if (true_end_align_nuc == len_align){
      break
    }
  }
  
  # if the region extension correspond to a previously run patch region, skip this region
  if (true_start_align_nuc < done_region | true_end_align_nuc < done_region){
    next
  }
  
  # get all the pos in the nucleotide alignment to use
  list_pos_align_nuc = as.character(true_start_align_nuc:true_end_align_nuc)
  
  # detect universal non-coding pos in these positions
  non_coding_pos = unlist(lapply(list_pos_align_nuc, function(pos_nuc) if (sum((df[which(df$pos_align_nuc == pos_nuc), "coding"]), na.rm = T) == 0) {pos_nuc} ))
  
  # store this information in a data frame
  pos_to_patch = data.frame( "coding" = unlist(lapply(list_pos_align_nuc, function(x) if (x %in% non_coding_pos) {F} else {T})),
                             "pos_align_nuc" = list_pos_align_nuc)
  
  # get the aa position of each coding nuc
  list_pos_align_aa = unlist(lapply(1:length(pos_to_patch$pos_align_nuc), function(i) if (pos_to_patch$coding[i]) {
    df$pos_align_aa[which(df$pos_align_nuc == pos_to_patch$pos_align_nuc[i])][ !is.na(df$pos_align_aa[which(df$pos_align_nuc == pos_to_patch$pos_align_nuc[i])])][1]
  } else {NA}))
  
  # update the data frame with this information
  pos_to_patch = data.frame(pos_to_patch, 
                            "pos_align_aa" = list_pos_align_aa)
  
  # detect indel in the non coding ie nucleotide position that are at least one time coding for a species and non coding for another one
  indel = c()
  for (pos in pos_to_patch$pos_align_nuc){
    coding_state = df$coding[which(df$pos_align_nuc == pos)]
    coding_state = coding_state[! is.na(coding_state)]
    nb_cand = length(coding_state)
    res = sum(coding_state)
    if (res != nb_cand & res != 0){
      indel = c(indel, pos)
    }
  }
  
  patch_nuc = data.frame("species" = unique(df$species)) # handler of the patched position
  done_pos_align_aa = c()
  # browse all position to patch for the given region
  for (i in 1:dim(pos_to_patch)[1]){
    temp_nuc = data.frame()
    for (species in unique(df$species)){
      # deal with non universal aligned coding position
      if (pos_to_patch$pos_align_nuc[i] %in% indel){
        # coding position for the given species: gap insertion
        gap_temp_handler = df$coding[which(df$pos_align_nuc == pos_to_patch$pos_align_nuc[i] & df$species == species)]
        if (is.na(gap_temp_handler)){
          temp_nuc = rbind(temp_nuc, "-")
        } else if(gap_temp_handler){
          temp_nuc = rbind(temp_nuc, "-")
        } else {
          temp_nuc = rbind(temp_nuc, df$nuc_species[which(df$pos_align_nuc == pos_to_patch$pos_align_nuc[i] & df$species == species)])
        }
        colnames(temp_nuc) = paste("indel_", pos_to_patch$pos_align_nuc[i], sep = "")
        next
      }
      ###
      # universal non coding position
      if (! pos_to_patch$coding[i]){
        patch = df$nuc_species[which(df$pos_align_nuc == pos_to_patch$pos_align_nuc[i] & df$species == species)]
        if (is.na(patch)){
          temp_nuc = rbind(temp_nuc, "-")
        } else {
          temp_nuc = rbind(temp_nuc, patch)
        }
        colnames(temp_nuc) = paste("nuc_", pos_to_patch$pos_align_nuc[i], sep = "")
      } else if (! pos_to_patch$pos_align_aa[i] %in% done_pos_align_aa) {
        patch = df$nuc_species[which(df$pos_align_aa == pos_to_patch$pos_align_aa[i] & df$species == species)]
        if (length(patch) == 1) {
          if (is.na(patch)){
            temp_nuc = rbind(temp_nuc, rep("-", 3))
          } 
        } else {
          temp_nuc = rbind(temp_nuc, patch)
        }
        colnames(temp_nuc) = paste("aa", pos_to_patch$pos_align_aa[i], 0:2, sep = "_")
      }
    }
    if (length(temp_nuc) > 0){
      done_pos_align_aa = c(done_pos_align_aa, pos_to_patch$pos_align_aa[i])
      patch_nuc = cbind(patch_nuc, temp_nuc)
    }
  } 
  
  # check if aa are well ordered, and reorder it if needed
  aa_pos = grep(pattern = "^aa_[0-9]*_[0-2]+$", x = colnames(patch_nuc))
  new_aa_order = aa_pos[order(colnames(patch_nuc)[aa_pos])]
  new_col_names = colnames(patch_nuc)
  new_col_names[aa_pos] = new_col_names[new_aa_order]
  patch_nuc[aa_pos] = patch_nuc[new_aa_order]
  colnames(patch_nuc) = new_col_names
  
  # first scan for missing aa
  aa_to_find = unique(pos_to_patch$pos_align_aa[!is.na(pos_to_patch$pos_align_aa)])
  aa_present = as.numeric(unique(gsub("_[0-9]+$", "", gsub("^aa_", "", colnames(patch_nuc)[grep("^aa_.*", colnames(patch_nuc))]))))
  
  aa_mising = sort(aa_to_find[!aa_to_find %in% aa_present])
  
  place_accord_indel = c()
  
  if (length(aa_mising) > 0){
    for (aa_pos in aa_mising){
      temp_nuc = c()
      for (species in unique(df$species)){
        patch = df$nuc_species[which(df$pos_align_aa == as.character(aa_pos) & df$species == species)]
        if (length(patch) == 0) {
          temp_nuc = rbind(temp_nuc, rep("-", 3))
        }
        if (length(patch) == 1) {
          if (is.na(patch)){
            temp_nuc = rbind(temp_nuc, rep("-", 3))
          } 
        } else {
          temp_nuc = rbind(temp_nuc, patch)
        }
        colnames(temp_nuc) = paste("aa", aa_pos, 0:2, sep = "_")
      }
      for (col_id in 1:ncol(temp_nuc)){
        # if no aa in the patch
        if (length(aa_present) == 0){
          # find indel region
          region_indel = split(as.numeric(indel), cumsum(c(1, diff(as.numeric(indel)) != 1)))
          # for a random species, find the closest indel
          pos_to_test = min(unlist(lapply(unique(df$species), function(species) df$pos_align_nuc[which(df$pos_align_aa == as.character(aa_pos) & df$species == species)])), na.rm = T)[1]
          region_to_keep = unlist(lapply(1:length(region_indel), function(x) if (pos_to_test %in% region_indel[[x]]) {return(x)}))
          r = which(colnames(patch_nuc) == paste("indel_", region_indel[[region_to_keep]][1], sep = ""))
        } else {
          
          #
          if (col_id == 1){
            prev_phase = 2
            prev_aa_correction = 1
          } else {
            prev_phase = col_id - 2
            prev_aa_correction = 0
          }
          
          prev_aa_pos = paste("aa_", as.numeric(sub("_.*$", "", sub("^aa_", "", colnames(temp_nuc)[1]))) - prev_aa_correction, "_", prev_phase, sep = "")
          temp_prev_aa_pos = which(colnames(patch_nuc) == prev_aa_pos)
          if (length(temp_prev_aa_pos) == 1){
            r = temp_prev_aa_pos + 1
            
            # check if the previous aa was place according to an indel, if yes ignore the indel and place this aa next to the precedent
            if (! prev_aa_pos %in% place_accord_indel){
              # check if the prev_aa_pos is followed by an indel if new aa (ie frame = 0)
              if (prev_phase == 2) {
                if (grepl("^indel_[0-9]*$", colnames(patch_nuc)[r])){
                  # find in wich indel region we are
                  region_indel = split(as.numeric(indel), cumsum(c(1, diff(as.numeric(indel)) != 1)))
                  # check if were in two successive indel region
                  if (length(region_indel) > 1){
                    region_to_keep = unlist(lapply(1:length(region_indel), function(x) if (as.numeric(sub("^indel_", "", colnames(patch_nuc)[r])) %in% region_indel[[x]]) {return(x)}))
                    if (region_to_keep != length(region_indel)){
                      if (which(colnames(patch_nuc) == paste("indel_", tail(region_indel[[region_to_keep]], 1), sep = "")) == (which(colnames(patch_nuc) == paste("indel_", region_indel[[region_to_keep + 1]][1], sep = "")) - 1)){
                        # get the position after the indel region
                        place_accord_indel = c(place_accord_indel, colnames(temp_nuc)[col_id])
                        r = which(colnames(patch_nuc) == paste("indel_", tail(region_indel[[1]], 1), sep = "")) + 1
                      }
                    }
                  }
                }
              } 
            } else {
              place_accord_indel = c(place_accord_indel, colnames(temp_nuc)[col_id])
            } 
          }
          #
          else {
            next_aa_pos = paste("aa_", as.numeric(sub("_.*$", "", sub("^aa_", "", colnames(temp_nuc)[1]))) + 1, "_0", sep = "")
            n = 2
            while(length(which(colnames(patch_nuc) == next_aa_pos)) == 0){
              next_aa_pos = paste("aa_", as.numeric(sub("_.*$", "", sub("^aa_", "", colnames(temp_nuc)[1]))) + n, "_0", sep = "")
              n = n + 1
            }
            r = which(colnames(patch_nuc) == next_aa_pos)
          }
        }
        col = colnames(patch_nuc)
        col = c(col[1:(r - 1)], colnames(temp_nuc)[col_id], col[(r):length(col)])
        patch_nuc = insertCol(patch_nuc, temp_nuc[, col_id], r)
        colnames(patch_nuc) = col
      }
    }
  }
  
  # second scan for missing aa pos
  col_handler = colnames(patch_nuc)
  col_handler = col_handler[grep("^aa_.*", col_handler)]
  unresolved = F
  if (length(col_handler) > 0){
    col_handler = as.numeric(sub("_.*$", "", sub("^aa_", "", col_handler)))[seq(1, length(col_handler), 3)]
    aa_region = split(col_handler, cumsum(c(1, diff(col_handler) != 1)))
    if (length(aa_region) > 1){
      for (i in 1:(length(aa_region) - 1)){
        born_min = tail(aa_region[[i]], 1)
        born_max = aa_region[[i + 1]][1]
        to_add = (born_min + 1) : (born_max - 1)
        if (to_add[1] > tail(to_add, 1)){
          to_add = rev(to_add)
        }
        for (aa_pos in to_add){
          temp_nuc = c()
          for (species in unique(df$species)){
            patch = df$nuc_species[which(df$pos_align_aa == as.character(aa_pos) & df$species == species)]
            if (length(patch) == 0) {
              temp_nuc = rbind(temp_nuc, rep("-", 3))
            }
            if (length(patch) == 1) {
              if (is.na(patch)){
                temp_nuc = rbind(temp_nuc, rep("-", 3))
              } 
            } else {
              temp_nuc = rbind(temp_nuc, patch)
            }
            colnames(temp_nuc) = paste("aa", aa_pos, 0:2, sep = "_")
          }
          for (col_id in 1:ncol(temp_nuc)){
            #
            if (col_id == 1){
              prev_phase = 2
              prev_aa_correction = 1
            } else {
              prev_phase = col_id - 2
              prev_aa_correction = 0
            }
            
            prev_aa_pos = paste("aa_", as.numeric(sub("_.*$", "", sub("^aa_", "", colnames(temp_nuc)[1]))) - prev_aa_correction, "_", prev_phase, sep = "")
            temp_prev_aa_pos = which(colnames(patch_nuc) == prev_aa_pos)
            if (length(temp_prev_aa_pos) == 1){
              r = temp_prev_aa_pos + 1
            } else {
              next_aa_pos = paste("aa_", as.numeric(sub("_.*$", "", sub("^aa_", "", colnames(temp_nuc)[1]))) + 1, "_0", sep = "")
              n = 2
              max_aa = max(as.numeric(str_sub(colnames(patch_nuc)[grep("^aa_", colnames(patch_nuc))], end = -3, start = 4)))
              while((length(which(colnames(patch_nuc) == next_aa_pos)) == 0) && (! unresolved)){
                next_aa_pos = paste("aa_", as.numeric(sub("_.*$", "", sub("^aa_", "", colnames(temp_nuc)[1]))) + n, "_0", sep = "")
                n = n + 1
                ###
                if (n > max_aa){
                  unresolved = T
                }
                ###
              }
              ###
              if (unresolved) next
              ###
              r = which(colnames(patch_nuc) == next_aa_pos)
            }
            ###
            if (unresolved) next
            ###
            col = colnames(patch_nuc)
            col = c(col[1:(r - 1)], colnames(temp_nuc)[col_id], col[(r):length(col)])
            patch_nuc = insertCol(patch_nuc, temp_nuc[, col_id], r)
            colnames(patch_nuc) = col
          }
          ###
          if (unresolved) next
          ###
        }
        ###
        if (unresolved) next
        ###
      }
    }
  }
  
  align_start = min(as.numeric(df$pos_align_nuc[which(df$pos_align_aa == min(list_pos_align_aa, na.rm = T))]), na.rm = T)
  align_end = max(as.numeric(df$pos_align_nuc[which(df$pos_align_aa == max(list_pos_align_aa, na.rm = T))]), na.rm = T)
  
  if(true_start_align_nuc < align_start){
    align_start = true_start_align_nuc
  }
  if(true_end_align_nuc > align_end){
    align_end = true_end_align_nuc
  }
  
  done_region = align_end
  
  ###
  if (unresolved) {
    unresolved_region = c(unresolved_region, 
                          paste0(align_start, "-", align_end))
    next
  }
  ###
  
  path_patching = paste(align_path, "/patching_", vers, "-vs-", vers_aa, "/", sep = "")
  dir.create(path_patching, showWarnings = F, recursive = T)
  path_patching_fasta = paste(path_patching, align_start, "-", align_end, ".fasta", sep = "")
  path_patching_table = paste(path_patching, align_start, "-", align_end, ".tab", sep = "")
  print(path_patching_fasta)
  t = as.data.frame(t(patch_nuc[, -1]))
  colnames(t) = patch_nuc$species
  write.fasta(sequences = t, file.out = path_patching_fasta, names = patch_nuc$species)
  
  write.table(patch_nuc, file = path_patching_table, col.names = T, row.names = F, quote = F, sep = '\t')
}

path_unresolved_region = paste(path_patching, "unresolved_region.csv", sep = "")
write.table(unresolved_region, file = path_unresolved_region, col.names = F, row.names = F, quote = F, sep = '\t')

