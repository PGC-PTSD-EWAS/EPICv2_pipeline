#!/usr/bin/env Rscript --vanilla

################################################################################

library(tidyverse)
library(data.table)
library(ggplot2)
library(stringr)

################################################################################
## - Functions.

`%ni%` <- Negate(`%in%`)

process_pidsleycsv_for_replicates = function(PIDSLEYCSV) {
  m = data.table::fread(PIDSLEYCSV, header=T, sep=',',fill=T)
  mreplicates = m[with(m, grepl("Y", paste(namerep,seqrep,posrep))),]
  ## - Set empty cells to NA in dataframe.
  mreplicates[mreplicates == ""] = NA 
  return(mreplicates)
}

get_flagged_replicates = function(FLAGGEDCSV, mreplicates){
  flagged_mapf = read.csv(FLAGGEDCSV, header=T, sep = ",")
  flagged_probes_replicate = mreplicates[mreplicates$IlmnID %in% flagged_mapf$IlmnID,]
  flagged_probes_replicate_nochr0 = flagged_probes_replicate[flagged_probes_replicate$CHR != "chr0",]
  return(flagged_probes_replicate_nochr0)
}

detect_if_rownames = function(BETAFILE) {
  firstfive = row.names(BETAFILE)[1:5]
  detectthis = stringr::str_detect(firstfive, "cg")
  detectthis = detectthis[!is.na(detectthis)]
  d = sum(detectthis)
  if (d > 0) {
    probes_as_rows = TRUE
  } else {
    probes_as_rows = FALSE
  }
  return(probes_as_rows)
}

process_replicates = function(BETAFILE, PIDSLEYCSV, FLAGGEDCSV) {

  mreplicates = process_pidsleycsv_for_replicates(PIDSLEYCSV = PIDSLEYCSV)
  position_replicates = mreplicates[mreplicates$posrep == "Y",]
  
  flagged_probes_replicate_nochr0 = get_flagged_replicates(FLAGGEDCSV, mreplicates)
  
  ## - Split the beta file to remove the non-replicate probes. Go by IlmnID.
  b_noreplicates = BETAFILE[rownames(BETAFILE) %ni% mreplicates$IlmnID,]
  b_replicates = BETAFILE[rownames(BETAFILE) %in% mreplicates$IlmnID,]
  
  ## - Flagged probe removal (sanity check. Already done earlier in pipeline):
  b_replicates = b_replicates[rownames(b_replicates) %ni% flagged_probes_replicate_nochr0$IlmnID,]
  b_noreplicates = b_noreplicates[rownames(b_noreplicates) %ni% flagged_probes_replicate_nochr0$IlmnID,]

  manifest_replicate_filtration = mreplicates[mreplicates$IlmnID %ni% flagged_probes_replicate_nochr0$IlmnID,]
  
  ## - Chr0 probe removal (sanity check. already done earlier in pipeline):
  chr0_replicate_probes = manifest_replicate_filtration[manifest_replicate_filtration$CHR == "chr0",]
  b_replicates = b_replicates[rownames(b_replicates) %ni% chr0_replicate_probes$IlmnID,]
  b_noreplicates = b_noreplicates[rownames(b_noreplicates) %ni% chr0_replicate_probes$IlmnID,]
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$CHR != "chr0",]
  
  ## - Superior probe handling (and removal of related probes in the sets):
  filters = process_for_best_probe(manifest_replicate_filtration, "Superior probe", "Rep_results_by_LOCATION")
  b_replicates = b_replicates[rownames(b_replicates) %ni% filters,]
  
  keeps = process_for_best_probe_keeps(manifest_replicate_filtration, "Superior probe", "Rep_results_by_LOCATION")
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% filters,]
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% keeps,]
  
  ## - Best precision probe handling (and removal of related probes in the sets):
  filters = process_for_best_probe(manifest_replicate_filtration, "Best precision", "Rep_results_by_LOCATION")
  b_replicates = b_replicates[rownames(b_replicates) %ni% filters,]
  
  keeps = process_for_best_probe_keeps(manifest_replicate_filtration, "Best precision", "Rep_results_by_LOCATION")
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% filters,]
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% keeps,]
  
  ## - Best precision by group mean.
  b_replicates = process_for_groupmean_probe(manifest_replicate_filtration, "Best precision by group mean", b_replicates, "Rep_results_by_LOCATION")
  
  repfilters = filter_for_groupmean_probes(manifest_replicate_filtration, "Best precision by group mean", "Rep_results_by_LOCATION")
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% repfilters,]
  
  ## - Superior group mean (by WGBS)
  b_replicates = process_for_groupmean_probe(manifest_replicate_filtration, "Superior group mean (by WGBS)", b_replicates, "Rep_results_by_LOCATION")
  
  repfilters = filter_for_groupmean_probes(manifest_replicate_filtration, "Superior group mean (by WGBS)", "Rep_results_by_LOCATION")
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% repfilters,]
  
  ## - Superior by WGBS probe handling (and removal of related probes in the sets):
  filters = process_for_best_probe(manifest_replicate_filtration, "Superior by WGBS", "Rep_results_by_LOCATION")
  b_replicates = b_replicates[rownames(b_replicates) %ni% filters,]
  
  keeps = process_for_best_probe_keeps(manifest_replicate_filtration, "Superior by WGBS", "Rep_results_by_LOCATION")
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% filters,]
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% keeps,]
  
  ## - Best sensitivity.
  filters = process_for_best_probe(manifest_replicate_filtration, "Best sensitivity", "Rep_results_by_LOCATION")
  b_replicates = b_replicates[rownames(b_replicates) %ni% filters,]
  
  keeps = process_for_best_probe_keeps(manifest_replicate_filtration, "Best sensitivity", "Rep_results_by_LOCATION")
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% filters,]
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% keeps,]
  
  
  ## - Processing INSUFFICIENT EVIDENCE probes.
  insufficient_evidence_manifest = manifest_replicate_filtration[manifest_replicate_filtration$Rep_results_by_LOCATION == "Insufficient evidence",]
  d = as.data.frame(table(insufficient_evidence_manifest$Name))
  ddups = d[d$Freq > 1,]
  singles_Names_insufficient_evidence = d[d$Freq == 1,]
  #KEEP THESE ONES IN B_REPLICATES:
  insufficient_evidence_singles_first_pass = insufficient_evidence_manifest[insufficient_evidence_manifest$Name %in% singles_Names_insufficient_evidence$Var1,] 
  #FURTHER PROCESS THESE ONES IN B_REPLICATES:
  insufficient_evidence_duplicates = insufficient_evidence_manifest[insufficient_evidence_manifest$Name %in% ddups$Var1,]
  
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% insufficient_evidence_singles_first_pass$IlmnID,]
  
  ## - INSUFFICIENT EVIDENCE PROBES - continued: Superior probe handling (and removal of related probes in the sets):
  filters = process_for_best_probe(manifest_replicate_filtration, "Superior probe", "Rep_results_by_SEQUENCE")
  b_replicates = b_replicates[rownames(b_replicates) %ni% filters,]
  
  keeps = process_for_best_probe_keeps(manifest_replicate_filtration, "Superior probe", "Rep_results_by_SEQUENCE")
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% filters,]
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% keeps,]
  
  ## - INSUFFICIENT EVIDENCE PROBES - continued: Best precision probe handling (and removal of related probes in the sets):
  filters = process_for_best_probe(manifest_replicate_filtration, "Best precision", "Rep_results_by_SEQUENCE")
  b_replicates = b_replicates[rownames(b_replicates) %ni% filters,]
  
  keeps = process_for_best_probe_keeps(manifest_replicate_filtration, "Best precision", "Rep_results_by_SEQUENCE")
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% filters,]
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% keeps,]
  
  ## - INSUFFICIENT EVIDENCE PROBES - continued: Best precision by group mean.
  b_replicates = process_for_groupmean_probe(manifest_replicate_filtration, "Best precision by group mean", b_replicates, "Rep_results_by_SEQUENCE")
  
  repfilters = filter_for_groupmean_probes(manifest_replicate_filtration, "Best precision by group mean", "Rep_results_by_SEQUENCE")
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% repfilters,]
  
  
  ## - INSUFFICIENT EVIDENCE PROBES - continued: Superior group mean (by WGBS)
  b_replicates = process_for_groupmean_probe(manifest_replicate_filtration, "Superior group mean (by WGBS)", b_replicates, "Rep_results_by_SEQUENCE")
  
  repfilters = filter_for_groupmean_probes(manifest_replicate_filtration, "Superior group mean (by WGBS)", "Rep_results_by_SEQUENCE")
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% repfilters,]
  
  ## - INSUFFICIENT EVIDENCE PROBES - continued: Superior by WGBS probe handling (and removal of related probes in the sets):
  filters = process_for_best_probe(manifest_replicate_filtration, "Superior by WGBS", "Rep_results_by_SEQUENCE")
  b_replicates = b_replicates[rownames(manifest_replicate_filtration) %ni% filters,]
  
  keeps = process_for_best_probe_keeps(manifest_replicate_filtration, "Superior by WGBS", "Rep_results_by_SEQUENCE")
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% filters,]
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% keeps,]
  
  ## - INSUFFICIENT EVIDENCE PROBES - continued: Best sensitivity.
  filters = process_for_best_probe(manifest_replicate_filtration, "Best sensitivity", "Rep_results_by_SEQUENCE")
  b_replicates = b_replicates[rownames(b_replicates) %ni% filters,]
  
  keeps = process_for_best_probe_keeps(manifest_replicate_filtration, "Best sensitivity", "Rep_results_by_SEQUENCE")
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% filters,]
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% keeps,]
  
  ## - Offtargets
  offtargets = manifest_replicate_filtration[manifest_replicate_filtration$Num_offtargets != 0,]
  
  no_offtargets = manifest_replicate_filtration[manifest_replicate_filtration$Num_offtargets == 0,] # Keep these
  
  b_replicates = b_replicates[rownames(b_replicates) %ni% offtargets$IlmnID,]
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% offtargets$IlmnID,]
  
  ## - INSUFFICIENT EVIDENCE PROBES - continued: RMSE values:
  use_b = b_replicates
  use_b = process_for_rmse(use_b, manifest_replicate_filtration)
  filters = process_for_rmse_filters(use_b, manifest_replicate_filtration)
  keeps = process_for_rmse_keeps(use_b, manifest_replicate_filtration)
  use_b = use_b[rownames(use_b) %ni% filters,]
  
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% filters,]
  manifest_replicate_filtration = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %ni% keeps,]
  
  ## - INSUFFICIENT EVIDENCE PROBES - continued: final group mean of remainder:
  use_b = obtain_groupmean_remainder_probes_final(manifest_replicate_filtration, use_b)
  
  ## - Concatenate the nonreplicate dataframe with the replicate dataframe.
  b_noreplicates = as.data.frame(b_noreplicates)
  beta_final = rbind(b_noreplicates, use_b)
  return(beta_final)
}

process_for_rmse = function(b_replicates, manifest_replicate_filtration) {
  beta_df = b_replicates
  j = 0
  keeps = c()
  filters = c()
  groupmean_filters = c()
  for (i in manifest_replicate_filtration$posrep_IlmnIDs){
    replicate_set = unlist(stringr::str_split(i, ";"))
    to_process_sup = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %in% replicate_set,]
    d = as.data.frame(table(to_process_sup$RMSE_with_WGBS))
    if (dim(d)[1] == 1) {
      
      dx = as.numeric(as.character(d[1,1]))
      keep_this_one = to_process_sup[to_process_sup$RMSE_with_WGBS == dx,]$IlmnID
      keeps = c(keeps, keep_this_one)
      filter_buddies = to_process_sup[to_process_sup$IlmnID %ni% keep_this_one,]
      filter_buddies_names = filter_buddies$IlmnID
      filters = c(filters, filter_buddies_names)
      
    } else {
      if (dim(to_process_sup)[1] > 1) {
        b_subset = beta_df[rownames(beta_df) %in% to_process_sup$IlmnID,]
        if (dim(b_subset)[1] != 0) {
          xtmp = as.data.frame(t(apply(b_subset, MARGIN=2, FUN=mean)))
          beta_df = beta_df[rownames(beta_df) %ni% to_process_sup$IlmnID,]
          prefix_name = unique(to_process_sup$Name)
          rownames(xtmp) = prefix_name
          beta_df = rbind(beta_df, xtmp)
          group_probes_unique = to_process_sup$IlmnID
          groupmean_filters = c(groupmean_filters, group_probes_unique)
        }
      }
    }
  }
  beta_df = beta_df[rownames(beta_df) %ni% filters,]
  return(beta_df)
}

process_for_rmse_filters = function(b_replicates, manifest_replicate_filtration) {
  beta_df = b_replicates
  j = 0
  keeps = c()
  filters = c()
  groupmean_filters = c()
  for (i in manifest_replicate_filtration$posrep_IlmnIDs){
    replicate_set = unlist(stringr::str_split(i, ";"))
    to_process_sup = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %in% replicate_set,]
    d = as.data.frame(table(to_process_sup$RMSE_with_WGBS))
    if (dim(d)[1] == 1) {
      
      dx = as.numeric(as.character(d[1,1]))
      keep_this_one = to_process_sup[to_process_sup$RMSE_with_WGBS == dx,]$IlmnID
      keeps = c(keeps, keep_this_one)
      filter_buddies = to_process_sup[to_process_sup$IlmnID %ni% keep_this_one,]
      filter_buddies_names = filter_buddies$IlmnID
      filters = c(filters, filter_buddies_names)
      
    } else {
      if (dim(to_process_sup)[1] > 1) {
        b_subset = beta_df[rownames(beta_df) %in% to_process_sup$IlmnID,]
        if (dim(b_subset)[1] != 0) {
          xtmp = as.data.frame(t(apply(b_subset, MARGIN=2, FUN=mean)))
          beta_df = beta_df[rownames(beta_df) %ni% to_process_sup$IlmnID,]
          prefix_name = unique(to_process_sup$Name)
          rownames(xtmp) = prefix_name
          beta_df = rbind(beta_df, xtmp)
          group_probes_unique = to_process_sup$IlmnID
          filters = c(filters, group_probes_unique)
        }
      }
    }
  }
  beta_df = beta_df[rownames(beta_df) %ni% filters,]
  return(filters)
}

process_for_rmse_keeps= function(b_replicates, manifest_replicate_filtration) {
  beta_df = b_replicates
  j = 0
  keeps = c()
  filters = c()
  groupmean_filters = c()
  for (i in manifest_replicate_filtration$posrep_IlmnIDs){
    replicate_set = unlist(stringr::str_split(i, ";"))
    to_process_sup = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %in% replicate_set,]
    d = as.data.frame(table(to_process_sup$RMSE_with_WGBS))
    if (dim(d)[1] == 1) {
      
      dx = as.numeric(as.character(d[1,1]))
      keep_this_one = to_process_sup[to_process_sup$RMSE_with_WGBS == dx,]$IlmnID
      keeps = c(keeps, keep_this_one)
      filter_buddies = to_process_sup[to_process_sup$IlmnID %ni% keep_this_one,]
      filter_buddies_names = filter_buddies$IlmnID
      filters = c(filters, filter_buddies_names)
      
    } else {
      if (dim(to_process_sup)[1] > 1) {
        b_subset = beta_df[rownames(beta_df) %in% to_process_sup$IlmnID,]
        if (dim(b_subset)[1] != 0) {
          xtmp = as.data.frame(t(apply(b_subset, MARGIN=2, FUN=mean)))
          beta_df = beta_df[rownames(beta_df) %ni% to_process_sup$IlmnID,]
          prefix_name = unique(to_process_sup$Name)
          rownames(xtmp) = prefix_name
          beta_df = rbind(beta_df, xtmp)
          group_probes_unique = to_process_sup$IlmnID
          filters = c(filters, group_probes_unique)
        }
      }
    }
  }
  beta_df = beta_df[rownames(beta_df) %ni% filters,]
  return(keeps)
}

process_for_best_probe = function(mreplicates, rep_target, column_name){
  # Used for superior probe.
  # Used for superior by WGBS probe.
  # Used for Best precision probe.
  position_replicates = mreplicates[mreplicates$posrep == "Y",]
  sup_probes = position_replicates[position_replicates[[column_name]] == rep_target,]
  j = 0
  keeps = c("character", length(sup_probes[[column_name]]))
  filters = c()
  for (i in sup_probes$posrep_IlmnIDs){
    j = j + 1
    replicate_set = unlist(stringr::str_split(i, ";"))
    to_process_sups = mreplicates[mreplicates$IlmnID %in% replicate_set,]
    ## - Check that only 1 superior per set.
    d = as.data.frame(table(to_process_sups[[column_name]]))
    if (d[d$Var1 == rep_target,]$Freq != 1) {
      stop("ERROR: Number of target probes in a set is > 1. Something is wrong. Exiting.")
    } else {
      keep_this_one = to_process_sups[to_process_sups[[column_name]] == rep_target,]$IlmnID
      keeps[j] = keep_this_one
      filter_buddies = to_process_sups[to_process_sups[[column_name]] != rep_target,]
      filter_buddies_names = filter_buddies$IlmnID
      filters = c(filters, filter_buddies_names)
    }
  }
  return(filters)
}

process_for_best_probe_keeps = function(mreplicates, rep_target, column_name){
  # Used for superior probe.
  # Used for superior by WGBS probe.
  # Used for Best precision probe.
  position_replicates = mreplicates[mreplicates$posrep == "Y",]
  sup_probes = position_replicates[position_replicates[[column_name]] == rep_target,]
  j = 0
  keeps = c("character", length(sup_probes[[column_name]]))
  filters = c()
  for (i in sup_probes$posrep_IlmnIDs){
    j = j + 1
    replicate_set = unlist(stringr::str_split(i, ";"))
    to_process_sups = mreplicates[mreplicates$IlmnID %in% replicate_set,]
    ## - Check that only 1 superior per set.
    d = as.data.frame(table(to_process_sups[[column_name]]))
    if (d[d$Var1 == rep_target,]$Freq != 1) {
      stop("ERROR: Number of target probes in a set is > 1. Something is wrong. Exiting.")
    } else {
      keep_this_one = to_process_sups[to_process_sups[[column_name]] == rep_target,]$IlmnID
      keeps[j] = keep_this_one
      filter_buddies = to_process_sups[to_process_sups[[column_name]] != rep_target,]
      filter_buddies_names = filter_buddies$IlmnID
      filters = c(filters, filter_buddies_names)
    }
  }
  return(keeps)
}


process_for_groupmean_probe = function(mreplicates, rep_target, beta_df, column_name){
  position_replicates = mreplicates[mreplicates$posrep == "Y",]
  group_probes = position_replicates[position_replicates[[column_name]] == rep_target,]
  
  j = 0
  group_probes_unique = unique(group_probes$posrep_IlmnIDs)
  for (i in group_probes_unique) {
    j = j+1
    replicate_set = unlist(stringr::str_split(i, ";"))
    to_process = mreplicates[mreplicates$IlmnID %in% replicate_set,]
    b_subset = beta_df[rownames(beta_df) %in% to_process$IlmnID,]
    
    if (dim(b_subset)[1] != 0) {
      xtmp = as.data.frame(t(apply((b_subset), MARGIN=2, FUN=mean)))
      beta_df = beta_df[rownames(beta_df) %ni% to_process$IlmnID,]
      prefix_name = unique(to_process$Name)
      rownames(xtmp) = prefix_name
      beta_df = rbind(beta_df, xtmp)
    }
  }
  return(beta_df)
}

filter_for_groupmean_probes = function(mreplicates, rep_target, column_name){
  position_replicates = mreplicates[mreplicates$posrep == "Y",]
  group_probes = position_replicates[position_replicates[[column_name]] == rep_target,]
  group_probes_unique = unique(group_probes$posrep_IlmnIDs)
  group_probes_all = unlist(stringr::str_split(group_probes_unique, ";"))
  return(group_probes_all)
}


obtain_groupmean_remainder_probes_final = function(manifest_replicate_filtration, beta_df){
  j = 0
  group_probes_unique = unique(manifest_replicate_filtration$posrep_IlmnIDs)
  for (i in group_probes_unique) {
    j = j+1
    replicate_set = unlist(stringr::str_split(i, ";"))
    to_process = manifest_replicate_filtration[manifest_replicate_filtration$IlmnID %in% replicate_set,]
    b_subset = beta_df[rownames(beta_df) %in% to_process$IlmnID,]
    
    if (dim(b_subset)[1] != 0) {
      xtmp = as.data.frame(t(apply((b_subset), MARGIN=2, FUN=mean)))
      beta_df = beta_df[rownames(beta_df) %ni% to_process$IlmnID,]
      prefix_name = unique(to_process$Name)
      rownames(xtmp) = prefix_name
      beta_df = rbind(beta_df, xtmp)
    }
  }
  return(beta_df)
}

## End. ##