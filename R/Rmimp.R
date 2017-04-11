BASE_DIR = system.file("extdata", "", package = "rmimp")
library(parallel)
library(mclust)
library(ROCR)

require(GenomicRanges)
require(data.table)
require(Biostrings)

if(F){
  setwd('~/Development/rmimp/')
  BASE_DIR = '~/Development/rmimp/inst/extdata/'
  setwd('~/Desktop/rmimp/')
  BASE_DIR = '~/Desktop/rmimp/inst/extdata/'
  source('R/display-functions.r')
  source('R/io-functions.r')
  source('R/pwm-functions.r')
  writeLines('Warning: remove base dir')
}

# Current version
.MIMP_VERSION = '1.2'

# Valid domains and corresonding descriptions
.VALID_DOMAINS <- list("central" = "phos", "not_central" = "sh3")
.DESCRIPTIONS <- list("phos" = "predicted kinase rewiring", "sh3" = "sh3 binding")

# This message appears on library or require call of package
.onAttach <- function(lib, pkg, ...) {
  packageStartupMessage(sprintf("MIMP v%s (%s)
                                Type '?mimp' for help or see the documentation 'help(package=rmimp)' for more details

                                If you use MIMP in your research, please cite:
                                Wagih O, Reimand J, Bader GD (2015). MIMP: predicting the impact of mutations on kinase-substrate phosphorylation. Nat. Methods 12(6):531-3. doi:10.1038/nmeth.3396", .MIMP_VERSION, format(Sys.Date(), '%Y')))
}

#' Converts all columns of a data frame of class factor to character
#'
#' @param string String to be manipulated
#' @keywords internal factor character
#' 
#' @examples
#' unfactor( data.frame(x=c('A', 'B')) )
unfactor <- function(df){
  i <- sapply(df, is.factor)
  df[i] <- lapply(df[i], as.character)
  df
}

#' Get mixture model parameters
#' @return means, sds, wts
.getParams = function(clust) {
  # Get weights
  nr = 1:length(clust[['parameters']]$pro)
  wts = clust[['parameters']]$pro

  # Get means
  means = clust[['parameters']]$mean

  # Get stdev, need to square root it to get standard deviation
  sds = clust$parameters$variance$sigmasq
  sds = sqrt(sds)

  # Return data frame
  data.frame(nr, wts, means, sds)
}

#' Get parameters for foreground and background distros
#' and return them as a list to be used later
.fgBgParams <- function(pos, neg){
  # Get gmms
  fg_gmm = suppressWarnings(Mclust(pos))
  bg_gmm = suppressWarnings(Mclust(neg))

  # and parameters for the models
  fg_params = .getParams(fg_gmm)
  bg_params = .getParams(bg_gmm)

  list(fg_params=fg_params, bg_params=bg_params)
}

#' Calculate AUC based on positive and negitive seqeunces of binding sites
#' @param posSeqPath Path to a directory containing all positive sequences of binding sites (one file per site)
#' @param negSeqFile Path to a file containing all negative sequences in a table format. See negSeq_example.txt in inst/extdata
#' @param cores Number of CPU cores that the function can use for calculation
#' @return AUCs calculated
.calculateAUC <- function(posSeqPath, negSeqFile, cores = 2) {
  # Get positive and negative seqs of binding sites
  posSeqsFiles <- list.files(posSeqPath, full.names = T)
  negSeqs <- read.table(negSeqFile, header = T, sep = "\t", check.names = F)
  
  # Score AUCs
  aucs <- mclapply(posSeqsFiles, function(fileName) {
    # Get positive sequences and binding site
    posSeqs <- unique(readLines(fileName))
    bindingSite <- unlist(strsplit(fileName, "/", fixed = T))
    bindingSite <- bindingSite[length(bindingSite)]
    
    # Find length of positive sequences
    seqLength <- unique(nchar(posSeqs))
    if (length(posSeqs) <= 1) {
      warning(sprintf("length of the positive sequences for binding site %s must be more than 1!", bindingSite))
      return(0)
    } else if (length(seqLength) != 1) {
      warning(sprintf("length of the positive sequences for binding site %s cannot be determined. There are %d lengths - %s",
                   bindingSite, length(seqLength), toString(seqLength)))
      return(0)
    }
    
    # Get negitive sequences with the same length as positive sequences
    if (!seqLength %in% as.integer(colnames(negSeqs))) {
      warning(sprintf("no matching negative sequences of binding site %s was found with the same length as of positive sequences.",
                      bindingSite))
      return(0)
    }
    selectedNegSeqs <- unique(as.character(negSeqs[[as.character(seqLength)]]))
    
    # Calcuate the average AUC
    auc <- mean(sapply(suppressWarnings(split(posSeqs, 1:ifelse(length(posSeqs) >= 10, 10, length(posSeqs)))), function(testset) {
      # Use the remaining sequences to generate PWM
      remaining <- setdiff(posSeqs, testset)
      pwm <- PWM(remaining, is.kinase.pwm = F)
      
      # Score the testset and negative seqs wit the PWM
      posScores <- mss(testset, pwm)
      negScores <- mss(selectedNegSeqs, pwm)
      
      # Calculate and return the AUC
      pred <- prediction(c(posScores, negScores), c(rep(1, length(posScores)), rep(0, length(negScores))))
      perf <- performance(pred, "auc")
      return(unlist(perf@y.values))
    }, simplify = T, USE.NAMES = F))

    return(auc)
  }, mc.cores = cores)
  
  # Get names of binding sites and assign to aucs
  names(aucs) <- list.files(posSeqPath)

  
  return(aucs)
}

#' Train GMM model
#' and return as a list to be used later. If file is passed,
#' the model will also be save to a .mimp file.
#' @param posScorePath the path to the directory contains positive scores
#' @param negScorePath the path to the directory contains negative scores
#' @param calculateAUC (optional) whether to calculate AUC for the model. The posSeqPath and negSeqPath are required for AUC calculation.
#' @param posSeqPath (optional) Path to a directory containing all positive sequences of binding sites (one file per site)
#' @param negSeqPath (optional) Path to a file containing all negative sequences in a table format. See negSeq_example.txt in inst/extdata
#' @param pwmPath the path to the directory contains PWMs
#' @param cores (optional) the number of CPU cores that can be used to train the model
#' @param file (optional) the path to save the model
#' @param threshold (optional) the minimum number of scores needed for each domain to train the model
#' @param minAUC (optional) the minimum number of AUC needed for each domain to train the model
#' @return a GMM model
#' @examples
#' No examples
#' @export
trainModel <- function(posScorePath, negScorePath, calculateAUC = F,
                       posSeqPath = NULL, negSeqPath = NULL, pwmPath, cores = 2, file = NULL, threshold = 10, minAUC = 0.65){
  # Get a list of files of scores and pwms from posScorePath, negScorePath and pwmPath
  # and check if the corresponding files exists
  fileNames <- intersect(list.files(path = c(posScorePath), full.names = F),
                         intersect(list.files(path = c(negScorePath), full.names = F), list.files(path = c(pwmPath), full.names = F)))
  allFiles <- unique(list.files(path = c(posScorePath, negScorePath, pwmPath), full.names = F))
  missingFiles <- which(!is.element(allFiles, fileNames))

  if (length(missingFiles) != 0) {
    warning(sprintf("A total of %d files are not found in all places(positive, negative and pwm): %s",
                    length(missingFiles), paste0(fileNames[missingFiles], collapse = ", ")))
  }

  cat("\rTraining model. \n")
  
  if(calculateAUC) {
    # Calculate AUCs
    aucs <- .calculateAUC(posSeqPath, negSeqPath, cores)
    names(aucs) <- unlist(strsplit(names(aucs), ".txt"))
  }
  
  # Get parameters for each pair of scores
  params <- mclapply(fileNames, function(fileName) {
    pos <- read.table(paste0(posScorePath, "/", fileName), header = F)
    neg <- read.table(paste0(negScorePath, "/", fileName), header = F)
    pwm <- read.table(paste0(pwmPath, "/", fileName), skip = 1)
    row.names(pwm) <- pwm[, 1]
    pwm <- pwm[, 2:ncol(pwm)]

    # Check if pwm is valid
    # Will still proceed if some rows are missing, but a warning message will be thrown.
    if (nrow(pwm) != length(AA)) {
      warning(sprintf("PWM of %s might be missing some AAs.", fileName))
    }

    # Check if both score tables are numeric
    if (!all(c(class(pos[1,]), class(neg[1,])) == "numeric")) {
      warning(sprintf("File %s is not entirely numeric. It is likely caused by the header in the score file. Please remove it and try again.", fileName))
      return(NULL)
    }

    # Check if both pos and neg have more scores than THRESHOLD
    if (any(c(nrow(pos), nrow(neg)) < threshold)) {
      warning(fileName, " skipped as at least ", threshold, " scores required for trainning.")
      return(NULL)
    }

    # Check if any prediction could be made,
    # and if not, this file will be skipped.
    if (all(is.na(mclustBIC(data = unlist(pos, use.names = F))))) {
      warning(fileName, " skipped as no prediction could be made.")
      return(NULL)
    } else {
      ret <- .fgBgParams(pos, neg)
      ret$pwm <- pwm
      
      # Check if any matching AUC
      # and attach the correct AUC to the model
      if (calculateAUC) {
        modelName <- unlist(strsplit(fileName, "--"))[1]
        if (!modelName %in% names(aucs)) {
          warning(fileName, "does not have corresponding AUC.")
          return(NULL)
        }
        if (as.double(aucs[[modelName]]) < minAUC) {
          warning(sprintf("AUC calculated for %s is %f - smaller than minAUC %f", fileName, as.double(aucs[[modelName]]), minAUC))
          return(NULL)
        }
        ret$auc <- aucs[[modelName]]
      }
      return(ret)
    }
  }, mc.cores = cores, mc.silent = F)
  seqNames <- unlist(strsplit(fileNames, ".txt", fixed = T))
  names(params) <- seqNames
  params <- params[!sapply(params, is.null)]
  
  cat("\rdone | Training model. \n")

  if(!is.null(file)) {
    saveRDS(params, file = file)
  }

  return(params)
}

#' Get flanking sequences of a position.
#'
#' This function obtains the flanking sequence at one or more position.
#' Out of bound indices are replaced by a blank character.
#'
#' @param seqs Character vector of sequences. If only one sequence is provided,
#'    indices from \code{inds} are assumed to all be from the same sequence.
#' @param inds Numerical vector of positions corresponding to the sequences provided in \code{seqs}.
#' @param flank Value indicating the number of characters to extract, before and after an index
#' @param empty_char Character used to replace out of bound flanking sequences
#' @keywords internal flank sequence
#' 
#' @examples
#' # One sequence and one index. Central character is 'B'
#' flankingSequence(seqs='ABC', inds=2, flank=1)
#' # An example showing the use of empty.char
#' flankingSequence(seqs='ABC', inds=2, flank=5)
#' # An example with multiple sequences and indices
#' flankingSequence(seqs=c('ABC', 'XYZ'), inds=c(2, 1), flank=1)
flankingSequence <- function(seqs, inds, flank=7, empty_char='-'){
  if(length(seqs) == 1 & length(inds) >= length(seqs)) seqs = rep(seqs, length(inds))
  if(length(seqs) != length(inds)) stop('Length of sequences must be equal to length of positions')
  
  border = paste0(rep(empty_char, flank), collapse="")
  seqs = paste0(border, unlist(seqs, F, F), border)
  substr(seqs, inds, inds+(flank*2))
}

#' Find sh3 related variants (sSNVs)
#'
#' Given mutation data, find variants that
#' exist in the flanking regions of the psite
#'
#' @param md Mutation data as data frame of two columns (1) name of gene or protein
#'    (2) mutation in the format X123Y, where X is the reference amino acid
#'    and Y is the alternative amino acid.
#' @param seqdata Phosphorylation data as a data frame of two columns (1) name of gene
#'    or protein (2) Position of the phosphorylated residue
#' @param flank Number of amino acids flanking the psite to be considered
#'
#' @keywords ssnv mutation snp
#' @export
#' @examples
#' # No examples
sSNVs <- function(md, seqdata, flank) {

  # Get start, end and mutation positions
  flank = flank - 1
  mutNames <- md$gene
  seqlengths <- nchar(seqdata[mutNames])
  names(seqlengths) <- NULL
  startPos <- pmax(1, md$mut_pos - flank)
  endPos <- pmin(seqlengths, md$mut_pos + flank)

  # Extract sequences based on positions
  seqs <- mapply(function(name, start, end, mut, mutAA) {
    wt.seq <- substr(seqdata[[name]], start, end)
    end <- end - start + 1
    mut <- mut - start + 1
    start <- 1
    mt.seq <- paste0(substr(wt.seq, start, mut - 1), mutAA, substr(wt.seq, mut + 1, end), collapse = "")
    return(c(wt.seq, mt.seq))
  }, mutNames, startPos, endPos, md$mut_pos, md$alt_aa)

  return(seqs)
}

#' Find phosphorylation related variants (pSNVs)
#'
#' Given mutation data and psites, find variants that
#' exist in the flanking regions of the psite
#'
#' @param muts Mutation data as data frame of two columns (1) name of gene or protein
#'    (2) mutation in the format X123Y, where X is the reference amino acid
#'    and Y is the alternative amino acid.
#' @param psites Phosphorylation data as a data frame of two columns (1) name of gene
#'    or protein (2) Position of the phosphorylated residue
#' @param seqs Sequence data as a name list. Names of the list correspond to the gene or
#'    protein name. Each entry contains the collapsed sequence.
#' @param flank Number of amino acids flanking the psite to be considered
#' 
#' @keywords internal psnv psites mutation snp
#' 
#' @import GenomicRanges
#' @importFrom IRanges IRanges
#' @examples
#' # No examples
pSNVs <- function(md, pd, seqdata, flank=7){
  
  md2 = md; pd2 = pd
  
  # Create psite and mutation Granges
  ps_gr = with(pd, GRanges(gene, IRanges(pos-flank, pos+flank))) 
  mt_gr = with(md, GRanges(gene, IRanges(mut_pos, mut_pos))) 
  
  # Find overlaps of mutations and psites
  ol = findOverlaps(mt_gr, ps_gr)
  
  if(!is.null(seqdata)){
    # Get flanking sequence from fasta file
    pd2$mt = pd2$wt = flankingSequence(seqdata[pd2$gene], pd2$pos, flank)
  }else{
    # We're processing psp data
    pd2$mt = pd2$wt = pd2$seq
  }
  
  df = cbind(pd2[subjectHits(ol),], md[queryHits(ol),])
  df$start = start(ps_gr[subjectHits(ol)])
  df$psite_pos = df$pos
  df$mut_dist = df$mut_pos - df$pos
  
  # Mutate mt sequence
  df$rel_ind = df$mut_pos - df$start + 1
  
  # Don't keep anything that doesnt match amino acid provided by user and amino acid in sequence
  seq_ref = base::substr(df$mt, df$rel_ind, df$rel_ind)
  df = subset(df, ref_aa == seq_ref)
  
  # Mutate
  base::substr(df$mt, df$rel_ind, df$rel_ind) = df$alt_aa
  
  
  if(is.null(seqdata)) df$gene = sprintf('%s (%s)', df$gene, df$symbol)
  
  df = df[,!names(df) %in% c('pos', 'start', 'end', 'gene.1', 'seq', 'symbol', 'rel_ind')]
  rownames(df) = NULL
  df
}

#' Compute likelihood using probability density function
#'
#' @param p scores
#' @param params parameters of GMM from getParams
#' @keywords internal
.probPoint = function(p, params) {
  # For all points, get likeli hood for each component
  component_vals = apply(params, 1, function(r){
    dnorm(p, mean = r['means'], sd = r['sds'], log=F) * r['wts']
  })

  if(length(p) == 1) return(sum(component_vals))
  # Summ accross all components
  apply(component_vals, 1, sum)
}

#' Computing posterior probability - ploss and pgain
#'
#' @param wt.scores Wild type score
#' @param mt.scores Mutant score
#' @param fg.params Distribution parameters of GMMs (foreground). This is precomputed and comes built into mimp.
#' @param bg.params Distribution parameters of GMMs (background). This is precomputed and comes built into mimp.
#' @param auc AUC of the model. This is precomputed and comes built into mimp.
#' @param intermediate If TRUE, intermediate likelihoods used to compute ploss and pgain is returned. Otherwise only ploss and pgain returned
#' @keywords internal
pRewiringPosterior <- function(wt.scores, mt.scores, fg.params, bg.params, auc=1, intermediate=F){

  fg_prior = auc/(1+auc)
  bg_prior = 1-fg_prior

  # Get probability of wt in fg/bg
  l.wt.fg = .probPoint(wt.scores, fg.params) * fg_prior
  l.wt.bg = .probPoint(wt.scores, bg.params) * bg_prior

  # Get wt posteriors
  post.wt.fg = (l.wt.fg) / (l.wt.bg + l.wt.fg)
  post.wt.bg = (l.wt.bg) / (l.wt.bg + l.wt.fg)

  if(missing(mt.scores))
    return(data.frame(l.wt.fg, l.wt.bg, post.wt.fg, post.wt.bg))


  # Get probability of mt in fg/bg
  l.mt.fg = .probPoint(mt.scores, fg.params) * fg_prior
  l.mt.bg = .probPoint(mt.scores, bg.params) * bg_prior

  # Get mt posteriors
  post.mt.fg = (l.mt.fg) / (l.mt.bg + l.mt.fg)
  post.mt.bg = (l.mt.bg) / (l.mt.bg + l.mt.fg)

  post.wt.fg[is.na(post.wt.fg)] = 0
  post.mt.bg[is.na(post.mt.bg)] = 1
  post.wt.bg[is.na(post.wt.bg)] = 1
  post.mt.fg[is.na(post.mt.fg)] = 0

  # Get probabilities of loss and gain
  ploss = post.wt.fg * post.mt.bg
  pgain = post.wt.bg * post.mt.fg

  if(intermediate){
    return(data.frame(l.wt.fg, l.mt.fg, l.wt.bg, l.mt.bg,
                      post.wt.fg, post.mt.fg, post.wt.bg, post.mt.bg,
                      ploss, pgain))
  }
  data.frame(ploss, pgain)
}


#' Score wt and mt sequences for a pwm
#'
#' @param obj MIMP kinase object containing PWM, auc, GMM parameters, family name, etc.
#' @param mut_ps psnvs data frame containing wt and mt sequences computed from pSNVs function
#' @param prob.thresh Probability threshold of gains and losses. This value should be between 0.5 and 1.
#' @param log2.thresh Threshold for the absolute value of log ratio between wild type and mutant scores. Anything less than this value is discarded (default: 1).
#' @param include.cent If TRUE, gains and losses caused by mutation in the central STY residue are kept
#' 
#' @keywords internal mut psites score snp snv
computeRewiring <- function(obj, mut_ps, prob.thresh=0.5, log2.thresh=1, include.cent=F, 
                            degenerate.pwms=F, .degenerate.groups=c('DE','KR','ILMV')){
  
  # Score wild type and mutant sequences
  pwm = obj$pwm

  # make PWM degenerate if needed
  if(degenerate.pwms) pwm = degeneratePWM(pwm, .degenerate.groups)
  mut_ps$score_wt = mss(mut_ps$wt, pwm)
  mut_ps$score_mt = mss(mut_ps$mt, pwm)
  
  # Remove cases where wt and mt are NA
  #mut_ps = mut_ps[!(is.na(mut_ps$score_wt) & is.na(mut_ps$score_mt)),]
  mut_ps = subset(mut_ps, !(is.na(score_wt) & is.na(score_mt)) )
  
  # Do we have any data left?
  if(nrow(mut_ps) == 0) return(NULL)
  
  # Compute the log2 ratio
  mut_ps$log_ratio = with(mut_ps, log2(score_mt/score_wt) )
  
  # Set the kinase name, family and number of sequences used to build the model
  mut_ps$pwm = obj$name
  mut_ps$pwm_fam = obj$family
  mut_ps$nseqs = attr(obj$pwm, 'nseqs')
  
  res = pRewiringPosterior(mut_ps$score_wt, mut_ps$score_mt, obj$fg.params, obj$bg.params, obj$auc)

  # Filter by probability threshold
  ind = with(res, ploss >= prob.thresh | pgain >= prob.thresh)
  
  # Nothing passes the prob threshold, return NULL
  if(sum(ind) == 0) return(NULL)

  # Join with main data frame
  mut_ps = cbind(mut_ps[ind,], res[ind,])

  # Filter by log2 ratio, things with NAs have no log ratio
  mut_ps = subset(mut_ps, is.na(log_ratio) | abs(log_ratio) >= log2.thresh)
  
  # Nothing passes the log2 ratio, return null
  if(nrow(mut_ps) == 0) return(NULL)

  # Set the effect
  effect = rep('gain', nrow(mut_ps))
  effect[mut_ps$ploss >= prob.thresh] = 'loss'
  
  # Set prob as maximum of ploss and pgain
  mut_ps$prob = pmax(mut_ps$ploss, mut_ps$pgain)
  mut_ps$effect = effect

  # Remove ploss pgain. Retain only 'prob'
  mut_ps = mut_ps[,!names(mut_ps) %in% c('ploss', 'pgain')]

  return(mut_ps)
}

#' Score wt and mt sequences for a pwm
#'
#' @param obj MIMP object containing PWM, GMM parameters, and etc.
#' @param mut_ss ssnvs data frame containing wt and mt sequences computed from sSNVs function
#' @param mut_location list of mutation locations
#' @param prob.thresh Probability threshold of gains and losses. This value should be between 0.5 and 1.
#' @param log2.thresh Threshold for the absolute value of log ratio between wild type and mutant scores. Anything less than this value is discarded (default: 1).
#'
#' @keywords mut sh3 score snp snv
computeSh3Binding <- function(obj, mut_ss, mut_location, prob.thresh = 0.5, log2.thresh = 1) {
  # Get wt and mut sequences and pwms
  pwm <- as.matrix(obj$pwm)
  pwm_name <- obj$pwm_name
  mutNames <- colnames(mut_ss)
  wtSeqs <- as.character(mut_ss[1,])
  mutSeqs <- as.character(mut_ss[2,])
  
  # Remove cases where wt and mt are damaged
  pwm_ncol <- ncol(pwm)
  invalidCases <- nchar(wtSeqs) < pwm_ncol | nchar(mutSeqs) < pwm_ncol
  mutNames <- mutNames[!invalidCases]
  mut_location <- mut_location[!invalidCases]
  wtSeqs <- wtSeqs[!invalidCases]
  mutSeqs <- mutSeqs[!invalidCases]
  
  # Check if there is no data left
  if(length(which(!invalidCases)) == 0) return(NULL)
  
  # Score wt and mut sequences
  wtScores <- mssSh3(wtSeqs, pwm)
  mutScores <- mssSh3(mutSeqs, pwm)
  
  # Predict
  predicted <- mapply(function(wt, mut, name, loc) {
    res = pRewiringPosterior(wt, mut, obj$fg_params, obj$bg_params, obj$auc)
    log2 <- log2(wt/mut)
    # Filter by probability threshold
    ind = res$ploss >= prob.thresh | res$pgain >= prob.thresh
    # Filter by log2 ratio, things with NAs have no log ratio
    ind = ind & (is.na(log2) | abs(log2) >= log2.thresh)
    # Return NULL if nothing passes the filters
    if(!any(ind)) {
      return(NULL)
    }
    # Find the most significant one
    log2[!ind] = 0
    max <- which.max(abs(log2))
    return(data.frame(name, loc, names(wt)[max], names(mut)[max], wt[[max]], mut[[max]], log2[[max]], res[max, "ploss"], res[max, "pgain"]))
  }, wtScores, mutScores, mutNames, mut_location, SIMPLIFY = T, USE.NAMES = F)
  
  # Return NULL if no prediction could be made
  if (all(sapply(predicted, is.null))) {
    return(NULL)
  }
  
  # Format the predicted
  predicted <- do.call("rbind", predicted[!sapply(predicted, is.null)])
  predicted <- cbind(predicted, pwm_name)
  colnames(predicted) <- c("gene", "mut", "wt", "mt", "score_wt", "score_mt", "log_ratio", "ploss", "pgain", "pwm")
  
  # Set the effect
  effect = rep('gain', nrow(predicted))
  effect[predicted$ploss >= prob.thresh] = 'loss'
  # Set prob as maximum of ploss and pgain
  predicted$prob = pmax(predicted$ploss, predicted$pgain)
  predicted$effect = effect
  
  # Remove ploss pgain. Retain only 'prob'
  predicted = predicted[,!names(predicted) %in% c('ploss', 'pgain')]
  return(predicted)
}

#' Predict the impact of single variants on phosphorylation.
#'
#' This function takes in mutation, sequence and phosphorylation data to predict the
#' impact the mutation has on phosphorylation.
#'
#' @param muts Mutation data file: a space delimited text file OR data frame containing two columns (1) gene and (1) mutation.
#' Example:
#' \tabular{ll}{
#'    TP53 \tab R282W\cr
#'    CTNNB1 \tab S33C\cr
#'    CTNNB1 \tab S37F\cr
#' }
#' @param seqs Sequence data file containing protein sequences in FASTA format OR
#'  named list of sequences where each list element is the uppercase sequence and the
#'  name of each element is that of the protein. Example: list(GENEA="ARNDGH", GENEB="YVRRHS")
#' @param central Whether the mutation site is at the central residue of the sequence
#' @param domain Which binding domain to run mimp for
#' @param psites Phosphorylation data file (optional): a space delimited text file OR data frame containing  two columns (1) gene and (1) positions of phosphorylation sites. Example:
#' \tabular{ll}{
#'    TP53 \tab 280\cr
#'    CTNNB1 \tab 29\cr
#'    CTNNB1 \tab 44\cr
#' }
#' @param pdz.range The number of amino acids used for predicting PDZ binding.
#' @param prob.thresh Probability threshold of gains and losses. This value should be between 0.5 and 1.
#' @param log2.thresh Threshold for the absolute value of log ratio between wild type and mutant scores. Anything less than this value is discarded (default: 1).
#' @param display.results If TRUE results are visualised in an html document after analysis is complete
#' @param include.cent If TRUE, gains and losses caused by mutation in the central STY residue are kept. Scores of peptides with a non-STY central residue is given a score of 0 (default: FALSE).
#' @param model.data Name of specificity model data to use, can be "hconf" : individual experimental kinase specificity models used to scan for rewiring events. For experimental kinase specificity models, grouped by family, set to "hconf-fam". Both are considered high confidence. For lower confidence predicted specificity models , set to "lconf". NOTE: Predicted models are purely speculative and should be used with caution
#'
#' @return
#' The data is returned in a \code{data.frame} with the following columns:
#' \item{gene}{Gene with the rewiring event}
#' \item{mut}{Mutation causing the rewiring event}
#' \item{psite_pos}{(Optional) Position of the phosphosite, if domain = "phos"}
#' \item{mut_dist}{(Optional) Distance of the mutation relative to the central residue, if domain = "phos"}
#' \item{wt}{Sequence of the wildtype phosphosite (before the mutation). Score is NA if the central residue is not S, T or Y}
#' \item{mt}{Sequence of the mutated phosphosite (after the mutation). Score is NA if the central residue is not S, T or Y}
#' \item{score_wt}{Matrix similarity score of the wildtype phosphosite}
#' \item{score_mt}{Matrix similarity score of the mutated phosphosite}
#' \item{log_ratio}{Log2 ratio between mutant and wildtype scores. A high positive log ratio represents a high confidence gain-of-phosphorylation event. A high negative log ratio represents a high confidence loss-of-phosphorylation event. This ratio is NA for mutations that affect the central phosphorylation sites}
#' \item{pwm}{Name of the kinase being rewired}
#' \item{pwm_fam}{(Optional, available only if domain = "phos") Family/subfamily of kinase being rewired. If a kinase subfamily is available the family and subfamily will be separated by an underscore e.g. "DMPK_ROCK". If no subfamily is available, only the family is shown e.g. "GSK"}
#' \item{nseqs}{(Optional, available only if domain = "phos") Number of sequences used to construct the PWM. PWMs constructed with a higher number of sequences are generally considered of better quality.}
#' \item{prob}{Joint probability of wild type sequence belonging to the foreground distribution and mutated sequence belonging to the background distribution, for loss and vice versa for gain.}
#' \item{effect}{Type of rewiring event, can be "loss" or "gain"}
#'
#' @keywords mimp psites mutation snp snv pwm rewiring phosphorylation kinase sh3 pdz
#' @export
#' @examples
#' # Get the path to example mutation data
#' mut.file = system.file("extdata", "mutation_data.txt", package = "rmimp")
#'
#' # Get the path to example FASTA sequence data
#' seq.file = system.file("extdata", "sequence_data.txt", package = "rmimp")
#'
#' # View the files in a text editor
#' browseURL(mut.file)
#' browseURL(seq.file)
#'
#' # Run rewiring analysis
#' results = mimp(mut.file, seq.file, display.results=TRUE)
#'
#' # Show head of results
#' head(results)
mimp <- function(muts, seqs, central=T, domain="phos", psites=NULL, pdz.range=5, prob.thresh=0.5, log2.thresh=1, display.results=T, include.cent=F, model.data='hconf'){
  # Constant, don't change
  flank = 7
  
  # Ensure valid domain
  if (!is.element(domain, ifelse(central, .VALID_DOMAINS$central, .VALID_DOMAINS$not_central))) {
    stop("Domain must be valid. Please check mimp documentation for a list of valid domains")
  }

  # Ensure valid thresholds
  if(!is.numeric(prob.thresh) | prob.thresh > 1 | prob.thresh < 0) stop('Probability threshold "prob.thresh" must be between 0.5 and 1')

  # Read in seq and mutation data
  seqdata = .readSequenceData(seqs)
  md = .readMutationData(muts, seqdata)

  if (central) {
    # Read in phos site data
    pd = .readPsiteData(psites, seqdata)
    if (nrow(pd) == 0) {
      warning('No phosphorylation data!')
      return(NULL)
    }

    # Keep only genes we have data for
    md = md[md$gene %in% pd$gene,]
    seqdata = seqdata[names(seqdata) %in% pd$gene]
  }
  
  # Validate that mutations mapp to the correct residue
  .validateMutationMapping(md, seqdata)

  # If we have nothing left, throw a warning and return NULL
  if(nrow(md) == 0 | length(seqdata) == 0){
    warning('No mutation or sequence data remaining after filtering!')
    return(NULL)
  }

  # Get data
  cat('\r.... | loading specificity models')
  mpath = .getModelDataPath(model.data, domain = domain)
  model.data = mpath$id
  mdata = readRDS(mpath$path)
  cat('\rdone\n')

  if (central) {
    # Ensure psite positions map to an S, T or Y
    .validatePsitesSTY(pd, seqdata)

    # Get mutations in flanking regions of psites
    mut_sites = pSNVs(md, pd, seqdata, flank)

    # Are we keeping central residues?
    if(!include.cent) mut_sites = mut_sites[mut_sites$mut_dist != 0,]

    # If no mutations map, throw warning and return NULL
    if(nrow(mut_sites) == 0){
      warning('No SNVs were found!')
      return(NULL)
    }
  } else if (domain == "sh3") {
    # Get wt and mut sequences
    pwms <- sapply(mdata, function(obj) return(obj$pwm))
    flank <- unique(unlist(sapply(pwms, ncol), use.names = F))
    mut_sites <- lapply(flank, sSNVs, md = md, seqdata = seqdata)
    names(mut_sites) <- flank

    # If no mutations map, throw warning and return NULL
    if(length(mut_sites) == 0){
      warning('No SNVs were found!')
      return(NULL)
    }
  }

  if (central) {
    scored = lapply(1:length(mdata), function(i){
      # Processing PWM i
      obj = mdata[[i]]

      # Print msg
      perc = round((i/length(mdata))*100)
      cat(sprintf('\r%3d%% | predicting %s events', perc, .DESCRIPTIONS[domain]))

      computeRewiring(obj, mut_sites, prob.thresh, log2.thresh, include.cent)
    })
  } else if (domain == "sh3") {
    # Score sequences
    scored = lapply(1:length(mdata), function(i){
      # Processing PWM i
      obj = mdata[[i]]
      obj$pwm_name <- names(mdata[i]) 
      pwm_ncol <- as.character(ncol(obj$pwm))

      # Print msg
      perc = round((i/length(mdata))*100)
      cat(sprintf('\r%3d%% | predicting %s events', perc, .DESCRIPTIONS[domain]))

      computeSh3Binding(obj, mut_sites[[pwm_ncol]], md[["mut"]], prob.thresh, log2.thresh)
    })
  }
  
  cat('\rdone | predicting', .DESCRIPTIONS[[domain]], 'events\n')

  # Check if we have any data left
  scored = scored[!sapply(scored, is.null)]
  if(length(scored) == 0){
    warning(sprintf("No %s event found, returning NULL", .DESCRIPTIONS[domain]))
    return(NULL)
  }

  # Merge all data frames
  scored_final = rbindlist(scored)
  scored_final = scored_final[order(prob, decreasing=T),]
  scored_final = as.data.frame(scored_final)
  
  attr(scored_final, 'model.data') = model.data
  
  if(display.results & model.data == "custom") {
    warning("cannot generate HTML when using custom model.")
  } else if (display.results) {
    results2html(scored_final, domain = domain)
  }

  # Remove useless columns, reset attributes and return
  cat(sprintf('\rdone | analysis complete with a total of %s %s events\a\n', nrow(scored_final), .DESCRIPTIONS[domain]))
  scored_final = scored_final[,!names(scored_final) %in% c('ref_aa', 'alt_aa', 'mut_pos')]
  attr(scored_final, 'model.data') = model.data

  return(scored_final)
}

#' Compute posterior probability of wild type phosphosites for kinases
#'
#' @param psites phosphorylation data, see \code{?mimp} for details
#' @param seqs sequence data, see \code{?mimp} for details
#' @param model.data MIMP model used, see \code{?mimp} for details
#' @param posterior_thresh posterior probability threshold that the score belongs to the foreground distribution of the kinase, probabilities below this value are discarded (default 0.8)
#' @param intermediate if TRUE intermediate MSS scores and likelihoods are reported (default FALSE)
#' @param kinases vector of kinases used for the scoring (e.g. c("AURKB", "CDK2")), if this isn't provided all kinases will be used .
#'
#' @export
#' @return
#' The data is returned in a \code{data.frame} with the following columns:
#' \item{gene}{Gene with the rewiring event}
#' \item{pos}{Position of the phosphosite}
#' \item{wt}{Sequence of the wildtype phosphosite}
#' \item{score_wt}{(intermediate value) matrix similarity score of sequence}
#' \item{l.wt.fg}{(intermediate value) likelihood of score given foreground distribution}
#' \item{l.wt.bg}{(intermediate value) likelihood of score given background distribution}
#' \item{post.wt.fg}{posterior probability of score in foreground distribution}
#' \item{post.wt.bg}{posterior probability of score in background distribution}
#' \item{pwm}{Name of the predicted kinase}
#' \item{pwm_fam}{Family/subfamily of the predicted kinase. If a kinase subfamily is available the family and subfamily will be seprated by an underscore e.g. "DMPK_ROCK". If no subfamily is available, only the family is shown e.g. "GSK"}
#' 
#' @importFrom data.table rbindlist 
#' 
#' If no predictions were made, function returns NULL
#' @examples
#' # Get the path to example phosphorylation data 
#' psite.file = system.file("extdata", "sample_phosphosites.tab", package = "rmimp")
#' 
#' # Get the path to example FASTA sequence data 
#' seq.file = system.file("extdata", "sample_seqs.fa", package = "rmimp")
#' 
#' # Run for all kinases
#' results_all = predictKinasePhosphosites(psite.file, seq.file)
#' 
#' # Run for select kinases
#' results_select = predictKinasePhosphosites(psite.file, seq.file, kinases=c("AURKB", "CDK2"))
predictKinasePhosphosites <- function(psites, seqs, model.data='hconf', posterior_thresh=0.8, intermediate=F, kinases){
  flank = 7
  # Read data
  seqdata = .readSequenceData(seqs)
  pd = .readPsiteData(psites, seqdata)

  # Get flanking sequence
  pd$wt = flankingSequence(seqdata[pd$gene], pd$pos, flank)
  seqdata = seqdata[names(seqdata) %in% pd$gene]

  # If we have nothing left, throw a warning and return NULL
  if(nrow(pd) == 0 | length(seqdata) == 0){
    warning('No phosphorylation or sequence data remaining after filtering!')
    return(NULL)
  }

  # Ensure psite positions map to an S, T or Y
  .validatePsitesSTY(pd, seqdata)

  # Get data
  cat('\r.... | loading specificity models')
  mdata = readRDS(.getModelDataPath(model.data)$path)
  cat('\rdone\n')

  # If we're missing kinase names
  if(missing(kinases)) kinases = names(mdata)

  # Throw error if invalid kinase names
  kinases = intersect(kinases, names(mdata))
  if(length(kinases) == 0)
    stop(sprintf('Invalid kinase names. Please choose from the following: %s',
                 paste(names(mdata), collapse=',')))

  # Score!
  scored = lapply(1:length(kinases), function(i){
    kin = kinases[i]
    # Processing PWM i
    obj = mdata[[kin]]
    # Print msg
    perc = round((i/length(kinases))*100)
    cat(sprintf('\r%3d%% | scoring sequences', perc))

    # PWM scores
    pwm_scores = mss(pd$wt, obj$pwm)

    # Get posteriors
    post = pRewiringPosterior(wt.scores = pwm_scores,
                              fg.params = obj$fg.params,
                              bg.params = obj$bg.params)

    # Attach mss score too
    post = cbind(mss=pwm_scores, post)
    # Set posteriors that have NA (ST vs Y) to -1 for filtering
    post$post.wt.fg[is.na(post$post.wt.fg)] = -1

    # Pick only those that pass the threshold
    ind = post$post.wt.fg >= posterior_thresh

    # Bind with phosphosite information
    z = cbind(pd[ind,], post[ind,])
    if(nrow(z) == 0) return(NULL)

    z$pwm = obj$name
    z$pwm_fam = obj$family
    z
  })

  cat('\rdone\n')

  # Remove null objects
  scored = scored[!sapply(scored, is.null)]
  if(length(scored) == 0) return(NULL)

  # Rbind everything
  final = as.data.frame( rbindlist(scored) )
  rownames(final) = NULL

  # Keep intermediates?
  if(!intermediate) final = final[,!names(final) %in% c('mss', 'l.wt.fg', 'l.wt.bg')]

  final
}

