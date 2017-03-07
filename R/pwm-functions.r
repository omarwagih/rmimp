
# Constants
AA = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
AA_PRIORS_HUMAN  =   c(A=0.070, R=0.056, N=0.036, D=0.048,C=0.023,Q=0.047,E=0.071,G=0.066,H=0.026,I=0.044,
                       L=0.100,K=0.058,M=0.021,F=0.037,P=0.063,S=0.083,T=0.053,W=0.012,Y=0.027,V=0.060)
AA_PRIORS_YEAST  =   c(A=0.055, R=0.045, N=0.061, D=0.058, C=0.013, Q=0.039, E=0.064, G=0.05, H=0.022, I=0.066,
                       L=0.096, K=0.073, M=0.021, F=0.045, P=0.044, S=0.091, T=0.059, W=0.01, Y=0.034, V=0.056)

#' Construct position weight matrix
#' 
#' Makes a position weight matrix given aligned sequences.
#'
#' @param seqs Aligned sequences all of the same length
#' @param pseudocount Pseudocount factor. Final pseudocount is background probability * this factor
#' @param relative.freq Set to TRUE if each column should be divided by the sum
#' @param is.kinase.pwm Set to TRUE if matrix is being built for a kinase
#' @param priors Named character vector containing priors of amino acids.
#' @param do.pseudocounts TRUE if we are to add pseudocounts
#' @keywords pwm construct
#' @export
#' @examples
#' # No examples
PWM <- function(seqs, pseudocount=0.01, relative.freq=T, is.kinase.pwm=T, priors=AA_PRIORS_HUMAN, do.pseudocounts=F){
  
  
  # Ensure same length characters 
  seq.len = sapply(seqs, nchar)
  num.pos = seq.len[1]
  if(! all(seq.len == num.pos)) stop('Unequal length of sequences')
  
  # List of valid amino acids, sorted
  namespace = AA
  
  # Match priors to aa 
  bg.prob = priors[match(namespace, names(priors))]
  
  # Make matrix of letters
  split = unlist( sapply(seqs, function(seq){strsplit(seq, '')}) )
  m = t( matrix(split, seq.len, length(split)/num.pos) )
  
  # Construct PWM
  pwm.matrix = apply(m, 2, function(pos.data){
    
    # Get frequencies 
    t = table(pos.data)
    
    # Match to aa
    ind = match(namespace, names(t))
    
    # Create column
    col = t[ind]
    col[is.na(col)] = 0
    names(col) = namespace
    
    # Do pseudocounts if were logging
    if(do.pseudocounts) col = col + (pseudocount * (20*bg.prob))
    
    # Do relative frequencies
    if(relative.freq) col = col / sum(col)
    
    col
  })
  
  # Information content for MATCH score
  ic2 = apply(pwm.matrix, 2, function(col) sum(col* log2(col/bg.prob), na.rm=T))
  
  attr(pwm.matrix, 'pseudocount') = pseudocount
  attr(pwm.matrix, 'match.ic') = ic2
  
  # Assign AA names to rows/pos col
  rownames(pwm.matrix) = namespace
  colnames(pwm.matrix) = 1:num.pos
  attr(pwm.matrix, 'is.kinase.pwm') = is.kinase.pwm
  
  attr(pwm.matrix, 'nseqs') = length(seqs)
  return(pwm.matrix)
}


#' Create a degenerate PWM
#' i.e. for each aa group, set weight to the best weight of the group at that position
#' e.g. R-2 has weight 0.7, K-2 has weight 0.1. Set both R-2 and K-2 to 0.7
#' 
#' @param pwm position weight matrix
#' @param dgroups groups of amino acids
degeneratePWM <- function(pwm, dgroups=c('DE','KR','ILMV','QN','ST')){
  .pwm = pwm
  sp = strsplit(dgroups, '')
  for(dgroup in sp){
    .pwm[dgroup,] = apply(.pwm[dgroup,], 2, function(z) (z * 0) + max(z))
  }
  .pwm
}

#' Get weight/probability for each amino acid in a sequence 
#' 
#' Gets weight/probability for the amino acid at each position of the sequence
#' as an array.
#'
#' @param seqs One or more sequences to be processed
#' @param pwm Position weight matrix
#'  
#' @keywords pwm mss match tfbs
#' @export
#' @examples
#' # No Examples
scoreArray <- function(seqs, pwm){
  # Split sequence
  sp = strsplit(seqs, '')
  # Number of positions
  seq.len = length(sp[[1]])
  # Iterate through sequences
  dat = lapply(sp, function(seq){
    # Match sequence to the PWM
    mat = matrix(c(match(seq, rownames(pwm)), 1:seq.len), seq.len, 2)
    pwm[mat]
  })
  names(dat) = seqs
  return(dat)
}

#' Get weight/probability for each amino acid in a sequence 
#' 
#' Gets weight/probability for the amino acid at each position of the sequence
#' as an array.
#'
#' @param seqs One or more sequences to be processed
#' @param pwm Position weight matrix
#'  
#' @keywords pwm mss match tfbs
#' @export
#' @examples
#' # No Examples
scoreArrayRolling <- function(seqs, pwm){
  # Split sequence
  sp = strsplit(seqs, '')
  
  # Number of positions
  seq.lens = sapply(sp, length)
  pwm_ncol = ncol(pwm)
  repeats = seq.lens - pwm_ncol + 1

  # Iterate through sequences
  dat = lapply(1:length(sp), function(index){
    seq <- sp[[index]]
    if(length(seq) < pwm_ncol) {
      return(NA)
    }
    # print(sprintf("original sequence: %s and total %d repeats", paste0(seq, collapse = ""), repeats[index]))
    # Generate sequences of ncol(pwm) length
    ret <- sapply(1:repeats[index], function(i) {
      currentSeq <- seq[i : (i + pwm_ncol - 1)]
      # print(sprintf("current sequence: %s", paste0(currentSeq, collapse = "")))
      # Match sequence to the PWM
      mat = matrix(c(match(currentSeq, rownames(pwm)), 1:pwm_ncol), pwm_ncol, 2)
      ret = list(pwm[mat])
      names(ret) = paste0(currentSeq, collapse = "")
      return(ret)
    })
    return(do.call("cbind", ret))
  })
  return(dat)
}

#' Compute matrix similarity score as described in MATCH algorithm
#' 
#' Computes matrix similarity score of a PWM with a k-mer.
#' Score ranges from 0-1, as described in [PMID: 12824369]
#'
#' @param seqs Sequences to be scored
#' @param pwm Position weight matrix
#' @param na.rm Remove NA scores?
#' @param ignore.central If TRUE, central residue is ignore from scoring.
#'  
#' @keywords pwm mss match tfbs
#' @export
#' @examples
#' # No Examples
mss <- function(seqs, pwm, na.rm=F, ignore.central=T){
  
  is.kinase.pwm = attr(pwm, 'is.kinase.pwm')
  if(is.null(is.kinase.pwm)) is.kinase.pwm = T
  
  central.res = '*'
  if(is.kinase.pwm){
    # Only score sequences which have a central residue S/T or Y depending on the PWM
    kinase.type = names(which.max(pwm[,ceiling(ncol(pwm)/2)]))
    kinase.type = ifelse(grepl('S|T', kinase.type), 'S|T', 'Y')
    central.res = kinase.type
  }
  
  # Central residue index
  central.ind = NA
  if(ignore.central)
    central.ind = ceiling(ncol(pwm)/2)
  
  
  # Best/worst sequence match
  oa = scoreArray(bestSequence(pwm), pwm)[[1]]
  wa = scoreArray(worstSequence(pwm), pwm)[[1]]
  
  # Which indicies do we keep
  keep = rep(T, ncol(pwm))
  keep[central.ind] = !ignore.central
  
  # Info content
  IC = attr(pwm, 'match.ic')[keep]
  
  # Best and worst scores
  opt.score   = sum( IC * (oa [keep]), na.rm=T )
  worst.score = sum( IC * (wa [keep]), na.rm=T )
  
  # Get array of scores
  keep.scores = rep(TRUE, length(seqs))
  if(is.kinase.pwm)
    keep.scores = grepl(central.res, substr(seqs, central.ind, central.ind))
  
  # Score only ones we're keeping
  score.arr = as.list(rep(NA, length(seqs)))
  names(score.arr) = seqs
  if(sum(keep.scores) != 0){
    score.arr[keep.scores] = scoreArray(seqs[keep.scores], pwm)
  }
  
  # Get current score
  scores = sapply(score.arr, function(sa) sum( IC * (sa [keep]), na.rm=T ))
  # Normalize
  scores = (scores - worst.score) / (opt.score - worst.score)
  
  scores[!keep.scores] = NA
  # Remove NA if requested
  if(na.rm) scores = scores[!is.na(scores)]
  
  return(scores)
}

#' Compute matrix similarity score as described in MATCH algorithm for sh3 domains.
#' 
#' Computes matrix similarity score of a PWM with a k-mer.
#' Score ranges from 0-1, as described in [PMID: 12824369]
#'
#' @param seqs Sequences to be scored
#' @param pwm Position weight matrix
#'  
#' @keywords pwm mss match tfbs
#' @export
#' @examples
#' # No Examples
mssSh3 <- function(seqs, pwm){
  # Best/worst sequence match
  oa = scoreArray(bestSequence(pwm), pwm)[[1]]
  wa = scoreArray(worstSequence(pwm), pwm)[[1]]

  # Info content
  if(!is.null(attr(pwm, 'match.ic'))) {
    # If info content is already in pwm matrix, get it directly from the variable
    IC = attr(pwm, 'match.ic')
  } else {
    # Otherwise, calculate the ic
    IC = apply(pwm, 2, function(col) sum(col * logb(length(AA) * col), na.rm=T))
  }
  
  # Best and worst scores
  opt.score   = sum( IC * (oa), na.rm=T )
  worst.score = sum( IC * (wa), na.rm=T )
  
  # Score
  score.arr = scoreArrayRolling(seqs, pwm)
  score.arr = score.arr[!is.na(score.arr)]

  # Get current score
  scores = lapply(score.arr, function(sa) {
    sum <- colSums(sa * IC, na.rm = T)
    return((sum - worst.score) / (opt.score - worst.score))
  })
  names(scores) = seqs
  
  return(scores)
}

#' Given a position weight matrix, find the best matching sequence
#' 
#' Finds the amino acid at each position of the PWM with the highest occurence.
#' Used in matrix similarity score calculation.  
#'
#' @param pwm Position weight matrix
#' @keywords pwm best
#' @export
#' @examples
#' # No Examples
bestSequence <- function(pwm){
  b = rownames(pwm)[apply(pwm, 2, which.max)]
  return(paste(b, collapse=''))
}

#' Given a position weight matrix, find the worst matching sequence
#' 
#' Finds the amino acid at each position of the PWM with the lowest occurence.
#' Used in matrix similarity score calculation.  
#'
#' @param pwm Position weight matrix
#' @keywords pwm worst
#' @export
#' @examples
#' # No Examples
worstSequence <- function(pwm){
  w = rownames(pwm)[apply(pwm, 2, which.min)]
  return(paste(w, collapse=''))
}