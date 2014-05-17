BASE_DIR = system.file("extdata", "", package = "MIMP")
#BASE_DIR = '~/Development/Rmimp/inst/extdata/'

AA = c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
DNA = c('A', 'T', 'G', 'C')
RNA = c('A', 'U', 'G', 'C')

AA_PRIORS_HUMAN  =   c(A=0.070, R=0.056, N=0.036, D=0.048,C=0.023,Q=0.047,E=0.071,G=0.066,H=0.026,I=0.044,
                       L=0.100,K=0.058,M=0.021,F=0.037,P=0.063,S=0.083,T=0.053,W=0.012,Y=0.027,V=0.060)
DNA_PRIORS_HUMAN =   c(A=0.25, C=0.25, G=0.25, T=0.25)

AA_PRIORS_YEAST  =   c(A=0.055, R=0.045, N=0.061, D=0.058, C=0.013, Q=0.039, E=0.064, G=0.05, H=0.022, I=0.066,
                       L=0.096, K=0.073, M=0.021, F=0.045, P=0.044, S=0.091, T=0.059, W=0.01, Y=0.034, V=0.056)
DNA_PRIORS_YEAST =   c(A=0.310, C=0.191, G=0.191, T=0.309)


#' Save a PWM matrix object to transfac format
#'
#' This saves an already generated PWM matrix object in R 
#' to transfac format, which can be read in by RWebLogo
#'
#' @param pwm PWM matrix object
#' @param file.out where the transfac matrix is written
#' @param type 'aa', 'dna' or 'rna' depending on the namespace
#' @keywords transfac
#' 
saveTransfac <- function(pwm, file.out=tempfile('transfac'), type='aa'){
  
  id = attr(pwm, 'id') 
  if(is.null(ac)) id = 'FOO'
  
  NAMESPACE = AA
  if(type=='dna') NAMESPACE = DNA
  if(type=='rna') NAMESPACE = RNA
  
  pwm = pwm[NAMESPACE,]
  
  out = c(paste0('AC', '\t', id), 
          'XX',
          paste0('ID', '\t', id),
          'XX', 
          paste0('P0', '\t', paste(rownames(pwm), collapse='\t')))
  
  for(i in 1:ncol(pwm)){
    s = as.character(i)
    if(i < 9) s = paste0('0', i)
    col.data = paste(pwm[,i], collapse='\t')
    max.aa = names( which.max(pwm[,i]) )
    out = c(out, paste0(s, '\t', col.data,  '\t', max.aa))
  }
  
  out = c(out, 'XX', '//')
  writeLines(out, file.out)
  return(file.out)
}


#' Replace charachters at certain positions of a string with another charachter.
#'
#' @param string String to be manipulated
#' @param pos One or more positions corresponding to charachters to be changed
#' @param char Replacement charachter
#' @keywords replace char
#' @export
#' @examples
#' replaceChar('ABC', 2, 'X')
replaceChar <- function(string, pos, char) { 
  for(i in pos)
    substr(string, i, i) <- char
  string 
} 

#' Extracts digits from a string and returns them in a numerical form
#'
#' @param string String to be manipulated
#' @keywords digits numerical
#' @export
#' @examples
#' extractDigits('A123F')
extractDigits <- function(string){
  as.numeric( gsub('[^\\d]', '', string, perl=T) )
}

#' Converts all columns of a data frame of class factor to character
#'
#' @param string String to be manipulated
#' @keywords factor character
#' @export
#' @examples
#' unfactor( data.frame(x=c('A', 'B')) )
unfactor <- function(df){
  i <- sapply(df, is.factor)
  df[i] <- lapply(df[i], as.character)
  df
}

#' Construct position weight matrix
#' 
#' Makes a position weight matrix given aligned sequences.
#'
#' @param seqs Aligned sequences all of the same length
#' @param pseudocount Pseudocount factor. Final pseudocount is background probability * this factor
#' @param relative.freq TRUE if each column should be divided by the sum
#' @param type Type of sequences 'AA' or 'DNA'
#' @param priors Named character vector containing priors of amino acids.
#' @keywords pwm construct
#' @export
#' @examples
#' # No examples
PWM <- function(seqs, pseudocount=0.001, relative.freq=T, type='AA', priors=AA_PRIORS_HUMAN){
  
  # Ensure same length characters 
  seq.len = sapply(seqs, nchar)
  num.pos = seq.len[1]
  if(! all(seq.len == num.pos)) stop('Unequal length of sequences')
  
  # Type validity
  if(!type %in% c('AA', 'DNA')) stop('Type must be AA or DNA')
  
  # List of valid amino acids, sorted
  namespace = AA
  if(type == 'DNA') namespace = DNA
  
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
    
    # Do pseudocounts
    col = col + (pseudocount * bg.prob) 
    
    # Do relative frequencies
    if(relative.freq) col = col / sum(col)
    col
  })
  
  attr(pwm.matrix, 'pseudocount') = pseudocount
  ic10 = apply(pwm.matrix, 2, function(col) sum(col * log10(col/bg.prob), na.rm=T) )
  ic2 = apply(pwm.matrix, 2, function(col) sum(col* log2(col/bg.prob), na.rm=T))
  attr(pwm.matrix, 'ic10') = ic10
  attr(pwm.matrix, 'ic2') = ic2
  
  # shannon entropy from http://en.wikipedia.org/wiki/Sequence_logo
  Hi = apply(pwm.matrix, 2, function(col) - sum(col * log2(col), na.rm=T) )
  en = (1/log(2)) * ( (20 - 1)/ (2*length(seqs)) )
  shan = log2(20) - Hi #+ en
  attr(pwm.matrix, 'shannon') = shan
  
  # Assign AA names to rows/pos col
  rownames(pwm.matrix) = namespace
  colnames(pwm.matrix) = 1:num.pos
  return(pwm.matrix)
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
#' @examples
#' # No Examples
scoreArray <- function(seqs, pwm){
  # Split sequence
  sp = strsplit(seqs, '')
  
  seq.lens = sapply(sp, length)
  seq.len = seq.lens[1]
  if(any(seq.lens != seq.len)) stop('Input sequences must be same length')
  
  # Iterate through sequences
  dat = lapply(sp, function(seq){
    # Match sequence to the PWM
    mat = matrix(c(match(seq, rownames(pwm)), 1:seq.len), seq.len, 2)
    prob.vector = pwm[ mat ]
    prob.vector
  })
  names(dat) = seqs
  return(dat)
}


#' Compute matrix similarity score as described in MATCH algorithm
#' 
#' Computes matrix similarity score of a PWM with a k-mer.
#' Score ranges from 0-1, as described in [PMID: 12824369]
#'
#' @param seqs Sequences to be scored
#' @param pwm Position weight matrix
#' @param is.kinase.pwm TRUE if PWM is that of a kinase
#' @param na.rm Remove NA scores?
#'  
#' @keywords pwm mss match tfbs
#' @examples
#' # No Examples
mss <- function(seqs, pwm, is.kinase.pwm=T, na.rm=F, ignore.ind=8){
  
  central.res = '*'
  if(is.kinase.pwm){
    # Only score sequences which have a central residue S/T or Y depending on the PWM
    kinase.type = names(which.max(pwm[,ceiling(ncol(pwm)/2)]))
    kinase.type = ifelse(grepl('S|T', kinase.type), 'S|T', 'Y')
    central.res = kinase.type
  }
  
  # Central residue index
  #central.ind = ceiling(ncol(pwm)/2)
  central.ind = ignore.ind
  
  # Best/worst sequence match
  oa = scoreArray(bestSequence(pwm), pwm)[[1]]
  wa = scoreArray(worstSequence(pwm), pwm)[[1]]
  
  I = attr(pwm, 'ic2')
  
  score.arr = scoreArray(seqs, pwm)
  scores = sapply(score.arr, function(sa){
    na = is.na(sa)
    na[central.ind] = is.kinase.pwm
    
    # Get information content of non-NA values
    IC = I[!na]
    
    # MSS method
    curr.score  = sum( IC * (sa [!na]), na.rm=T ) 
    opt.score   = sum( IC * (oa [!na]), na.rm=T )
    worst.score = sum( IC * (wa [!na]), na.rm=T )
    score.final = ( (curr.score - worst.score) / (opt.score - worst.score) )
    score.final
  })
  
  if(central.res != '*'){
    keep = grepl(central.res, substr(seqs, central.ind, central.ind))
    scores[!keep] = NA
  }
  
  # Remove NA if requested
  if(na.rm) scores = scores[!is.na(scores)]
  
  return(scores)
}


#' Given a position weight matrix, find the best matching sequence
#' 
#' Finds the amino acid at each position of the PWM with the highest occurence.
#' Used in matrix similarity score calculation.  
#'
#' @param pwm Position weight matrix
#' @keywords pwm best
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
#' @examples
#' # No Examples
worstSequence <- function(pwm){
  w = rownames(pwm)[apply(pwm, 2, which.min)]
  return(paste(w, collapse=''))
}

#' Get flanking sequences of a position.
#'
#' This function obtains the flanking sequence at one or more position. 
#' Out of bound indices are replaced by a blank charachter. 
#'
#' @param seqs Charachter vector of sequences. If only one sequence is provided, 
#'    indices from \code{inds} are assumed to all be from the same sequence. 
#' @param inds Numerical vector of positions corresponding to the sequences provided in \code{seqs}. 
#' @param flank Value indicating the number of charachters to extract, before and after an index
#' @param empty.char Charachter used to replace out of bound flanking sequences
#' @keywords flank sequence
#' @examples
#' # One sequence and one index. Central charachter is 'B'
#' flankingSequence(seqs='ABC', inds=2, flank=1)
#' # An example showing the use of empty.char 
#' flankingSequence(seqs='ABC', inds=2, flank=5)
#' # An example with muliple sequences and indicies
#' flankingSequence(seqs=c('ABC', 'XYZ'), inds=c(2, 1), flank=1)
flankingSequence <- function(seqs, inds, flank=7, empty.char='-'){
  if(length(seqs) == 1 & length(inds) >= length(seqs)) seqs = rep(seqs, length(inds))
  if(length(seqs) != length(inds)) stop('Length of sequences must be equal to length of positions')
  
  border = paste0(rep(empty.char, flank), collapse="")
  seqs = sapply(seqs, function(s) paste0(border,s,border))
  
  ret = sapply(1:length(seqs), function(i){
    p = flank + inds[i]
    substr(seqs[i], p-flank, p+flank)
  })
  return(ret)
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
#' @param multicore If true, will use mclapply to speed things up!
#' 
#' @keywords psnv psites mutation snp
#' @examples
#' # No examples
pSNVs <- function(muts, psites, seqs, flank=7, multicore=F){
  # Remove any factors
  muts = unfactor(muts)
  psites = unfactor(psites)
  
  names(psites)[1:2] = c('gene', 'psite_pos')
  names(muts)[1:2] = c('gene', 'mut')
  # If we have a char class instead of numerical, extract numerical values
  if(is.character( psites$psite_pos ))
    psites$psite_pos = extractDigits(psites$psite_pos)
  
  # Split psites by gene
  psites = split(psites, psites$gene)
  
  # Extract ref/alt amino acid and mutation position
  tt = t( sapply(strsplit(muts$mut, ''), function(s){
    c(s[1], s[length(s)], paste0(s[2:(length(s)-1)],collapse=''))
  }) )
  tt = as.data.frame(tt, stringsAsFactors=F)
  names(tt) = c('ref_aa', 'alt_aa', 'mut_pos')
  tt$mut_pos = as.numeric(tt$mut_pos)
  muts = cbind(muts, tt)
  
  mutAA = sapply(1:nrow(muts), function(i){
    j = muts$mut_pos[i]
    substr( seqs[[ muts$gene[i] ]], j, j )
  })
  wrongMut = mutAA != '' & mutAA != muts$ref_aa
  
  if(sum(wrongMut) > 0){
    wr = head(muts[wrongMut,])
    wr = sprintf('%s: expected %s at %s found %s', wr$gene, wr$ref_aa, wr$mut_pos, mutAA[wrongMut])
    warning(sprintf('The reference amino acid for %s mutation(s) do not correspond to the amino acid in the sequence:\n%s',
                    sum(wrongMut), paste0(wr, collapse='\n'),
                    ifelse(sum(wrongMut) > 6, '\n...', '') ))
  }
  
  # Are we using mclapply?
  myapply = lapply
  if(multicore) myapply = mclapply
  
  # Go through each mutation
  mut_ps = myapply(1:nrow(muts), function(i){
    prot = muts$gene[i]
    mut_pos = muts$mut_pos[i]
    
    # Get psites for the mutated gene
    prot_psites = psites[[prot]]
    
    # Check if the mutation exists in any of the psites 
    t = ( mut_pos <= (prot_psites$psite_pos + flank) 
          & mut_pos >= (prot_psites$psite_pos - flank) )
    if(sum(t) == 0){
      NULL
    }else{
      # Mutate sequence
      mut_seq = replaceChar(seqs[[prot]], mut_pos, muts[i,]$alt_aa)
      # Append mutation information
      prot_psites = prot_psites[t,]
      # Append mutated flanking sequence
      prot_psites$mt = flankingSequence(mut_seq, prot_psites$psite_pos, flank)
      # Bind everything
      z = prot_psites[-1]
      suppressWarnings( ret <- cbind(muts[i,], z) )
      ret$mut_dist = mut_pos - prot_psites$psite_pos
      ret
    }
  })
  
  mut_ps = mut_ps[!sapply(mut_ps, is.null)]
  mut_ps = do.call(rbind, mut_ps)
  mut_ps$wt = flankingSequence(seqs[mut_ps$gene], mut_ps$psite_pos, flank)
  mut_ps = mut_ps[,c('gene', 'mut', 'ref_aa', 'alt_aa', 'mut_pos', 
                     'psite_pos', 'mut_dist', 'wt', 'mt')]
  return(mut_ps)
}


#' Score wt and mt sequences for a pwm
#'
#' @param pwm Position weight matrix of interest
#' @param mut_ps psnvs data frame containing wt and mt sequences computed from pSNVs function
#' @param is.kinase.pwm TRUE if pwm is that of a kinase
#' @param thresh.bg Anything below this threshold is considered a negative hit
#' @param thresh.fg Anything above this threshold is considered a positive hit
#' @param thresh.log2 Threshold for the absolute value of log ratio. Anything less than this value is discarded.
#' 
#' @keywords mut psites score snp
scoreWtMt <- function(pwm, mut_ps, is.kinase.pwm=T, thresh.bg=1, thresh.fg=0, thresh.log2=0){
  
  mut_ps$score_wt = mss(mut_ps$wt, pwm, is.kinase.pwm)
  mut_ps$score_mt = mss(mut_ps$mt, pwm, is.kinase.pwm)
  mut_ps$log_ratio = log2(mut_ps$score_mt/mut_ps$score_wt)
  mut_ps$pwm = ''
  
  mut_ps = mut_ps[!is.na(mut_ps$score_wt) & !is.na(mut_ps$score_mt),]
  
  a = mut_ps$score_wt > thresh.fg & mut_ps$score_mt < thresh.bg
  b = mut_ps$score_mt > thresh.fg & mut_ps$score_wt < thresh.bg
    
  mut_ps = mut_ps[a|b,]
  mut_ps = mut_ps[abs(mut_ps$log_ratio) >= thresh.log2, ]
  
  return(mut_ps)
}

#' Predict the impact of single variants on phosphorylation.
#'
#' This function takes in mutation, sequence and phosphorylation data to predict the 
#' impact the mutation has on phosphorylation.
#' 
#' @param muts Mutation data file: a space delimited text file or data frame containing two columns (1) gene and (1) mutation. 
#' Example:
#' \tabular{ll}{
#'    TP53 \tab R282W\cr
#'    CTNNB1 \tab S33C\cr
#'    CTNNB1 \tab S37F\cr
#' }
#' @param seqs Sequence data file containing protein sequences in FASTA format OR 
#'  named list of sequences where each list element is the uppercase sequence and the 
#'  name of each element is that of the protein. Example: list(TP53="ABCXYZ", CDK2="HJKEWR")
#' @param psites Phosphorylation data file (optional): a space delimited text file containing positions of phosphorylation sites. Example:
#' \tabular{ll}{
#'    TP53 \tab 280\cr
#'    CTNNB1 \tab 29\cr
#'    CTNNB1 \tab 44\cr
#' }
#' @param perc.bg Percentile value between 0 - 100. This value is used to compute a threshold, β from the negative (background) distribution of scores. 
#'  By default this is the 90th percentile of the background distribution of scores. Anything below the threshold is considered a negative hit.
#' @param perc.fg Percentile value between 0 - 100. This value is used to compute a threshold, α from the positive (foreground) distribution of scores. 
#'  By default this is the 10th percentile of the foreground distribution of scores. Anything above the threshold is considered a positive hit.
#' @param thresh.log2 Threshold for the absolute value of log ratio. Anything less than this value is discarded (default: 0).
#' 
#' @return 
#' The data is returned in a \code{data.frame} with the following columns:
#' \item{gene}{gene with the rewiring event}
#' \item{mut}{mutation causing the rewiring event}
#' \item{psite_pos}{position of the central residue of the phosphosite}
#' \item{mut_dist}{distance of the mutation from the central phosphosite}
#' \item{wt}{sequence of the wildtype phosphosite (before the mutation) }
#' \item{mt}{sequence of the mutated phosphosite (after the mutation)}
#' \item{score_wt}{matrix similarity score of the wildtype phosphosite }
#' \item{score_mt}{matrix similarity score of the mutated phosphosite }
#' \item{log_ratio}{Log2 ratio between mutant and wildtype scores. A high positive log ratio represents a high confidence gain-of-signaling event. A high negative log ratio represents a high confidence loss-of-signaling event.}
#' \item{pwm}{name of the kinase being rewiried}
#' \item{perc_wt}{Percentile rank of the wt score}
#' \item{perc_mt}{Percentile rank of the mutant score}
#' 
#' 
#' @keywords mimp psites mutation snp snv pwm rewiring phosphorylation kinase
#' @export
#' @examples
#' # Get the path to example mutation data 
#' mut.file = system.file("extdata", "mutation_data.txt", package = "MIMP")
#' 
#' # Get the path to example FASTA sequence data 
#' seq.file = system.file("extdata", "sequence_data.txt", package = "MIMP")
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
mimp <- function(muts, seqs, psites, perc.bg=90, perc.fg=10, thresh.log2=0, display.results=T){
  flank=7
  MUT_REGEX = '^[A-Z]\\d+[A-Z]$'
  DIG_REGEX = '^\\d+$'
  
  perc.bg = as.integer(perc.bg)
  perc.fg = as.integer(perc.fg)
  if(!is.numeric(perc.bg) | perc.bg > 100 | perc.bg < 0) stop('perc.bg must be an integer between 0-100')
  if(!is.numeric(perc.fg) | perc.fg > 100 | perc.fg < 0) stop('perc.fg must be an integer between 0-100')
  
  # 1. Read sequence data
  writeLines('Reading fasta data from file ...')
  seqdata = seqs
  if( !is.list(seqs) ) 
    seqdata = .read.fasta(seqs, seqtype='AA', forceDNAtolower=F, as.string=T)
  
  # 2. Read mutation data
  writeLines('Reading mutation data from file ...')
  md = muts
  if( !is.data.frame(muts) ) 
    md = read.table(muts, header=F, stringsAsFactors=F)
  names(md)[1:2] = c('gene', 'mut')
  
  # 2. Validate mutation data
  # Ensure valid regex for mutations
  if(!(all(grepl(MUT_REGEX, md$mut) ))) 
    stop('Mutations must follow the following format X123Y for example: A78R!')
  
  # Ensure gene in mutation data is matched to fasta file, if not ignore these rows
  z = setdiff(md$gene, names(seqdata))
  if(length(z) > 0){
    warning(sprintf('%s genes found in mutation data could not be matched to headers of the fasta file and will be ignored. Genes: %s', 
                    length(z), paste0(z, collapse=', ')))
    md = md[!md$gene %in% z,]
  }
  
  # 3. Read p-site data
  have.ps = T
  if(missing(psites)){
    have.ps = F
    psites = NULL
  }
  if(! is.null(psites)){
    writeLines('Reading phosphosites data from file ...')
    pd = psites
    if( !is.data.frame(psites) ) 
      pd = read.table(psites, header=F, stringsAsFactors=F)
    names(pd)[1:2] = c('gene', 'pos')
  }else{
    # No psite file, generate psites using all STYs
    writeLines('No phosphosite data found, generating potential phosphosites using all STYs ...')
    sd_sp = strsplit(unlist(seqdata), '')
    pd = lapply(names(sd_sp), function(n){
      ind = grep('[STY]', sd_sp[[n]])
      data.frame(gene=n, pos=ind, stringsAsFactors=F)
    })
    pd = do.call(rbind, pd)
  }
  
  # 3. Validate p-site data
  # Ensure numerical values in p-site position data
  if(!(all(grepl(DIG_REGEX, pd$pos) ))) 
    stop('All positions of p-sites must be non-negative and non-zero!')
  
  # Ensure gene in p-site data is matched to fasta file, if not ignore these rows
  z = setdiff(pd$gene, names(seqdata))
  if(length(z) > 0){
    warning(sprintf('%s genes found in p-site data could not be matched to headers of the fasta file and will be ignored. Genes: %s', 
                    length(z), paste0(z, collapse=', ')))
    pd = pd[!pd$gene %in% z,]
  }
  
  psiteAA = sapply(1:nrow(pd), function(i){
    j = pd$pos[i]
    substr( seqdata[[pd$gene[i]]], j,j ) 
  })
  
  wrongPsite = psiteAA != '' & !grepl('S|T|Y', psiteAA)
  
  if(sum(wrongPsite) > 0){
    wr = head(pd[wrongPsite,])
    wr = sprintf('%s: expected STY at %s found %s', wr$gene, wr$pos, psiteAA[wrongPsite])
    warning(sprintf('There are %s p-site(s) that do not correspond to an S, T or Y. These will be ignored:\n%s%s', 
                    sum(wrongPsite), paste(wr, collapse='\n'), 
                    ifelse(sum(wrongPsite) > 6, '\n...', '') ))
  }
  
  mut_psites = pSNVs(md, pd, seqdata, flank)
  
  
  writeLines('Loading kinase specificity models and cutoffs ...')
  mdata = readRDS( file.path(BASE_DIR, 'mimp_data.rds'))
  pwms = mdata$pwms
  perc = mdata$perc
  cutoffs = mdata$cutoffs
  
  writeLines('Predicting kinase rewiring events ...')
  pb <- txtProgressBar(min = 0, max = length(pwms), style = 3)
  
  ct = t(sapply(names(pwms), function(name){
    c_bg = cutoffs[[name]]$bg[perc.bg+1]
    c_fg = cutoffs[[name]]$fg[perc.fg+1]
    c(c_bg, c_fg)
  }))
  ct = as.data.frame(ct, stringsAsFactors=F)
  ct$pwm = names(pwms)
  colnames(ct) = c('bg', 'fg', 'pwm')
  
  scored = lapply(1:length(pwms), function(i){
    setTxtProgressBar(pb, i)
    pwm = pwms[[i]]
    name = names(pwms)[i]
    thresh.bg = cutoffs[[name]]$bg[perc.bg+1]
    thresh.fg = cutoffs[[name]]$fg[perc.fg+1]
    ss = scoreWtMt(pwm, mut_psites, T, thresh.bg, thresh.fg, thresh.log2)
    if(nrow(ss)>0) ss$pwm = name
    ss
  })
  
  close(pb)
  
  z = scored
  scored = scored[sapply(scored, function(x) nrow(x) > 0)]
  if(length(scored) == 0){
    return(z[[1]])
  }
  
  
  # Merge all data frames
  s = do.call('rbind', scored)
  
  s$effect = 'Gain'
  s$effect[s$score_wt > s$score_mt] = 'Loss'
  
  # Percentile rank
  pfg = unlist( perc$fg[s$pwm] )
  pbg = unlist( perc$bg[s$pwm] )
  
  
  r = sapply(1:nrow(s), function(i){
    e = s$effect[i]
    fg = pfg[[i]]
    bg = pbg[[i]]
    wt = s$score_wt[i]
    mt = s$score_mt[i]
    if(e == "Loss"){
      c( fg(wt), bg(mt))
    }else{
      c( bg(wt), fg(mt))
    }
  })
  r = as.data.frame(t(r))
  names(r) = c('perc_wt', 'perc_mt')
  s$perc_wt = r$perc_wt
  s$perc_mt = r$perc_mt
  
  if(display.results) results2html(s)
  
  writeLines('Analysis complete!')
  rownames(s) = NULL
  
  s = s[,!names(s) %in% c('ref_aa', 'alt_aa', 'mut_pos')]
  attr(s, 'cutoffs') = ct
  return(s)
}

#' Helper function for \code{dohtml}
#'
#' Adds colors for mutants at distance
#'
#' @param s Data frame resulting from mimp call.
#' @param dist Distance of mutation.
#' 
#' @keywords helper mimp
.htmlSeq <- function(s, dist){
  s = strsplit(s,'')[[1]]
  s[8] = sprintf('<a class="psite">%s</a>', s[8])
  p = 8 + dist
  s[p] = sprintf('<a class="mut">%s</a>', s[p])
  paste0(s, collapse='')
}

#' Helper function for \code{results2html}
#'
#' @param x Data frame resulting from mimp call.
#' @param LOGO_DIR Directory containing sequence logo images.
#' 
#' @keywords display mimp
#' @export
dohtml <- function(x, LOGO_DIR){
  x = unfactor(x)
  x$score_wt = signif(x$score_wt, 3)
  x$score_mt = signif(x$score_mt, 3)
  x$log_ratio = signif(x$log_ratio, 3)
  x$perc_wt = signif(x$perc_wt, 3)
  x$perc_mt = signif(x$perc_mt, 3)
  
  x$pwm = gsub('-', '_', x$pwm)
  
  cnt = as.list( table(x$effect) )
  n_gain = ifelse( is.null( cnt[['Gain']] ), 0, cnt[['Gain']])
  n_loss = ifelse( is.null( cnt[['Loss']] ), 0, cnt[['Loss']])
  n_mut = nrow(unique(x[,c('gene', 'mut')]))
  lines = sapply(1:nrow(x), function(i){
    r = x[i,]
    d = r$mut_dist
    seq = sprintf('%s<br>%s', .htmlSeq(r$wt, d), .htmlSeq(r$mt, d))
    scr = sprintf('%s<br>%s', r$score_wt, r$score_mt)
    prc = sprintf('%s<br>%s', r$perc_wt, r$perc_mt)
    logo = file.path(LOGO_DIR, paste0(r$pwm, '.png'))
    t = sprintf('<a class="hide name">%s</a><img src="%s" class="logo" alt=%s/>', r$pwm, normalizePath(logo), r$pwm)
    
    if(r$effect == 'Loss'){
      eff = '<a class="loss">Loss</a>'
    }else{
      eff = '<a class="gain">Gain</a>'
    }
    
    gene = sprintf('<a target="_blank" class="gene-link" href="http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s">%s</a>', r$gene, r$gene)
    
    # sub - _ logos
    sprintf('<tr> <td class="gene-name">%s</td>
            <td class="psite-pos">%s</td>
            <td class="mut-abbr">%s</td>
            <td class="mut-dist">%s</td>
            <td class="sequence">%s</td>
            <td class="wt-score">%s</td>
            <td class="mt-score">%s</td>
            <td class="wt-perc">%s</td>
            <td class="mt-perc">%s</td>
            <td class="log-ratio">%s</td>
            <td class="effect">%s</td>
            <td class="seq-logo">%s</td>
            </tr>', 
            r$gene, r$psite_pos, r$mut, r$mut_dist, seq, r$score_wt, r$score_mt, r$perc_wt, r$perc_mt, r$log_ratio, eff, t)
    
  })
  
  tt = '<div id="%s" style="display:none">%s</div>'
  
  lines = append(lines,
                 c(sprintf(fmt=tt, 'n_mut', n_mut),
                   sprintf(fmt=tt, 'n_gain', n_gain),
                   sprintf(fmt=tt, 'n_loss', n_loss)))
  return(lines)
}

#' Display MIMP results interactively in browser
#'
#' @param x Data frame resulting from mimp call.
#' @param max.rows If data contains more rows than this value, results won't be displayed.
#' 
#' @keywords display mimp
#' @export
#' 
results2html <- function(x, max.rows=5000){
  if(nrow(x) > max.rows){
    warning(sprintf('Rows of resulting data exceeds %s and cannot be displayed', max.rows))
    return(1)
  }
  #BASE_DIR = system.file("extdata", "", package = "MIMP")
  LOGO_DIR = file.path(BASE_DIR, 'html', 'images', 'logos')
  lines = dohtml(x, LOGO_DIR)
  save = file.path(BASE_DIR, 'html', 'MIMP_results.html')
  zz <- file(save,"w")
  tt = readLines( file.path(BASE_DIR, 'html', 'index.html'))
  ind = grep('<DATA>', tt)
  writeLines(tt[1: (ind-1)],con=zz,sep="\n")
  writeLines(lines,con=zz,sep="\n")
  writeLines(tt[(ind+1): length(tt)],con=zz,sep="\n")
  close(zz)
  
  browseURL(save)
}


#' Read fasta file as list (from seqinr package)
#'
#' @param file file containing fasta sequences
#' 
#' @keywords fasta sequence
#' 
.read.fasta = function (file, 
                        seqtype = c("DNA", "AA"), as.string = FALSE, forceDNAtolower = TRUE, 
                        set.attributes = TRUE, legacy.mode = TRUE, seqonly = FALSE, 
                        strip.desc = FALSE, bfa = FALSE, sizeof.longlong = .Machine$sizeof.longlong, 
                        endian = .Platform$endian, apply.mask = TRUE) 
{
  seqtype <- match.arg(seqtype)
  if (!bfa) {
    lines <- readLines(file)
    if (legacy.mode) {
      comments <- grep("^;", lines)
      if (length(comments) > 0) 
        lines <- lines[-comments]
    }
    ind <- which(substr(lines, 1L, 1L) == ">")
    nseq <- length(ind)
    if (nseq == 0) {
      stop("no line starting with a > character found")
    }
    start <- ind + 1
    end <- ind - 1
    end <- c(end[-1], length(lines))
    sequences <- lapply(seq_len(nseq), function(i) paste(lines[start[i]:end[i]], 
                                                         collapse = ""))
    if (seqonly) 
      return(sequences)
    nomseq <- lapply(seq_len(nseq), function(i) {
      firstword <- strsplit(lines[ind[i]], " ")[[1]][1]
      substr(firstword, 2, nchar(firstword))
    })
    if (seqtype == "DNA") {
      if (forceDNAtolower) {
        sequences <- as.list(tolower(sequences))
      }
    }
    if (as.string == FALSE) 
      sequences <- lapply(sequences, s2c)
    if (set.attributes) {
      for (i in seq_len(nseq)) {
        Annot <- lines[ind[i]]
        if (strip.desc) 
          Annot <- substr(Annot, 2L, nchar(Annot))
        attributes(sequences[[i]]) <- list(name = nomseq[[i]], 
                                           Annot = Annot, class = switch(seqtype, AA = "SeqFastaAA", 
                                                                         DNA = "SeqFastadna"))
      }
    }
    names(sequences) <- nomseq
    return(sequences)
  }
  if (bfa) {
    if (seqtype != "DNA") 
      stop("binary fasta file available for DNA sequences only")
    mycon <- file(file, open = "rb")
    r2s <- words(4)
    readOneBFARecord <- function(con, sizeof.longlong, endian, 
                                 apply.mask) {
      len <- readBin(con, n = 1, what = "int", endian = endian)
      if (length(len) == 0) 
        return(NULL)
      name <- readBin(con, n = 1, what = "character", endian = endian)
      ori_len <- readBin(con, n = 1, what = "int", endian = endian)
      len <- readBin(con, n = 1, what = "int", endian = endian)
      seq <- readBin(con, n = len * sizeof.longlong, what = "raw", 
                     size = 1, endian = endian)
      mask <- readBin(con, n = len * sizeof.longlong, what = "raw", 
                      size = 1, endian = endian)
      if (endian == "little") {
        neword <- sizeof.longlong:1 + rep(seq(0, (len - 
                                                    1) * sizeof.longlong, by = sizeof.longlong), 
                                          each = sizeof.longlong)
        seq <- seq[neword]
        mask <- mask[neword]
      }
      seq4 <- c2s(r2s[as.integer(seq) + 1])
      seq4 <- substr(seq4, 1, ori_len)
      if (apply.mask) {
        mask4 <- c2s(r2s[as.integer(mask) + 1])
        mask4 <- substr(mask4, 1, ori_len)
        npos <- gregexpr("a", mask4, fixed = TRUE)[[1]]
        for (i in npos) substr(seq4, i, i + 1) <- "n"
      }
      return(list(seq = seq4, name = name))
    }
    sequences <- vector(mode = "list")
    nomseq <- vector(mode = "list")
    i <- 1
    repeat {
      res <- readOneBFARecord(mycon, sizeof.longlong, endian, 
                              apply.mask)
      if (is.null(res)) 
        break
      sequences[[i]] <- res$seq
      nomseq[[i]] <- res$name
      i <- i + 1
    }
    close(mycon)
    nseq <- length(sequences)
    if (seqonly) 
      return(sequences)
    if (as.string == FALSE) 
      sequences <- lapply(sequences, s2c)
    if (set.attributes) {
      for (i in seq_len(nseq)) {
        if (!strip.desc) 
          Annot <- c2s(c(">", nomseq[[i]]))
        attributes(sequences[[i]]) <- list(name = nomseq[[i]], 
                                           Annot = Annot, class = "SeqFastadna")
      }
    }
    names(sequences) <- nomseq
    return(sequences)
  }
}
