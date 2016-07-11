#' Read FASTA data from file or use directly if list
#'
#' @param seqs Path to fasta sequences or a name list, with collapsed sequences
.readSequenceData <- function(seqs){
  
  if( is.character(seqs) ){
    msg = sprintf('reading fasta file from %s', normalizePath(seqs))
    cat(sprintf('\r.... | %s', msg))
    
    seqdata = .read.fasta(seqs, seqtype='AA', forceDNAtolower=F, as.string=T)
  }else if( is.list(seqs) ){
    msg = 'processing fasta file from list'
    cat(sprintf('\r.... | %s', msg))
    
    seqdata = seqs
  }else{
    stop('Input sequences must be a path to the fasta file, or a list of sequences')
  }
  cat('\rdone\n')
  return(seqdata)
}

#' Read mutation data from file or use directly if data.frame
#'
#' @param muts Path to mutation data or a two column data.frame [gene,mut]
#' @param seqdata Named list of sequences
.readMutationData <- function(muts, seqdata){
  msg = 'reading mutation data'
  MUT_REGEX = '^[A-Z]\\d+[A-Z]$'
  
  md = muts
  
  if(!is.data.frame(md) & !is.character(md)){
    stop('Mutation input should be a path to the mutation data file or a data.frame')
  }
  
  if( is.character(muts) ){
    msg = sprintf('reading mutation data file from %s', normalizePath(muts))
    cat(sprintf('\r.... | %s', msg))
    md = read.table(muts, header=F, stringsAsFactors=F)
  }else{
    msg = 'processing mutation data from data frame'
    cat(sprintf('\r.... | %s', msg))
  }
  if(ncol(md) < 2){
    cat('\n')
    stop('Mutation data should have at least two columns, gene and mutation')
  }
  names(md)[1:2] = c('gene', 'mut')
  
  # Remove all factors
  md = unfactor(md)
  
  # 2. Validate mutation data
  # Ensure valid regex for mutations
  if(!(all(grepl(MUT_REGEX, md$mut) ))){
    cat('\n')
    stop('Mutations must follow the following format X123Y for example: A78R!')
  }
  
  # Ensure gene in mutation data is matched to fasta file, if not ignore these rows
  z = setdiff(md$gene, names(seqdata))
  if(length(z) > 0){
    warning(sprintf('%s gene(s) found in mutation data could not be matched to headers of the fasta file and will be ignored: %s', 
                    length(z), paste0(z, collapse=', ')))
    md = md[!md$gene %in% z,]
  }
  
  
  # Get mutation position/ref/alt
  end = nchar(md$mut)
  md$ref_aa = substr(md$mut, 1,1)
  md$alt_aa = substr(md$mut, end,end)
  md$mut_pos = as.integer(substr(md$mut, 2,end-1))
  
  
  cat('\rdone\n')
  return(md)
}



#' Read phosphorylation data from file or use directly if data.frame (if any)
#' If no phosphorylation data is provided (psites = NULL), then use all STYs
#'
#' @param muts Path to phosphorylation data or a two column data.frame [gene,psite]
#' @param seqdata Named list of sequences
.readPsiteData <- function(psites, seqdata){
  
  
  DIG_REGEX = '^\\d+$'
  
  # Do we have phosphorylation sites?
  have.ps = !is.null(psites)
  
  if(have.ps){
    # Yes, read from file
    pd = psites
    if(!is.data.frame(pd) & !is.character(pd)){
      stop('Phosphorylation input should be a path to the phosphorylation data file or a data.frame')
    }
    
    # Read in from file if necessary
    if( is.character(psites) ){
      msg = sprintf('reading phosphorylation sites file from %s', normalizePath(psites))
      cat(sprintf('\r.... | %s', msg))
      pd = read.table(psites, header=F, stringsAsFactors=F)
    }else{
      msg = 'processing phosphorylation sites from data frame'
      cat(sprintf('\r.... | %s', msg))
    }
    if(ncol(pd) < 2){
      cat('\n')
      stop('Phosphorylation data should have at least two columns, gene and position')
    }
    names(pd)[1:2] = c('gene', 'pos')
    
  }else{
    
    msg = 'no phosphorylation data found, generating potential phosphosites using all STY residues'
    cat(sprintf('\r.... | %s', msg))
    # No psite file, generate psites using all STYs
    sd_sp = strsplit(unlist(seqdata), '')
    
    pd = lapply(names(sd_sp), function(n){
      ind = grep('[STY]', sd_sp[[n]])
      data.frame(gene=n, pos=ind, stringsAsFactors=F)
    })
    pd = do.call(rbind, pd)
    
  }
  
  # Validate p-site data
  # Ensure numerical values in p-site position data
  if(!(all(grepl(DIG_REGEX, pd$pos) ))){
    cat('\n'); stop('All positions of psites must be non-negative and non-zero!')
  }
  
  # Ensure gene in p-site data is matched to fasta file, if not ignore these rows
  z = setdiff(pd$gene, names(seqdata))
  if(length(z) > 0){
    warning(sprintf('%s gene(s) found in psite data could not be matched to headers of the fasta file and will be ignored: %s', 
                    length(z), paste0(z, collapse=', ')))
    pd = pd[!pd$gene %in% z,]
  }
  
  pd$pos = as.integer(pd$pos)
  cat('\rdone\n')
  return(pd)
}


#' Ensure phosphosites map to an STY residue in the sequences
#'
#' @param muts Phosphorylation data.frame [gene,psite]
#' @param seqdata Named list of sequences
.validatePsitesSTY <- function(pd, seqdata){
  
  # Check that psite positions matches an S, T or Y
  psite.aa = substr( unlist(seqdata[pd$gene]), pd$pos, pd$pos )
  wrong.psite = which( psite.aa != '' & !grepl('S|T|Y', psite.aa) )
  n.wrong = length(wrong.psite)
  
  # If we have any wrong, throw warning
  if(n.wrong > 0){
    n = pmin(10, n.wrong)
    wrong.psite = wrong.psite[1:n]
    wr = pd[wrong.psite,]
    wr = sprintf('%s: expected STY at %s found %s', wr$gene, wr$pos, psite.aa[wrong.psite])
    
    warning(sprintf('There are %s psite(s) that do not correspond to an S, T or Y. These will be ignored:\n%s%s', 
                    n.wrong, paste(wr, collapse='\n'), 
                    ifelse(n.wrong > 10, '\n...', '') ))
  }
}


#' Ensure mutations map to the reference amino acid in sequences
#'
#' @param muts Mutation data.frame [gene,mut,ref_aa,alt_aa,mut_pos]
#' @param seqdata Named list of sequences
.validateMutationMapping <- function(muts, seqdata){
  
  mut.aa = substr( unlist(seqdata[muts$gene]), muts$mut_pos, muts$mut_pos )
  wrong.mut = which( mut.aa != '' & mut.aa != muts$ref_aa )
  n.wrong = length(wrong.mut)
  
  if(n.wrong > 0){
    n = pmin(10, n.wrong)
    wrong.mut = wrong.mut[1:n]
    wr = muts[wrong.mut,]
    wr = sprintf('%s: expected %s at %s found %s', wr$gene, wr$ref_aa, wr$mut_pos, mut.aa[wrong.mut])
    warning(sprintf('The reference amino acid for %s mutation(s) do not correspond to the amino acid in the sequence:\n%s',
                    n.wrong, paste0(wr, collapse='\n'),
                    ifelse(n.wrong > 10, '\n...', '') ))
  }
  
}



#' Get path to precomputed models in MIMP
#'
#' @param model.data name of models to retrieve. This can be hconf, hconf-fam, or lconf. It can also be a path to an RDS file containing custom models
#' @export
.getModelDataPath <- function(model.data){
  mdata = c('hconf'     ='kinase_individual_human_experimental.mimp', 
            'hconf-fam' ='kinase_family_human_experimental.mimp', 
            'lconf'     ='kinase_individual_human_predicted.mimp')
  
  if(length(model.data) != 1 | !is.character(model.data)){ 
    cat('\n'); stop('model.data must be a character of length one')
  }
  
  custom.models = file.exists(model.data)
  if(!model.data %in% names(mdata) & !custom.models){
    cat('\n'); stop('model.data must be one of the following: "hconf", "hconf-fam", or "lconf"')
  }
  
  model.data.path = model.data
  id = model.data
  if(!custom.models){
    model.data.path = file.path(BASE_DIR, mdata[model.data])
  }else{
    id = 'custom'
  }
  
  list(path=model.data.path, id=id)
}

#' Read fasta file as list (from seqinr package)
#'
#' @param file file containing fasta sequences
#' 
#' @keywords fasta sequence
#' @export
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


