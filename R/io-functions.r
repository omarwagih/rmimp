#' Read FASTA data from file or use directly if list
#'
#' @param seqs Path to fasta sequences or a name list, with collapsed sequences
#' 
#' @importFrom Biostrings readAAStringSet
#' @keywords internal
.readSequenceData <- function(seqs){
  
  if( is.character(seqs) ){
    msg = sprintf('reading fasta file from %s', normalizePath(seqs))
    cat(sprintf('\r.... | %s', msg))
    
    tmp = Biostrings::readAAStringSet(seqs)
    seqdata = as.list( as.character( tmp ) )
  }else if( is.list(seqs) ){
    msg = 'processing fasta file from list'
    cat(sprintf('\r.... | %s', msg))
    
    seqdata = seqs
  }else{
    stop('Input sequences must be a path to the fasta file, or a named list of sequences')
  }
  cat('\rdone\n')
  return(seqdata)
}

#' Read mutation data from file or use directly if data.frame
#'
#' @param muts Path to mutation data or a two column data.frame [gene,mut]
#' @param seqdata Named list of sequences
#' 
#' @importFrom data.table fread
#' @keywords internal
.readMutationData <- function(muts, seqdata){
  msg = 'reading mutation data'
  MUT_REGEX = '^[A-Z]\\d+[A-Z]$'
  
  md = muts
  
  # Ensure mutation data is path or data frame
  if(!is.data.frame(md) & !is.character(md))
    stop('Mutation input should be a path to the mutation data file or a data.frame')
  
  if( is.character(muts) ){
    # We have a path to mutation data - read in using data table
    msg = sprintf('reading mutation data file from %s', normalizePath(muts))
    cat(sprintf('\r.... | %s', msg))
    md = fread(muts, header=F, stringsAsFactors = F, data.table = F)
  }else{
    msg = 'processing mutation data from data frame'
    cat(sprintf('\r.... | %s', msg))
  }
  
  if(ncol(md) < 2){
    cat('\n')
    stop('Mutation data should have at least two columns, gene and mutation')
  }
  
  # Set column names and remove all factors - if any
  names(md)[1:2] = c('gene', 'mut')
  md = rapply(md, as.character, classes="factor", how="replace")
  
  # 2. Validate mutation data
  # Ensure valid regex for mutations
  if(!(all(grepl(MUT_REGEX, md$mut) ))){
    cat('\n')
    stop('All mutations must follow the following format X123Y for example: A78R!')
  }
  
  # If we have sequence data
  if(!is.null(seqdata)){
    # Ensure gene in mutation data is matched to fasta file, if not ignore these rows
    z = setdiff(md$gene, names(seqdata))
    if(length(z) > 0){
      warning(sprintf('%s gene(s) found in mutation data could not be matched to headers of the fasta file and will be ignored: %s', 
                      length(z), paste0(z, collapse=', ')))
      md = md[!md$gene %in% z,]
    }
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
#' 
#' @importFrom data.table fread rbindlist
#' @keywords internal
.readPsiteData <- function(psites, seqdata){
  
  DIG_REGEX = '^\\d+$'
  
  if(is.null(seqdata)){
    psp_path = file.path(BASE_DIR, 'psp_phosphosites.rds')
    psp = readRDS(psp_path)
    return(psp)
  }
  
  # Do we have phosphorylation sites?
  have.ps = !is.null(psites)
  
  if(have.ps){
    # Yes, read from file
    pd = psites
    if(!is.data.frame(pd) & !is.character(pd)){
      stop('psites should be a path to the phosphorylation data file or a data.frame')
    }
    
    if( is.character(psites) ){
      # Read in from path
      msg = sprintf('reading phosphorylation sites file from %s', normalizePath(psites))
      cat(sprintf('\r.... | %s', msg))
      pd = fread(psites, header=F, stringsAsFactors=F, data.table = F)
    }else{
      # Assume we have data frame passed in
      msg = 'processing phosphorylation sites from data frame'
      cat(sprintf('\r.... | %s', msg))
    }
    
    # Cannot have less than two columns
    if(ncol(pd) < 2){
      cat('\n')
      stop('Phosphorylation data should have at least two columns, gene and position')
    }
    
    # Set first two column names
    names(pd)[1:2] = c('gene', 'pos')
    
  }else{
    
    msg = 'no phosphorylation data found, generating potential phosphosites using all STY residues'
    cat(sprintf('\r.... | %s', msg))
    # No psite file, generate psites using all STYs
    sd_sp = strsplit(unlist(seqdata), '')
    
    pd = lapply(names(sd_sp), function(n){
      ind = grep('[STY]', sd_sp[[n]])
      
      # Skip if phosphorylation sites are not found
      if (identical(ind, integer(0))) {
        warning(sprintf("%s gene is skipped because it has no phosphorylation sites.", n))
        return(NULL)
      }
      
      data.frame(gene=n, pos=ind, stringsAsFactors=F)
    })
    
    pd = as.data.frame( rbindlist(pd) )
  }
  
  # Validate p-site data
  # Ensure numerical values in p-site position data
  if(!(all(grepl(DIG_REGEX, pd$pos) ))){
    cat('\n')
    stop('All positions of psites must be non-negative and non-zero!')
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
#' 
#' @keywords internal
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
#' 
#' @keywords internal
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
#' @param central Whether the mutation site is at the central residue of the sequence
#' @param domain Which binding domain to run mimp for
#' @export
.getModelDataPath <- function(model.data, domain="phos", species="human"){
  if (domain == "phos") {
    mdata = c('hconf'     = sprintf('kinase_individual_%s_experimental.mimp', species), 
              'hconf-fam' = sprintf('kinase_family_human_experimental.mimp', species), 
              'lconf'     = sprintf('kinase_individual_human_predicted.mimp', species))
  } else if (domain %in% c("sh3", "sh2", "pdf")) {
    mdata = c('hconf'     = sprintf('%s_individual_human_experimental.mimp', domain, species))
  }
  
  if(length(model.data) != 1 | !is.character(model.data)){ 
    cat('\n'); stop('model.data must be a character of length one')
  }
  
  custom.models = file.exists(model.data)
  if(!model.data %in% names(mdata) & !custom.models){
    # Automatically retrieve the possible valid options of model.data
    cat('\n'); stop(sprintf("model.data must be one of the following: %s, or a path to custom model.", paste0(names(mdata), collapse = ", ")))
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


