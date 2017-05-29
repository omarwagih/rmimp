# Valid domains
.VALID_DOMAINS <- c("phos", "sh3", "sh2", "pdz")

#BASE_DIR = system.file("extdata", "", package = "rmimp")

#' Helper function for \code{dohtml}
#'
#' Adds colors for mutants at distance
#'
#' @param s Data frame resulting from mimp call.
#' @param dist Distance of mutation.
#'
#' @keywords helper mimp
.htmlSeq <- function(s, dist) {
  s = strsplit(s, '')[[1]]
  if (dist != 0)
    s[8] = sprintf('<a class="psite">%s</a>', s[8])
  p = 8 + dist
  s[p] = sprintf('<a class="mut">%s</a>', s[p])
  paste0(s, collapse = '')
}

#' Helper function for \code{dohtml} with sh3 binding.
#'
#' Adds colors for mutants at distance
#'
#' @param s Data frame resulting from mimp call.
#' @param dist Distance of mutation.
#'
#' @keywords helper mimp
.htmlSeqSh3 <- function(wt, mut) {
  wtSplit <- strsplit(wt, '')[[1]]
  mutSplit <- strsplit(mut, '')[[1]]
  p <- matrix(c(wtSplit, mutSplit), ncol = 2)
  p <- which(p[, 1] != p[, 2])
  wtSplit[p] = sprintf('<a class="psite">%s</a>', wtSplit[p])
  mutSplit[p] = sprintf('<a class="mut">%s</a>', mutSplit[p])
  wt <- paste0(wtSplit, collapse = '')
  mut <- paste0(mutSplit, collapse = '')
  return(c("wt" = wt, "mut" = mut, "loc" = p))
}


#' Helper function for \code{results2html}
#'
#' @param x Data frame resulting from mimp call.
#' @param LOGO_DIR Directory containing sequence logo images.
#' @param HL_DIR Directory containing overlays
#' @param logoExt Extension of logo files
#' @param .webserver Request coming from webserver?
#'
#' @keywords display mimp
dohtml <- function(x, LOGO_DIR, HL_DIR, logoExt = ".svg", .webserver = F) {
  x = unfactor(x)
  x$score_wt = signif(x$score_wt, 3)
  x$score_mt = signif(x$score_mt, 3)
  x$log_ratio = signif(x$log_ratio, 3)
  x$prob = signif(x$prob, 3)
  
  
  x$log_ratio[is.na(x$log_ratio)] = '-'
  x$score_wt[is.na(x$score_wt)] = '-'
  x$score_mt[is.na(x$score_mt)] = '-'
  
  #x$pwm = gsub('-', '_', x$pwm)
  
  cnt = as.list(table(x$effect))
  n_gain = ifelse(is.null(cnt[['gain']]), 0, cnt[['gain']])
  n_loss = ifelse(is.null(cnt[['loss']]), 0, cnt[['loss']])
  n_mut = nrow(unique(x[, c('gene', 'mut')]))
  
  have.logos = file.exists(LOGO_DIR) | .webserver
  
  lines = sapply(1:nrow(x), function(i) {
    r = x[i,]
    d = r$mut_dist
    seq = sprintf('%s<br>%s', .htmlSeq(r$wt, d), .htmlSeq(r$mt, d))
    logo = file.path(LOGO_DIR, paste0(r$pwm, logoExt))
    hl = file.path(HL_DIR, paste0(r$mut_dist, logoExt))
    
    if (have.logos) {
      t = sprintf(
        '<a class="hide name">%s</a>
        <img style="background: url(%s) no-repeat;"
        src="%s" class="logo" alt="%s" />',
        r$pwm,
        normalizePath(hl),
        normalizePath(logo),
        r$pwm
      )
    } else{
      t = sprintf('<p class="name">%s</p><img src="#" class="logo hide"/>',
                  r$pwm)
    }
    
    if (r$effect == 'loss') {
      eff = '<a class="loss">Loss</a>'
    } else{
      eff = '<a class="gain">Gain</a>'
    }
    
    gene = sprintf(
      '<a target="_blank" class="gene-link" href="http://useast.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=%s">%s</a>',
      r$gene,
      r$gene
    )
    
    # sub - _ logos
    sprintf(
      '<tr> <td class="gene-name">%s</td>
      <td class="psite-pos">%s</td>
      <td class="mut-abbr">%s</td>
      <td class="mut-dist">%s</td>
      <td class="sequence">%s</td>
      <td class="wt-score">%s</td>
      <td class="mt-score">%s</td>
      <td class="prob">%s</td>
      <td class="log-ratio">%s</td>
      <td class="effect">%s</td>
      <td class="seq-logo">%s</td>
      </tr>',
      gene,
      r$psite_pos,
      r$mut,
      r$mut_dist,
      seq,
      r$score_wt,
      r$score_mt,
      r$prob,
      r$log_ratio,
      eff,
      t
    )
    
  })
  
  tt = '<div id="%s" style="display:none">%s</div>'
  lines = append(lines,
                 c(
                   sprintf(fmt = tt, 'n_mut', n_mut),
                   sprintf(fmt = tt, 'n_gain', n_gain),
                   sprintf(fmt = tt, 'n_loss', n_loss)
                 ))
  return(lines)
}

#' Helper function for \code{results2html} and sh3 binding
#'
#' @param x Data frame resulting from mimp call.
#' @param LOGO_DIR Directory containing sequence logo images.
#' @param HL_DIR Directory containing overlays
#' @param logoExt Extension of logo files
#' @param .webserver Request coming from webserver?
#'
#' @keywords display mimp
dohtmlSh3 <- function(x, LOGO_DIR, HL_DIR, logoExt = ".svg", .webserver = F) {
  x = unfactor(x)
  x$score_wt = signif(x$score_wt, 3)
  x$score_mt = signif(x$score_mt, 3)
  x$log_ratio = signif(x$log_ratio, 3)
  x$prob = signif(x$prob, 3)
  
  x$log_ratio[is.na(x$log_ratio)] = '-'
  x$score_wt[is.na(x$score_wt)] = '-'
  x$score_mt[is.na(x$score_mt)] = '-'

  cnt = as.list(table(x$effect))
  n_gain = ifelse(is.null(cnt[['gain']]), 0, cnt[['gain']])
  n_loss = ifelse(is.null(cnt[['loss']]), 0, cnt[['loss']])
  n_mut = nrow(unique(x[, c('gene', 'mut')]))
  
  have.logos = file.exists(LOGO_DIR) | .webserver
  
  lines = sapply(1:nrow(x), function(i) {
    r = x[i,]
    seq = .htmlSeqSh3(r$wt, r$mt)
    mut_loc = seq[["loc"]]
    seq = sprintf('%s<br>%s', seq[["wt"]], seq[["mut"]])
    pwm_length = nchar(r$wt)
    if (is.element(pwm_length, 12:14)) {
      pwm_length = "12-14"
    }
    logo = file.path(LOGO_DIR, paste0(r$pwm, logoExt))
    hl = file.path(HL_DIR, pwm_length, paste0(mut_loc, ".svg"))
    
    if (have.logos) {
      t = sprintf( '<a class="hide name">%s</a><img style="background: url(%s) no-repeat;" src="%s" class="logo" alt="%s" />', 
                   r$pwm, normalizePath(hl), normalizePath(logo), r$pwm)
    } else{
      t = sprintf('<p class="name">%s</p><img src="#" class="logo hide"/>', r$pwm)
    }
    
    if (r$effect == 'loss') {
      eff = '<a class="loss">Loss</a>'
    } else{
      eff = '<a class="gain">Gain</a>'
    }
    
    gene = sprintf('<a target="_blank" class="gene-link" href="http://useast.ensembl.org/Homo_sapiens/Transcript/Summary?db=core;t=%s">%s</a>', r$gene, r$gene)
    
    # sub - _ logos
    sprintf(
      '<tr>
      <td class="gene-name">%s</td>
      <td class="mut-abbr">%s</td>
      <td class="sequence">%s</td>
      <td class="wt-score">%s</td>
      <td class="mt-score">%s</td>
      <td class="prob">%s</td>
      <td class="log-ratio">%s</td>
      <td class="effect">%s</td>
      <td class="seq-logo">%s</td>
      </tr>', gene, r$mut, seq, r$score_wt, r$score_mt, r$prob, r$log_ratio, eff, t)
  })
  
  tt = '<div id="%s" style="display:none">%s</div>'
  lines = append(lines, c(sprintf(fmt = tt, 'n_mut', n_mut), 
                          sprintf(fmt = tt, 'n_gain', n_gain), sprintf(fmt = tt, 'n_loss', n_loss)))
  return(lines)
}

#' Display MIMP results interactively in browser
#'
#' @param x Data frame resulting from mimp call.
#' @param domain Which binding domain to run mimp for
#' @param max.rows If data contains more rows than this value, results won't be displayed.
#'
#' @keywords display mimp
#' @export
#'
results2html <- function(x, domain = "phos", max.rows = 5000) {
  # Ensure valid domain
  if (!is.element(domain, .VALID_DOMAINS)) {
    stop("domain must be valid. Please check mimp documentation for a list of valid domains.")
  }
  
  # Ensure proper size of rows
  if (nrow(x) > max.rows) {
    stop(sprintf('Rows of resulting data exceeds %s and cannot be displayed',max.rows))
  }
  
  # Ensure the input data is generated by mimp
  model.data = attr(x, 'model.data')
  if (is.null(model.data)) {
    stop('Input data must be generated using the mimp function. Please try again.')
  }
  
  # Set possible names of logo directories
  if (domain == "phos") {
    logoDirNames = c(
      'hconf' = 'logos',
      'hconf-fam' = 'logos_fam',
      'lconf' = 'logos_newman')
  } else if (domain == "sh3") {
    logoDirNames = c('hconf' = 'logos_sh3')
  }

  # Check if the given model is valid
  if (!is.element(model.data, names(logoDirNames))) {
    stop('Input data must be generated using a valid model. Please check the model.data parameter of your mimp function call.')
  }
  
  # Get file paths
  LOGO_DIR = file.path(BASE_DIR, 'html', 'images', logoDirNames[model.data])
  
  # Generate result html from template
  if (domain == "phos") {
    HL_DIR = file.path(BASE_DIR, 'html', 'images', 'highlight')
    lines = dohtml(x, LOGO_DIR, HL_DIR)
    template = readLines(file.path(BASE_DIR, 'html', 'index.html'))
  } else if (domain == "sh3") {
    HL_DIR = file.path(BASE_DIR, 'html', 'images', 'highlight_sh3')
    lines = dohtmlSh3(x, LOGO_DIR, HL_DIR)
    template = readLines(file.path(BASE_DIR, 'html', 'index_sh3.html'))
  }
  
  predictionType = sprintf('<div id="predictionType" style="display:none">%s</div>', domain)
  outputPath = file.path(BASE_DIR, 'html', 'MIMP_results.html')
  resultHtml <- file(outputPath, "w")
  ind = grep('<DATA>', template)
  
  # Write result html to file
  writeLines(template[1:(ind - 1)], con = resultHtml, sep = "\n")
  writeLines(lines, con = resultHtml, sep = "\n")
  writeLines(predictionType, con = resultHtml, sep = "\n")
  writeLines(template[(ind + 1):length(template)], con = resultHtml, sep = "\n")
  close(resultHtml)
  
  # Open the result html with web browser
  browseURL(outputPath)
}