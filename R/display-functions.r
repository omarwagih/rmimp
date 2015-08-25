
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
  if(dist != 0) s[8] = sprintf('<a class="psite">%s</a>', s[8])
  p = 8 + dist
  s[p] = sprintf('<a class="mut">%s</a>', s[p])
  paste0(s, collapse='')
}


#' Helper function for \code{results2html}
#'
#' @param x Data frame resulting from mimp call.
#' @param LOGO_DIR Directory containing sequence logo images.
#' @param HL_DIR Directory containing overlays
#' @param .webserver Request coming from webserver?
#' 
#' @keywords display mimp
#' @export
dohtml <- function(x, LOGO_DIR, HL_DIR, .webserver=F){
  x = unfactor(x)
  x$score_wt = signif(x$score_wt, 3)
  x$score_mt = signif(x$score_mt, 3)
  x$log_ratio = signif(x$log_ratio, 3)
  x$prob = signif(x$prob, 3)
  
  
  x$log_ratio[is.na(x$log_ratio)] = '-'
  x$score_wt[is.na(x$score_wt)] = '-'
  x$score_mt[is.na(x$score_mt)] = '-'
  
  #x$pwm = gsub('-', '_', x$pwm)
  
  cnt = as.list( table(x$effect) )
  n_gain = ifelse( is.null( cnt[['gain']] ), 0, cnt[['gain']])
  n_loss = ifelse( is.null( cnt[['loss']] ), 0, cnt[['loss']])
  n_mut = nrow(unique(x[,c('gene', 'mut')]))
  
  have.logos = file.exists(LOGO_DIR) | .webserver
  
  lines = sapply(1:nrow(x), function(i){
    r = x[i,]
    d = r$mut_dist
    seq = sprintf('%s<br>%s', .htmlSeq(r$wt, d), .htmlSeq(r$mt, d))
    scr = sprintf('%s<br>%s', r$score_wt, r$score_mt)
    
    logo = file.path(LOGO_DIR, paste0(r$pwm, '.svg'))
    hl = file.path(HL_DIR, paste0(r$mut_dist, '.svg'))
    
    
    if(have.logos){
      t = sprintf('<a class="hide name">%s</a>
                <img style="background: url(%s) no-repeat;"
                src="%s" class="logo" alt="%s" />', 
                  r$pwm, normalizePath(hl),  normalizePath(logo), r$pwm)
    }else{
      t = sprintf('<p class="name">%s</p><img src="#" class="logo hide"/>', r$pwm)
    }
    
    if(r$effect == 'loss'){
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
            <td class="prob">%s</td>
            <td class="log-ratio">%s</td>
            <td class="effect">%s</td>
            <td class="seq-logo">%s</td>
            </tr>', 
            r$gene, r$psite_pos, r$mut, r$mut_dist, seq, r$score_wt, r$score_mt, r$prob, r$log_ratio, eff, t)
    
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
  #BASE_DIR = system.file("extdata", "", package = "rmimp")
  
  model.data = attr(x, 'model.data')
  if(is.null(model.data)) stop('Input data must be generated using the mimp function. Please try again.')
  z = c('hconf'='logos', 'hconf-fam'='logos_fam', 'lconf'='logos_newman')
  
  LOGO_DIR = file.path(BASE_DIR, 'html', 'images', z[model.data])
  HL_DIR = file.path(BASE_DIR, 'html', 'images', 'highlight')
  lines = dohtml(x, LOGO_DIR, HL_DIR)
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