require(mclust)
require(parallel)
require(ROCR)

#' Get parameters for foreground and background distros
#' and return them as a list to be used later
fgBgParams <- function(pos, neg){
  # Get gmms
  fg.gmm = suppressWarnings(Mclust(pos))
  bg.gmm = suppressWarnings(Mclust(neg))
  
  # and parameters for the models
  fg.params = getParams(fg.gmm)
  bg.params = getParams(bg.gmm)
  
  list(fg.params=fg.params, bg.params=bg.params)
}


rocCurveHelper <- function(pos, neg, pwm, return.auc=T){
  set.seed(123)
  if(length(pos) == 0) return(NA)
  
  # Score positive and negative
  sc.pos = mss(pos, pwm, na_rm = T)
  sc.neg = mss(neg, pwm, na_rm = T)
  
  mypred = prediction( c(sc.neg, sc.pos), c(rep(0, length(sc.neg)), rep(1, length(sc.pos))) )
  auc = performance(mypred, "auc")
  auc = unlist(slot(auc, "y.values"))
  
  if(return.auc) return(auc)
  
  # Else return roc + auc
  roc = performance(mypred, "tpr","fpr")
  return(list(roc=roc, auc=auc))
}

# Splits data into k sections, use 1/k of the data for testing, k-1/k for training 
# @param data - list of data
# @param k - split number
k.fold.split <- function(data, k=5){
  set.seed(111)
  # Randomize 
  data = sample(data)
  
  # Split
  max <- length(data)/k
  x <- seq_along(data)
  d <- split(data, ceiling(x/max))
  
  d = lapply(1:length(d), function(i){
    set = list(train=as.vector(unlist(d[-i])), test=d[[i]])
  })
  
  return(d)
}

cv <- function(pos, neg, k=10, avg=T, avg.fun=mean){
  sp = k.fold.split(pos, k = k)
  
  aucs = sapply(sp, function(x){
    rocCurveHelper(pos = x$test, neg = neg, pwm = PWM(x$train), return.auc = T)
  })
  
  if(avg) return(avg.fun(aucs, na.rm=T))
  return(aucs)
}


setwd('~/Development/rmimp/train/')

# From https://raw.githubusercontent.com/omarwagih/rmimp/master/R/pwm-functions.r
source('pwm-functions.r')

# Get path rmimp data path
BASE_DIR = system.file("extdata", "", package = "rmimp")

# Get path to hconf object
hconf_path = file.path(BASE_DIR, 'kinase_individual_human_experimental.mimp')

# Create a backup
file.copy(hconf_path, file.path(BASE_DIR, 'kinase_individual_human_experimental.mimp.backup'), overwrite = F)

# Pre-trained models
pretrained_models = readRDS(paste0(BASE_DIR, 'kinase_individual_human_experimental.mimp'))


set.seed(10)
models = list()

# Randomly generate 100 sequences with S/T in the center (position 8) and P/R at position +1
pos = sapply(1:100, function(i){
  sp = sample(AA, 15, replace = T)
  sp[8] = sample(c('S', 'T'), 1)
  sp[9] = sample(c('P', 'R'), 1)
  paste0(sp, collapse='')
})

# Randomly generate 1000 sequences with S/T in the center 
neg = lapply(1:100, function(i){
  sp = sample(AA, 15, replace = T)
  sp[8] = sample(c('S', 'T'), 1)
  paste0(sp, collapse='')
})

# Test kinase
kinase_name = 'KINASE_X'
kinase_family = 'FAMILY_X'
pwm = PWM(pos) # construct PWM

pos.score = mss(pos, pwm, na_rm = T)
neg.score = mss(neg, pwm, na_rm = T)
params = fgBgParams(pos.score, neg.score)

# Compute 5-fold cross validation AUC and average 
auc = cv(pos, neg, k=5, avg=T, avg.fun = mean)

# Append to models
models[[kinase_name]] = list(
  name = kinase_name, 
  family = kinase_family, 
  pwm = pwm, 
  fg.params = params$fg.params,
  bg.params = params$bg.params,
  auc = auc
)

# Repeat above for each kinase

# Replace the hconf file - you can always restore this from the backup
saveRDS(models, hconf_path)


