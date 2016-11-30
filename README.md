<img src="https://cdn.rawgit.com/omarwagih/rmimp/master/inst/extdata/html/images/mimp_logo.svg" alt="rmimp logo" width="250px"><br> Predicting the impact of mutations on kinase-substrate phosphorylation in R
===============================================================

## Installation

To install MIMP, first make sure your R version is at least R 3.0. You can check by typing the following into your R console:

```r
R.Version()$major
```

Next, install and load `devtools` package:

```r
install.packages("devtools")
library("devtools")
```

Then download and install the `rmimp` package from github (**this may take a few minutes, please be patient**):

```r
install_github("omarwagih/rmimp")
```

Load the `rmimp` package into R, and you're ready to go!

```r
library("rmimp")
```

## Running MIMP on sample data:

To start using MIMP, try loading paths to the sample data, which come with the package:

```r
# Get the path to example mutation data 
mut.file = system.file("extdata", "sample_muts.tab", package = "rmimp")

# Get the path to example FASTA sequence file 
seq.file = system.file("extdata", "sample_seqs.fa", package = "rmimp")

# Get the path to example FASTA sequence file 
psite.file = system.file("extdata", "sample_phosphosites.tab", package = "rmimp")
```


The mutation file contains the following lines:

```
TP53 R282W
TP53 R248P
TP53 W146S
CTNNB1 S33C
CTNNB1 S37F
```


The phophosite file contains the following lines:

```
TP53	284
TP53	215
CTNNB1	33
```


To start the analysis, simply run the `mimp` function on these files:

```r
# Run rewiring analysis
results = mimp(mut.file, seq.file, psite.file, display.results=TRUE)
```

The output is stored in the `results` variable and should show up in your browser. To suppress browser display, set `display.results=FALSE`

If you'd like to redisplay the results in your browser at a later time, run the following:

```r
results2html(results)
```

## Running MIMP without phosphosite data:

If you don't pass phosphosite data to the function call, MIMP will use positions of all S, T and Y residues in your FASTA file as potential phosphosites.

```r
# Run rewiring analysis
results2 = mimp(mut.file, seq.file, display.results=TRUE)
```


## Running MIMP without phosphosite or sequence data:
If you pass only a mutation file to the function call, MIMP will use positions of experimentally identified phosphosites from [PhosphoSitePlus](phosphosite.org), and their corresponding sequences. PhosphoSitePlus uses UniProt accessions as identifiers, so for this please make sure the IDs used in your mutation file are uniprot accessions:

Here's an example using a file with five mutations:

```
P04637 R282W
P04637 R248P
P04637 W146S
P35222 S33C
P35222 S37F
```


```r
# Get the path to example mutation data 
mut.file.up = system.file("extdata", "sample_muts_uniprot.tab", package = "rmimp")

# run mimp
results3 = mimp(mut.file.up, display.results = TRUE)
```


## Documentation

For a full list of all options, take a look at the documentation by typing the following in your R console:

```r
?mimp
```

## Contact
If you have any feedback, suggestions or questions, please drop me a line at (wagih(at)ebi.ac.uk) or open an issue on github.