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

Then install and load `MIMP` package from github (**this may take a few minutes, please be patient**):

```r
# Download and install rmimp package
install_github("omarwagih/rmimp")
# Load the package
library("rmimp")
```

## Running MIMP on sample data:

To start using MIMP, try loading paths to the sample data, which comes with the package:
```r
# Get the path to example mutation data 
mut.file = system.file("extdata", "mutation_data.txt", package = "rmimp")

# Get the path to example FASTA sequence data 
seq.file = system.file("extdata", "sequence_data", package = "rmimp")

```

Take a look at the sample data. Try running the following:

```r
browseURL(mut.file)
browseURL(seq.file)
```

Now to start the analysis, simply run the `mimp` function on these files:

```r
# Run rewiring analysis
results = mimp(mut.file, seq.file, display.results=TRUE)
```

The output is stored in the `results` variable and should show up in your browser. To suppress browser display, set `display.results=FALSE`

If you'd like to redisplay the results in your browser at a later time, run the following:
```r
results2html(results)
```

## Documentation

For a full list of all options, take a look at the documentation by typing the following in your R console:

```r
?mimp
```

## Contact
If you have any feedback, suggestions or questions, please drop me a line at (wagih(at)ebi.ac.uk) or open an issue on github.