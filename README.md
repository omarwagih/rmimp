
<img align="center" src="https://github.com/omarwagih/mimp/blob/master/inst/extdata/html/images/mimp_logo.png?raw=trueg" alt="MIMP">

## Installation

To install MIMP, first make sure your R version is at least R 3.0. You can check by typing the following into your R console:

<pre>
R.Version()$major
</pre>

Next, install and load `devtools` package:

<pre>
install.packages("devtools")
library("devtools")
</pre>

Then install and load `MIMP` package from github:

<pre>
install_github("MIMP", "omarwagih")
library("MIMP")
</pre>

## Running MIMP on sample data:

To start using MIMP, try loading paths to the sample data, which comes with the package:
<pre>
# Get the path to example mutation data 
mut.file = system.file("extdata", "mutation_data.txt", package = "MIMP")

# Get the path to example FASTA sequence data 
seq.file = system.file("extdata", "sequence_data", package = "MIMP")

</pre>

To take a look at the sample data, try running the following:

<pre>
browseURL(mut.file)
browseURL(seq.file)
</pre>

Now to start the analysis, simply run the `mimp` function on these files:

<pre>
# Run rewiring analysis
results = mimp(mut.file, seq.file, display.results=TRUE)
</pre>

The output is stored in the `results` variable and should show up in your browser. To suppress browser display, set `display.results=FALSE`

If you'd like to redisplay the results in your browser at a later time, run the following:
<pre>
results2html(results)
</pre>

## Documentation

For a full list of all options, take a look at the documentation by typing the following in your R console:

<pre>
?mimp
</pre>