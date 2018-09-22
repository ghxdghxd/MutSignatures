# [deconstructSigs](!https://github.com/raerose01/deconstructSigs)

## depend on packages

```R
library(optparse)

Usage: deconstructSigs [options]

Options:
    -i FILE, --input=FILE           the input file with header: sample chr pos ref alt
    -o STR, --outname=STR           the output name
    -p INT, --processNumber=INT     the number of process
    -r FILE, --refSignature=FILE    the file of reference signature
    -c FLOAT, --cutoff=FLOAT        Discard any signature contributions with a weight less than this amount[0.06]
    -h, --help                      Show this help message and exit
```