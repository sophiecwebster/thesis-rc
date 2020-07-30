import ipyrad.analysis as ipa
import pandas as pd
import toyplot
data = '/n/holyscratch01/hopkins_lab/webster/radseq/ipyRAD_assembly02/assembly2_outfiles/assembly2.snps.hdf5'

# group individuals into populations
imap = {
    "arb": ["1011L001_", "1011L002_", "1019L001_", "1019L002_", "1027L001_", "1027L002_", "1035L001_", "1035L002_"],
    "hf": ["1042L001_", "1042L002_", "1051L001_", "1051L002_", "1059L001_", "1059L002_", "1067L001_", "1067L002_", "1075L001_", "1075L002_"]
}

# require that 50% of samples have data in each group
minmap = {i: 0.5 for i in imap}

# initialize pca object with input data and optional parameters options
pca = ipa.pca(
    data=data,
    imap=imap,
    minmap=minmap,
    mincov=0.75,
    impute_method="sample"
)

# run the PCA analysis
pca.run()

# store
