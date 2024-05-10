### Load packages
import numpy as np
import pandas as pd

### Load file
pheno = snakemake.wildcards["trait"]
df = pd.read_csv(snakemake.input[0], compression='gzip', sep = '\t')
print(df.head(2))

### Rename columns accordingly
df = df.rename(columns = {'rsID':'SNP','Effect_allele':'A1', 'Other_allele':'A2', 'EAF_HRC':'freq', 'P':'p'}) 

### Add N 
df['N'] = 765283

### Effect size standardization
df['b'] = df['Beta']/df['SE']/np.sqrt(df['N'])
df['se'] = 1/np.sqrt(df['N'])

### Export to .ma
df[['SNP', 'A1', 'A2', 'freq', 'b', 'se', 'p', 'N']].to_csv(snakemake.output[0], index = False, sep = '\t')