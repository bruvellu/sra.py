#!/usr/bin/env python
'''Search, fetch, and filter SRA packages.

Example script that queries SRA records for RNA-Seq data from Gallus gallus
(chicken) sequenced with Illumina and published after 2006. Results are then
filtered to select only paired end libraries with read length of at least 70 bp.

Querying and fetching done with Biopython and filtering with Pandas.
'''

from datetime import date
from sra import SRASearch, SRAPackage, FilterPackages

# Today:
today = date.today()

# Email required by NCBI.
email = 'your@email.com'

# Define basename for output files.
basename = 'paired_gte70bp_%s' % today.strftime('%Y%m%d')

# Search query. Use http://www.ncbi.nlm.nih.gov/sra/advanced to build yours.
query = '''(((strategy rna seq[Properties]) AND platform illumina[Properties]) AND chicken[Organism]) AND ("2006/01/01"[Modification Date] : "3000"[Modification Date])'''

# Maximum number of returned results.
retmax = 100

# Instantiate search object.
sra_search = SRASearch(query=query, retmax=retmax, email=email)

# Execute search itself.
sra_search.esearch()

# Fetch metadata from packages.
packages = [SRAPackage(sra_id) for sra_id in sra_search.idlist]

# Store packages in data frame for filtering.
package_filter = FilterPackages(packages)

# Copy working data frame.
df = package_filter.data_frame

# Filter booleans.
filtered_df = df[df.library_layout == 'PAIRED'][df.nreads > 1][df.read_average >= 70]

# Sort data buy lineage.
sorted_df = filtered_df.sort('lineage')

# Write CSV out.
package_filter.filtered_data_frame = sorted_df
package_filter.write_csv(basename)

# Write unique list of taxa.
unique = package_filter.filtered_data_frame.lineage.unique()
unique.tofile(basename + '_unique' + '.txt', sep='\n')
