#!/usr/bin/env python
'''Search, fetch, and filter SRA packages for chicken.'''

from datetime import date
from sra import SRASearch, SRAPackage, FilterPackages
from email import email_address

# Today:
today = date.today()
# Define basename for output files.
basename = 'paired_gte100bp_%s' % today.strftime('%Y%m%d')

# Search for RNAseq by Illumina later than year 2006.
query = '''(((strategy rna seq[Properties]) AND platform illumina[Properties]) AND chicken[Organism]) AND ("2006/01/01"[Modification Date] : "3000"[Modification Date])'''

# Maximum number of returned results.
retmax = 500

# Instantiate search object.
sra_search = SRASearch(query=query, retmax=retmax, email=email_address)

# Execute search itself.
sra_search.esearch()

# Fetch metadata from packages.
packages = [SRAPackage(sra_id) for sra_id in sra_search.idlist]

# Store packages in data frame for filtering.
package_filter = FilterPackages(packages)

# copy working data frame.
df = package_filter.data_frame

# Filter booleans.
filtered_df = df[df.library_layout == 'PAIRED'][df.nreads > 1][df.read_average >= 100]

# Sort data buy lineage.
sorted_df = filtered_df.sort('lineage')

# Write CSV out.
package_filter.filtered_data_frame = sorted_df
package_filter.write_csv(basename)

# Write unique list of taxa.
unique = package_filter.filtered_data_frame.lineage.unique()
unique.tofile(basename + '_unique' + '.txt', sep='\n')
