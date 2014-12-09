sra.py
======

Python module to search, fetch and filter records from NCBI's [Sequence Read
Archive](http://www.ncbi.nlm.nih.gov/sra/). In addition to the search options
available on [SRA](http://www.ncbi.nlm.nih.gov/sra/advanced) it makes it
possible to filter by read length or library layout. Outputs a CSV file.

Requirements:

* [Biopython](http://biopython.org/wiki/Main_Page)
* [Pandas](http://pandas.pydata.org/)

Quick and rough:

    ./sra.py -s "agalma[Organism]" -m 3 -o sra_output -e your@email.com

For a more refined filtering check out [fetch_sra.py](fetch_sra.py).

Parsed fields of an SRA record that are exposed to filtering: accession, title,
study title, library strategy, library layout, instrument model, taxon id,
scientific name, taxonomic lineage, run accession, total spots, total bases,
size, published, nreads, read average.
