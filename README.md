sra.py
======

Python module to search, fetch and filter records from NCBI's [Sequence Read
Archive](http://www.ncbi.nlm.nih.gov/sra/). In addition to the search options
available on [SRA](http://www.ncbi.nlm.nih.gov/sra/advanced) it is possible to
filter by read length and library layout. Output to CSV file.

Requirements:

* [Biopython](http://biopython.org/wiki/Main_Page)
* [Pandas](http://pandas.pydata.org/)

Quick and rough:

    ./sra.py -s 'agalma[Organism]' -m 3 -o sra_output -e your@email.com

For a more refined filtering check out [fetch_sra.py](fetch_sra.py).
