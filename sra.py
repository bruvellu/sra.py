#!/usr/bin/env python
'''Search & Fetch records from NCBI's Sequence Read Archive.

Example:

./fetch_sra.py -s 'agalma[Organism]' -m 3 -o sra_output.csv -e your@email.com

Right now filtering needs to be done on a separate script. Either importing the
unfiltered CSV file elsewhere or importing this program as a module and
hard-coding the filter booleans. An example of how to filter is below the CSV
writing line in the main() function.
'''

import argparse
import os
import pandas as pd
import re

from Bio import Entrez


class SRADatabase:
    '''General information about SRA Database'''

    def __init__(self):
        einfo_handle = Entrez.einfo(db='sra')
        einfo = Entrez.read(einfo_handle, validate=False)

        # Define attributes.
        self.count = einfo['DbInfo']['Count']
        self.last_update = einfo['DbInfo']['LastUpdate']
        self.menu_name = einfo['DbInfo']['MenuName']
        self.description = einfo['DbInfo']['Description']
        self.link_list = einfo['DbInfo']['LinkList']
        self.field_list = einfo['DbInfo']['FieldList']
        self.einfo = einfo

        print('\n' + self.description + ' (%s)' % self.menu_name)
        print(self.count + ' entries, last updated in ' + self.last_update)
        print('Available search fields: ')
        for field in self.field_list:
            print('\t[%s] %s | %s (%s terms)' % (field['Name'], field['FullName'], field['Description'], field['TermCount']))

    # TODO Use def __unicode__ or __str__ to identify class objects.


class SRASearch:
    '''Perform search and keep IDs of SRA packages.

    Example of query:

        ((((("strategy rna seq"[Properties]) AND "platform illumina"[Properties])
        AND metazoa[Organism]) NOT vertebrata[Organism]) NOT insects[Organism]) AND
        ("2000/01/01"[Modification Date] : "3000"[Modification Date])
    '''

    def __init__(self, query, retmax, email):
        # Required arguments.
        # TODO Cache search based on query. Store list of IDs!
        self.query = query
        if int(retmax) > 100000:
            # Limit defined by Entrez.
            self.retmax = 100000
        else:
            self.retmax = retmax
        Entrez.email = email

        # Search metadata.
        self.count = None
        self.retstart = None
        self.query_translation = None
        self.idlist = None

        # Additional attributes.
        self.results = None
        self.database = SRADatabase()
        # TODO Add timestamp.

    def esearch(self):
        '''Search SRA packages with Entrez using query.'''
        handle = Entrez.esearch(db='sra', term=self.query, retmax=self.retmax)
        self.results = Entrez.read(handle)
        self.parse_results()
        print('\nSuccess! %s results found.' % self.count)
        print('Your query was: %s' % self.query_translation)
        print('Returned IDs (max=%s): %s' % (self.retmax, ', '.join(self.idlist)))
        return self.results

    def parse_results(self):
        '''Populate class attributes by parsing results.'''
        self.count = self.results['Count']
        self.retstart = self.results['RetStart']
        self.query_translation = self.results['QueryTranslation']
        self.idlist = self.results['IdList']
        print('\nFetched %d package IDs.' % len(self.idlist))


class SRAPackage:
    '''Fetch and store metadata from a SRA package.'''

    def __init__(self, sra_id):
        self.id = sra_id
        self.record = None
        self.cached_filepath = os.path.join('.cache', self.id)

        self.accession = None
        self.title = None
        self.study_title = None
        self.library_strategy = None
        self.library_layout = None
        self.instrument_model = None
        self.taxon_id = None
        self.scientific_name = None
        self.lineage = None
        self.run_accession = None
        self.nreads = None
        self.read_average = None
        self.total_spots = None
        self.total_bases = None
        self.size = None
        self.published = None

        # Define header section for CSV. Must match self.metadata.
        self.header = ['id', 'accession', 'title', 'lineage',
                       'taxon_id', 'scientific_name',
                       'library_strategy', 'library_layout',
                       'instrument_model', 'run_accession',
                       'nreads', 'read_average', 'total_spots',
                       'total_bases', 'size', 'published']

        # Do the actual metadata fetching.
        self.efetch()

        # Retrieve whole lineage by taxon ID.
        self.get_lineage()

        # Fill metadata set for later processing.
        self.metadata = (self.id, self.accession, self.title, self.lineage,
                         self.taxon_id, self.scientific_name,
                         self.library_strategy, self.library_layout,
                         self.instrument_model, self.run_accession,
                         self.nreads, self.read_average, self.total_spots,
                         self.total_bases, self.size, self.published,)

        print('Done!')

    def efetch(self):
        '''Fetch package metadata from Entrez'''
        print('\nProcessing ID=%s' % self.id)
        cached_file = self.cache()
        if cached_file:
            self.record = cached_file.read()
            cached_file.close()
        else:
            print('Record not in cache. Fetching...')
            handle = Entrez.efetch(db='sra', id=self.id)
            self.record = handle.read()
            # Write cache file.
            new_cache = open(self.cached_filepath, 'w')
            new_cache.write(self.record)
            new_cache.close()
        self.extract()

    def cache(self):
        '''Write and read cache files.'''
        # Make sure folder exists.
        cache_folder = '.cache'
        if not os.path.isdir(cache_folder):
            os.mkdir(cache_folder)

        # Try to get cache file.
        try:
            cached = open(self.cached_filepath)
            return cached
        except:
            return None

    def extract(self):
        '''Extract relevant fields from summary.'''

        # Fields with attributes.
        fields = {}

        # Fields to be parsed.
        regexes = {
            'accession': '<EXPERIMENT\s+.*?accession="(?P<accession>.*?)".*?>',
            'title': '<EXPERIMENT\s+.*?>.*?<TITLE>(?P<title>.*?)<\/TITLE>',
            'study_title': '<STUDY_TITLE>(?P<study_title>.*?)<\/STUDY_TITLE>',
            'library_strategy': '<LIBRARY_STRATEGY>(?P<library_strategy>.*?)<\/LIBRARY_STRATEGY>',
            'library_layout': '<LIBRARY_LAYOUT>\s*<(?P<library_layout>SINGLE|PAIRED)',
            'instrument_model': '<INSTRUMENT_MODEL>(?P<instrument_model>.*?)<\/INSTRUMENT_MODEL>',
            'taxon_id': '<TAXON_ID>(?P<taxon_id>.*?)<\/TAXON_ID>',
            'scientific_name': '<SCIENTIFIC_NAME>(?P<scientific_name>.*?)<\/SCIENTIFIC_NAME>',
            'run_accession': '<RUN\s+.*?accession="(?P<run_accession>.*?)"\s+.*?total_spots="(?P<total_spots>.*?)"\s+.*?total_bases="(?P<total_bases>.*?)"\s+.*?size="(?P<size>.*?)"\s+.*?published="(?P<published>.*?)"\s+.*?>',
            'nreads': '<Statistics\s+.*?nreads="(?P<nreads>.*?)"\s+.*?>',
            'read_average': '<Read\s+.*?average="(?P<read_average>.*?)"\s+.*?\/>',
        }

        # Iterate over regexes to parse attributes.
        # TODO handle multiple matches like "runs", "nreads", and "average"?
        # Right now it only gets the first run accession, nreads and
        # read_average. This is OK for now, since it is only a primary filter.
        for field, regex in regexes.iteritems():
            re_search = re.search(regex, self.record)
            if re_search:
                re_groups = re_search.groupdict()
                if re_groups:
                    for k, v in re_groups.iteritems():
                        fields[k] = v
                else:
                    if field in ['taxon_id', 'nreads', 'read_average',
                                 'total_spots', 'total_bases', 'size']:
                        fields[field] = 0
                    else:
                        fields[field] = ''
            else:
                if field in ['taxon_id', 'nreads', 'read_average',
                                'total_spots', 'total_bases', 'size']:
                    fields[field] = 0
                else:
                    fields[field] = ''

        self.accession = fields['accession']
        self.title = fields['title']
        if not self.title:
            self.title = fields['study_title']
        self.library_strategy = fields['library_strategy']
        self.library_layout = fields['library_layout']
        self.instrument_model = fields['instrument_model']
        self.taxon_id = int(fields['taxon_id'])
        self.scientific_name = fields['scientific_name']
        self.nreads = int(fields['nreads'])
        self.read_average = int(float(fields['read_average']))
        self.run_accession = fields['run_accession']
        if self.run_accession:
            self.total_spots = int(fields['total_spots'])
            self.total_bases = int(fields['total_bases'])
            self.size = int(fields['size'])
            self.published = fields['published']
        else:
            self.total_spots = 0
            self.total_bases = 0
            self.size = 0
            self.published = ''

    def get_lineage(self):
        '''Fetch hierarchy from NCBI's Taxonomy database.'''
        # Open taxa cache file.
        try:
            cached_taxa = pd.DataFrame.from_csv('.cache/taxa.csv')
        except:
            cached_taxa = pd.DataFrame(columns=['taxon_id', 'lineage', 'scientific_name'])

        # Fetch row with taxon_id.
        taxon_row = cached_taxa[cached_taxa.taxon_id == self.taxon_id]
        if not taxon_row.empty:
            self.scientific_name = taxon_row.scientific_name.values[0]
            self.lineage = taxon_row.lineage.values[0]
        else:
            print('Taxon %d not in cache. Fetching...' % self.taxon_id)
            handle = Entrez.efetch(db='taxonomy', id=str(self.taxon_id))
            taxon = Entrez.read(handle)
            self.scientific_name = taxon[0]['ScientificName']
            self.lineage = taxon[0]['Lineage'] + '; ' + self.scientific_name
            new_cache = cached_taxa.append([{'taxon_id': self.taxon_id,
                                             'lineage': self.lineage,
                                             'scientific_name': self.scientific_name}])
            new_cache.to_csv('.cache/taxa.csv')


class FilterPackages:
    '''Build data frame with package metadata for filtering.'''

    def __init__(self, packages, filter=None):
        self.packages = packages
        self.data_frame = None
        self.build_data_frame()
        self.filtered_data_frame = self.data_frame

    def build_data_frame(self):
        '''Get metadata from each package and save to data frame'''
        data = []
        index_ids = []
        header = []
        for package in self.packages:
            data.append(package.metadata)
            index_ids.append(package.metadata[0])
            header = package.header
        self.data_frame = pd.DataFrame(data, index=index_ids, columns=header)

    def write_csv(self, basename):
        '''Write CSV file from data frame.'''
        self.data_frame.to_csv(basename + '_unfiltered' + '.csv', index=False)
        self.filtered_data_frame.to_csv(basename + '.csv', index=False)
        print('\n%d of %d packages written to "%s" after filtering.\n' % (self.filtered_data_frame.index.size, self.data_frame.index.size, basename + '.csv'))


def main():
    '''Parse arguments and call SRA search.

    Main function simply parses arguments from command line input and assures
    everything is ok to instantiate the SRA search class.
    '''

    # Parse arguments.
    parser = argparse.ArgumentParser(description='Search & Fetch records from NCBI\'s Sequence Read Archive.',
                                     epilog='Work out those reads, dude.')
    parser.add_argument('-s', '--search',
                        help='put search terms between "quotes"',
                        type=str, required=True)
    parser.add_argument('-m', '--maximum',
                        help='maximum number of records to be retrieved',
                        default='20')
    parser.add_argument('-o', '--output',
                        help='indicate output CSV file',
                        required=True)
    parser.add_argument('-e', '--email',
                        help='an email address is required for Entrez',
                        required=True)
    args = parser.parse_args()

    # Instantiate search object.
    sra_search = SRASearch(query=args.search, retmax=args.maximum,
                           email=args.email)

    # Execute search itself.
    sra_search.esearch()

    # Fetch metadata from packages.
    packages = [SRAPackage(sra_id) for sra_id in sra_search.idlist]

    # Store packages in data frame for filtering.
    packages_to_filter = FilterPackages(packages)

    # Example of filter booleans using Pandas data frame.
    # Only get rows whose 'library_layout' column equals 'PAIRED', 'nreads' is
    # greater than 1, and 'read_average' is greater or equal than 70.
    #filtered_df = df[df.library_layout == 'PAIRED'][df.nreads > 1][df.read_average >= 70]

    # Write CSV out.
    packages_to_filter.write_csv(args.output)

if __name__ == '__main__':
    main()
