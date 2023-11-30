import pandas as pd
import pkg_resources
from pkg_resources import resource_filename
from logging import Logger, INFO, CRITICAL, WARNING

LOG = Logger('HGNC_converter')

# todo use INFO logger for things like missing symbols
# todo, provide the full table '/Users/johnc.thomas/Dropbox/crispr/database_stuff/reference_tables/HGNC_noFilter.tsv'
#   but if it's use for symbol matching, privilage tab.Locus_group == 'protein-coding gene'

class MultipleMatchesError(Exception):
    def __init__(self, bad_query, hits):
        self.bq = bad_query
        self.hits = hits

    def __str__(self):
        hits = ', '.join([str(h) for h in self.hits])
        return "Multiple hits for query: {}\nMatches symbols: {}".format(hits, str(self.bq))

class NotFoundError(Exception):
    def __init__(self, message):
        print(message)

class HGNC_Converter:

    def __init__(self, null_value=pd.np.nan, table_path=None):
        """null_value returned when a search fails to find a hit.

        Some previous approved symbols have been used to refer to multiple
        genes in the past, and these are forever ambigous because of this. If
        a symbol matches a current symbol it is assumed to be this, otherwise
        ambiguous symbols return null.

        Use .hgnc_symbol or .hgnc_id with the query sequence.

        prefer_plus_one leads to incorrect matches for at least one common
        gene name, so don't use it.
        """
        if table_path is None:
            fn = resource_filename(__name__, "data/HGNC_table.tsv")
        else:
            fn = table_path
        self.table = None
        # most of these
        self.null = null_value
        self.prev_alias = dict()
        self.ambiguous = set()
        self.symbol_id = dict()
        self.approved = dict()
        self.approved_set = set()

        self.load_table(fn)


    def load_table(self, fn):
        """Load a table with HGNC data. Is expected to be tab separated with
        columns: HGNC ID, Previous symbol, Alias symbol

        Column names will have spaces replaced with underscores"""

        self.table = pd.read_csv(fn, '\t', dtype=str)
        self.table.columns = [c.replace(' ', '_') for c in self.table.columns]
        self.table = self.table.set_index('HGNC_ID', drop=False, )

        for k in 'Previous_symbol', 'Alias_symbol':
            new_ambiguous = set()
            subtab = self.table[['Approved_symbol', k]].dropna()

            # ambigious names will be dealt with later
            for symb, inds in subtab.groupby(k).groups.items():
                if len(inds) > 1:
                    if not all(inds[0] == inds):
                        new_ambiguous.add(symb)

            # mapping alt name to new name
            subtab = subtab.drop_duplicates(k)
            d = dict(zip(subtab[k].str.lower().values, subtab.Approved_symbol.values))

            # deal with ambiguous
            for gn in new_ambiguous:
                d[gn] = pd.np.nan

            self.ambiguous.update(new_ambiguous)

            self.prev_alias[k.split('_')[0]] = d

        symbolId = self.table.set_index('Approved_symbol')['HGNC_ID']
        self.symbol_id = symbolId.drop_duplicates().to_dict()

        # switching to using dicts
        self.approved = {s.lower():s for s in self.table.Approved_symbol.unique()}
        self.approved_set = set(self.approved.values())


    def get_hgnc_symbol(self, query:str, ignore_nonstring=True, return_query_if_missing=False):
        """null_value returned when a search fails to find a hit, unless
        return_query_if_missing is True, in which case the query will be returned instead.

         symbols."""

        origQuery = query
        res = None
        try:
            query = query.lower()
        except AttributeError:
            if not ignore_nonstring:
                raise ValueError('Non-string query: '+str(query))
            else:
                return self.null
        # look for (unique) match in the symbols first
        try:
            return self.approved[query]
        except KeyError:
            pass

        if query in self.ambiguous:
            LOG.info("Multiple hits for query: {}\nMatches symbols: {}".format(query, res))
            return self.null
        else:
            for k in 'Previous', 'Alias':
                try:
                    return self.prev_alias[k][query]
                except KeyError:
                    pass

        if res is None:
            LOG.info("Can't find match for " + str(origQuery))
            if return_query_if_missing:
                return origQuery
            return self.null

    def get_hgnc_id(self, symbol, ignore_nonstring=True):
        if symbol not in self.approved_set:
            symbol = self.get_hgnc_symbol(symbol, ignore_nonstring)
            if symbol is self.null:
                return symbol
        return self.symbol_id[symbol]


def run_test():
    #todo test library, should include all caps matches
    converter = HGNC_Converter()
    gene_list = pd.read_csv('/Users/johnc.thomas/Dropbox/crispr/papers/genetic map of DDR/list_of_connected_genes.txt',
                           header=None)
    print(converter.get_hgnc_symbol('ADPRS'))
    # gene_list.head()
    mapped = gene_list.applymap(converter.get_hgnc_symbol)
    print(mapped)
    mapped_ids = gene_list.applymap(converter.get_hgnc_id)
    print(mapped_ids.head(3))


if __name__ == '__main__':
    LOG.setLevel(INFO)
    run_test()
