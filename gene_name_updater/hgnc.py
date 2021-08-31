import numpy as np
import pandas as pd
import pkg_resources
from pkg_resources import resource_filename
from logging import Logger, INFO, CRITICAL, WARNING
import pickle
import copy

LOG = Logger(__name__)

def _load_data():
    symbol_ids_table = pd.read_csv(resource_filename(__name__, "data/symbol_ids_table.csv"))
    _objs = {}
    for fn in ['alt_symbols.dict', 'approved.set']:
        with open(resource_filename(__name__, f"data/{fn}", ), 'rb') as f:
            _objs[fn] = pickle.load(f)
    alt_symbols = _objs['alt_symbols.dict']
    approved = _objs['approved.set']
    return symbol_ids_table, alt_symbols, approved

symbol_ids_table, alt_symbols, approved = _load_data()

def hgnc_approved_symbol(g, null=np.nan):
    """"""
    if g in approved:
        return g
    try:
        found = alt_symbols[g]
        if type(found) is np.ndarray:
            LOG.info(f"{g} is ambiguous, maps to {found}.")
        return found
    except KeyError:
        LOG.info(f"{g} not found.")
        return null


def update_lookup_lists(hgnc_table_path):

    hgnctab = pd.read_csv(hgnc_table_path, '\t', )
    hgnctab.columns = hgnctab.columns.map(lambda x: x.replace(' ', '_'))

    approved_hgnc = hgnctab['Approved_symbol'].dropna().unique()
    alias_hgnc = hgnctab['Alias_symbol'].dropna().unique()
    prev_hgnc = hgnctab['Previous_symbol'].dropna().unique()

    # unique symbol->IDs mapping
    hgnc_approved_ids = hgnctab.drop_duplicates('Approved_symbol').set_index('Approved_symbol')
    hgnc_approved_ids.loc[:, 'NCBI_gene_ID'] = hgnc_approved_ids.NCBI_gene_ID.astype(pd.Int64Dtype())
    hgnc_approved_ids.head()

    alt_symbol_sets = {'Previous_symbol':prev_hgnc, 'Alias_symbol':alias_hgnc}

    # get reference that will allow idenfication of symbols that appear as alt to multiple approved symbols
    ambig_groups = {}
    for k in 'Previous_symbol', 'Alias_symbol':
        subtab = hgnctab[['Approved_symbol', k]].dropna().drop_duplicates()

        # ambigious names will be dealt with later
        ambig_groups[k] =  subtab.groupby(k).groups

    def map_symbols_get_ambiguous(gset, alt_symbol_column):
        """Give mapping for a gene set from old to Approved_symbol, identify some ambiguous
        genes (those that appear multiple times in Alias/Previous)"""
        mapping = pd.Series(index=gset)
        ambig = set()
        for g in gset:

            if len(ambig_groups[alt_symbol_column][g]) > 1:
                ambig.add(g)

            approved = hgnctab.loc[hgnctab[alt_symbol_column] == g, 'Approved_symbol'].unique()
            # len > 1
            if len(approved) == 1:
                approved = approved[0]
            mapping[g] = approved
        return mapping, ambig

    confirmed_ambig = set()

    # go through previous and alias symbol lists, get mappings to Approved and identify ambiguous symbols
    alt_approved_map = {}
    for col, gset in alt_symbol_sets.items():
        mappings, ambig = map_symbols_get_ambiguous(gset, col)
        confirmed_ambig.update(ambig)
        alt_approved_map[col] = mappings

    am, pm = alt_approved_map['Alias_symbol'], alt_approved_map['Previous_symbol']
    # identify symbols that appear in both mappings, "crossed"
    m = am.index.isin(pm.index)
    crossed_genes = alt_approved_map['Alias_symbol'].index[m]

    for g in crossed_genes:
        ag, pg = am[g], pm[g]
        if type(ag) != np.ndarray:
            ag = np.array([ag])
        if type(pg) != np.ndarray:
            pg = np.array([pg])
        comb = np.concatenate((ag, pg))
        # crossed genes added to a list and the mappings both go to the list
        am.loc[g] = comb
        pm.loc[g] = comb
        confirmed_ambig.add(g)

    # get objects for writing
    # combine the previous/alias lists as there should be no conflict now
    both = {**am, **pm}
    approved_set = set(approved_hgnc)

    # get ID table, mapping approved symbol to IDs provided by HGNC, inc NCBI and Ensembl
    symbol_ids_table = hgnctab.drop_duplicates('Approved_symbol').set_index('Approved_symbol', drop=False)
    symbol_ids_table = symbol_ids_table.drop(['Approved_name','Previous_symbol','Alias_symbol'], 1)
    # cast the NCBI ids as int then strings to remove decimal point (they'll still get read as floats on the other end)
    symbol_ids_table.loc[:, 'NCBI_gene_ID'] = symbol_ids_table.NCBI_gene_ID[~symbol_ids_table.NCBI_gene_ID.isna()].apply(lambda x: str(int(x)))

    # write the outputs
    symbol_ids_table.to_csv(resource_filename(__name__, "data/symbol_ids_table.csv"))

    for obj, fn in [(both, 'alt_symbols.dict'), (approved_set, 'approved.set')]:
        with open(resource_filename(__name__, f"data/{fn}"), 'wb') as f:
            pickle.dump(obj, f)

def test():
    LOG.setLevel('INFO')
    # ambiguous, alias, correct, previous, missing
    genes = ['ASP', 'APOBEC1CF', 'XRCC1', 'MRE11A', 'not a gene']
    print([hgnc_approved_symbol(g) for g in genes])

# update_lookup_lists('data/hgnc_table.20210702.tsv')
#test()