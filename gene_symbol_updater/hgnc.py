import numpy as np
import pandas as pd
from datetime import datetime
import requests
from pkg_resources import resource_filename
#from logging import Logger, INFO, CRITICAL, WARNING
import logging
import pickle
#import copy

LOG = logging.getLogger(__name__)

def load_data():
    symbol_ids_table = pd.read_csv(resource_filename(__name__, "data/symbol_ids_table.csv"))
    # convert the NCBI ids to strings to avoid float issue
    s = symbol_ids_table.NCBI_gene_ID.fillna(-1).astype(int).astype(str)
    s[s=='-1'] = ''
    symbol_ids_table['NCBI_gene_ID'] = s

    # load pickles
    objs = {}
    for fn in ['alt_symbols.dict', 'approved.set', 'hg19_ambiguous_mapping.dict']:
        with open(resource_filename(__name__, f"data/{fn}", ), 'rb') as f:
            objs[fn] = pickle.load(f)

    return (symbol_ids_table, objs['alt_symbols.dict'],  objs['approved.set'],
            objs['hg19_ambiguous_mapping.dict'])


symbol_ids_table, alt_symbols, approved, hg19map = load_data()

def hgnc_approved_symbol(g, null=np.nan, map_ambig_with_hg19=True):
    """Return HGNC approved symbol for query, searching previous or alias
    lists. Return null if not found.

    If a symbol matches a current approved symbol, it is returned unchanged.

    With map_ambig_with_hg19=True, symbols that map to multiple previous
    symbols are searched against approved symbols retreved from GRCh37,
    Ensembl release 75. If found, the recent approved symbol is returned.
    """
    if g in approved:
        return g
    try:
        found = alt_symbols[g]
        # ambiguous results return ndarray
        if type(found) is np.ndarray:
            if map_ambig_with_hg19:
                if g in hg19map:
                    found = hg19map[g]
            else:
                LOG.info(f"{g} is ambiguous, maps to {found}.")
        return found
    except KeyError:
        LOG.info(f"{g} not found.")
        return null


def update_hgnc_table(run_update_lookup_lists=True):
    """Download a table from biomart.genenames.org and update file
    used in mapping gene symbols.

    The table is written to data directory by default. Table written
    to data/hgnc_table.{YYYYMMDD}.tsv. update_lookup_lists then
    used to update the symbol mapping files."""

    # download from biomart
    restxml = """<!DOCTYPE Query><Query client="biomartclient" processor="TSV" limit="-1" header="1"><Dataset name="hgnc_gene_mart" config="hgnc_gene_config"><Filter name="hgnc_gene__status_1010" value="Approved" filter_list=""/><Attribute name="hgnc_gene__hgnc_gene_id_1010"/><Attribute name="hgnc_gene__approved_symbol_1010"/><Attribute name="hgnc_gene__approved_name_1010"/><Attribute name="hgnc_gene__hgnc_alias_symbol__alias_symbol_108"/><Attribute name="hgnc_gene__hgnc_previous_symbol__previous_symbol_1012"/><Attribute name="hgnc_gene__chromosome_1010"/><Attribute name="hgnc_gene__locus_group_1010"/><Attribute name="hgnc_gene__locus_type_1010"/><Attribute name="hgnc_gene__hgnc_family__hgnc_family_name_109"/><Attribute name="hgnc_gene__date_symbol_changed_1010"/><Attribute name="hgnc_gene__ensembl_gene__ensembl_gene_id_104"/><Attribute name="hgnc_gene__ncbi_gene__gene_id_1026"/><Attribute name="hgnc_gene__uniprot__uniprot_accession_1036"/></Dataset></Query>"""
    res = requests.get(f'http://biomart.genenames.org/martservice/results?query={restxml}')
    if res.status_code != 200:
        raise RuntimeError(f"Download from biomart.genenames.org failed with code {res.status_code}.")

    # write the file
    out_fn =  f"data/hgnc_table.{datetime.today().strftime('%Y%m%d')}.tsv"
    out_fn = resource_filename(__name__, out_fn)
    with open(out_fn, 'w') as f:
        f.write(res.text)

    if run_update_lookup_lists:
        update_lookup_lists(out_fn)


def update_lookup_lists(hgnc_table_path):
    """Using a table downloaded from HGNC, update the look up lists
    used to quickly map queries to symbols.

    Called by update_hgnc_table."""
    global symbol_ids_table, alt_symbols, approved

    hgnctab = pd.read_csv(hgnc_table_path, sep='\t', )
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
    symbol_ids_table = symbol_ids_table.drop(['Approved_name','Previous_symbol','Alias_symbol'], axis=1)
    # cast the NCBI ids as int then strings to remove decimal point (they'll still get read as floats on the other end)
    symbol_ids_table.loc[:, 'NCBI_gene_ID'] = symbol_ids_table.NCBI_gene_ID[~symbol_ids_table.NCBI_gene_ID.isna()].apply(lambda x: str(int(x)))

    # write the outputs
    symbol_ids_table.to_csv(resource_filename(__name__, "data/symbol_ids_table.csv"))

    for obj, fn in [(both, 'alt_symbols.dict'), (approved_set, 'approved.set')]:
        with open(resource_filename(__name__, f"data/{fn}"), 'wb') as f:
            pickle.dump(obj, f)

    symbol_ids_table, alt_symbols, approved = load_data()

def test_lookup():
    LOG.setLevel('INFO')
    # ambiguous, alias, correct, previous, missing
    genes = ['ASP', 'APOBEC1CF', 'XRCC1', 'MRE11A', 'FAM155B', 'not a gene']
    print([hgnc_approved_symbol(g) for g in genes])

# update_lookup_lists('data/hgnc_table.20210702.tsv')
# test_lookup()
if __name__ == '__main__':
    input("Press enter to update tables, or ctrl-C to cancel.")
    update_hgnc_table()