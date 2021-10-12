

from .hgnc import hgnc_approved_symbol, symbol_ids_table
from .ncbi import entrez_name_id, set_Entrez_email
import pandas as pd
import numpy as np

# this should just return a single DF, with 3rd col indicating status, i.e. MISSING/ENTREZ/HGNC/NOCHANGE
def update_gene_symbols(gset:np.ndarray, email=None) -> dict:
    """Returns a map of original names to udpated, array of ambiguous names,
    and array of genes not found in HGNC or Entrez databases.

    A dictionary is returned with keys:
        'gene_series': pd.Series. indexed by gset with updated names (if found)
            as values. Ambiguous queries will have a list of possible names.
        'ambiguous': np.ndarray. Genes with multiple possible values according
            to HGNC database.
        'no_hits': Genes not found in either the HGNC or NCBI databases.

    Final Series has original name when no other is found

    Args:
        gset: the genes one wishes to be updated
        email: your email to be used when querying NCBI

    """

    if email:
        set_Entrez_email(email)

    found = pd.Series(gset, index=gset).apply(hgnc_approved_symbol)

    ambiguous = found.loc[found.apply(lambda x: (type(x) is list) or (type(x) is np.ndarray))].index
    print('Ambiguous:', len(ambiguous))
    # try NCBI
    nfm = found.isna()
    print('Searching NCBI for:', nfm.sum())
    found.loc[nfm] = found.loc[nfm].index.map(lambda g: entrez_name_id(g, null_value=np.nan)[0])
    # take the raw array so we can apply it to the index directly
    nfm = found.isna().values
    print("No result for :", sum(nfm))
    no_hits = found.index[nfm].values
    found.loc[nfm] = found.index[nfm].values
    return {'genes':found, 'ambiguous':ambiguous, 'no_hits':no_hits}

