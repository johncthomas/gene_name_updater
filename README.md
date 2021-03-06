# gene_name_updater
For updating old human gene symbols to the latest approved symbols.

`update_gene_symbols()` is the primary function. It searches first the HGNC table, then queries NCBI for missing symbols. Returns a dictionary with genes mapped to updated names, a list of ambiguous symbols, and a list of symbols that didn't hit anything.
  
  `hgnc_approved_symbol()`. Approved symbols sourced from table downloaded from https://genenames.org, packaged as efficient data structures to quickly look up the approved symbols for known alias or previous approved symbols. Returns a null value (default np.nan) when query is not found, and a list of symbols when the query is ambiguous.
  
  `entrez_name_id()`. Queries the Entrez database to find the symbol and ID used by NCBI.
  
  `symbol_ids_table`. Pandas DataFrame giving mapping of approved symbol to various IDs (NCBI, Entrez, HGNC). Symbols other than HGNC approved are not included.
  
  There are also function for updating reference tables, documented in the code.
