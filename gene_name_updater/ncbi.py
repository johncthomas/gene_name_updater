import typing
from io import BytesIO
from urllib.request import urlopen
from pkg_resources import resource_filename

import pandas as pd
from Bio import Entrez

import gzip
from time import sleep
import logging, shutil
LOG = logging.getLogger('entrez_converter')
LOG.setLevel('WARNING')

# todo WITHDRAWN is a thing, different from discontinued?
# todo do something better with multiple hits? Currently returns nothing.

"""
efetch for 101362076 includes:

Entrezgene_comments:
[
    Gene-commentary_type:  254
Gene-commentary_heading:  RefSeq Status
Gene-commentary_label:  WITHDRAWN
]"""
# todo sort out how discontinued is dealt with

setUrEmail = 'Use gene_name_updater.ncbi.set_Entrez_email("yourmail@here.com") before searching'
print(setUrEmail)

fn_oldId = resource_filename(__name__, "data/NCBI_oldId_to_newId.csv")
ncbiOldIdTable = None
def load_oldIdTable(path):
    global ncbiOldIdTable
    ncbiOldIdTable = pd.read_csv(path, index_col='OldId', dtype=str)
load_oldIdTable(fn_oldId)


def update_ncbiOldIdTable():
    """

    Downloads file https://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz, filters
    for human genes and writes the output. This

    Overwrites file to package data folder as NCBI_oldId_to_newId.csv, old one
    is backed up as to NCBI_oldId_to_newId.csv.old. If function fails the new
    file will be incomplete.

    Retired IDs with no new ID are written as '{oldId},DISCONTINUED'
    set HGNC_converter.entrez.ncbiOldIdTable = """
    # remember to check any changes here against all times ncbiOldIdTable is used
    global ncbiOldIdTable

    try:
        shutil.copy(fn_oldId, fn_oldId+'.old')
    except FileNotFoundError:
        pass
    print('Writing to', fn_oldId)

    response = urlopen('https://ftp.ncbi.nih.gov/gene/DATA/gene_history.gz')
    reader = gzip.open(BytesIO(response.read()), 'r')

    with reader as f:
        header = next(f)
        header = header.decode('utf-8')
        if header != '#tax_id	GeneID	Discontinued_GeneID	Discontinued_Symbol	Discontinue_Date\n':
            print(header)
            raise RuntimeError("They've changed the header, check it. New file Not created.")
        with open(fn_oldId, 'w') as out_f:
            out_f.write('OldId,GeneId\n')
            for line in f:

                line = line.decode('utf-8')
                spline = line.split('\t')
                if spline[0] == '9606':
                    if spline[1] == '-':
                        out_f.write(f"{spline[2]},DISCONTINUED\n")
                    else:
                        out_f.write(f"{spline[2]},{spline[1]}\n")
    load_oldIdTable(fn_oldId)


def set_Entrez_email(email):
    assert "@" in email
    Entrez.email = email

#todo pep8 these names
def entrez_name_id(query, fullResultsOnFail=False, null_value ='', extra_sleep=0,
                   taxid='9606') -> (typing.Any, str):
    """Returns tuple of (name, NCBI-ID) if a single alive record is found.
    The ID field can have values 'discontinued', 'multiple hits' & 'no hits'
    indicating failed searches. The name will have `null_value` in these cases.

    """
    LOG.debug('Query = '+str(query))
    if not Entrez.email:
        raise RuntimeError(setUrEmail)

    ids = []

    # deal with entrez ID derived names, pull the record with this ID
    if query.startswith('LOC'):
        # check it ends in a number
        try:
            idn = int(query[3:])
            ids = [str(idn)]
        # if not it's some other alias probably
        except ValueError:
            pass

    # do a search if it's not a LOC id
    if not ids:
        fullQuery = f'{query}[Gene Name] AND {taxid}[Taxonomy ID]'
        # print(fullQuery)
        sleep(extra_sleep)

        handle = Entrez.esearch('gene', fullQuery)
        res = Entrez.read(handle)
        handle.close()
        ids = res['IdList']

    # deal with readthrough names
    if not ids:
        # it might not recognise a particular readthrough of outdated queries
        if '-' in query:
            bits = []
            for subq in query.split('-'):
                # get valid entrez name if one
                bits.append(entrez_name_id(subq, fullResultsOnFail=False)[0])
            if all(bits):
                return entrez_name_id('-'.join(bits))

    # if still nothing...
    if not ids:
        LOG.info(query + ' was not found')
        return null_value, 'no hits'

    names = []
    reses = []
    goodIds = []
    discontinued = False
    LOG.debug('IDs = {}'.format(ids) )
    for idd in ids:
        handle = Entrez.efetch(db='gene', id=idd, retmode='xml')
        sleep(extra_sleep)
        res = Entrez.read(handle)
        reses.append(res)
        # if the id is discontinued, get the new one
        if 'Gene-track_discontinue-date' in res[0]['Entrezgene_track-info']['Gene-track']:
            discontinued = True
            #newId = newIdWhenDiscontinued(res[0])
            try:
                newId = ncbiOldIdTable.loc[int(idd), 'GeneId']
            # if we get a new id, add it to the list
            except KeyError:
                LOG.warning(f'Discontinued ID {idd} not in data/NCBI_oldId_to_newId.csv, it needs updating.')
                continue
            if newId in ids:
                continue
            if newId == 'DISCONTINUED':
                LOG.info(idd+'discontinued without replacement')
            else:
                ids.append(newId)
                continue
        try:
            names.append(res[0]['Entrezgene_gene']['Gene-ref']['Gene-ref_locus'])
            goodIds.append(idd)
        except KeyError:
            print(idd, "does not have ['Entrezgene_gene']['Gene-ref']['Gene-ref_locus']")
            pass
    if (len(names) != 0) and all([names[0] == n for n in names]):
        # entrezNames[query] = names[0]
        LOG.info('ID for '+query+' was found, returning '+names[0])
        return names[0], goodIds[0]
    if len(names) == 0:
        if len(reses) == 1 and discontinued:

            LOG.info(query+' discontinued with no replacement')
            if fullResultsOnFail:
                return null_value, reses
            else:
                return null_value, 'discontinued'
        else:
            # I don't think this will ever be reached
            LOG.info(query + ' was not found')
            return null_value, 'no hits'
    else:
        LOG.info(query + ' returned multiple valid hits')
        return null_value, 'multiple hits'


def keys_crawler(d, level=0):
    """Accepts a dictionary, gets the type of values for each key,
    if it's a list the type of the FIRST item obtained. Lists and dicts
    encountered are explored recursively.

    Output gives the structure of the object as a string, levels of
    indentation show the hierarchy. Lists indicated with [...]"""

    dictLike = Entrez.Parser.DictionaryElement
    listLike = Entrez.Parser.ListElement
    out = ''
    for k in d.keys():
        out += '  '*level+k+':\n'
        if type(d[k]) is dict or type(d[k]) is dictLike:
            out += keys_crawler(d[k], level=level+1)
        elif type(d[k]) is listLike:
            out += '  '*(level+1)+'[\n'
            if type(d[k][0]) is dictLike:
                out += keys_crawler(d[k][0], level=level+2)
            out += '  '*(level+1)+']\n'
        else:# type(d[k]) is str:
            out = out[:-1] + f'  {d[k]}\n'
    return out

def print_efetch_results(res):
    for i, r in enumerate(res):
        print(f'Record {i}')
        print(keys_crawler(r))

def run_test():

    gns = ['LOC653602', 'LOC100289561', 'LOC554223', 'WI2-2373I1.2', 'LOC403312', 'LOC101928841', 'TARP',
           'LOC101928093', 'LOC113230', 'LOC100130357', 'MAGEA10-MAGEA5', 'PHOSPHO2-KLHL23', 'TNFAIP8L2-SCNM1',
           'LOC728392', 'LOC100996634', 'LOC388436', 'ZNF664-FAM101A', 'LOC401052', '-09sdf0qnroif-sod']

    results = []

    for gn in gns:
        results.append(
            entrez_name_id(gn)
        )

    return results

if __name__ == '__main__':
    input("Press enter to update NCBI discontinued table, or ctrl-C to cancel.")
    update_ncbiOldIdTable()

