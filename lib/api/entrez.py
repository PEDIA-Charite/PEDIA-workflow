'''
Bindings to query NCBI, Pubmed using Entrez Etools
---
This library is very incomplete, also take caution, that after April 2018
strict ratelimiting will be implemented in etools api calls.
'''
from ratelimit import rate_limited
from Bio import Entrez


class RefSeq:
    '''Extends entrez function of biocommons to expose usable functions.
    '''

    organism = '(Homo sapiens[orgn])'

    database = 'gene'

    def __init__(self, email: str):
        Entrez.email = email

    def _build_query(self, *args):
        '''AND concatenation of all provided terms
        '''
        return " AND ".join(args)

    # limit rate to one request per second
    @rate_limited(1)
    def get_refseq(self, gene: str) -> [str]:
        '''Returns a complicated json, which does contain the reference
        sequences somewhere.
        '''
        response = self.search_gene(gene)
        if response['IdList']:
            # try to get the first entry
            rec = self._get_entry(response['IdList'][0])
            if rec:
                # TODO: actually make returned reference sequences usable
                return rec

    # limit our number of requests to three per second
    def search_gene(self, gene: str) -> dict:
        '''Query the gene database. This returns a gene identifier which can
        be fetched.
        '''
        # limit our search to human genes
        gene = "({}[gene])".format(gene)
        query = self._build_query(gene, self.organism)
        handle = Entrez.esearch(self.database, query)
        record = Entrez.read(handle)
        handle.close()
        return record

    def _get_entry(self, gene_id: str) -> dict:
        '''Get entry from gene database by Gene ID.
        '''
        # TODO fetches should be batched across gene_symbol queries
        handle = Entrez.efetch(db=self.database, id=gene_id, retmode='xml')
        record = Entrez.read(handle)
        handle.close()
        return record


def main():
    '''Only for testing purposes.
    '''
    rs_num = RefSeq('foo@bar.lol')
    response = rs_num.search_gene('PGAP1')
    response = rs_num.get_refseq('PGAP1')
    print(response)


if __name__ == '__main__':
    main()
