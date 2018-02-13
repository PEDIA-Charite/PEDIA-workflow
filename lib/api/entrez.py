from Bio import Entrez

class RefSeq:

    organism = '(Homo sapiens[orgn])'

    database = 'gene'

    def __init__(self, email):
        Entrez.email = email

    def _build_query(self, *args):
        '''AND concatenation of all provided terms
        '''
        return " AND ".join(args)

    def get_refseq(self, gene):
        '''Returns a complicated json, which does contain the reference sequences somewhere.
        '''
        r = self.search_gene(gene)
        if r['IdList']:
            # try to get the first entry
            rec = self._get_entry(r['IdList'][0])
            if rec:
                # TODO: actually make returned reference sequences usable
                return rec

    def search_gene(self, gene):
        # limit our search to human genes
        gene = "({}[gene])".format(gene)
        query = self._build_query(gene, self.organism)
        r = Entrez.esearch(self.database, query)
        record = Entrez.read(r)
        r.close()
        return record

    def _get_entry(self, gene_id):
        '''Get entry from gene database by Gene ID.
        '''
        handle = Entrez.efetch(db=self.database, id=gene_id, retmode='xml')
        record = Entrez.read(handle)
        handle.close()
        return record


def main():
    rs = RefSeq('foo@bar.lol')
    r = rs.search_gene('PGAP1')
    r = rs.get_refseq('PGAP1')
    print(r)

if __name__ == '__main__':
    main()
