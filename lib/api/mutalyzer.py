'''
Bindings for Mutalyzer. A service which can be used to look up RS numbers and
check the validity of hgvs strings.
'''
import requests


class Mutalyzer(requests.Session):
    '''Implements API bindings for the Mutalyzer.
    '''

    base_url = 'https://mutalyzer.nl/json/'

    def __init__(self):
        super().__init__()

    def _call_service(self, method: str, params: dict) -> list:
        '''Call the Mutalyzer API service.
        '''
        url = self.base_url + method
        response = self.get(url, params=params)
        response.raise_for_status()
        return response.json()

    def getdbSNPDescriptions(self, rs_id: str) -> [str]:
        '''Return a list of possible RS numbers for the given RS code.
        '''
        params = {'rs_id': rs_id}
        return self._call_service('getdbSNPDescriptions', params)


def main():
    '''Only for testing purposes.
    '''
    mut = Mutalyzer()
    hgvs_list = mut.getdbSNPDescriptions("rs386834107")
    print(hgvs_list)


if __name__ == '__main__':
    main()
