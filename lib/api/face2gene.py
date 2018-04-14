'''
Use the Face2Gene library to get syndrome information.

'''
import logging
from typing import Union

import hashlib
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

LOGGER = logging.getLogger(__name__)


class Face2Gene(requests.Session):
    '''Implements login into face2gene to access f2g services.
    '''

    base_url = 'https://app.face2gene.com/'

    def __init__(self, user: str='', password: str='',
                 config: 'ConfigParser'=None):
        super().__init__()
        if config:
            user = config.face2gene['user']
            password = config.face2gene['password']
        retries = Retry(total=5, backoff_factor=1, status_forcelist=[429, 500])
        self.mount('http://', HTTPAdapter(max_retries=retries))
        # explicitly set user-agent to identify as a known scraper
        self.headers.update({'User-Agent': 'pediaScript'})
        self._login(user, password)

        # cache previous searches by keywords used
        self._cache = {}

    def goto_library(self):
        '''Go to specific library page. This is not needed for search reqests.
        Since XCSRF is not enforced for these.
        '''
        url = self.base_url + 'library/search'
        response = self.get(url)
        return response

    def search_library(self, query: str) -> dict:
        '''Search for a string in Face2Gene and return the resulting
        search dictionary.
        '''
        if query in self._cache:
            return self._cache[query]

        url = self.base_url + 'library/search_api'
        params = {'keywords': query}
        search_results = self.get(url, params=params).json()
        self._cache[query] = search_results
        return search_results

    def browse_library(self, field: str, key: str) -> list:
        '''Browse entries in the face2gene library for given key,
        returning all entries by increasing the offset.'''
        url = self.base_url + "library/browse"
        params = {
            "field": field,
            "key": key,
            "from": 0
        }
        received = []
        total = 1
        while len(received) < total:
            while True:
                try:
                    response = self.get(url, params=params)
                    data = response.json()
                except:
                    LOGGER.warning("F2G error occurred %s", response.text)
                else:
                    break

            total = data["total"]
            received += data["data"]
            # increase the offset by the number of received elements
            params["from"] += len(data["data"])
        return received

    def browse_syndromes(self, key: str) -> list:
        '''Browse syndromes.'''
        syndromes = self.browse_library(field="syndromes", key=key)
        return syndromes

    def browse_all_syndromes(self) -> list:
        '''Browse all syndromes from A-Z.'''
        # crawl alphabet and get the empty key entries
        syndrome_keys = "ABCDEFGHIJKLMNOPQRSTUVWXYZ#"
        syndromes = []
        for letter in syndrome_keys:
            syndromes += self.browse_syndromes(key=letter)
        return syndromes

    def search_syndrome(self, query: str, omim_list: list=[],
                        return_first: bool=True) -> Union[str, list]:
        '''Search for syndromes on Face2Gene library.
        Args:
            query: Query string, such as the syndrome name.
            omim_list: A list of given omim ids, against which search results
                       will be matched.
            return_first: The first matched result will be returned. Otherwise
                          a list of all filtered results will be returned.
        '''
        result = self.search_library(query)
        syndromes = result['syndromes']['data']
        # return the first entry with a omim id
        for syndrome in syndromes:
            if not syndrome['omim_id']:
                continue
            if str(syndrome['omim_id']) != '0':
                if omim_list:
                    if int(syndrome['omim_id']) in omim_list:
                        return syndrome['omim_id']
                else:
                    return syndrome['omim_id']
        return ''

    def _login(self, user: str, password: str):
        '''Login to Face2Gene using md5 hashed password.
        '''
        payload = {
            "email": user,
            "password": hashlib.md5(password.encode('utf-8')).hexdigest(),
            "screenHeight": 875,
            "screenWidth": 1167,
            "userAgent": (
                "Mozilla/5.0 (Windows NT 6.1; Win64; x64) "
                "AppleWebKit/537.36 (KHTML, like Gecko) "
                "Chrome/40.0.2214.85 Safari/537.36"
            )
        }
        r = self.post('https://app.face2gene.com/access/login', data=payload)
