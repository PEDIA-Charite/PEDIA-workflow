import re
import hashlib
import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry

'''
Use the Face2Gene library to get syndrome information.

'''

class Face2Gene(requests.Session):
    base_url = 'https://app.face2gene.com/'

    def __init__(self, user='', password='', config=None):
        super().__init__()
        if config:
            user = config.face2gene['user']
            password = config.face2gene['password']
        retries = Retry(total=5, backoff_factor=1, status_forcelist=[ 429, 500 ])
        self.mount('http://', HTTPAdapter(max_retries=retries))
        self.headers.update({'User-Agent' : 'pediaScript'})
        self._login(user, password)

        # cache previous searches by keywords used
        self._cache = {}

    def goto_library(self):
        url = self.base_url + 'library/search'
        r = self.get(url)
        return r

    def search_library(self, string):
        if string in self._cache:
            return self._cache[string]

        url = self.base_url + 'library/search_api'
        params = {'keywords' : string}
        r = self.get(url, params=params).json()
        self._cache[string] = r
        return r

    def search_syndrome(self, string, omim_list=[], return_first=True):
        '''Search for syndromes on Face2Gene library.
        '''
        r = self.search_library(string)
        syndromes = r['syndromes']['data']
        # return the first entry with a omim id
        for s in syndromes:
            if not s['omim_id']:
                continue
            if str(s['omim_id']) != '0':
                if omim_list:
                    if int(s['omim_id']) in omim_list:
                        return s['omim_id']
                else:
                    return s['omim_id']
        return ''

    def _login(self, user, password):
        payload = { "email" : user,
                "password" : hashlib.md5(password.encode('utf-8')).hexdigest() }
        self.post('https://app.face2gene.com/access/login', data=payload)

