import requests

class Mutalyzer(requests.Session):

    base_url = 'https://mutalyzer.nl/json/'

    def __init__(self):
        super().__init__()

    def _call_service(self, method, params):
        url = self.base_url + method
        r = self.get(url, params=params)
        r.raise_for_status()
        return r.json()

    def getdbSNPDescriptions(self, rs_id):
        params = {'rs_id' : rs_id}
        return self._call_service('getdbSNPDescriptions', params)

def main():
    m = Mutalyzer()
    j = m.getdbSNPDescriptions("rs386834107")


if __name__ == '__main__':
    main()
