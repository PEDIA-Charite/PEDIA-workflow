'''
LAB Download Function
---
Downloading from the F2G LAB API requires public and secret AWS Keys. They can
be saved to environment variables AWS_ACCESS_KEY_ID and AWS_ACCESS_SECRET_KEY
'''

import os
import logging
import sys

import time
import datetime
import boto3
import requests
import json

from lib.singleton import LazyConfigure
from lib import visual


LOGGER = logging.getLogger(__name__)


def load_timestamp(path: str) -> datetime.datetime:
    '''Get timestamp from file.'''
    if os.path.exists(path):
        with open(path, "r") as tfile:
            timestamp = float(tfile.read())
    else:
        timestamp = 0.0
    return datetime.datetime.fromtimestamp(timestamp, datetime.timezone.utc)


def save_current_time(path: str) -> None:
    '''Save current time as UNIX timestamp to the specified file.'''
    with open(path, "w") as tfile:
        tfile.write(str(time.time()))

class Lab(LazyConfigure):

    def __init__(self):
        super().__init__()
        self.access_key = None
        self.secret_key = None
        self.lab_id = None

    def configure(
            self,
            lab_id: str = "",
            key: str = "",
            secret: str = "",
    ):
        super().configure()
        self.access_key = key
        self.secret_key = secret
        self.lab_id = lab_id
        self.connect_key = {"api_key": key, "secret": secret}

    def download_lab_case(
            self,
            download_location: str = '',
            lab_case_id: str = ''
        ):
        '''Save entire AWS bucket defined by AWS_BUCKET_NAME to download
        location.  '''
        os.makedirs(download_location, exist_ok=True)

        lab_connect = "https://app.face2gene.com/api/labs/auth/connect"
        lab_get_case_list = "https://app.face2gene.com/api/lab-cases?lab_id=" + self.lab_id
        lab_get_case = "https://app.face2gene.com/api/lab-cases/" + lab_case_id
        response = requests.post(
                lab_connect,
                data=json.dumps(self.connect_key),
                headers={"Content-Type": "application/json", "Accept": "application/json"})

        token = ''
        if response.status_code == 200:
            token = json.loads(response.content.decode('utf-8'))['jwt']
        else:
            err_str = "There is error while connecting to lab " + self.lab_id + ". Error code: " + str(response.status_code)
            sys.exit(err_str)
        if token != '':
            auth = "Bearer " + token
            response = requests.get(lab_get_case, headers={"Accept":"application/json", "Authorization":auth})
            if response.status_code == 200:
                case_content = json.loads(response.content.decode('utf-8'))
                case_id = case_content["case_data"]["case_id"]
                output = case_content
                file_path = os.path.join(download_location, case_id + ".json")
                out_file = open(file_path, 'w')
                json.dump(output, out_file)
                out_file.close()
            else:
                err_str = "There is error while downloading lab case " + lab_case_id + " from lab " + self.lab_id + ". Error code: " + str(response.status_code)
                sys.exit(err_str)
        return file_path

