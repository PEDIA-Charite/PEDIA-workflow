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
from datetime import datetime
import boto3
import requests
import json

from lib.singleton import LazyConfigure
from lib import visual


LOGGER = logging.getLogger(__name__)


def load_last_updated(path: str) -> datetime:
    '''load the last update time'''
    last_update = None
    if os.path.exists(path):
        with open(path, "r") as tfile:
            content = json.load(tfile)
            last_update = datetime.strptime(content['updated_at'], '%Y-%m-%d %H:%M:%S')
    return last_update

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
        self.header = {"Accept":"application/json", "Authorization":""}

    def connect(self):
        '''connect to LAB and update header'''
        lab_connect = "https://app.face2gene.com/api/labs/auth/connect"
        response = requests.post(
                lab_connect,
                data=json.dumps(self.connect_key),
                headers={"Content-Type": "application/json", "Accept": "application/json"})

        if response.status_code == 200:
            token = json.loads(response.content.decode('utf-8'))['jwt']
            self.header["Authorization"] = "Bearer " + token
        else:
            err_str = "There is error while connecting to lab %s. Error code: %d" % \
                    (self.lab_id, response.status_code)
            LOGGER.error(err_str)

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

    def connect_case(self, lab_case_id):
        '''connect to LAB and get list of cases on specific page'''
        retry = 3
        content = {}
        while retry > 0:
            lab_get_case = "https://app.face2gene.com/api/lab-cases/%d" % lab_case_id
            response = requests.get(
                    lab_get_case,
                    headers = self.header
                    )
            if response.status_code == 200:
                content = json.loads(response.content.decode('utf-8'))
                break
            elif retry > 0:
                self.connect()
                err_str = "Retry There is error while downloading lab case %d from lab %s. Error code: %d" % \
                        (lab_case_id, self.lab_id, response.status_code)
                LOGGER.debug(err_str)
            else:
                err_str = "There is error while downloading lab cases " + self.lab_id + ". Error code: " + str(response.status_code)
                LOGGER.error(err_str)
                break
            retry -= 1
        return content

    def download_case(
            self,
            download_location,
            lab_case_id,
            count
        ):
        file_path = None
        case_content = self.connect_case(lab_case_id)
        if case_content:
            case_id = case_content["case_data"]["case_id"]
            file_path = os.path.join(download_location, "cases", case_id + ".json")
            last_update = load_last_updated(file_path)
            new_update = datetime.strptime(case_content['updated_at'], '%Y-%m-%d %H:%M:%S')
            if not load_last_updated(file_path):
                out_file = open(file_path, 'w')
                json.dump(case_content, out_file)
                out_file.close()
                status = 'Download case %s' % case_id
                count['download'] += 1
            elif new_update > last_update:
                update_count += 1
                out_file = open(file_path, 'w')
                json.dump(case_content, out_file)
                out_file.close()
                status = 'Update case %s' % case_id
                count['update'] += 1
            else:
                status = 'Skip case %s, file already exists' % case_id

            LOGGER.debug(status)
        else:
            err_str = "There is error while downloading lab case %d from lab %s. Error code: %d" % \
                    (lab_case_id, self.lab_id, response.status_code)
            count['error'] += 1
            LOGGER.error(err_str)
        return file_path, count

    def connect_all_cases(self, page):
        '''connect to LAB and get list of cases on specific page'''
        retry = 3
        content = {}
        self.connect()
        while retry > 0:
            lab_get_all_case = "https://app.face2gene.com/api/labs/v2/%s/lab-cases?page=%d" % \
                    (self.lab_id, page)
            response = requests.get(
                    lab_get_all_case,
                    headers = self.header
                    )
            if response.status_code == 200:
                content = json.loads(response.content.decode('utf-8'))
            elif retry > 0:
                self.connect()
                err_str = "Retry There is error while downloading lab cases on page %d from lab %s. Error code: %d" % \
                        (page, self.lab_id, response.status_code)
                LOGGER.debug(err_str)
            else:
                err_str = "There is error while downloading lab cases on page %d from lab %s. Error code: %d" % \
                        (page, self.lab_id, response.status_code)
                LOGGER.error(err_str)
                break
            retry -= 1
        return content

    def download_all_lab_case(
            self,
            download_location: str = '',
            lab_case_id: str = ''
        ):
        '''Save entire AWS bucket defined by AWS_BUCKET_NAME to download
        location.  '''
        os.makedirs(download_location, exist_ok=True)

        page = 1
        last_page = 1
        current_num = 0
        file_path = []
        count = {'download': 0, 'update': 0, 'error': 0}
        # while not the last page, download next page
        while not (page > last_page):
            content = self.connect_all_cases(page)
            if content:
                last_page = content['last_page']
                total = content['total']
                case_list = content['labCasesList']
                case_in_page = len(case_list)
                for i, case_content in enumerate(case_list):
                    lab_case_id = case_content['lab_case_id']
                    case_path, count = self.download_case(
                        download_location,
                        lab_case_id,
                        count
                        )
                    if case_path:
                        file_path.append(case_path)
                    current_num += 1
                    visual.print_status("Checking LAB", 20, current_num, total)
                page += 1
            else:
                err_str = "There is error while downloading lab cases " + self.lab_id + ". Error code: " + str(response.status_code)
                LOGGER.error(err_str)
        print("")
        LOGGER.info(
            "LAB Stats: Downloaded %d new %d updated of %d total files.",
            count['download'], count['update'], current_num
        )
        return file_path
