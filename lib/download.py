#!/usr/bin/env python
'''
AWS Download Function
---
Get all files in specified AWS_BUCKET and download them to the directory
specified as the first argument or as a parameter to download_location.

Downloading from the AWS Bucket requires public and secret AWS Keys. They can
be saved to environment variables AWS_ACCESS_KEY_ID and AWS_ACCESS_SECRET_KEY,
when calling the function as a program or directly to the backup_s3_folder
function.
'''

import os
import logging
from argparse import ArgumentParser

import boto3

from lib.constants import AWS_BUCKET_NAME


def backup_s3_folder(
        aws_access_key: str='',
        aws_secret_key: str='',
        download_location: str='',
        config: 'ConfigParser'=None):
    '''Save entire AWS bucket defined by AWS_BUCKET_NAME to download location.
    '''
    if config:
        aws_access_key = config.aws['access_key']
        aws_secret_key = config.aws['secret_key']
        download_location = config.aws['download_location']
    # use Amazon S3 resource
    s3_resource = boto3.resource('s3',
                                 aws_access_key_id=aws_access_key,
                                 aws_secret_access_key=aws_secret_key)

    if not os.path.exists(download_location):
        logging.info("Making download directory")
        os.makedirs(download_location)

    bucket = s3_resource.Bucket(AWS_BUCKET_NAME)

    downloaded = 0
    # iterate over all items and check whether they already exist on the drive
    # take note, that this does not check for the integrity, so corrupted
    # files will need to be removed manually
    for i, key in enumerate(bucket.objects.all()):
        # strip leading slashes for path joining
        path, filename = os.path.split(key.key)
        filename = filename.strip('/')
        path = path.strip('/')
        print(filename)
        dlpath = os.path.join(download_location, path, filename)
        all_files = i
        if os.path.exists(dlpath):
            logging.info("File {} already exists. Skipping...".format(dlpath))
            continue
        # download file to specified location this process might be unreliable
        # some exceptions might need to be caught for the process to run
        # reliably
        if not os.path.exists(os.path.join(download_location, path)):
            os.makedirs(os.path.join(download_location, path))
        logging.info("Downloading {}".format(filename))
        bucket.download_file(key.key, dlpath)
        downloaded += 1

    logging.info(
        "S3 Stats: Downloaded {} of {} files.".format(downloaded, all_files))


def main():
    '''The main should only be used for testing purposes.
    '''
    logging.basicConfig(filename='aws_download.log', level=logging.INFO)
    parser = ArgumentParser()
    parser.add_argument('dl_dir')
    download_location = parser.parse_args().dl_dir
    aws_access_key = os.getenv("AWS_ACCESS_KEY_ID")
    aws_secret_key = os.getenv("AWS_ACCESS_SECRET_KEY")
    if aws_secret_key is None or aws_access_key is None:
        logging.error(
            'Please set env vars: AWS_ACCESS_KEY_ID and AWS_ACCESS_SECRET_KEY')
    backup_s3_folder(aws_access_key, aws_secret_key, download_location)


if __name__ == '__main__':
    main()
