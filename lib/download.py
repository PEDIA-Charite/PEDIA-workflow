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

import boto3

from lib.constants import AWS_BUCKET_NAME
from lib import visual


LOGGER = logging.getLogger(__name__)


def backup_s3_folder(
        aws_access_key: str = '',
        aws_secret_key: str = '',
        download_location: str = '',
        config: 'ConfigParser' = None):
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
    LOGGER.info("Downloading files from AWS Bucket %s", AWS_BUCKET_NAME)
    bucket_objs = list(bucket.objects.all())
    num_objs = len(bucket_objs)
    for i, key in enumerate(bucket_objs):
        visual.print_status("Checking AWS", 20, i+1, num_objs)
        # strip leading slashes for path joining
        path, filename = os.path.split(key.key)
        # skip directories
        if not filename:
            continue
        filename = filename.strip('/')
        path = path.strip('/')
        dlpath = os.path.join(download_location, path, filename)
        all_files = i
        if os.path.exists(dlpath):
            LOGGER.debug("File %s already exists. Skipping...", dlpath)
            continue
        # download file to specified location this process might be unreliable
        # some exceptions might need to be caught for the process to run
        # reliably
        if not os.path.exists(os.path.join(download_location, path)):
            os.makedirs(os.path.join(download_location, path))
        LOGGER.debug("Downloading %s", filename)
        bucket.download_file(key.key, dlpath)
        downloaded += 1
    # stop progressbar
    print("")

    LOGGER.info("S3 Stats: Downloaded %d new of %d total files.",
                downloaded, all_files)
