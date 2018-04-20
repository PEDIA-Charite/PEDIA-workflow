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

import time
import datetime
import boto3

from lib.constants import AWS_BUCKET_NAME
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

    os.makedirs(download_location, exist_ok=True)

    bucket = s3_resource.Bucket(AWS_BUCKET_NAME)

    # time information
    time_path = os.path.join(download_location, ".LAST_UPDATE")
    last_update = load_timestamp(time_path)

    downloaded = 0
    updated = 0
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
        edit_time = key.last_modified
        if os.path.exists(dlpath) and edit_time <= last_update:
            LOGGER.debug("File %s already exists. Skipping...", dlpath)
            continue
        elif edit_time > last_update:
            LOGGER.debug(
                "File %s change time %s will be updated.",
                dlpath, str(edit_time)
            )
            updated += 1
        else:
            downloaded += 1
        # download file to specified location this process might be unreliable
        # some exceptions might need to be caught for the process to run
        # reliably
        os.makedirs(os.path.join(download_location, path), exist_ok=True)
        LOGGER.debug("Downloading %s", filename)
        bucket.download_file(key.key, dlpath)
    # stop progressbar
    print("")

    save_current_time(time_path)

    LOGGER.info("S3 Stats: Downloaded %d new %d updated of %d total files.",
                downloaded, updated, num_objs)
