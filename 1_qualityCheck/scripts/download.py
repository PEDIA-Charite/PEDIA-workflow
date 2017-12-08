#!/usr/bin/env python

import boto3
import sys, os
from argparse import ArgumentParser
import logging

AWS_BUCKET_NAME = "fdna-pedia-dump"

def backup_s3_folder(aws_access_key, aws_secret_key, download_location):
    # use Amazon S3 resource
    s3 = boto3.resource('s3', aws_access_key_id=aws_access_key,
    aws_secret_access_key=aws_secret_key)

    if not os.path.exists(download_location):
        logging.info("Making download directory")
        os.mkdir(download_location)

    bucket = s3.Bucket(AWS_BUCKET_NAME)

    for key in bucket.objects.all():
        # strip leading slashes for path joining
        path, filename = os.path.split(key.key)
        filename = filename.strip('/')
        path = path.strip('/')
        dlpath = os.path.join(download_location, path, filename)
        try:
            logging.info("Downloading ", filename)
            bucket.download_file(key.key, dlpath)
        except (OSError,S3ResponseError) as e:
            logging.error(e)

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('dl_dir')
    download_location = parser.parse_args().dl_dir
    aws_access_key = os.getenv("AWS_ACCESS_KEY_ID")
    aws_secret_key = os.getenv("AWS_ACCESS_SECRET_KEY")
    if aws_secret_key is None or aws_access_key is None:
        logging.error('Please set env vars: AWS_ACCESS_KEY_ID and AWS_ACCESS_SECRET_KEY')
    backup_s3_folder(aws_access_key, aws_secret_key, download_location)
