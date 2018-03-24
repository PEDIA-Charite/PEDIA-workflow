'''
VCF operations to move all vcfs to identical vcf.gz format.
'''

import os
import gzip
import zipfile
from typing import Callable


def compress_gz(instr: bytes, outfile: str) -> None:
    '''Gzip binary data to outfile.'''
    os.makedirs(os.path.split(outfile)[0], exist_ok=True)
    with gzip.open(outfile, "wb") as gzipped:
        gzipped.write(instr)


def handle_vcf(openfunc: Callable, infile: str, outfile: str) -> None:
    '''Handle uncompressed vcf files, by validating basic vcf properties
    and gzipping them.'''
    with openfunc(infile) as vcf_file:
        first_byte = vcf_file.readline()
        first = first_byte.decode("utf-8")
        if "VCF" not in first:
            raise TypeError("Uncompressed text file is not vcf format.")
        # reset pointer for gzipping to destination
        rawdata = first_byte + vcf_file.read()

        # create output folder
    compress_gz(rawdata, outfile)


def handle_uncompressed(infile: str, outfile: str) -> None:
    '''Handle uncompressed files. Such as with ending vcf.  These do not
    have a matchable mimetype and are mapped to the text metaclass.
    '''
    def fopen(filepath):
        return open(filepath, "rb")
    handle_vcf(fopen, infile, outfile)


def handle_zip(infile: str, outfile: str) -> None:
    '''Handle zipped files by unzipping and checking vcf.'''
    with zipfile.ZipFile(infile, "r") as inzip:
        filename = inzip.namelist()
        assert len(filename) == 1

        def fopen(filepath):
            return inzip.open(filepath)

        handle_vcf(fopen, filename[0], outfile)


def handle_gzip(infile: str, outfile: str) -> None:
    '''Handle gzipped files by directly moving them to the destination.'''
    def fopen(filepath):
        return gzip.open(filepath, "rb")
    return handle_vcf(fopen, infile, outfile)


MIMETYPES = {
    'application/zip': handle_zip,
    'application/gzip': handle_gzip,
    'text': handle_uncompressed
}


def move_vcf(orig_path: str, new_path: str, mimetype: str) -> None:
    '''Convert vcf file to vcf.gz and move to the new directory.'''
    if mimetype not in MIMETYPES:
        raise TypeError("Not supported mime {}".format(mimetype))
    MIMETYPES[mimetype](orig_path, new_path)
