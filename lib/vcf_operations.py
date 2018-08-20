'''
VCF operations to move all vcfs to identical vcf.gz format.
'''

import os
import gzip
import zipfile
from typing import Callable

import filetype


def read_bytes(openfunc: Callable, infile: str) -> None:
    '''Handle uncompressed vcf files, by validating basic vcf properties
    and gzipping them.'''
    with openfunc(infile) as vcf_file:
        first_line = vcf_file.readline()
        first = first_line.decode("utf-8")
        if "VCF" not in first:
            print(first)
            raise TypeError("Uncompressed text file is not vcf format.")
        # not every archive format supports seeking correctly
        rawdata = first_line + vcf_file.read()
    return rawdata


def read_text(infile: str) -> bytes:
    '''Handle uncompressed files. Such as with ending vcf.  These do not
    have a matchable mimetype and are mapped to the text metaclass.
    '''
    def fopen(filepath):
        return open(filepath, "rb")
    return read_bytes(fopen, infile)


def read_zip(infile: str) -> bytes:
    '''Handle zipped files by unzipping and checking vcf.'''
    with zipfile.ZipFile(infile, "r") as inzip:
        filename = inzip.namelist()
        assert len(filename) == 1

        def fopen(filepath):
            return inzip.open(filepath)

        data = read_bytes(fopen, filename[0])
    return data


def read_gzip(infile: str) -> bytes:
    '''Handle gzipped files by directly moving them to the destination.'''
    def fopen(filepath):
        return gzip.open(filepath, "rb")
    return read_bytes(fopen, infile)


def read_vcf(path: str) -> bytes:
    '''Read vcf to raw bytestring.'''
    # get mimetype
    kind = filetype.guess(path)
    mimetype = kind.mime if kind is not None else "text"
    if mimetype == "application/zip":
        data = read_zip(path)
    elif mimetype == "application/gzip":
        data = read_gzip(path)
    elif mimetype == "text":
        data = read_text(path)
    else:
        raise TypeError("Not supported mime {}".format(mimetype))

    return data


def compress_gz(instr: bytes, outfile: str) -> None:
    '''Gzip binary data to outfile.'''
    os.makedirs(os.path.split(outfile)[0], exist_ok=True)
    with gzip.open(outfile, "wb") as gzipped:
        gzipped.write(instr)


def write_vcf(data: bytes, path: str) -> None:
    '''Write VCF data to output path.'''
    compress_gz(data, path)


def move_vcf(orig_path: str, new_path: str) -> None:
    '''Convert vcf file to vcf.gz and move to the new directory.'''
    data = read_vcf(orig_path)
    write_vcf(data, new_path)
