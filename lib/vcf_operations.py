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
        fixed_data = first_line + fix_header(vcf_file)
    return fixed_data


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

def fix_header(vcf_file):
    header = \
        '##contig=<ID=1,assembly=b37,length=249250621>\n' + \
        '##contig=<ID=2,assembly=b37,length=243199373>\n' + \
        '##contig=<ID=3,assembly=b37,length=198022430>\n' + \
        '##contig=<ID=4,assembly=b37,length=191154276>\n' + \
        '##contig=<ID=5,assembly=b37,length=180915260>\n' + \
        '##contig=<ID=6,assembly=b37,length=171115067>\n' + \
        '##contig=<ID=7,assembly=b37,length=159138663>\n' + \
        '##contig=<ID=8,assembly=b37,length=146364022>\n' + \
        '##contig=<ID=9,assembly=b37,length=141213431>\n' + \
        '##contig=<ID=10,assembly=b37,length=135534747>\n' + \
        '##contig=<ID=11,assembly=b37,length=135006516>\n' + \
        '##contig=<ID=12,assembly=b37,length=133851895>\n' + \
        '##contig=<ID=13,assembly=b37,length=115169878>\n' + \
        '##contig=<ID=14,assembly=b37,length=107349540>\n' + \
        '##contig=<ID=15,assembly=b37,length=102531392>\n' + \
        '##contig=<ID=16,assembly=b37,length=90354753>\n' + \
        '##contig=<ID=17,assembly=b37,length=81195210>\n' + \
        '##contig=<ID=18,assembly=b37,length=78077248>\n' + \
        '##contig=<ID=19,assembly=b37,length=59128983>\n' + \
        '##contig=<ID=20,assembly=b37,length=63025520>\n' + \
        '##contig=<ID=21,assembly=b37,length=48129895>\n' + \
        '##contig=<ID=22,assembly=b37,length=51304566>\n' + \
        '##contig=<ID=X,assembly=b37,length=155270560>\n' + \
        '##contig=<ID=Y,assembly=b37,length=59373566>'
    out = bytes('', 'utf-8')
    stop = False
    line = vcf_file.readline()
    while line:
        line = line.decode('utf-8')
        if '##contig=<ID=' not in line:
            if '#CHROM' in line:
                line = header + '\n' + line
                stop = True
            out += bytes(line, 'utf-8')
        if stop == True:
            break
        line = vcf_file.readline()
    fixed_content = out + vcf_file.read()
    return fixed_content

def move_vcf(orig_path: str, new_path: str) -> None:
    '''Convert vcf file to vcf.gz and move to the new directory.'''
    move_flag = True
    data = read_vcf(orig_path)
    if os.path.exists(new_path):
        new_data = read_vcf(new_path)
        move_flag = False if new_data == data else True
    if move_flag:
        print("\nCopy VCF file to data/PEDIA/vcfs/original")
        write_vcf(data, new_path)
    else:
        print("\nVCF file is existed and identical.")
