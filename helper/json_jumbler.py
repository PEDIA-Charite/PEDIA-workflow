'''
Merge json files with same id from two different directories and
save them to another.
'''
import json
import os


def main():
    '''Merge old and new files and save them to another directory.
    '''
    directory_a = '/data/users/Max/PEDIA-workflow/process/aws_dir/cases'
    directory_b = '/data/users/Max/PEDIA-workflow/convert'
    files_a = os.listdir(directory_a)
    files_b = os.listdir(directory_b)

    output_dir = '/data/users/Max/PEDIA-workflow/merged'
    os.makedirs(output_dir, exist_ok=True)

    jsons_a = {v: json.load(open(os.path.join(directory_a, v), 'r'))
               for v in files_a}
    jsons_b = {v: json.load(open(os.path.join(directory_b, v), 'r'))
               for v in files_b}

    for filename, data in jsons_b.items():
        data_a = jsons_a[filename]
        merged = dict(data_a, **data)
        json.dump(merged, open(os.path.join(output_dir, filename), 'w'))


if __name__ == '__main__':
    main()
