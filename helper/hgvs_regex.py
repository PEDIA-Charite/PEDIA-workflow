import hgvs.parser
import re

test_strings = [
        'NM_005221:exon2:c.A533C:p.Q178P:c.',
        'NM_020732.3:c.3323_3324AAdel'
        ]

ACCN_RE = re.compile("([A-Za-z]\w+(_\w+)?(.\d+)?)")

# consume the string step by step
def find_parts(string):
    i = 0
    while True:
        m  = cur_re.search(string)
        if m:
            accn = m.group(1)
            string = string[m.end():]
        else:
            cur_re = hgvs_re.next()
        i += 1
        if not string:
            break


parser = hgvs.parser.Parser()
