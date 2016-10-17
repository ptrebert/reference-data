# coding=utf-8

import csv as csv
import gzip as gz


def make_5p_window(inputfile, outputfile, adjust):
    """
    :param inputfile:
    :param outputfile:
    :return:
    """

    rowbuffer = []
    with gz.open(inputfile, 'rt', newline='') as inf:
        rows = csv.DictReader(inf, delimiter='\t')
        fields = rows.fieldnames
        for r in rows:
            if r['strand'] == '+' or r['strand'] == '.':
                r['start'] = str(int(r['start']) - adjust)
                r['end'] = str(int(r['start']) + adjust)
            else:
                r['start'] = str(int(r['end']) - adjust)
                r['end'] = str(int(r['end']) + adjust)
            rowbuffer.append(r)
    with gz.open(outputfile, 'wt') as outf:
        writer = csv.DictWriter(outf, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rowbuffer)
    return outputfile
