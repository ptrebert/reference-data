# coding=utf-8

import csv as csv
import gzip as gz


def make_5p_window(inputfile, outputfile, upstream, downstream):
    """
    :param inputfile:
    :param outputfile:
    :return:
    """
    rowbuffer = []
    with gz.open(inputfile, 'rt', newline='') as inf:
        rows = csv.DictReader(inf, delimiter='\t')
        fields = rows.fieldnames
        fields.append('5p_end')
        for r in rows:
            if r['strand'] == '+' or r['strand'] == '.':
                tss_5p = int(r['start'])
                r['5p_end'] = tss_5p
                r['start'] = tss_5p - upstream
                r['end'] = tss_5p + downstream
            else:
                tss_5p = int(r['end'])
                r['5p_end'] = tss_5p
                r['start'] = tss_5p - downstream
                r['end'] = tss_5p + upstream
            assert r['end'] - r['start'] == upstream + downstream, 'Length mismatch: {}'.format(r)
            rowbuffer.append(r)
    with gz.open(outputfile, 'wt') as outf:
        writer = csv.DictWriter(outf, fieldnames=fields, delimiter='\t')
        writer.writeheader()
        writer.writerows(rowbuffer)
    return outputfile
