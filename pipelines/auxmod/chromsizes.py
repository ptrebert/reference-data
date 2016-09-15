# coding=utf-8

import re as re


def filter_chromosomes(inputfile, outputfiles, selector):
    """
    :param inputfile:
    :param outputfiles:
    :param selector:
    :return:
    """
    if isinstance(selector, str):
        selector = re.compile(selector)
    selected = []
    with open(inputfile, 'r') as inf:
        for line in inf:
            if not line.strip():
                continue
            name, size = line.strip().split()[:2]
            if selector.match(name) is None:
                continue
            selected.append((name, size))
    assert selected, 'No chromosome names extracted from {} with pattern {}'.format(inputfile, selector.pattern)
    selected = sorted(selected, key=lambda x: int(x[1]), reverse=True)
    # first output file is usual 2-column text
    with open(outputfiles[0], 'w') as outf:
        _ = outf.write('\n'.join('\t'.join(t) for t in selected) + '\n')
    # second output file is BED file
    with open(outputfiles[1], 'w') as outf:
        for name, size in selected:
            _ = outf.write('{}\t{}\t{}\n'.format(name, 0, size))
    return outputfiles
