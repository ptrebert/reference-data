#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import io as io
import traceback as trb
import argparse as argp
import gzip as gz
import operator as op
import functools as fnt


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--target', '-t', type=str, dest='target')
    parser.add_argument('--query', '-q', type=str, dest='query')
    parser.add_argument('--output', '-o', type=str, dest='output')
    parser.add_argument('--switch', '-s', action='store_true', default=False, dest='switch',
                        help='Switch target and query in the output')
    parser.add_argument('--filter', '-f', type=int, dest='filter', default=0,
                        help='Skip blocks smaller than this size. Default: 0')
    args = parser.parse_args()
    return args


def join_parts(switch, *args):
    """
    :param switch:
    :param args: tc, ts, te, tstr, bc, qc, qs, qe, qstr
    :return:
    """
    # had an annoying bug here - passed "(switch,)" instead of "switch"
    # which always evaluated to True; but did not affect the one file
    # where I used switch, so maybe introduced the error later...?
    # anyway, just to be sure here, some manual type checking...
    assert isinstance(switch, bool), 'Received wrong type for switch: {}'.format(switch)
    if switch:
        items = op.itemgetter(*(5, 6, 7, 3,  # the new target / original query region
                                4,  # block ID
                                0, 1, 2, 8))  # the new query / original target region
    else:
        items = op.itemgetter(*(0, 1, 2, 3,  # the target region
                                4,  # block ID
                                5, 6, 7, 8))  # the query region
    joined = '\t'.join(items(args))
    return joined


def main():
    """
    :return:
    """
    args = parse_command_line()
    outbuffer = io.StringIO()
    bufsize = 0
    block_count = 0
    block_ids = set()
    build_block = fnt.partial(join_parts, *(args.switch, ))
    with open(args.target, 'r') as trgfile:
        with open(args.query, 'r') as qryfile:
            while 1:
                tb = trgfile.readline().strip()
                qb = qryfile.readline().strip()
                try:
                    tc, ts, te, tid = tb.split()
                    qc, qs, qe, qid = qb.split()
                    assert tid == qid,\
                        'Block mismatch for files {} and {}\nLines {} and {}'.format(os.path.basename(args.target),
                                                                                     os.path.basename(args.query),
                                                                                     tb, qb)
                    assert tid not in block_ids,\
                        'Block ID duplicate in files {} and {}\nLines {} and {}'.format(os.path.basename(args.target),
                                                                                        os.path.basename(args.query),
                                                                                        tb, qb)
                    tl = int(te) - int(ts)
                    ql = int(qe) - int(qs)
                    assert tl == ql,\
                        'Coverage mismatch for files {} and {}\nLines {} and {}'.format(os.path.basename(args.target),
                                                                                        os.path.basename(args.query),
                                                                                        tb, qb)
                    if tl < args.filter:
                        continue
                    block_count += 1
                    qstrand = qid.split('_')[-1]
                    #blockline = '\t'.join([tc, ts, te, '+', str(block_count),
                    #                       qc, qs, qe, qstrand])
                    blockline = build_block(tc, ts, te, '+',
                                            str(block_count),
                                            qc, qs, qe, qstrand)
                    bufsize += len(blockline)
                    outbuffer.write(blockline + '\n')
                    if bufsize > 100000:
                        with gz.open(args.output, 'at') as outfile:
                            _ = outfile.write(outbuffer.getvalue())
                            outfile.flush()
                        outbuffer = io.StringIO()
                        bufsize = 0
                except ValueError:
                    break
    with gz.open(args.output, 'at') as outfile:
        _ = outfile.write(outbuffer.getvalue())
        # head a corrupted gzip once - not sure about the cause... I/O interrupted???
        outfile.flush()
    return


if __name__ == '__main__':
    try:
        main()
    except Exception as err:
        trb.print_exc(file=sys.stderr)
        sys.stderr.write('\nError: {}\n'.format(str(err)))
        sys.exit(1)
    else:
        sys.exit(0)
