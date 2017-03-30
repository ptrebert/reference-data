#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import io as io
import traceback as trb
import argparse as argp
import gzip as gz


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--target', '-t', type=str, dest='target')
    parser.add_argument('--query', '-q', type=str, dest='query')
    parser.add_argument('--output', '-o', type=str, dest='output')
    args = parser.parse_args()
    return args


def main():
    """
    :return:
    """
    args = parse_command_line()
    outbuffer = io.StringIO()
    bufsize = 0
    block_count = 0
    block_ids = set()
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
                    block_count += 1
                    qstrand = qid.split('_')[-1]
                    blockline = '\t'.join([tc, ts, te, '+', str(block_count),
                                           qc, qs, qe, qstrand])
                    bufsize += len(blockline)
                    outbuffer.write(blockline + '\n')
                    if bufsize > 100000:
                        with gz.open(args.output, 'at') as outfile:
                            _ = outfile.write(outbuffer.getvalue())
                        outbuffer = io.StringIO()
                        bufsize = 0
                except ValueError:
                    break
    with gz.open(args.output, 'at') as outfile:
        _ = outfile.write(outbuffer.getvalue())

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
