#!/usr/bin/env python3
# coding=utf-8

import os as os
import sys as sys
import traceback as trb
import argparse as argp
import collections as col
import multiprocessing as mp
import functools as fnt
import gzip as gz

import pandas as pd


CrossBlock = col.namedtuple('CrossBlock', ['chrom', 'start', 'end', 'block'])

SymmBlock = col.namedtuple('SymmBlock', ['trg_chrom', 'trg_start', 'trg_end', 'trg_strand',
                                         'block',
                                         'qry_chrom', 'qry_start', 'qry_end', 'qry_strand'])


def _parse_crossmap_line(num, explode, line):
    """
    :param line:
    :return:
    """
    chrom, start, end, block = explode(line)
    cb = CrossBlock(chrom, num(start), num(end), num(block))
    return cb


def _parse_symm_line(num, explode, line):
    """
    :param line:
    :return:
    """
    tchrom, tstart, tend, tstrand, block, qchrom, qstart, qend, qstrand = explode(line)
    sb = SymmBlock(tchrom, num(tstart), num(tend), tstrand,
                   num(block),
                   qchrom, num(qstart), num(qend), qstrand)
    return sb


def _compare_blocks(symblock, crossblock):
    """
    :param symblock:
    :param crossblock:
    :return:
    """
    v = symblock.qry_chrom == crossblock.chrom
    v &= symblock.qry_start == crossblock.start
    v &= symblock.qry_end == crossblock.end
    return v


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser()
    parser.add_argument('--symm-map', '-s', type=str, dest='symmmap')
    parser.add_argument('--cross-map', '-c', type=str, dest='crossmap')
    parser.add_argument('--chain-lut', '-lu', type=str, dest='chainlut')
    parser.add_argument('--output', '-o', type=str, dest='output')
    args = parser.parse_args()
    return args


def confirm_mapping(lut, symb, cmb):
    """
    :param lut:
    :param symb:
    :param cmb:
    :return:
    """
    chains = lut.query('qchrom == @cmb.chrom and qstart < @cmb.end and '
                       'qend > @cmb.start and score > 1000 and '
                       'tlength > 100 and qlength > 100', inplace=False)
    assert not chains.empty, 'Blocks {} and {} resulted in empty lookup'.format(symb, cmb)
    if chains.shape[0] == 1:
        assert symb.trg_chrom == chains['tchrom'].item(), 'Chromosome mismatch: {} vs {}'.format(symb, chains)
        assert symb.trg_start < chains['tend'].item() and symb.trg_end > chains['tstart'].item(),\
            'Coordinates incompatible: {} vs {}'.format(symb, chains)
        assert symb.qry_strand == chains['qstrand'].item(), 'Strand mismatch: {} vs {}'.format(symb, chains)
        # only one chain matches - optimal case
        new_symb = SymmBlock(symb.trg_chrom, symb.trg_start, symb.trg_end, symb.trg_strand,
                             symb.block,
                             symb.qry_chrom, cmb.start, cmb.end, symb.qry_strand)
        return new_symb
    elif chains.shape[0] == 2:
        min_val = chains['score'].min()
        chains = chains.query('score == @min_val', inplace=False)
        assert symb.trg_chrom == chains['tchrom'].item(), 'Chromosome mismatch: {} vs {}'.format(symb, chains)
        assert symb.trg_start < chains['tend'].item() and symb.trg_end > chains['tstart'].item(),\
            'Coordinates incompatible: {} vs {}'.format(symb, chains)
        new_symb = SymmBlock(symb.trg_chrom, symb.trg_start, symb.trg_end, symb.trg_strand,
                             symb.block,
                             symb.qry_chrom, cmb.start, cmb.end, chains['qstrand'].item())
        return new_symb
    else:
        raise ValueError('Query blocks\n{}\nand\n{}\nreturned several matches:\n\n{}'.format(symb, cmb, chains))


def lookup_chain(chainlut, recv_lines, send_checked):
    """
    :param chainlut:
    :param recv_lines:
    :param send_checked:
    :return:
    """
    with pd.HDFStore(chainlut, 'r') as hdf:
        lut = hdf['chainlut']
    check_map = fnt.partial(confirm_mapping, *(lut, ))
    while 1:
        if recv_lines.poll():
            item = recv_lines.recv()
            if item is None:
                break
            symb, cmb = item
            try:
                new_symb = check_map(symb, cmb)
                send_checked.send(new_symb)
            except Exception as err:
                send_checked.send(None)
                raise err
    send_checked.send(None)
    return


def read_map_files(symmap, crossmap, send_lines, recv_checked, send_blocks):
    """
    :param symmap:
    :param crossmap:
    :param send_lines:
    :param recv_checked:
    :param send_blocks:
    :return:
    """
    read_symm = fnt.partial(_parse_symm_line, *(int, str.split))
    read_cross = fnt.partial(_parse_crossmap_line, *(int, str.split))
    comp = _compare_blocks
    clean = str.strip
    read_lines = 0
    try:
        with gz.open(symmap, 'rt') as symmfile:
            with open(crossmap, 'r') as crossfile:
                while 1:
                    try:
                        crossblock = read_cross(crossfile.readline())
                    except ValueError:
                        break
                    symblock = read_symm(clean(symmfile.readline()))
                    read_lines += 1
                    assert symblock.block == crossblock.block,\
                        'Block ID mismatch - files not sorted? {} vs {}'.format(symblock, crossblock)
                    if comp(symblock, crossblock):
                        send_blocks.send(symblock)
                    else:
                        send_lines.send((symblock, crossblock))
                    while recv_checked.poll():
                        item = recv_checked.recv()
                        if item is None:
                            raise RuntimeError('Reader received None before completion, read {} lines'.format(read_lines))
                        send_blocks.send(item)
    except Exception as err:
        send_lines.send(None)
        send_blocks.send(None)
        raise err
    else:
        send_lines.send(None)
        while 1:
            if recv_checked.poll():
                item = recv_checked.recv()
                if item is None:
                    break
                send_blocks.send(item)
        return


def dump_buffer(outbuffer):
    """
    :param outbuffer:
    :return:
    """
    conv = str
    for tup in outbuffer:
        _ = sys.stdout.write(tup[0] + '\t' + conv(tup[1]) + '\t' + conv(tup[2]) + '\t' +
                             tup[3] + '\t' + conv(tup[4]) + '\t' + tup[5] + '\t' +
                             conv(tup[6]) + '\t' + conv(tup[7]) + '\t' + tup[8] + '\n')
    return


def main():
    """
    :return:
    """
    mp.set_start_method('forkserver')
    args = parse_command_line()
    outbuffer = []
    children = []
    recv_lines, send_lines = mp.Pipe(duplex=False)
    recv_checked, send_checked = mp.Pipe(duplex=False)
    recv_blocks, send_blocks = mp.Pipe(duplex=False)
    deadcount = 0
    try:

        reader = mp.Process(target=read_map_files, args=(args.symmmap, args.crossmap, send_lines,
                                                         recv_checked, send_blocks),
                            name='read')
        reader.daemon = True
        checker = mp.Process(target=lookup_chain, args=(args.chainlut, recv_lines, send_checked),
                             name='check')
        checker.daemon = True
        children = [reader, checker]
        for p in children:
            p.start()
            p.join(0.01)
        while 1:
            if not (reader.is_alive() and checker.is_alive()):
                deadcount += 1
                if deadcount > 9:
                    raise RuntimeError('I think he is dead, Jim.\n'
                                       'Reader process is alive: {}\n'
                                       'Checker process is alive {}'.format(reader.is_alive(), checker.is_alive()))
            if recv_blocks.poll():
                item = recv_blocks.recv()
                if item is None:
                    break
                outbuffer.append(item)
            if len(outbuffer) > 100000:
                dump_buffer(outbuffer)
                outbuffer = []
    except Exception as err:
        for p in children:
            if p.name == 'check':
                try:
                    send_lines.send(None)
                except:
                    pass
            elif p.name == 'read':
                try:
                    send_checked.send(None)
                except:
                    pass
        raise err
    else:
        dump_buffer(outbuffer)
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
