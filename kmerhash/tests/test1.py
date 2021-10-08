
import numpy as np
import random
import time
import os
import datetime
import sys

from kmerhash import kmerhasher, same_as_str_split, char_to_nibble_array, hashes2seq

random.seed(0)


basemap = np.array(['.', 'A', 'C', '.', 'G', '.', '.', '.', 'T', '.', '.', '.', '.', '.', 'N'])
basemap_2_int = [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 13,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  2,
        0,  0,  0,  4,  0,  0,  0,  0,  0,  0, 14,  0,  0,  0,  0,  0,  8]


assert basemap_2_int[ord("A")] == 1
assert basemap_2_int[ord("C")] == 2
assert basemap_2_int[ord("G")] == 4
assert basemap_2_int[ord("T")] == 8
assert basemap_2_int[ord("N")] == 14

import resource

basemap_2_int = {k: i for i, k in enumerate(basemap)}

pth = os.path.dirname(os.path.realpath(__file__))

# in 2bit format:
# A = 0 (bits = 00), C = 1 (01), G = 2 (10), T = 3 (11), N = 0 (00)
# if input is ATCG
# (00) (11) (01) (10) = 54

# in nibble format
# A = 0001, C = 0010, G = 0100, T = 1000
# ATCG = first byte = (0001) (1000)  second byte = (0010) (0100) = [24, 36]


assert all(i == j for i, j in zip(char_to_nibble_array([basemap_2_int[i] for i in "ATCG"]), [24, 36]))


def test():
    k = 21
    seq_lengths = np.array([k, 100, 1000, 10000, 100000])
    times = []
    times_normal = []
    for sl in seq_lengths:
        avg = []
        avg_n = []
        inputseq = "".join([random.choice("ATCG") for i in range(sl)])
        for n in range(1):

            t0 = time.time()
            kmer_hashes = [hash(inputseq[i:i+k]) for i in range(len(inputseq) - k + 1)]
            avg_n.append(time.time() - t0)

            t0 = time.time()
            # inputseq = "TCTTGCCGG"
            # k = 5
            # print(inputseq, k, sl)
            # print([inputseq[i:i+k] for i in range(len(inputseq) - k + 1)])
            extracted_hashes = kmerhasher(inputseq, kmer_length=k)
            avg.append(time.time() - t0)
            # print(k, sl)

            same_as_str_split(extracted_hashes, inputseq, k)
            # quit()

        # print("Mean speedup", np.mean(avg_n) / np.mean(avg))
        times.append(np.mean(avg))
        times_normal.append(np.mean(avg_n))

    print("time normal", times_normal)
    print((seq_lengths / times_normal) / 1e6)

    print("time kmerhash", times)
    print((seq_lengths / times) / 1e6)


test()


def demo():

    seq = "TTCGGACCGGATT"
    k = 11
    a = kmerhasher(seq, k)
    print(a)
    print(hashes2seq(a, k))

# demo()