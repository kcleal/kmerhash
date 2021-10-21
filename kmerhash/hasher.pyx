#cython: language_level=3, boundscheck=False, c_string_encoding=utf8, infer_types=True, wraparound=False
import numpy as np
cimport numpy as np
from libc.stdint cimport uint64_t, uint8_t
import array

__all__ = ["kmerhasher", "char_to_nibble_array", "hash2seq", "same_as_str_split",
           "seq_2_nibbles", "hashes2seq"]


basemap = np.array(['.', 'A', 'C', '.', 'G', '.', '.', '.', 'T', '.', '.', '.', '.', '.', 'N'])

cdef int[85] basemap_2_int = [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 13,  0,  0,  0,  0,
        0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  2,
        0,  0,  0,  4,  0,  0,  0,  0,  0,  0, 14,  0,  0,  0,  0,  0,  8]



# in 2bit format:
# A = 0 (bits = 00), C = 1 (01), G = 2 (10), T = 3 (11), N = 0 (00)
# if input is ATCG
# (00) (11) (01) (10) = 54

# in nibble format
# A = 0001, C = 0010, G = 0100, T = 1000
# ATCG = first byte = (0001) (1000)  second byte = (0010) (0100) = [24, 36]


def nibble_2bit():
    return np.array([0, 0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0], dtype="uint8")


def twobit_2_base():
    return {0: "A", 1: "C", 2: "G", 3: "T"}


def char_to_nibble_array(t):
    # encode as nibble
    test_nibble = array.array("B", [])

    t_itr = iter(t)
    for base1 in t_itr:
        v = 0
        v = v << 4 | base1
        try:
            base2 = next(t_itr)
        except StopIteration:
            pass
        else:
            v = v << 4 | base2
        test_nibble.append(v)

    return np.array(test_nibble, dtype="uint8")


def seq_2_nibbles(seq):
    input_int_nibbles = [basemap_2_int[i] for i in seq]
    return char_to_nibble_array(input_int_nibbles)


cpdef hash2seq(uint64_t v, int k):
    twobit_2_base_d = ["A", "C", "G", "T"]
    seq = ""
    cdef int i
    cdef bint is_odd = k % 2 != 0
    for i in range(k):
        base = twobit_2_base_d[v & 3]
        seq += base
        v = v >> 2
    return seq[::-1]


def same_as_str_split(extracted, seq, k):
    kmers = [seq[i:i + k] for i in range(len(seq) - k + 1)]
    rolling_kmers = [hash2seq(i, k) for i in extracted]
    cdef int index = 0
    for i, j in zip(kmers, rolling_kmers):
        if i != j:
            raise ValueError("i != j: {}, {}, occurred at index {}".format(i, j, index))
        index += 1
    return True


cpdef hashes2seq(np.ndarray[np.uint64_t, ndim=1] t, int k):
    ret = []
    cdef uint64_t v
    for v in t:
        ret.append(hash2seq(v, k))
    return ret


cdef uint8_t[16] nib = [0, 0, 1, 0, 2, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0]


cpdef np.ndarray[np.uint64_t, ndim=1] kmerhasher(str seq, int kmer_length):
# def kmerhasher(seq, kmer_length):

    cdef int seq_len = len(seq)
    if seq_len < kmer_length:
        raise ValueError("Input sequence must be >= kmer length")
    if kmer_length > 32:
        raise ValueError("Max kmer length is 32")

    cdef bytes seq_b = seq.encode("ascii")
    cdef uint8_t[::1] t  # t is a nibble array
    cdef int i = 0
    cdef int block = 0

    if seq_len % 2 == 0:
        t = np.zeros(int(seq_len / 2), dtype="uint8")
    else:
        t = np.zeros(int(seq_len / 2) + 1, dtype="uint8")

    cdef uint8_t base1, base2, v

    while i < seq_len:

        base1 = basemap_2_int[seq_b[i]]
        v = 0
        v = v << 4 | base1
        i += 1
        if i == seq_len:
            t[block] = v << 4
            break

        base2 = basemap_2_int[seq_b[i]]
        i += 1
        v = v << 4 | base2
        t[block] = v
        block += 1

    cdef np.ndarray[np.uint64_t, ndim=1] a

    if seq_len % 2 == 0:
        a = np.zeros(len(t) * 2 - kmer_length + 1, dtype="uint64")
    else:
        a = np.zeros(len(t) * 2 - kmer_length, dtype="uint64")

    # cdef np.ndarray[np.uint64_t, ndim=1] a = np.zeros(len(t) * 2 - kmer_length + 1, dtype="uint64")
    cdef int bases_remaining = seq_len
    cdef uint64_t h = 0
    # https://stackoverflow.com/questions/15816927/bit-manipulation-clearing-range-of-bits
    cdef uint64_t mask = ~(h & 0) >> (64 - (kmer_length * 2))
    cdef uint8_t first_2bit, second_2bit
    cdef int index = 1
    cdef int end = <int>(kmer_length / 2)
    cdef int len_t
    bb = twobit_2_base()
    i = 0
    while i < end:  # fill the first block
        v = t[i]

        second_2bit = nib[v & 15]
        first_2bit = nib[(v >> 4) & 15]

        h = h << 2 | first_2bit
        h = h << 2 | second_2bit

        i += 1
        bases_remaining -= 2

    cdef bint is_odd = kmer_length % 2 != 0
    # if odd length kmer push last nibble
    if is_odd:
        first_2bit = nib[(t[i] >> 4) & 15]
        h = h << 2 | first_2bit


    a[0] = h

    if bases_remaining == 0:
        return a

    len_t = len(t)
    while i < len_t:

        v = t[i]

        # these are encoded as nibbles
        # second = v & 15
        # first = (v >> 4) & 15
        # this converts from nibble to 2bit representation

        first_2bit = nib[(v >> 4) & 15]
        second_2bit = nib[v & 15]

        if bases_remaining == 1:
            if is_odd:
                h = (h << 2 | second_2bit) & mask
            else:
                h = (h << 2 | first_2bit) & mask
            a[index] = h
            break

        # break this into three steps:
        # h << 2 makes space for the next 2bit
        # h | first_2bit adds the base into the hash
        # h & mask  drops any bits outside of the kmer length, prevents overflow of the int

        if is_odd:
            if i == end:  # bases remaining > 0
                # just add one bit first bit was already shifted
                h = (h << 2 | second_2bit) & mask
                a[index] = h
                index += 1
                bases_remaining -= 1
                i += 1

            else:
                h = (h << 2 | first_2bit) & mask
                a[index] = h
                index += 1

                h = (h << 2 | second_2bit) & mask
                a[index] = h
                index += 1
                bases_remaining -=2
                i += 1

        else:
            h = (h << 2 | first_2bit) & mask
            a[index] = h
            index += 1

            h = (h << 2 | second_2bit) & mask
            a[index] = h
            index += 1
            bases_remaining -=2
            i += 1

    return a
