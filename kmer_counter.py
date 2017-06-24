#!/usr/bin/env python

# Main idea is from https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-12-333

from bitarray import bitarray

# According to this post, best hash is murmur3 https://stackoverflow.com/questions/11954086/which-hash-functions-to-use-in-a-bloom-filter
import mmh3

# Getting the arguments from command line
import sys

# Basic bloomfilter implementation using murmur and bitarray
class BloomFilter:
    def __init__(self, size, hash_count):
        self.bit_array = bitarray(size)
        self.bit_array.setall(0)
        self.size = size
        self.hash_count = hash_count

    def add(self, item):
        for ii in xrange(self.hash_count):
            index = mmh3.hash(item, ii) % self.size
            self.bit_array[index] = 1
        return self

    def lookup(self, item):
        out = True
        for ii in xrange(self.hash_count):
            index = mmh3.hash(item, ii) % self.size
            if self.bit_array[index] == 0:
                out = False
        return out

    def printer(self):
        print self.bit_array

# One million bits and 7 different hashes are enough ?!
bloom_filter = BloomFilter(10000000, 7)
hashmap = {}

# Part1 of the algorithm, which is checking firstly bloom filter then hashmap
def initializer(seq):
    global bloom_filter
    global hashmap

    if bloom_filter.lookup(seq) == False:
        bloom_filter.add(seq)
    else:
        if (seq in hashmap) == False:
            hashmap[seq] = 0

if __name__ == '__main__':
    filename = sys.argv[1]
    kmersize = int(sys.argv[2])
    count = int(sys.argv[3])
    final_seq_list = []

    kmer_start = 0
    kmer_end = kmersize

    # First iteration
    with open(filename) as f:
        line_count = 0
        for line in f:
            if line_count == 1:
                seq = line[0 : -1]
                # We are only interested in this line which is actual seq
                for i in xrange(len(seq) - kmersize + 1):
                    current_seq = seq[kmer_start + i: kmer_end + i]
                    initializer(current_seq)
            elif line_count == 3:
                line_count = -1
            line_count += 1

    # Second iteration
    with open(filename) as f:
        line_count = 0
        for line in f:
            if line_count == 1:
                seq = line[0 : -1]
                # We are only interested in this line which is actual seq
                for i in xrange(len(seq) - kmersize + 1):
                    current_seq = seq[kmer_start + i: kmer_end + i]
                    # Part2 of the algorithm
                    if (current_seq in hashmap) == True:
                        hashmap[current_seq] += 1
            elif line_count == 3:
                line_count = -1
            line_count += 1

    # Part3 of the algorithm which is deleting the unique aka one time only items
    for key, value in hashmap.items():
        if value == 1:
            del hashmap[key]

    # Sort the dict using values and get most n frequent seqs
    seq_counter = 0
    for w in sorted(hashmap, key = hashmap.get, reverse = True):
        if seq_counter == count:
            break
        final_seq_list.append((w, hashmap[w]))
        seq_counter += 1

    # Printing the final sequence list with their freqs.
    print final_seq_list
