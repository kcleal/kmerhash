========
kmerhash
========

Fast generation of hashed kmers from a python string

Installation
------------

To install, wheels are available for Python >= 3.7 < 3.10 on linux/MacOSX::

    pip install kmerhash

Or build from source::

    git clone https://github.com/kcleal/kmerhash
    cd kmerhash
    pip install .

Overview
--------

kmers are packed into an int64 using a rolling hash. The maximum kmer size is 32 bp and
currently an alphabet of {A, T, C, G} is supported (uppercase only):

.. code-block:: python

    import numpy as np
    from kmerhash import kmerhasher, hashes2seq

    seq = "TTCGGACCGGATT"
    k = 11
    a = kmerhasher(seq, k)  # numpy array returned

    print(a)
    # [4039016 3573155 1709711]


kmer integers can be converted back to string format:

.. code-block:: python

    print(hashes2seq(a, k))
    # ['TTCGGACCGGA', 'TCGGACCGGAT', 'CGGACCGGATT']


kmers can be counted efficiently using e.g. `Counter`:

.. code-block:: python

    from collections import Counter

    counter = Counter(a)



Benchmarks
----------

kmerhash was compared to a pure python implementation with k=21:

.. code-block:: python

    # pure python
    [hash(seq[i:i+k]) for i in range(len(seq) - k + 1)]

    # kmerhash
    kmerhasher(seq, k)

.. list-table::
   :widths: 25 25 50
   :header-rows: 1

   * - Sequence length
     - Throughput python (Mbp / s)
     - Throughput kmerhash (Mbp / s)
   * - 21
     - 9.8
     - 1.1
   * - 100
     - 5.6
     - 9.1
   * - 1000
     - 5.2
     - 50.5
   * - 10,000
     - 4.8
     - 90.79
   * - 100,000
     - 4.9
     - 102.9


Limitations
-----------

- Maximum kmer length is 32
- Only {A, T, C, G} bases are supported
- Uppercase bases only
- N bases will be intepreted as A's


Acknowledgements
----------------
Alexander Kearsey https://github.com/kearseya
