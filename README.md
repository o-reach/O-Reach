[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/o-reach/O-Reach/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/o-reach/O-Reach)
[![SWH](https://archive.softwareheritage.org/badge/swh:1:dir:55d23a5b940f1ead285729c8dbd82c71e28d504a/)](https://archive.softwareheritage.org/swh:1:dir:55d23a5b940f1ead285729c8dbd82c71e28d504a)

O'Reach
=====
O'Reach is a fast, deterministic algorithm for answering reachability queries in a directed graph, i.e.,
it spends some time to preprocess a given directed (acyclic) graph and then quickly answers for arbitrary pairs of vertices s and t whether the graph has a directed path from s to t.

O'Reach can either be run standalone or use another algorithm as fallback in
case that it cannot answer a reachability query in constant time.

This repository contains the source code accompanying our paper
[O'Reach: Even Faster Reachability in Large Graphs](https://doi.org/10.4230/LIPIcs.SEA.2021.13).
Further information as well as instances are available on our [paper website](https://oreach.taa.univie.ac.at/).

Installation Notes
=====

Before you can start you need to install the following software packages:

- [CMake](https://cmake.org/) (at least version 3.10)
- [Boost](https://www.boost.org/) (at least version 1.61)

Both are also available as packages on many Linux distributions, e.g.:
- Debian/Ubuntu: ``apt install cmake libboost-dev``
- Fedora: ``dnf install cmake boost``


Afterwards, in the main project directory, type `./compile_withcmake.sh`. Once
you did that you can try to run the following command:
```
./build/reachability ./examples/arXiv.metis
```

External Projects
======
We use code from the following projects (shipped in extern/):

- [KaHIP v3.10](https://github.com/KaHIP/KaHIP/)
- [argtable3](https://github.com/argtable/argtable3)
- [MurmurHash3](https://github.com/aappleby/smhasher/blob/92cf3702fcfaadc84eb7bef59825a23e0cd84f56/src/MurmurHash3.cpp)
- [ska flat hashmap](https://github.com/skarupke/flat_hash_map/blob/2c4687431f978f02a3780e24b8b701d22aa32d9c/flat_hash_map.hpp)


License
=====
Our code is licensed under MIT license.
If you publish results using our algorithms, please acknowledge our work by
citing the [corresponding paper](https://doi.org/10.4230/LIPIcs.SEA.2021.13):

```
@inproceedings{hst-oreach-2021,
  author    = {Hanauer, Kathrin and Schulz, Christian and Trummer, Jonathan},
  editor    = {Coudert, David and Natale, Emanuele},
  title     = {O'Reach: Even Faster Reachability in Large Graphs},
  booktitle = {19th International Symposium on Experimental Algorithms, {SEA} 2021,
               June 7-9, 2021, Valrose, France},
  series    = {LIPIcs},
  volume    = {190},
  pages     = {13:1--13:24},
  publisher = {Schloss Dagstuhl - Leibniz-Zentrum f{\"{u}}r Informatik},
  year      = {2021},
  url       = {https://doi.org/10.4230/LIPIcs.SEA.2021.13},
  doi       = {10.4230/LIPIcs.SEA.2021.13},
}
```

Project Contributors (sorted by last name)
=====
- Kathrin Hanauer
- Christian Schulz
- Jonathan Trummer
