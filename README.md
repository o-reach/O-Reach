O'Reach
=====
O'Reach is a fast, deterministic algorithm for answering reachability queries in a directed graph, i.e.,
it spends some time to preprocess a given directed (acyclic) graph and then quickly answers for arbitrary pairs of vertices s and t whether the graph has a directed path from s to t.

O'Reach can either be run standalone or use another algorithm as fallback in
case that it cannot answer a reachability query in constant time.

This repository contains the source code accompanying our paper
[O'Reach: Even Faster Reachability in Large Graphs](https://arxiv.org/abs/2008.10932).
Further information and instances are available on our [paper website](https://oreach.taa.univie.ac.at/).

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
The program is licensed under MIT licence.
If you publish results using our algorithms, please acknowledge our work by citing the following paper:

```
@article{hst-oreachfull-2021,
  author        = {Hanauer, Kathrin and Schulz, Christian and Trummer, Jonathan},
  title         = {O'Reach: Even Faster Reachability in Large Graphs},
  journal       = {CoRR},
  volume        = {abs/2008.10932},
  year          = {2021},
  url           = {https://arxiv.org/abs/2008.10932},
  archivePrefix = {arXiv},
  eprint        = {2008.10932},
  primaryClass  = {cs.DS},
}
```

Project Contributors (sorted by last name)
=====
Kathrin Hanauer

Christian Schulz

Jonathan Trummer
