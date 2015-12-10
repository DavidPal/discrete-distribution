# discrete-distribution

Fast algorithm for sampling from discrete distributions.

Generating a sample takes O(1) time. This is in contrast with the naive
algorithm that takes O(log N) time to generate a sample, where N is the size of
support of the distributions. The naive algorithm is commonly used in many
implementations of C++ standard library (clang, GCC).

The details of the algorithm are described in the `doc/` directory.

The ultimate goal of the project is to make implementation conform to C++
ISO standard and have it accepted to major open source implementations (clang,
GCC).
