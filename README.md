# Parallelised-Temperature-Simulation
Coursework for masters degree in which OpenMP is used to speed up a temperature simulation algorithm. Contains a report on findings.

# File Structure
The modified code and build is under the "03_coursework1_code" folder. The original is under the "03_coursework1_original" folder.

# System Specifications for Benchmarking
The Warwick supercomputer 'Kudu' has 6 NUMA nodes. At each node there is a single Intel Xeon core with 40 threads and a clock speed of 2.4Ghz. For this project only a single node was used at once and the program speed was increased through splitting the work between threads.

# Modifying Program Variables
There is a single modifyable variable defined in the 'DIVISIONS.h' header. This control how many times the physics mesh is divided along the x and y axis. After editing the program must be recompiled using 'make clean', 'make' commands.

2 divisions
  |
-----
  |

3 divisions
  | |
-------
  | |
-------
  | |

etc...