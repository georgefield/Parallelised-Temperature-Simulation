## Parallelised Temperature Simulation
This is coursework for my masters degree in which OpenMP is used to speed up a temperature simulation algorithm called DEQN. This repository contains a 4 page report detailing the findings. The report is named 'CS402_coursework_1.pdf'.

## File Structure
The modified code and build is under the "03_coursework1_code" folder. The original is under the "03_coursework1_original" folder.

## System Specifications for Benchmarking
The Warwick supercomputer 'Kudu' has 6 NUMA nodes. At each node there is a single Intel Xeon core with 40 threads and a clock speed of 2.4Ghz. For this project only a single node was used and the program speed was increased through splitting the work between threads.

## Modifying Program Variables
There is a single modifyable variable defined in the 'DIVISIONS.h' header. This control how many times the physics mesh is divided along the x and y axis. After editing the program must be recompiled using 'make clean', 'make' commands.

2 divisions (4 physics cells)
| Cell 0 | Cell 1 |
|--------|--------|
| Cell 2 | Cell 3 |

3 divisions (9 physics cells)
| Cell 0 | Cell 1 | Cell 2 |
|--------|--------|--------|
| Cell 3 | Cell 4 | Cell 5 |
| Cell 6 | Cell 7 | Cell 8 |

etc...
