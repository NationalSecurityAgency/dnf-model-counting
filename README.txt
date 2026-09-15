DNF Model Counter

A high‑performance polynomial-time DNF model counter. The estimate returned is
probably approximately correct, with arbitrary precision eps and certainty 
delta.


Authors
Paul Burkhardt, Research Directorate, National Security Agency, 
pburkha@nsa.gov
David G. Harris, Department of Computer Science, University of Maryland, 
davidgharris29@gmail.com
Kevin T. Schmitt, Research Directorate, National Security Agency, 
ktschm2@nsa.gov

Dependencies: GNU Multiple Precision Arithmetic Library (GMP)

Installation

1. Install GMP library (if not already present)

2. C++ Compilation

    g++ -lgmp -I include/ -O3 -o count count.cpp rand.cpp counthelp.cpp

    counthelp.cpp defines many subroutines useful for logical formulations, 
    such as Lazy Assignments for a sampled clause and SAT checking.

    rand.cpp handles random number generation for sampling

    count.cpp contains the main algorithm

Basic Example

Inputs to the executable accept a DNF in the below format consisting of a 
header and rows. The header only specifies the number of variables n and 
number of clauses m. The number of variables n may exceed the number of
variables which appear in the DNF, the same is not true for the number of
clauses. Each row after the header represents a clause, which consists of a 
variable amount of literals represented by non-zero integers; the last entry in 
each row is 0 to signal the end of the clause. 

The example below specifies that there are 12 variables and 4 clauses, but only 
5 variables appear in the clauses. The first clause has two literals, 1 and -4, 
meaning it is true if 1 and -4 are  true.

p DNF 12 4
1 -4 0
1 2 3 0
2 -4 10 0
3 10 0

For the DNF above, we run the executable with eps=delta=0.05 and obtain the 
estimate mu=0.501296.