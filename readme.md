# Distributed matrix multiplication

The project was created for the sake of the High-Performance Computing course. BTW IMO that's the best course I have during my Bachelor's study at MIMUW. Anyway, I've implemented two similar distributed algorithms using MPI:
1. 1.5D Blocked Column A
2. 1.5D Blocked Inner ABC

These algorithms solve the problem of multiplying sparse matrix times dense. They replicate redundantly the data over the nodes to reduce transmission time. You can find out more [here](https://people.eecs.berkeley.edu/~yelick/papers/spdmmm16.pdf)




