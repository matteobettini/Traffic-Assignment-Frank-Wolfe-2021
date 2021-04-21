# Traffic-Assignment-Frank-Wolfe-2021

This simple script computes the traffic assignment using the **Frank-Wolfe algorithm (FW)** or the **Method of succesive averages (MSA)**.

It can compute the **User Equilibrium (UE)** assignment or the **System Optimal (SO)** assignment.

The travel time cost function that models the effect of congestion on travel time is pluggable and definable by the users.

Currently, we provide three cost function implementations:
* BPR cost function ([see more](https://rdrr.io/rforge/travelr/man/bpr.function.html))
* Greenshields cost function (see Greenshields, B. D., et al. "A study of traffic capacity." Highway research board proceedings. Vol. 1935. National Research Council (USA), Highway Research Board, 1935.)
* Constant cost function (no congestion effects)

# How to use

To use the script, simply call the `computeAssingment` method in the `assignment.py` script.

The documentation of the method provides a through description of all the available parameters and their meaning.

# Importing networks
 Networks and demand files must be specified in the TNTP data format.
 
 A through description of the TNTP format and a wide raange of real transportation networks to test the algorithm on is avaialble at [TransportationNetworks](https://github.com/bstabler/TransportationNetworks)
 
 # Acknowledgments
 
* This work is based on [Traffic-Assignment](https://github.com/prameshk/Traffic-Assignment). We focused on fixing this implementation and extanding it to pluggable cost functions and user optimal flows.
* All the networks we used for testing the correctness of the algorithm are available at [TransportationNetworks](https://github.com/bstabler/TransportationNetworks).
