At the core of the sequential playHard algorithm one most find a travel with stops in
the airports A, B, C and D. Let us call the routes AB, BC, CD "subtravels". A first observation to make is that the 3 subtravels can be computed individually. In fact, if we assumed no additional discount was possible when merging the subtravels, one could apply a Dijkstra-like algorithm and find for each subtravel one single cheapest
subtravel, then merge them by concatenation.

However, due to possible additional discount when merging, one cannot stop when the best subtravel has been found. One has to look for additional subtravels, until eventually one finds a subtravel which is so expensive, that no possible discount after merging could make it interesting. For example, any subtravel A->B can be pruned if:

   cost(subtravelAB) > cost(bestAB) + 0.3 * cost(maxToB) + 0.3 * cost(maxFromB)

where bestAB is the best subtravel for AB and maxToB/maxFromB is the priciest flight to/from B. It can easily be proven, that any subtravel on AB more expensive than this cannot lead to a cheaper travel after merger than if using the best subtravel on AB. Similarly, pruning can efficiently be done on subtravels C->D.

For subtravel B->C, one can do "less" pruning, as discounts can come both from flights to/from B and to/from C. Therefore, we first compute subtravels for A->B and C->D, then prune on B->C as follows:

   cost(subtravelBC) > cost(bestBC) + 0.3 * cost(maxFoundToB) + 0.3 * cost(maxFromB)
                                    + 0.3 * cost(maxToC) + 0.3 * cost(maxFoundFromC)

where maxFoundFromC is the most expensive flight from C among the subtravels C->D (after pruning) and maxFoundToB is the most expensive flight from B among the subtravels A->B (again, after pruning).

Additionally, similar pruning can be done while merging. These observations greatly reduce the search space and improve time-to-solution.

