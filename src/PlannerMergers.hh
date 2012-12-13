#ifndef PLANNERMERGERS_HH
#define PLANNERMERGERS_HH

#include "PlannerHelpers.hh"

namespace Planner {

float computeCostAfterMerger(const Alliances &alliances, const Travel &travelAB, const Travel &travelBC, const Travel &travelCD = Travel())
{
	float totalCost = travelAB.totalCost + travelBC.totalCost + travelCD.totalCost;

	float discountInB = getDiscount(alliances, *travelAB.flights.back(), *travelBC.flights.front());
	float prevDiscountAB = travelAB.discounts.back();
	float prevDiscountBC = travelBC.discounts.front();
	totalCost -= (prevDiscountAB - min(prevDiscountAB, discountInB)) * travelAB.flights.back()->cost;
	totalCost -= (prevDiscountBC - min(prevDiscountBC, discountInB)) * travelBC.flights.front()->cost;

	if (travelCD.flights.size()) {
		float discountInC = getDiscount(alliances, *travelBC.flights.back(), *travelCD.flights.front());
		float prevDiscountBC = travelBC.discounts.back();
		if (travelBC.discounts.size() == 1) /* if travelBC only have one single flight */
			prevDiscountBC = min(prevDiscountBC, discountInB); /* make sure we don't apply the discount on it twice */
		float prevDiscountCD = travelCD.discounts.front();
		totalCost -= (prevDiscountBC - min(prevDiscountBC, discountInC)) * travelBC.flights.back()->cost;
		totalCost -= (prevDiscountCD - min(prevDiscountCD, discountInC)) * travelCD.flights.front()->cost;
	}

	return totalCost;
}

Travel mergeTravels(const Alliances &alliances, const Travel &travelAB, const Travel &travelBC, const Travel &travelCD = Travel())
{
	/* Create a new travel with all the flights and discounts */
	Travel mergedTravel = travelAB;

	mergedTravel.flights.insert  (mergedTravel.flights.end()  , travelBC.flights.begin()  , travelBC.flights.end());
	mergedTravel.discounts.insert(mergedTravel.discounts.end(), travelBC.discounts.begin(), travelBC.discounts.end());

	if (travelCD.flights.size()) {
		mergedTravel.flights.insert  (mergedTravel.flights.end()  , travelCD.flights.begin()  , travelCD.flights.end());
		mergedTravel.discounts.insert(mergedTravel.discounts.end(), travelCD.discounts.begin(), travelCD.discounts.end());
	}

	/* Update discounts */
	mergedTravel.discounts[0] = 1;
	for (size_t i = 1; i < mergedTravel.flights.size(); i++) {
		float discount = getDiscount(alliances, *mergedTravel.flights[i-1], *mergedTravel.flights[i]);
		mergedTravel.discounts[i-1] = min(mergedTravel.discounts[i-1], discount);
		mergedTravel.discounts[i] = discount;
	}

	/* Recompute total */
	float totalCost = 0;
	for (size_t i = 0; i < mergedTravel.flights.size(); i++)
		totalCost += mergedTravel.discounts[i] * mergedTravel.flights[i]->cost;
	mergedTravel.totalCost = totalCost;

	return mergedTravel;
}

Travel findCheapestAndMerge(const Alliances &alliances, vector<Travel> &travelsAB, vector<Travel> &travelsBC, vector<Travel> &travelsCD)
{
	/* Variables to store best solution */
	float bestCost = INFINITY;
	const Travel *bestTravelAB = NULL, *bestTravelBC = NULL, *bestTravelCD = NULL;

	/* Variables to store information for "fine" pruning */
	float maxToB, maxFromB, maxToC, maxFromC, minAB, minBC, minCD;

	/* "Coarse pruning */
	while (true) {
		int prevSize = travelsAB.size() + travelsBC.size() + travelsCD.size();

		maxToB   = priciestLastFlight(travelsAB);
		maxFromB = priciestFirstFlight(travelsBC);
		maxToC   = priciestLastFlight(travelsBC);
		maxFromC = priciestFirstFlight(travelsCD);
		minAB = cheapestTravel(travelsAB);
		minBC = cheapestTravel(travelsBC);
		minCD = cheapestTravel(travelsCD);

		float maxABToExplore = minAB + 0.3 * maxToB + 0.3 * maxFromB;
		float maxCDToExplore = minCD + 0.3 * maxToC + 0.3 * maxFromC;
		float maxBCToExplore = minBC + 0.3 * maxToB + 0.3 * maxFromB + 0.3 * maxToC + 0.3 * maxFromC;

		vector<Travel>::iterator it;

		/* We can use upper_bound here, because computePath sorts travels by increasing totalCost */
		it = upper_bound(travelsAB.begin(), travelsAB.end(), maxABToExplore);
		travelsAB.erase(it, travelsAB.end());

		it = upper_bound(travelsBC.begin(), travelsBC.end(), maxBCToExplore);
		travelsBC.erase(it, travelsBC.end());

		it = upper_bound(travelsCD.begin(), travelsCD.end(), maxCDToExplore);
		travelsCD.erase(it, travelsCD.end());

		int newSize = travelsAB.size() + travelsBC.size() + travelsCD.size();
		if (newSize == prevSize)
			break;
	}

	/* Do cartezian product, do pruning, find best */
#pragma omp parallel for schedule(dynamic, 1)
	for (size_t i = 0; i < travelsAB.size(); i++) {
		/* To avoid data dependency between threads, we first compute a
		 * thread-local best, then merge it into the global best in a critical section
		 */
		const Travel *localBestTravelAB = NULL, *localBestTravelBC = NULL, *localBestTravelCD = NULL;
		float localBestCost = INFINITY;

		const Travel &travelAB = travelsAB[i];

		/* Prune this travel, if, assuming best discounts, it cannot be better than the cheapest choice */
		const Flight &lastFlight = *travelAB.flights.back();
		if (travelAB.totalCost - 0.3 * lastFlight.cost - maxFromB * 0.3 > minAB)
			continue;

		for (const Travel &travelCD : travelsCD) {
			/* Prune this travel, if, assuming best discounts, it cannot be better than the cheapest choice */
			const Flight &firstFlight = *travelCD.flights.front();
			if (travelCD.totalCost - 0.3 * firstFlight.cost - maxToC * 0.3 > minCD)
				continue;

			for (const Travel &travelBC : travelsBC) {
				/* Prune this travel, if, assuming best discounts, it cannot be better than the cheapest choice */
				const Flight &lastFlight = *travelBC.flights.back();
				const Flight &firstFlight = *travelBC.flights.front();

				if (travelBC.totalCost
					- 0.3 * firstFlight.cost - 0.3 * lastFlight.cost
					- maxToB * 0.3 - maxFromC > minBC)
					continue;

				float newCost = computeCostAfterMerger(alliances, travelAB, travelBC, travelCD);
				if (newCost < localBestCost) {
					localBestCost = newCost;
					localBestTravelAB = &travelAB;
					localBestTravelBC = &travelBC;
					localBestTravelCD = &travelCD;
				}
			}

		}
		/* Reduce */
		if (localBestCost < bestCost) {
			bestCost = localBestCost;
			bestTravelAB = localBestTravelAB;
			bestTravelBC = localBestTravelBC;
			bestTravelCD = localBestTravelCD;
		}
	}

	if (isfinite(bestCost))
		return mergeTravels(alliances, *bestTravelAB, *bestTravelBC, *bestTravelCD);
	else {
		Travel travel = Travel();
		travel.totalCost = INFINITY;
		return travel;
	}
}

Travel findCheapestAndMerge(const Alliances &alliances, const vector<Travel> &travelsAB, const vector<Travel> &travelsBC)
{
	/* This function is very similar to the one above
	 * one day, they should probabily be merge together, but they have pretty different pruning rules
	 * so for now we have not bothered studying their merger
	 */
	const Travel *bestTravelAB = NULL, *bestTravelBC = NULL;
	float bestCost = INFINITY;

	/* Compute valid prunings */
	float maxToB   = priciestLastFlight(travelsAB);
	float maxFromB = priciestFirstFlight(travelsBC);
	float minAB = cheapestTravel(travelsAB);
	float minBC = cheapestTravel(travelsBC);

	/* Do cartezian product, do pruning, find best */
	for (const Travel &travelAB : travelsAB) {
		const Flight &lastFlight = *travelAB.flights.back();

		/* Prune this travel, if, assuming best discounts, it cannot be better than the cheapest choice */
		if (travelAB.totalCost - 0.3 * lastFlight.cost - maxFromB * 0.3 > minAB)
			continue;

		for (const Travel &travelBC : travelsBC) {
			const Flight &firstFlight = *travelBC.flights.front();

			/* Prune this travel, if, assuming best discounts, it cannot be better than the cheapest choice */
			if (travelBC.totalCost - 0.3 * firstFlight.cost - maxToB * 0.3 > minBC)
				continue;

			/* Compute cost after merger */
			float newCost = computeCostAfterMerger(alliances, travelAB, travelBC);
			if (newCost < bestCost) {
				bestCost = newCost;
				bestTravelAB = &travelAB;
				bestTravelBC = &travelBC;
			}

		}
	}

	if (bestTravelAB)
		return mergeTravels(alliances, *bestTravelAB, *bestTravelBC);
	else {
		Travel travel = Travel();
		travel.totalCost = INFINITY;
		return travel;
	}
}

} /* namespace Planner */

#endif /* PLANNERMERGERS_HH */
