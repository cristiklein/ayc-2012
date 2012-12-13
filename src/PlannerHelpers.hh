#ifndef PLANNERHELPERS_HH
#define PLANNERHELPERS_HH

#include "PlannerDataTypes.hh"

namespace Planner {

float priciestFirstFlight(const vector<Travel> &travels)
{
	float result = 0;
	for (const Travel &travel : travels) {
		float cost = travel.flights.front()->cost;
		if (result < cost)
			result = cost;
	}
	return result;
}

float priciestLastFlight(const vector<Travel> &travels)
{
	float result = 0;
	for (const Travel &travel : travels) {
		float cost = travel.flights.back()->cost;
		if (result < cost)
			result = cost;
	}
	return result;
}

float cheapestTravel(const vector<Travel> &travels)
{
	float result = INFINITY;
	for (const Travel &travel : travels) {
		float cost = travel.totalCost;
		if (result > cost)
			result = cost;
	}
	return result;
}

float priciestFlight(const Flights::Range &flights)
{
	float result = 0;
	for (const Flight &flight : flights) {
		if (result < flight.cost)
			result = flight.cost;
	}
	return result;
}

float getDiscount(const Alliances &alliances, const Flight &f1, const Flight &f2)
{
	if (f1.company == f2.company)
		return 0.7;
	if (alliances.areAllied(f1.company, f2.company))
		return 0.8;
	return 1;
}

} /* namespace Planner */

#endif /* PLANNERHELPERS_HH */
