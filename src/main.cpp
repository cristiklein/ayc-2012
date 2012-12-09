/*!
 * \file main.cpp
 * \brief This file contains source code that solves the Work Hard - Play Hard problem for the Acceler8 contest
 */
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <string>
#include <list>
#include <map>
#include <math.h>
#include <set>
#include <vector>
#include <fstream>
#include <sys/time.h>
#include <string.h>
#include <time.h>
#include <unordered_map>
#include <queue>
#include <vector>

#include "Alliances.hh"
#include "Flights.hh"
#include "TimeMe.hh"
#include "UniqueId.hh"

namespace Planner {

using namespace std;

/* Instantiate templates */
typedef ::Alliances<> Alliances;
typedef ::Flight<> Flight;
typedef ::Flights<> Flights;
typedef ::UniqueId<>::Id Id;
typedef ::UniqueId<>::Name Name;
typedef unsigned long Time;
typedef float Cost;

/**
 * \struct Parameters
 * \brief Store the program's parameters.
 * This structure don't need to be modified but feel free to change it if you want.
 */
struct Parameters{
	Id from;/*!< The city where the travel begins */
	Id to;/*!< The city where the conference takes place */
	unsigned long dep_time_min;/*!< The minimum departure time for the conference (epoch). No flight towards the conference's city must be scheduled before this time. */
	unsigned long dep_time_max;/*!< The maximum departure time for the conference (epoch). No flight towards the conference's city must be scheduled after this time.  */
	unsigned long ar_time_min;/*!< The minimum arrival time after the conference (epoch). No flight from the conference's city must be scheduled before this time.  */
	unsigned long ar_time_max;/*!< The maximum arrival time after the conference (epoch). No flight from the conference's city must be scheduled after this time.  */
	unsigned long max_layover_time;/*!< You don't want to wait more than this amount of time at the airport between 2 flights (in seconds) */
	unsigned long vacation_time_min;/*!< Your minimum vacation time (in seconds). You can't be in a plane during this time. */
	unsigned long vacation_time_max;/*!< Your maximum vacation time (in seconds). You can't be in a plane during this time. */
	list<Id> airports_of_interest;/*!< The list of cities you are interested in. */
	string flights_file;/*!< The name of the file containing the flights. */
	string alliances_file;/*!< The name of the file containing the company alliances. */
	string work_hard_file;/*!< The file used to output the work hard result. */
	string play_hard_file;/*!< The file used to output the play hard result. */
	int nb_threads;/*!< The maximum number of worker threads */
};

/*!
 * \brief Store a travel
 */
struct Travel {
	vector<const Flight *> flights; //!< List of flights
	vector<float> discounts; //!< Applied discounts
	float totalCost; //!< Total cost of the travel
};

time_t convert_string_to_timestamp(const string &s);
void print_params(Parameters &parameters);
void print_flight(const Flight& flight, float discount, ostream& output);
void read_parameters(Parameters& parameters, int argc, char **argv);
void split_string(vector<string>& result, string line, char separator);
void parse_flight(Flight& flight, string& line);
float parse_flights(Flights& flights, string filename);
void parse_alliance(Alliances &alliance, string line);
void parse_alliances(Alliances &alliances, string filename);
bool has_just_traveled_with_company(const Flight& flight_before, const Flight& current_flight);
bool has_just_traveled_with_alliance(const Flight& flight_before, const Flight& current_flight, const Alliances& alliances);
vector<float> apply_discount(const Travel & travel, const Alliances&alliances);
float compute_cost(const Travel & travel, const Alliances&alliances);
void print_alliances(const Alliances &alliances);
void print_flights(const vector<Flight>& flights, const vector<float>& discounts, ostream& output);
bool never_traveled_to(Travel travel, Id city);
void print_travel(const Travel& travel, const Alliances&alliances, ostream& output);
void compute_path(const Flights& flights, Id to, vector<Travel>& travels, unsigned long t_min, unsigned long t_max, const Parameters &parameters, const Alliances &alliances);
Travel find_cheapest(const vector<Travel>& flights, const Alliances& alliances);
Travel find_cheapest(const vector<Travel>& inbounds, const vector<Travel>& outbounds, const Alliances&alliances);
Travel find_cheapest(const vector<Travel>& inbounds, const vector<Travel>& vias, const vector<Travel>& outbounds, const Alliances&alliances);
void fill_travel(vector<Travel>& travels, const Flights& flights, Id starting_point, unsigned long t_min, unsigned long t_max);
Travel work_hard(const Flights& flights, Parameters& parameters, const Alliances& alliances);
vector<Travel> play_hard(const Flights& flights, Parameters& parameters, const Alliances& alliances);
void output_play_hard(const Flights& flights, Parameters& parameters, const Alliances& alliances);
void output_work_hard(const Flights& flights, Parameters& parameters, const Alliances& alliances);

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
	float discount = 1;
	if (f1.company == f2.company)
		discount = 0.7;
	else if (alliances.areAllied(f1.company, f2.company))
		discount = 0.8;
	return discount;
}

vector<Travel> computePath(
	const Alliances& alliances,
	const Flights& flights,
	Id from, Id to,
	Time tMin, Time tMax,
	Time maxLayover,
	Cost maxDiscount)
{
	struct Segment {
		const struct Segment *prev; //!< pointer to previous segment
		Id airport; //!< airport in which we are now
		Time landingTime; //!< time from which we can take next plan
		const Flight *flight; //!< flight taken to reach this city
		float prevTotalCost; //!< total cost - minus price (cost * discount) of last flight
		float totalCost; //!< total cost
		float discount; //!< discount applied to last flight
	};

	vector<const Segment *> allSegments; //!< Store all segments for garbage collection
	priority_queue<pair<float, const Segment *>> queue; //!< Queue of segments to process

	/* Add initial segments to queue */
	auto range = flights.takeoffs(from,
		tMin,
		numeric_limits<Time>::max());
	for(const Flight &newFlight : range) {
		if (newFlight.landingTime > tMax)
			continue;

		Segment *seg = new Segment();
		seg->prev = NULL;
		seg->airport = newFlight.to;
		seg->landingTime = newFlight.landingTime;
		seg->flight = &newFlight;
		seg->prevTotalCost = 0;
		seg->totalCost = newFlight.cost;
		seg->discount = 1;
		queue.push(pair<float, const Segment *>(-seg->totalCost, seg));
	}

	float cheapestTravel = INFINITY;
	vector<const Segment *> finalSegments;
	while (!queue.empty()) {
		const Segment *seg = queue.top().second;
		queue.pop();

		/* Pruning */
		if (seg->totalCost - maxDiscount > cheapestTravel)
			break;

		/* Check if we reached our destination */
		if (seg->airport == to) {
			finalSegments.push_back(seg);
			cheapestTravel = min(cheapestTravel, seg->totalCost);
			continue;
		}

		/* Where can we go to from here? */
		auto range = flights.takeoffs(seg->airport,
			seg->landingTime,
			seg->landingTime + maxLayover);
		for(const Flight &newFlight : range) {
			/* Out of time? */
			if (newFlight.landingTime > tMax)
				continue;

			/* Avoid cycles */
			bool wasHereBefore = false;
			for (const Segment *s = seg; s != NULL; s = s->prev) {
				if (s->airport == newFlight.to) {
					wasHereBefore = true;
					break;
				}
			}
			if (wasHereBefore)
				continue;

			/* Everything seems okey, prepare new segment to be added to the queue */
			Segment *newSeg = new Segment();
			newSeg->prev = seg;
			newSeg->airport = newFlight.to;
			newSeg->landingTime = newFlight.landingTime;
			newSeg->flight = &newFlight;

			/* Compute cost (quite tricky, but works) */
			float discount = getDiscount(alliances, newFlight, *seg->flight);
			newSeg->prevTotalCost = seg->prevTotalCost + seg->flight->cost * min(discount, seg->discount);
			newSeg->totalCost = seg->prevTotalCost + seg->flight->cost * discount;
			newSeg->discount = discount;

			/* Done preparing new segment, add to queue */
			allSegments.push_back(newSeg);
			queue.push(make_pair<float, const Segment *>(-newSeg->totalCost, newSeg));
		}
	}

	vector<Travel> travels;
	for (const Segment *seg : finalSegments) {
		Travel travel;
		for (const Segment *s = seg; s != NULL && s->flight != NULL; s = s->prev) {
			travel.flights.push_back(s->flight);
			travel.discounts.push_back(s->discount);
		}
		travel.totalCost = seg->totalCost;
		reverse(travel.flights.begin(), travel.flights.end());
		travels.push_back(travel);
	}

	for (const Segment *s : allSegments)
		delete s;
	return travels;
}

float computeCostAfterMerger(const Alliances &alliances, const Travel &travelAB, const Travel &travelBC, const Travel &travelCD = Travel())
{
	float totalCost = travelAB.totalCost + travelBC.totalCost + travelCD.totalCost;

	float discountInB = getDiscount(alliances, *travelAB.flights.back(), *travelBC.flights.front());
	totalCost -= (travelAB.discounts.back()  - discountInB) * travelAB.flights.back()->cost;
	totalCost -= (travelBC.discounts.front() - discountInB) * travelBC.flights.front()->cost;

	if (travelCD.flights.size()) {
		float discountInC = getDiscount(alliances, *travelBC.flights.back(), *travelCD.flights.front());
		totalCost -= (travelBC.discounts.back()  - discountInC) * travelBC.flights.back()->cost;
		totalCost -= (travelCD.discounts.front() - discountInC) * travelCD.flights.front()->cost;
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

Travel findCheapestAndMerge(const Alliances &alliances, const vector<Travel> &inbounds, const vector<Travel> &outbounds)
{
	Travel bestTravel;
	bestTravel.totalCost = INFINITY;

	/* Compute valid prunings */
	float highestInboundLast = 0, highestOutboundFirst = 0;
	float cheapestInbound = INFINITY, cheapestOutbound = INFINITY;
	for (const Travel &inbound : inbounds) {
		highestInboundLast = max(highestInboundLast, inbound.flights.back()->cost);
		cheapestInbound = min(cheapestInbound, inbound.totalCost);
	}
	for (const Travel &outbound : outbounds) {
		highestOutboundFirst = max(highestOutboundFirst, outbound.flights.front()->cost);
		cheapestOutbound = min(cheapestOutbound, outbound.totalCost);
	}

	/* Do cartezian product, do pruning, find best */
	for (const Travel &inbound : inbounds) {
		const Flight &lastInbound = *inbound.flights.back();

		/* Prune this travel, if, assuming best discounts, it cannot be better than the cheapest inbound */
		if (inbound.totalCost - 0.3 * lastInbound.cost - highestOutboundFirst * 0.3 > cheapestInbound)
			continue;

		for (const Travel &outbound : outbounds) {
			const Flight &firstOutbound = *outbound.flights.front();

			/* Prune this travel, if, assuming best discounts, it cannot be better than the cheapest outbound */
			if (outbound.totalCost - 0.3 * firstOutbound.cost - highestInboundLast * 0.3 > cheapestOutbound)
				continue;

			/* Compute cost after merger */
			float discount = getDiscount(alliances, lastInbound, firstOutbound);
			float cost =
				inbound.totalCost  - (inbound.discounts.back()  - discount) * lastInbound.cost +
				outbound.totalCost - (outbound.discounts.back() - discount) * firstOutbound.cost;
			
			if (cost < bestTravel.totalCost) {
				bestTravel = inbound;
				bestTravel.flights.insert(bestTravel.flights.end(), outbound.flights.begin(), outbound.flights.end());
				bestTravel.discounts.insert(bestTravel.discounts.end(), outbound.discounts.begin(), outbound.discounts.end());
				bestTravel.discounts[inbound.discounts.size()-1] = discount;
				bestTravel.discounts[inbound.discounts.size()  ] = discount;
				bestTravel.totalCost = cost;
			}

		}
	}
	return bestTravel;
}

Travel findCheapestAndMerge(const Alliances &alliances, const vector<Travel> &travelsAB, const vector<Travel> &travelsBC, const vector<Travel> &travelsCD)
{
	/* Notation: the airports are named as follows A -> B -> C -> D
	 * Hence, the name of the flights are AB, BC and CD
	 */

	const Travel *bestTravelAB = NULL, *bestTravelBC = NULL, *bestTravelCD = NULL;
	float bestCost = INFINITY;

	/* Compute valid prunings */
	float maxToB   = priciestLastFlight(travelsAB);
	float maxFromB = priciestFirstFlight(travelsBC);
	float maxToC   = priciestLastFlight(travelsBC);
	float maxFromC = priciestFirstFlight(travelsCD);
	float minAB = cheapestTravel(travelsAB);
	float minBC = cheapestTravel(travelsBC);
	float minCD = cheapestTravel(travelsCD);

	/* Do cartezian product, do pruning, find best */
	for (const Travel &travelAB : travelsAB) {
		const Flight &lastFlight = *travelAB.flights.back();

		/* Prune this travel, if, assuming best discounts, it cannot be better than the cheapest choice */
		if (travelAB.totalCost - 0.3 * lastFlight.cost - maxFromB * 0.3 > minAB)
			continue;

		for (const Travel &travelCD : travelsCD) {
			const Flight &firstFlight = *travelCD.flights.front();

			/* Prune this travel, if, assuming best discounts, it cannot be better than the cheapest choice */
			if (travelCD.totalCost - 0.3 * firstFlight.cost - maxToC * 0.3 > minCD)
				continue;

			for (const Travel &travelBC : travelsBC) {
				const Flight &lastFlight = *travelBC.flights.back();
				const Flight &firstFlight = *travelBC.flights.front();

				/* Prune this travel, if, assuming best discounts, it cannot be better than the cheapest choice */
				/* Probably useless, since we compute this later anyway */
				if (travelBC.totalCost
					- 0.3 * firstFlight.cost - 0.3 * lastFlight.cost
					- maxToB * 0.3 - maxFromC > minBC)
					continue;

				float newCost = computeCostAfterMerger(alliances, travelAB, travelBC, travelCD);
				/* check */
				Travel mergedTravel = mergeTravels(alliances, travelAB, travelBC, travelCD);
				assert(newCost == mergedTravel.totalCost);

				if (newCost < bestCost) {
					bestCost = newCost;
					bestTravelAB = &travelAB;
					bestTravelBC = &travelBC;
					bestTravelCD = &travelCD;
				}
			}

		}
	}

	return mergeTravels(alliances, *bestTravelAB, *bestTravelBC, *bestTravelCD);
}

Travel workHard(const Alliances& alliances, const Flights& flights, const Parameters& parameters)
{
	/* Get most expensive flights from/to conference */
	float maxFromConf = priciestFlight(flights.takeoffs(parameters.to, parameters.ar_time_min, parameters.ar_time_max));
	float maxToConf   = priciestFlight(flights.landings(parameters.to, parameters.dep_time_min, parameters.dep_time_max));

	vector<Travel> inbounds = computePath(
		alliances, flights, /* description about the world */
		parameters.from, parameters.to, /* source, destination airport */
		parameters.dep_time_min, parameters.dep_time_max, /* interval of time during which to fly */
		parameters.max_layover_time, /* other trip parameters */
		(maxToConf + maxFromConf) * 0.3 /* pruning parameter */);
	vector<Travel> outbounds = computePath(
		alliances, flights,
		parameters.to, parameters.from,
		parameters.ar_time_min, parameters.ar_time_max,
		parameters.max_layover_time,
		(maxToConf + maxFromConf) * 0.3 /* pruning parameter */);

	return findCheapestAndMerge(alliances, inbounds, outbounds);
}

/**
 * \fn vector<Travel> play_hard(vector<Flight>& flights, Parameters& parameters, const Alliances& alliances)
 * \brief Solve the "Play Hard" problem.
 * This problem can be considered as the hard one. The goal is to find the cheapest way to join a point B from a point A regarding some parameters and for each city in the vacation destination list.
 * \param flights The list of available flights.
 * \param parameters The parameters.
 * \param alliances The alliances between companies.
 * \return The cheapest trips found ordered by vacation destination (only the best result for each vacation destination).
 */
#if 0
vector<Travel> play_hard(Flights& flights, Parameters& parameters, const Alliances& alliances){
	vector<Travel> results;
	list<Id>::iterator it = parameters.airports_of_interest.begin();
	for(; it != parameters.airports_of_interest.end(); it++){
		Id current_airport_of_interest = *it;
		vector<Travel> all_travels;
		/*
		 * The first part compute a travel from home -> vacation -> conference -> home
		 */
		vector<Travel> home_to_vacation, vacation_to_conference, conference_to_home;
		//compute the paths from home to vacation
		fill_travel(home_to_vacation, flights, parameters.from, parameters.dep_time_min-parameters.vacation_time_max, parameters.dep_time_min-parameters.vacation_time_min);
		compute_path(flights, current_airport_of_interest, home_to_vacation, parameters.dep_time_min-parameters.vacation_time_max, parameters.dep_time_min-parameters.vacation_time_min, parameters, alliances);
		//compute the paths from vacation to conference
		fill_travel(vacation_to_conference, flights, current_airport_of_interest, parameters.dep_time_min, parameters.dep_time_max);
		compute_path(flights, parameters.to, vacation_to_conference, parameters.dep_time_min, parameters.dep_time_max, parameters, alliances);
		//compute the paths from conference to home
		fill_travel(conference_to_home, flights, parameters.to, parameters.ar_time_min, parameters.ar_time_max);
		compute_path(flights, parameters.from, conference_to_home, parameters.ar_time_min, parameters.ar_time_max, parameters, alliances);
		Travel best1 = find_cheapest(home_to_vacation, vacation_to_conference, conference_to_home, alliances);

		/*
		 * The second part compute a travel from home -> conference -> vacation -> home
		 */
		vector<Travel> home_to_conference, conference_to_vacation, vacation_to_home;
		//compute the paths from home to conference
		fill_travel(home_to_conference, flights, parameters.from, parameters.dep_time_min, parameters.dep_time_max);
		compute_path(flights, parameters.to, home_to_conference, parameters.dep_time_min, parameters.dep_time_max, parameters, alliances);
		//compute the paths from conference to vacation
		fill_travel(conference_to_vacation, flights, parameters.to, parameters.ar_time_min, parameters.ar_time_max);
		compute_path(flights, current_airport_of_interest, conference_to_vacation, parameters.ar_time_min, parameters.ar_time_max, parameters, alliances);
		//compute paths from vacation to home
		fill_travel(vacation_to_home, flights, current_airport_of_interest, parameters.ar_time_max+parameters.vacation_time_min, parameters.ar_time_max+parameters.vacation_time_max);
		compute_path(flights, parameters.from, vacation_to_home, parameters.ar_time_max+parameters.vacation_time_min, parameters.ar_time_max+parameters.vacation_time_max, parameters, alliances);
		Travel best2 = find_cheapest(home_to_conference, conference_to_vacation, vacation_to_home, alliances);

		float cost1 = best1.flights.size() ? compute_cost(best1, alliances) : INFINITY;
		float cost2 = best2.flights.size() ? compute_cost(best2, alliances) : INFINITY;
		if (cost1 < cost2)
			results.push_back(best1);
		else
			results.push_back(best2);
	}
	return results;
}
#endif

/**
 * \fn void apply_discount(Travel & travel, const Alliances&alliances)
 * \brief Apply a discount when possible to the flights of a travel.
 * \param travel A travel (it will be modified to apply the discounts).
 * \param alliances The alliances.
 */
vector<float> apply_discount(const Travel & travel, const Alliances &alliances){
	vector<float> discounts(travel.flights.size(), 1);
	for(unsigned int i=1; i<travel.flights.size(); i++){
		const Flight& flight_before = *travel.flights[i-1];
		const Flight& current_flight = *travel.flights[i];
		if (flight_before.company == current_flight.company){
			discounts[i-1] = 0.7;
			discounts[i  ] = 0.7;
		}else if(alliances.areAllied(flight_before.company, current_flight.company)){
			if(discounts[i-1] > 0.8)
				discounts[i-1] = 0.8;
			discounts[i] = 0.8;
		}
	}
	return discounts;
}

/**
 * \fn float compute_cost(Travel & travel, const Alliances&alliances)
 * \brief Compute the cost of a travel and uses the discounts when possible.
 * \param travel The travel.
 * \param alliances The alliances.
 */
float compute_cost(const Travel & travel, const Alliances&alliances){
	vector<float> discounts = apply_discount(travel, alliances);
	float result = 0;
	for(unsigned int i=0; i<travel.flights.size(); i++){
		result += (travel.flights[i]->cost * discounts[i]);
	}
	return result;
}


/**
 * \fn Travel find_cheapest(vector<Travel>& travels, const Alliances&alliances)
 * \brief Finds the cheapest travel amongst the travels's vector.
 * \param travels A vector of acceptable travels
 * \param alliances The alliances
 * \return The cheapest travel found.
 */
Travel find_cheapest(const vector<Travel>& travels, const Alliances&alliances){
	const Travel *bestTravel = NULL;
	float bestCost = INFINITY;

	for (const Travel &travel : travels) {
		float currentCost = compute_cost(travel, alliances);
		if (currentCost < bestCost) {
			bestCost = currentCost;
			bestTravel = &travel;
		}
	}
	return *bestTravel;
}

/**
 * \fn void fill_travel(vector<Travel>& travels, vector<Flight>& flights, Id starting_point, unsigned long t_min, unsigned long t_max)
 * \brief Fills the travels's vector with flights that take off from the starting_point.
 * This function might probably be improved.
 * \param travels A vector of travels under construction
 * \param flights All the flights that are available.
 * \param starting_point The starting point.
 * \param travels The list of possible travels that we are building.
 * \param t_min You must not be in a plane before this value (epoch).
 * \param t_max You must not be in a plane after this value (epoch).
 */
void fill_travel(vector<Travel>& travels, const Flights& flights, Id from, unsigned long t_min, unsigned long t_max){
	for(const Flight &flight : flights.takeoffs(from, t_min, t_max)) {
		if(flight.landingTime <= t_max){
			Travel t;
			t.flights.push_back(&flight);
			travels.push_back(t);
		}
	}
}

Travel find_cheapest(const vector<Travel> &inbounds, const vector<Travel> &outbounds, const Alliances &alliances) {
	Travel bestTravel;
	float bestCost = INFINITY;

	/* Compute valid cutoffs */
	float highestInboundLast = 0, highestOutboundFirst = 0;
	for (const Travel &inbound : inbounds)
		highestInboundLast = max(highestInboundLast, inbound.flights.back()->cost);
	for (const Travel &outbound : outbounds)
		highestOutboundFirst = max(highestOutboundFirst, outbound.flights.front()->cost);

	float cheapestInbound  = compute_cost(find_cheapest(inbounds, alliances), alliances);
	float cheapestOutbound = compute_cost(find_cheapest(outbounds, alliances), alliances);

	for (const Travel &inbound : inbounds) {
		if (compute_cost(inbound, alliances) - 0.3 * inbound.flights.back()->cost - highestOutboundFirst * 0.3 > cheapestInbound) continue;
		for (const Travel &outbound : outbounds) {
			if (compute_cost(outbound, alliances) - 0.3 * outbound.flights.front()->cost - highestInboundLast * 0.3 > cheapestOutbound)
				continue;
			Travel travel;
			travel.flights = inbound.flights;
			travel.flights.insert(travel.flights.end(), outbound.flights.begin(), outbound.flights.end());
			float cost = compute_cost(travel, alliances);
			if (cost < bestCost) {
				bestTravel = travel;
				bestCost = cost;
			}
		}
	}
	return bestTravel;
}

Travel find_cheapest(const vector<Travel> &inbounds, const vector<Travel> &vias, const vector<Travel> &outbounds, const Alliances &alliances) {
	Travel bestTravel;
	float bestCost = INFINITY;

	if (inbounds.size() == 0 ||
		vias.size() == 0 ||
		outbounds.size() == 0)
		return bestTravel;
	
	/* Compute valid cutoffs */
	float highestInboundLast = 0, highestViaFirst = 0, highestViaLast = 0, highestOutboundFirst = 0;
	for (const Travel &inbound : inbounds)
		highestInboundLast = max(highestInboundLast, inbound.flights.back()->cost);
	for (const Travel &outbound : outbounds)
		highestOutboundFirst = max(highestOutboundFirst, outbound.flights.front()->cost);
	for (const Travel &via : vias) {
		highestViaFirst = max(highestViaFirst, via.flights.front()->cost);
		highestViaLast  = max(highestViaLast , via.flights.back()->cost);
	}

	float cheapestIn  = compute_cost(find_cheapest(inbounds, alliances), alliances);
	float cheapestVia = compute_cost(find_cheapest(vias, alliances), alliances);
	float cheapestOut = compute_cost(find_cheapest(outbounds, alliances), alliances);

	for (const Travel &inbound : inbounds) {
		if (compute_cost(inbound, alliances) - 0.3 * inbound.flights.back()->cost - 0.3 * highestViaFirst > cheapestIn) continue;
		for (const Travel &via : vias) {
			if (compute_cost(via, alliances)
				- 0.3 * via.flights.front()->cost
				- 0.3 * via.flights.back()->cost
				- 0.3 * highestInboundLast
				- 0.3 * highestOutboundFirst > cheapestVia)
				continue;
			for (const Travel &outbound : outbounds) {
				if (compute_cost(outbound, alliances) - 0.3 * outbound.flights.front()->cost - 0.3 * highestViaLast > cheapestOut) continue;
				Travel travel;
				travel.flights = inbound.flights;
				travel.flights.insert(travel.flights.end(), via.flights.begin(), via.flights.end());
				travel.flights.insert(travel.flights.end(), outbound.flights.begin(), outbound.flights.end());
				float cost = compute_cost(travel, alliances);
				if (cost < bestCost) {
					bestTravel = travel;
					bestCost = cost;
				}
			}
		}
	}
	return bestTravel;
}


/**
 * \fn time_t convert_string_to_timestamp(const string &s)
 * \brief Parses the string s and returns a timestamp (epoch)
 * \param s A string that represents a date with the following format MMDDYYYYhhmmss with
 * M = Month number
 * D = Day number
 * Y = Year number
 * h = hour number
 * m = minute number
 * s = second number
 * You shouldn't modify this part of the code unless you know what you are doing.
 * \return a timestamp (epoch) corresponding to the given parameters.
 */
time_t convert_string_to_timestamp(const string &s){
	if(s.size() != 14){
		cerr<<"The given string is not a valid timestamp"<<endl;
		exit(0);
	}else{
		const char *c = s.c_str();
		tm time;
		time.tm_year = (c[4]-'0') * 1000 + (c[5]-'0') * 100 + (c[6]-'0') * 10 + (c[7]-'0') - 1900;
		time.tm_mon  = (c[ 0]-'0') * 10 + c[ 1]-'0' - 1;
		time.tm_mday = (c[ 2]-'0') * 10 + c[ 3]-'0';
		time.tm_hour = (c[ 8]-'0') * 10 + c[ 9]-'0';
		time.tm_min  = (c[10]-'0') * 10 + c[11]-'0';
		time.tm_sec  = (c[12]-'0') * 10 + c[13]-'0';
		return timegm(&time);
	}
}

/**
 * \fn void print_params(Parameters &parameters)
 * \brief You can use this function to display the parameters
 */
void print_params(Parameters &parameters){
	cout<<"From : "					<<parameters.from					<<endl;
	cout<<"To : "					<<parameters.to						<<endl;
	cout<<"dep_time_min : "			<<parameters.dep_time_min			<<endl;
	cout<<"dep_time_max : "			<<parameters.dep_time_max			<<endl;
	cout<<"ar_time_min : "			<<parameters.ar_time_min			<<endl;
	cout<<"ar_time_max : "			<<parameters.ar_time_max			<<endl;
	cout<<"max_layover_time : "		<<parameters.max_layover_time		<<endl;
	cout<<"vacation_time_min : "	<<parameters.vacation_time_min		<<endl;
	cout<<"vacation_time_max : "	<<parameters.vacation_time_max		<<endl;
	cout<<"flights_file : "			<<parameters.flights_file			<<endl;
	cout<<"alliances_file : "		<<parameters.alliances_file			<<endl;
	cout<<"work_hard_file : "		<<parameters.work_hard_file			<<endl;
	cout<<"play_hard_file : "		<<parameters.play_hard_file			<<endl;
	list<Id>::iterator it = parameters.airports_of_interest.begin();
	for(; it != parameters.airports_of_interest.end(); it++)
		cout<<"airports_of_interest : "	<<*it	<<endl;
	cout<<"flights : "				<<parameters.flights_file			<<endl;
	cout<<"alliances : "			<<parameters.alliances_file			<<endl;
	cout<<"nb_threads : "			<<parameters.nb_threads				<<endl;
}

/**
 * \fn void print_flight(const Flight& flight, ostream& output)
 * \brief You can use this function to display a flight
 */
void print_flight(const UniqueId<> &uniqueId, const Flight& flight, float discount, ostream& output){
	struct tm * take_off_t, *land_t;
	take_off_t = gmtime(((const time_t*)&(flight.takeoffTime)));
	output<<uniqueId.getName(flight.company)<<"-";
	output<<""<<uniqueId.getName(flight.id)<<"-";
	output<<uniqueId.getName(flight.from)<<" ("<<(take_off_t->tm_mon+1)<<"/"<<take_off_t->tm_mday<<" "<<take_off_t->tm_hour<<"h"<<take_off_t->tm_min<<"min"<<")"<<"/";
	land_t = gmtime(((const time_t*)&(flight.landingTime)));
	output<<uniqueId.getName(flight.to)<<" ("<<(land_t->tm_mon+1)<<"/"<<land_t->tm_mday<<" "<<land_t->tm_hour<<"h"<<land_t->tm_min<<"min"<<")-";
	output<<flight.cost<<"$"<<"-"<<discount*100<<"%"<<endl;

}

/**
 * \fn void read_parameters(Parameters& parameters, int argc, char **argv)
 * \brief This function is used to read the parameters
 * \param parameters Represents the structure that will be filled with the parameters.
 */
void read_parameters(UniqueId<> &uniqueId, Parameters& parameters, int argc, char **argv){
	for(int i=0; i<argc; i++){
		string current_parameter = argv[i];
		if(current_parameter == "-from"){
			parameters.from = uniqueId.getId(argv[++i]);
		}else if(current_parameter == "-arrival_time_min"){
			parameters.ar_time_min = convert_string_to_timestamp(argv[++i]);
		}else if(current_parameter == "-arrival_time_max"){
			parameters.ar_time_max = convert_string_to_timestamp(argv[++i]);
		}else if(current_parameter == "-to"){
			parameters.to = uniqueId.getId(argv[++i]);
		}else if(current_parameter == "-departure_time_min"){
			parameters.dep_time_min = convert_string_to_timestamp(argv[++i]);
		}else if(current_parameter == "-departure_time_max"){
			parameters.dep_time_max = convert_string_to_timestamp(argv[++i]);
		}else if(current_parameter == "-max_layover"){
			parameters.max_layover_time = atol(argv[++i]);
		}else if(current_parameter == "-vacation_time_min"){
			parameters.vacation_time_min = atol(argv[++i]);
		}else if(current_parameter == "-vacation_time_max"){
			parameters.vacation_time_max = atol(argv[++i]);
		}else if(current_parameter == "-vacation_airports"){
			while(i+1 < argc && argv[i+1][0] != '-'){
				parameters.airports_of_interest.push_back(uniqueId.getId(argv[++i]));
			}
		}else if(current_parameter == "-flights"){
			parameters.flights_file = argv[++i];
		}else if(current_parameter == "-alliances"){
			parameters.alliances_file = argv[++i];
		}else if(current_parameter == "-work_hard_file"){
			parameters.work_hard_file = argv[++i];
		}else if(current_parameter == "-play_hard_file"){
			parameters.play_hard_file = argv[++i];
		}else if(current_parameter == "-nb_threads"){
			parameters.nb_threads = atoi(argv[++i]);
		}

	}
}

/**
 * \fn void split_string(vector<string>& result, string line, char separator)
 * \brief This function split a string into a vector of strings regarding the separator.
 * \param result The vector of separated strings
 * \param line The line that must be split.
 * \param separator The separator character.
 */
void split_string(vector<string>& result, string line, char separator){
	while(line.find(separator) != string::npos){
		size_t pos = line.find(separator);
		result.push_back(line.substr(0, pos));
		line = line.substr(pos+1);
	}
	result.push_back(line);
}

/**
 * \fn void parse_flight(vector<Flight>& flights, string& line)
 * \brief This function parses a line containing a flight description.
 * \param flights The vector of flights.
 * \param line The line that must be parsed.
 */
void parse_flight(UniqueId<> &uniqueId, Flight& flight, string& line){
	vector<string> splittedLine;
	split_string(splittedLine, line, ';');
	if(splittedLine.size() == 7){
		flight.id = uniqueId.getId(splittedLine[0]);
		flight.from = uniqueId.getId(splittedLine[1]);
		flight.takeoffTime = convert_string_to_timestamp(splittedLine[2]);
		flight.to = uniqueId.getId(splittedLine[3]);
		flight.landingTime = convert_string_to_timestamp(splittedLine[4]);
		flight.cost = atof(splittedLine[5].c_str());
		flight.company = uniqueId.getId(splittedLine[6]);
	}
}

/**
 * \fn void parse_flights(vector<Flight>& flights, string filename)
 * \brief This function parses the flights from a file.
 * \param flights The vector of flights.
 * \param filename The name of the file containing the flights.
 */
float parse_flights(UniqueId<> &uniqueId, Flights& flights, string filename){
	string line = "";
	ifstream file;
	file.open(filename.c_str());
	if(!file.is_open()){
		cerr<<"Problem while opening the file "<<filename<<endl;
		exit(0);
	}

	float highestCost = 0;
	while (!file.eof())
	{
		Flight flight;
		getline(file, line);
		parse_flight(uniqueId, flight, line);
		flights.add(flight);
		highestCost = max(highestCost, flight.cost);
	}
	return highestCost;
}

/**
 * \fn void parse_alliance(Alliances &alliance, string line)
 * \brief This function parses a line containing alliances between companies.
 * \param alliance A vector of companies sharing a same alliance.
 * \param line A line that contains the name of companies in the same alliance.
 */
void parse_alliance(UniqueId<> &uniqueId, Alliances &alliances, string line){
	vector<string> splittedLine;
	split_string(splittedLine, line, ';');
	for(unsigned int i=0; i<splittedLine.size(); i++)
	for(unsigned int j=i+1; j<splittedLine.size(); j++){
		Id c1 = uniqueId.getId(splittedLine[i]);
		Id c2 = uniqueId.getId(splittedLine[j]);
		if (c1 == c2)
			continue;

		alliances.add(c1, c2);
	}
}

/**
 * \fn void parse_alliances(Alliances &alliances, string filename)
 * \brief This function parses a line containing alliances between companies.
 * \param alliances A 2D vector representing the alliances. Companies on the same line are in the same alliance.
 * \param filename The name of the file containing the alliances description.
 */
void parse_alliances(UniqueId<> &uniqueId, Alliances &alliances, string filename){
	string line = "";
	ifstream file;

	file.open(filename.c_str());
	if(!file.is_open()){
		cerr<<"Problem while opening the file "<<filename<<endl;
		exit(0);
	}
	while (!file.eof())
	{
		getline(file, line);
		parse_alliance(uniqueId, alliances, line);
	}
}

/**
 * \fn void print_flights(vector<Flight>& flights, ostream& output)
 * \brief Display the flights on the standard output.
 * \param flights The flights.
 */
void print_flights(const UniqueId<> &uniqueId, const vector<const Flight *>& flights, const vector<float> &discounts, ostream& output){
	for(unsigned int i=0; i<flights.size(); i++)
		print_flight(uniqueId, *flights[i], discounts[i], output);
}

/**
 * \fn void print_travel(Travel& travel, const Alliances&alliances)
 * \brief Display a travel on the standard output.
 * \param travel The travel.
 * \param alliances The alliances (used to compute the price).
 */
void print_travel(const UniqueId<> &uniqueId, const Travel& travel, const Alliances&alliances, ostream& output){
	vector<float> discounts = apply_discount(travel, alliances);
	output<<"Price : "<<compute_cost(travel, alliances)<<endl;
	print_flights(uniqueId, travel.flights, discounts, output);
	output<<endl;
}

/**
 * \fn void output_play_hard(vector<Flight>& flights, Parameters& parameters, const Alliances& alliances)
 * \brief Display the solution of the "Play Hard" problem by solving it first.
 * \param flights The list of available flights.
 * \param parameters The parameters.
 * \param alliances The alliances between companies.
 */
void output_play_hard(const UniqueId<> &uniqueId, Flights& flights, Parameters& parameters, const Alliances& alliances){
	ofstream output;
	output.open(parameters.play_hard_file.c_str());
	vector<Travel> travels;// = play_hard(flights, parameters, alliances);
	list<Id> cities = parameters.airports_of_interest;
	for(unsigned int i=0; i<travels.size(); i++){
		output<<"“Play Hard” Proposition "<<(i+1)<<" : "<<uniqueId.getName(cities.front())<<endl;
		print_travel(uniqueId, travels[i], alliances, output);
		cities.pop_front();
		output<<endl;
	}
	output.close();
}

/**
 * \fn void output_work_hard(vector<Flight>& flights, Parameters& parameters, const Alliances& alliances)
 * \brief Display the solution of the "Work Hard" problem by solving it first.
 * \param flights The list of available flights.
 * \param parameters The parameters.
 * \param alliances The alliances between companies.
 */
void output_work_hard(const UniqueId<> &uniqueId, Flights& flights, Parameters& parameters, const Alliances& alliances){
	ofstream output;
	output.open(parameters.work_hard_file.c_str());
	Travel travel = workHard(alliances, flights, parameters);
	output<<"“Work Hard” Proposition :"<<endl;
	print_travel(uniqueId, travel, alliances, output);
	output.close();
}

} /* namespace Planner */

int main(int argc, char **argv) {
	using Planner::Alliances;
	using Planner::Flights;
	using Planner::Parameters;
	using Planner::parse_alliances;
	using Planner::read_parameters;
	using Planner::parse_flights;

	UniqueId<> uniqueId;

	timeMe("start");
	//Declare variables and read the args
	Parameters parameters;
	Alliances alliances;
	read_parameters(uniqueId, parameters, argc, argv);
//	cout<<"Printing parameters..."<<endl;
//	print_params(parameters);
	Flights flights;
	timeMe("parse params");
	/*parameters.highestCost = */parse_flights(uniqueId, flights, parameters.flights_file);
	timeMe("parse flights");
//	cout<<"Printing flights..."<<endl;
//	print_flights(flights);
//	cout<<"flights printed "<<endl;
	parse_alliances(uniqueId, alliances, parameters.alliances_file);
//	cout<<"Printing alliances..."<<endl;
//	print_alliances(alliances);
	timeMe("parse alliances");
	output_play_hard(uniqueId, flights, parameters, alliances);
	timeMe("play hard");
	output_work_hard(uniqueId, flights, parameters, alliances);
	timeMe("work hard");

	return 0;
}

//./run -from Paris -to Los\ Ang/eles -departure_time_min 11152012000000 -departure_time_max 11172012000000 -arrival_time_min 11222012000000 -arrival_time_max 11252012000000 -max_layover 100000 -vacation_time_min 432000 -vacation_time_max 604800 -vacation_airports Rio London Chicago -flights flights.txt -alliances alliances.txt
