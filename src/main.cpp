/*!
 * \file main.cpp
 * \brief This file contains source code that solves the Work Hard - Play Hard problem for the Acceler8 contest
 */
#include <algorithm>
#include <assert.h>
#include <iostream>
#include <cstdlib>
#include <omp.h>
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
	vector<Id> airports_of_interest;/*!< The list of cities you are interested in. */
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

/* Used for bisection */
inline bool operator<(float cost, const Travel &travel) { return cost < travel.totalCost; }
inline bool operator<(const Travel &travel, float cost) { return travel.totalCost < cost; }

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
		const Flight *flight; //!< flight taken to reach this city
		float prevTotalCost; //!< total cost - minus price (cost * discount) of last flight (may be different from prev->totalCost, if a lower discount was applied to prev->flight)
		float totalCost; //!< total cost
		float discount; //!< discount applied to last flight (i.e., the one in flight)
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

		/* Faster retrieval */
		Id airport = seg->flight->to;
		Time landingTime = seg->flight->landingTime;

		/* Pruning */
		if (seg->totalCost - maxDiscount > cheapestTravel)
			break;

		/* Check if we reached our destination */
		if (airport == to) {
			finalSegments.push_back(seg);
			cheapestTravel = min(cheapestTravel, seg->totalCost);
			continue;
		}

		/* Where can we go to from here? */
		auto range = flights.takeoffs(airport,
			landingTime,
			landingTime + maxLayover);
		for(const Flight &newFlight : range) {
			/* Out of time? */
			if (newFlight.landingTime > tMax)
				continue;

			/* Avoid cycles */
			/* Implementation note: we tried accelerating this with an unordered_set,
			 * but we would actually observe a 2x slowdown */
			if (newFlight.to == from)
				continue;

			bool wasHereBefore = false;
			for (const Segment *s = seg; s != NULL && !wasHereBefore; s = s->prev)
				if (s->flight->to == newFlight.to)
					wasHereBefore = true;
			if (wasHereBefore)
				continue;

			/* Everything seems okey, prepare new segment to be added to the queue */
			Segment *newSeg = new Segment();
			newSeg->prev = seg;
			newSeg->flight = &newFlight;

			/* Compute cost (quite tricky, but works) */
			float discount = getDiscount(alliances, newFlight, *seg->flight);
			newSeg->prevTotalCost = seg->prevTotalCost + seg->flight->cost * min(discount, seg->discount);
			newSeg->totalCost = newSeg->prevTotalCost + newSeg->flight->cost * discount;
			newSeg->discount = discount;

			/* Done preparing new segment, add to queue */
			allSegments.push_back(newSeg);
			queue.push(make_pair<float, const Segment *>(-newSeg->totalCost, newSeg));
		}
	}

	/* Store segments in a nice travel structure
	 * There is quite complex machinery for fast recomputation of discounts
	 */
	vector<Travel> travels;
	for (const Segment *seg : finalSegments) {
		Travel travel;
		float lastDiscount = 1;
		for (const Segment *s = seg; s != NULL && s->flight != NULL; s = s->prev) {
			travel.flights.push_back(s->flight);
			float discount = 1;
			if (s->prev && s->prev->flight)
				discount = getDiscount(alliances, *s->flight, *s->prev->flight);
			travel.discounts.push_back(min(discount, lastDiscount));
			lastDiscount = discount;
		}
		travel.totalCost = seg->totalCost;
		reverse(travel.flights.begin(), travel.flights.end());
		reverse(travel.discounts.begin(), travel.discounts.end());
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
	/* We need this because OpenMP does not have a construct to customize reduction
	 * Should we use TBB next time?
	 */
	int maxThreads = omp_get_max_threads();
	const Travel *localBestTravelAB[maxThreads], *localBestTravelBC[maxThreads], *localBestTravelCD[maxThreads];
	float localBestCost[maxThreads];
	
	for (int i = 0; i < maxThreads; i++) {
		localBestCost[i] = INFINITY;
		localBestTravelAB[i] = NULL;
		localBestTravelBC[i] = NULL;
		localBestTravelCD[i] = NULL;
	}

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
				int tid = omp_get_thread_num();
				if (newCost < localBestCost[tid]) {
					localBestCost[tid] = newCost;
					localBestTravelAB[tid] = &travelAB;
					localBestTravelBC[tid] = &travelBC;
					localBestTravelCD[tid] = &travelCD;
				}
			}

		}
	}

	/* Reduce */
	float bestCost = INFINITY;
	const Travel *bestTravelAB = NULL, *bestTravelBC = NULL, *bestTravelCD = NULL;
	for (int i = 0; i < maxThreads; i++) {
		if (localBestCost[i] < bestCost) {
			bestCost = localBestCost[i];
			bestTravelAB = localBestTravelAB[i];
			bestTravelBC = localBestTravelBC[i];
			bestTravelCD = localBestTravelCD[i];
		}
	}

	if (bestTravelAB)
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

Travel workHard(const Alliances& alliances, const Flights& flights, const Parameters& parameters, vector<Travel> &homeToConf, vector<Travel> &confToHome)
{
	/* Get most expensive flights from/to conference */
	float maxToConf   = priciestFlight(flights.landings(parameters.to, parameters.dep_time_min, parameters.dep_time_max));
	float maxFromConf = priciestFlight(flights.takeoffs(parameters.to, parameters.ar_time_min, parameters.ar_time_max));

	homeToConf = computePath(
		alliances, flights, /* description about the world */
		parameters.from, parameters.to, /* source, destination airport */
		parameters.dep_time_min, parameters.dep_time_max, /* interval of time during which to fly */
		parameters.max_layover_time, /* other trip parameters */
		(maxToConf + maxFromConf) * 0.3 /* pruning parameter */);
	confToHome = computePath(
		alliances, flights,
		parameters.to, parameters.from,
		parameters.ar_time_min, parameters.ar_time_max,
		parameters.max_layover_time,
		(maxToConf + maxFromConf) * 0.3 /* pruning parameter */);

	return findCheapestAndMerge(alliances, homeToConf, confToHome);
}

struct ABCDTravelParameters {
	Id a, b, c, d;
	Time tminAB, tmaxAB, tminBC, tmaxBC, tminCD, tmaxCD;
	Time maxLayover;
	const vector<Travel> *aToB; //!< Potentially pre-computed values
	const vector<Travel> *cToD; //!< Potentially pre-computed values
};
typedef tuple<vector<Travel>, vector<Travel>, vector<Travel>> ABCDTravelResult;

ABCDTravelResult computeABCDTravel(const Alliances &alliances, const Flights &flights, const ABCDTravelParameters &p)
{
	/* Get most expensive flights from/to vacation (for pruning) */
	float maxToB   = priciestFlight(flights.landings(p.b, p.tminAB, p.tmaxAB));
	float maxFromB = priciestFlight(flights.takeoffs(p.b, p.tminBC, p.tmaxBC));
	float maxToC   = priciestFlight(flights.landings(p.c, p.tminBC, p.tmaxBC));
	float maxFromC = priciestFlight(flights.takeoffs(p.c, p.tminCD, p.tmaxCD));

	vector<Travel> aToB, bToC, cToD;

	/* A -> B */
	if (p.aToB == NULL) {
		aToB = computePath(
			alliances, flights, /* description about the world */
			p.a, p.b, /* source, destination airport */
			p.tminAB, p.tmaxAB, /* interval of time during which to fly */
			p.maxLayover, /* other trip parameters */
			(maxToB + maxFromB) * 0.3 /* pruning parameter */);
	}
	else {
		aToB = *p.aToB;
	}

	/* Update (hopefully reducing) maxB */
	maxToB = 0;
	for (const Travel &travel : aToB)
		maxToB = max(maxToB, travel.flights.back()->cost);
	
	/* C -> D */
	if (p.cToD == NULL) {
		cToD = computePath(
			alliances, flights, /* description about the world */
			p.c, p.d, /* source, destination airport */
			p.tminCD, p.tmaxCD, /* interval of time during which to fly */
			p.maxLayover, /* other trip parameters */
			(maxToC + maxFromC) * 0.3 /* pruning parameter */);
	}
	else {
		cToD = *p.cToD;
	}

	/* Update (hopefully reducing) maxC */
	maxFromC = 0;
	for (const Travel &travel : cToD)
		maxFromC = max(maxFromC, travel.flights.front()->cost);

	/* B -> C */
	bToC = computePath(
		alliances, flights, /* description about the world */
		p.b, p.c, /* source, destination airport */
		p.tminBC, p.tmaxBC, /* interval of time during which to fly */
		p.maxLayover, /* other trip parameters */
		(maxToB + maxFromB + maxToC + maxFromC) * 0.3 /* pruning parameter */);
	
	return ABCDTravelResult(aToB, bToC, cToD);
}

map<Id, Travel> playHard(const Alliances &alliances, const Flights& flights, Parameters& parameters, const vector<Travel> &homeToConf, const vector<Travel> &confToHome)
{
	map<Id, Travel> results;
	int n = parameters.airports_of_interest.size();

	map<Id, ABCDTravelResult> blobs1, blobs2;
	map<Id, Travel> bests1, bests2;

	/* Compute intermediate results */
#pragma omp parallel
#pragma omp single
	for (int i = 0; i < n; i++) {
		Id vacation = parameters.airports_of_interest[i];

		/*
		 * The first part compute a travel from home -> vacation -> conference -> home
		 * We'll use the terminology A -> B -> C -> D and AB BC CD for travels
		 */
		ABCDTravelParameters p1;

		/* Cities */
		p1.a = parameters.from;
		p1.b = vacation;
		p1.c = parameters.to;
		p1.d = parameters.from;

		/* Compute valid travel times */
		p1.tminAB = parameters.dep_time_min - parameters.vacation_time_max;
		p1.tmaxAB = parameters.dep_time_min - parameters.vacation_time_min;
		p1.tminBC = parameters.dep_time_min;
		p1.tmaxBC = parameters.dep_time_max;
		p1.tminCD = parameters.ar_time_min;
		p1.tmaxCD = parameters.ar_time_max;
		p1.maxLayover = parameters.max_layover_time;

		/* Add cached values */
		p1.aToB = NULL;
		p1.cToD = &confToHome;

#pragma omp task shared(alliances, flights) firstprivate(p1) untied
		{
			auto blob1 = computeABCDTravel(alliances, flights, p1);
#pragma omp critical
			blobs1[vacation] = blob1;
		}

		/*
		 * The second part compute a travel from home -> conference -> vacation -> home
		 */
		ABCDTravelParameters p2;

		/* Cities */
		p2.a = parameters.from;
		p2.b = parameters.to;
		p2.c = vacation;
		p2.d = parameters.from;

		/* Compute valid travel times */
		p2.tminAB = parameters.dep_time_min;
		p2.tmaxAB = parameters.dep_time_max;
		p2.tminBC = parameters.ar_time_min;
		p2.tmaxBC = parameters.ar_time_max;
		p2.tminCD = parameters.ar_time_max + parameters.vacation_time_min;
		p2.tmaxCD = parameters.ar_time_max + parameters.vacation_time_max;
		p2.maxLayover = parameters.max_layover_time;

		/* Add cached values */
		p2.aToB = &homeToConf;
		p2.cToD = NULL;

#pragma omp task shared(alliances, flights) firstprivate(p2) untied
		{
			auto blob2 = computeABCDTravel(alliances, flights, p2);
#pragma omp critical
			blobs2[vacation] = blob2;
		}
	}
#pragma omp taskwait

	timeMe("computePath");

	/*
	 * Run time-consuming findCheapest
	 */
	for (size_t i = 0; i < parameters.airports_of_interest.size(); i++) {
		Id vacation = parameters.airports_of_interest[i];
		auto &thisBlob1 = blobs1[vacation];
		bests1[vacation] = findCheapestAndMerge(alliances, get<0>(thisBlob1), get<1>(thisBlob1), get<2>(thisBlob1));
		auto &thisBlob2 = blobs2[vacation];
		bests2[vacation] = findCheapestAndMerge(alliances, get<0>(thisBlob2), get<1>(thisBlob2), get<2>(thisBlob2));

		timeMe("findCheapest");
	}
	timeMe("done findCheapest");

	/*
	 * Reduce results
	 */
	for (Id vacation : parameters.airports_of_interest) {
		if (bests1[vacation].totalCost > bests2[vacation].totalCost)
			results[vacation] = bests2[vacation];
		else
			results[vacation] = bests1[vacation];
	}
	timeMe("done");

	return results;
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

void printTravel(const UniqueId<> &uniqueId, const Travel& travel, ostream& output)
{
	float priceToDisplay = isfinite(travel.totalCost) ? travel.totalCost : 0;
	output << "Price : " << priceToDisplay << endl;

	for (size_t i = 0; i < travel.flights.size(); i++) {
		float discount = travel.discounts[i];
		const Flight &flight = *travel.flights[i];

		struct tm *take_off_t, *land_t;
		take_off_t = gmtime(((const time_t*)&(flight.takeoffTime)));
		output << uniqueId.getName(flight.company) << "-";
		output << "" << uniqueId.getName(flight.id) << "-";
		output << uniqueId.getName(flight.from)<<" ("<<(take_off_t->tm_mon+1)<<"/"<<take_off_t->tm_mday<<" "<<take_off_t->tm_hour<<"h"<<take_off_t->tm_min<<"min"<<")"<<"/";
		land_t = gmtime(((const time_t*)&(flight.landingTime)));
		output<<uniqueId.getName(flight.to)<<" ("<<(land_t->tm_mon+1)<<"/"<<land_t->tm_mday<<" "<<land_t->tm_hour<<"h"<<land_t->tm_min<<"min"<<")-";
		output<<flight.cost<<"$"<<"-"<<discount*100<<"%"<<endl;
	}

	output<<endl;
}

void outputPlayHard(const UniqueId<> &uniqueId, const Parameters &parameters, const std::map<Id, Travel> &travels)
{
	ofstream output;
	output.open(parameters.play_hard_file.c_str());
	vector<Id> cities = parameters.airports_of_interest;
	for(unsigned int i=0; i<travels.size(); i++){
		output<<"“Play Hard” Proposition "<<(i+1)<<" : "<<uniqueId.getName(cities[i])<<endl;
		printTravel(uniqueId, travels.at(cities[i]), output);
		output<<endl;
	}
	output.close();
}

void outputWorkHard(const UniqueId<> &uniqueId, const Parameters& parameters, const Travel &travel)
{
	ofstream output;
	output.open(parameters.work_hard_file.c_str());
	output << "“Work Hard” Proposition :" << endl;
	printTravel(uniqueId, travel, output);
	output.close();
}

} /* namespace Planner */

int main(int argc, char **argv) {
	using Planner::Alliances;
	using Planner::Flights;
	using Planner::Parameters;
	using Planner::Travel;
	using Planner::parse_alliances;
	using Planner::read_parameters;
	using Planner::parse_flights;

	UniqueId<> uniqueId;

	timeMe("start");
	//Declare variables and read the args
	Parameters parameters;
	Alliances alliances;
	read_parameters(uniqueId, parameters, argc, argv);
	omp_set_num_threads(parameters.nb_threads);
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

	std::vector<Travel> homeToConf, confToHome; /* cache some results from workHard for playHard */
	Travel workHardTravel = workHard(alliances, flights, parameters, homeToConf, confToHome);
	outputWorkHard(uniqueId, parameters, workHardTravel);
	timeMe("work hard");
	
	auto travels = playHard(alliances, flights, parameters, homeToConf, confToHome);
	outputPlayHard(uniqueId, parameters, travels);
	timeMe("play hard");

	return 0;
}

//./run -from Paris -to Los\ Ang/eles -departure_time_min 11152012000000 -departure_time_max 11172012000000 -arrival_time_min 11222012000000 -arrival_time_max 11252012000000 -max_layover 100000 -vacation_time_min 432000 -vacation_time_max 604800 -vacation_airports Rio London Chicago -flights flights.txt -alliances alliances.txt
