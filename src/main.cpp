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
		Id airport; //!< airport in which we are now
		Time landingTime; //!< time from which we can take next plan
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
		float discountInC = getDiscount(alliances, *travelAB.flights.back(), *travelBC.flights.front());
		float prevDiscountBC = travelBC.discounts.back();
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

Travel findCheapestAndMerge(const Alliances &alliances, const vector<Travel> &travelsAB, const vector<Travel> &travelsBC, const vector<Travel> &travelsCD)
{
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
				if (newCost < bestCost) {
					bestCost = newCost;
					bestTravelAB = &travelAB;
					bestTravelBC = &travelBC;
					bestTravelCD = &travelCD;
				}
			}

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

Travel workHard(const Alliances& alliances, const Flights& flights, const Parameters& parameters)
{
	/* Get most expensive flights from/to conference */
	float maxToConf   = priciestFlight(flights.landings(parameters.to, parameters.dep_time_min, parameters.dep_time_max));
	float maxFromConf = priciestFlight(flights.takeoffs(parameters.to, parameters.ar_time_min, parameters.ar_time_max));

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

struct ParametersABCDTravel {
	Id a, b, c, d;
	Time tminAB, tmaxAB, tminBC, tmaxBC, tminCD, tmaxCD;
	Time maxLayover;
};

Travel computeBestABCDTravel(const Alliances &alliances, const Flights &flights, const ParametersABCDTravel &p)
{
	/* Get most expensive flights from/to vacation (for pruning) */
	float maxToB   = priciestFlight(flights.landings(p.b, p.tminAB, p.tmaxAB));
	float maxFromB = priciestFlight(flights.takeoffs(p.b, p.tminBC, p.tmaxBC));
	float maxToC   = priciestFlight(flights.landings(p.c, p.tminBC, p.tmaxBC));
	float maxFromC = priciestFlight(flights.takeoffs(p.c, p.tminCD, p.tmaxCD));

	vector<Travel> aToB = computePath(
		alliances, flights, /* description about the world */
		p.a, p.b, /* source, destination airport */
		p.tminAB, p.tmaxAB, /* interval of time during which to fly */
		p.maxLayover, /* other trip parameters */
		(maxToB + maxFromB) * 0.3 /* pruning parameter */);
	vector<Travel> bToC = computePath(
		alliances, flights, /* description about the world */
		p.b, p.c, /* source, destination airport */
		p.tminBC, p.tmaxBC, /* interval of time during which to fly */
		p.maxLayover, /* other trip parameters */
		(maxToB + maxFromB + maxToC + maxFromC) * 0.3 /* pruning parameter */);
	vector<Travel> cToD = computePath(
		alliances, flights, /* description about the world */
		p.c, p.d, /* source, destination airport */
		p.tminCD, p.tmaxCD, /* interval of time during which to fly */
		p.maxLayover, /* other trip parameters */
		(maxToC + maxFromC) * 0.3 /* pruning parameter */);

	return findCheapestAndMerge(alliances, aToB, bToC, cToD);
}

vector<Travel> playHard(const Alliances &alliances, const Flights& flights, Parameters& parameters)
{
	vector<Travel> results;
	for (Id vacation : parameters.airports_of_interest) {
		/*
		 * The first part compute a travel from home -> vacation -> conference -> home
		 * We'll use the terminology A -> B -> C -> D and AB BC CD for travels
		 */
		ParametersABCDTravel p1;

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

		Travel best1 = computeBestABCDTravel(alliances, flights, p1);

		/*
		 * The second part compute a travel from home -> conference -> vacation -> home
		 */
		ParametersABCDTravel p2;

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

		Travel best2 = computeBestABCDTravel(alliances, flights, p2);

		/*
		 * Compare the two solutions
		 */
		if (best1.totalCost > best2.totalCost)
			results.push_back(best2);
		else
			results.push_back(best1);
	}
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
	output << "Price : " << travel.totalCost << endl;

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
	vector<Travel> travels = playHard(alliances, flights, parameters);
	list<Id> cities = parameters.airports_of_interest;
	for(unsigned int i=0; i<travels.size(); i++){
		output<<"“Play Hard” Proposition "<<(i+1)<<" : "<<uniqueId.getName(cities.front())<<endl;
		printTravel(uniqueId, travels[i], output);
		cities.pop_front();
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

	Travel workHardTravel = workHard(alliances, flights, parameters);
	outputWorkHard(uniqueId, parameters, workHardTravel);
	timeMe("work hard");

	return 0;
}

//./run -from Paris -to Los\ Ang/eles -departure_time_min 11152012000000 -departure_time_max 11172012000000 -arrival_time_min 11222012000000 -arrival_time_max 11252012000000 -max_layover 100000 -vacation_time_min 432000 -vacation_time_max 604800 -vacation_airports Rio London Chicago -flights flights.txt -alliances alliances.txt
