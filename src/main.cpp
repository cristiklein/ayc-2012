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

#include "UniqueId.hh"

#define TIMEME \
	timeme(__FILE__, __LINE__, __func__);

void timeme(const char *file, int line, const char *func)
{
	static struct timeval last;
	struct timeval now, diff;

	gettimeofday(&now, NULL);
	if (last.tv_sec != 0)
		timersub(&now, &last, &diff);
	else
		diff = { 0, 0 };
	last = now;

	fprintf(stderr, "+%2ld.%06ld %s:%d %s\n", diff.tv_sec, diff.tv_usec, file, line, func);
}

using namespace std;

UniqueId<> g_ids;
typedef UniqueId<>::Id Id;
typedef UniqueId<>::Name Name;
Id inline getId(const Name &s) { return g_ids.getId(s); }
Name inline getName(Id id) { return g_ids.getName(id); }

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

	float highestCost;/*!< The cost of the most expensive flight. Used for pruning */
};

/**
 * \struct Flight
 * \brief Store a single flight data.
 *This structure don't need to be modified but feel free to change it if you want.
 */
struct Flight{
	Id id;/*!< Unique id of the flight. */
	Id from;/*!< City where you take off. */
	Id to;/*!< City where you land. */
	unsigned long take_off_time;/*!< Take off time (epoch). */
	unsigned long land_time;/*!< Land time (epoch). */
	Id company;/*!< The company's name. */
	float cost;/*!< The cost of the flight. */
};

typedef unordered_map<Id, multimap<unsigned long, Flight>> Flights; /* ordered by take_off_time */

/**
 * \struct Travel
 * \brief Store a travel.
 * This structure don't need to be modified but feel free to change it if you want.
 */
struct Travel{
	vector<const Flight *> flights;/*!< A travel is just a list of Flight(s). */
};

/**
 * \struct Alliances
 * \brief Store whether two companies are in the same alliance
 */
struct Alliances{
	set<pair<Id, Id>> areAllied;
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
bool inline company_are_in_a_common_alliance(Id c1, Id c2, const Alliances& alliances);
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
void merge_path(vector<Travel>& travel1, vector<Travel>& travel2);
Travel work_hard(const Flights& flights, Parameters& parameters, const Alliances& alliances);
vector<Travel> play_hard(const Flights& flights, Parameters& parameters, const Alliances& alliances);
void output_play_hard(const Flights& flights, Parameters& parameters, const Alliances& alliances);
void output_work_hard(const Flights& flights, Parameters& parameters, const Alliances& alliances);

/**
 * \fn Travel work_hard(vector<Flight>& flights, Parameters& parameters, const Alliances& alliances)
 * \brief Solve the "Work Hard" problem.
 * This problem can be considered as the easy one. The goal is to find the cheapest way to join a point B from a point A regarding some parameters.
 * \param flights The list of available flights.
 * \param parameters The parameters.
 * \param alliances The alliances between companies.
 * \return The cheapest trip found.
 */
Travel work_hard(Flights& flights, Parameters& parameters, const Alliances& alliances){
	vector<Travel> travels;
	//First, we need to create as much travels as it as the number of flights that take off from the
	//first city
	fill_travel(travels, flights, parameters.from, parameters.dep_time_min, parameters.dep_time_max);
	compute_path(flights, parameters.to, travels, parameters.dep_time_min, parameters.dep_time_max, parameters, alliances);
	vector<Travel> travels_back;
	//Then we need to travel back
	fill_travel(travels_back, flights, parameters.to, parameters.ar_time_min, parameters.ar_time_max);
	compute_path(flights, parameters.from, travels_back, parameters.ar_time_min, parameters.ar_time_max, parameters, alliances);

	return find_cheapest(travels, travels_back, alliances);
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
		if (has_just_traveled_with_company(flight_before, current_flight)){
			discounts[i-1] = 0.7;
			discounts[i  ] = 0.7;
		}else if(has_just_traveled_with_alliance(flight_before, current_flight, alliances)){
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
 * \fn void compute_path(vector<Flight>& flights, Id to, vector<Travel>& travels, unsigned long t_min, unsigned long t_max, Parameters parameters)
 * \brief Computes a path from a point A to a point B. The flights must be scheduled between t_min and t_max. It is also important to take the layover in consideration.
 * 	You should try to improve and parallelize this function. A lot of stuff can probably done here. You can do almost what you want but the program output must not be modified.
 * \param flights All the flights that are available.
 * \param to The destination.
 * \param travels The list of possible travels that we are building.
 * \param t_min You must not be in a plane before this value (epoch)
 * \param t_max You must not be in a plane after this value (epoch)
 * \param parameters The program parameters
 */
void compute_path(const Flights& flights, Id to, vector<Travel>& final_travels, unsigned long t_min, unsigned long t_max, const Parameters &parameters, const Alliances &alliances){
	struct Segment {
		const struct Segment *prev;
		const Flight *flight;
		float prevTotalCost;
		float totalCost;
		float discount;
	};

	vector<const Segment *> allSegments; /* our way of doing malloc() and cleaning up */

	priority_queue<pair<float, const Segment *>> queue;

	for (const Travel &travel : final_travels) {
		Segment *s = new Segment();
		s->prev = NULL;
		s->flight = travel.flights.front();
		s->prevTotalCost = 0;
		s->totalCost = s->flight->cost;
		s->discount = 1;
		allSegments.push_back(s);
		queue.push(pair<float, const Segment *>(-s->totalCost, s));
	}

	float cheapest = INFINITY;
	vector<const Segment *> finalSegments;
	while (!queue.empty()) {
		const Segment *currentSegment = queue.top().second;
		queue.pop();

		const Flight &currentFlight = *currentSegment->flight;

		/* Pruning */
		if (currentSegment->totalCost - 0.3 * currentFlight.cost - 0.3 * parameters.highestCost > cheapest)
			break;

		if (currentFlight.to == to) {
			finalSegments.push_back(currentSegment);
			cheapest = min(cheapest, currentSegment->totalCost);
			continue;
		}

		auto itFlightsFromCurrentCity = flights.find(currentFlight.to);
		if (itFlightsFromCurrentCity == flights.end())
			continue;
		const multimap<unsigned long, Flight> &flightFromCurrentCity = itFlightsFromCurrentCity->second;
		auto itlo = flightFromCurrentCity.lower_bound(currentFlight.land_time);
		auto itup = flightFromCurrentCity.upper_bound(currentFlight.land_time + parameters.max_layover_time);
		for(auto it = itlo; it != itup; it++){
			const Flight &flight = it->second;
			if (flight.land_time > t_max) continue;

			bool wasHereBefore = false;
			for (const Segment *s = currentSegment; s != NULL; s = s->prev)
				if (s->flight->from == flight.to)
				{
					wasHereBefore = true;
					break;
				}
			if (wasHereBefore) continue;

			Segment *s = new Segment();
			s->prev = currentSegment;
			s->flight = &flight;

			/* Compute cost */
			float discount = 1;
			if (flight.company == currentFlight.company)
				discount = 0.7;
			else if (company_are_in_a_common_alliance(flight.company, currentFlight.company, alliances))
				discount = 0.8;

			s->prevTotalCost = currentSegment->prevTotalCost + currentFlight.cost * min(discount, currentSegment->discount);
			s->totalCost = s->prevTotalCost + flight.cost * discount;
			s->discount = discount;

			allSegments.push_back(s);
			queue.push(make_pair<float, const Segment *>(-s->totalCost, s));
		}
	}

	final_travels.clear();
	for (const Segment *s : finalSegments) {
		Travel travel;
		const Segment *currentSegment = s;
		while (currentSegment) {
			travel.flights.push_back(currentSegment->flight);
			currentSegment = currentSegment->prev;
		}
		reverse(travel.flights.begin(), travel.flights.end());
		final_travels.push_back(travel);
		float cost = compute_cost(travel, alliances);
		if (cost != s->totalCost)
		{
			fprintf(stderr, "Cost error: compute_cost %f, s->cost %f\n", cost, s->totalCost);
			fprintf(stderr, "Cost error: compute_cost %X, s->cost %X\n", *(int*)&cost, *(int*)&s->totalCost);
			print_travel(travel, alliances, cerr);
			assert(0);
		}
	}

	for (const Segment *s : allSegments) {
		delete s;
	}
	allSegments.clear();
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
void fill_travel(vector<Travel>& travels, const Flights& flights, Id starting_point, unsigned long t_min, unsigned long t_max){
	auto itFlightsFromCurrentCity = flights.find(starting_point);
	if (itFlightsFromCurrentCity == flights.end())
		return;
	const multimap<unsigned long, Flight> &flightFromCurrentCity = itFlightsFromCurrentCity->second;
	auto itlo = flightFromCurrentCity.lower_bound(t_min);
	auto itup = flightFromCurrentCity.end();
	for(auto it = itlo; it != itup; it++){
		const Flight &flight = it->second;
		if(flight.land_time <= t_max){
			Travel t;
			t.flights.push_back(&flight);
			travels.push_back(t);
		}
	}
}

/**
 * \fn void merge_path(vector<Travel>& travel1, vector<Travel>& travel2)
 * \brief Merge the travel1 with the travel2 and put the result in the travel1.
 * \param travel1 The first part of the trip.
 * \param travel2 The second part of the trip.
 */
void merge_path(vector<Travel>& travel1, vector<Travel>& travel2){
	vector<Travel> result;
	for(unsigned int i=0; i<travel1.size(); i++){
		Travel t1 = travel1[i];
		for(unsigned j=0; j<travel2.size(); j++){
			Travel t2 = travel2[j];
			const Flight &last_flight_t1 = *t1.flights.back();
			const Flight &first_flight_t2 = *t2.flights[0];
			if(last_flight_t1.land_time < first_flight_t2.take_off_time){
				Travel new_travel = t1;
				new_travel.flights.insert(new_travel.flights.end(), t2.flights.begin(), t2.flights.end());
				result.push_back(new_travel);
			}
		}
	}
	travel1 = result;
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
void print_flight(const Flight& flight, float discount, ostream& output){
	struct tm * take_off_t, *land_t;
	take_off_t = gmtime(((const time_t*)&(flight.take_off_time)));
	output<<getName(flight.company)<<"-";
	output<<""<<getName(flight.id)<<"-";
	output<<getName(flight.from)<<" ("<<(take_off_t->tm_mon+1)<<"/"<<take_off_t->tm_mday<<" "<<take_off_t->tm_hour<<"h"<<take_off_t->tm_min<<"min"<<")"<<"/";
	land_t = gmtime(((const time_t*)&(flight.land_time)));
	output<<getName(flight.to)<<" ("<<(land_t->tm_mon+1)<<"/"<<land_t->tm_mday<<" "<<land_t->tm_hour<<"h"<<land_t->tm_min<<"min"<<")-";
	output<<flight.cost<<"$"<<"-"<<discount*100<<"%"<<endl;

}

/**
 * \fn void read_parameters(Parameters& parameters, int argc, char **argv)
 * \brief This function is used to read the parameters
 * \param parameters Represents the structure that will be filled with the parameters.
 */
void read_parameters(Parameters& parameters, int argc, char **argv){
	for(int i=0; i<argc; i++){
		string current_parameter = argv[i];
		if(current_parameter == "-from"){
			parameters.from = getId(argv[++i]);
		}else if(current_parameter == "-arrival_time_min"){
			parameters.ar_time_min = convert_string_to_timestamp(argv[++i]);
		}else if(current_parameter == "-arrival_time_max"){
			parameters.ar_time_max = convert_string_to_timestamp(argv[++i]);
		}else if(current_parameter == "-to"){
			parameters.to = getId(argv[++i]);
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
				parameters.airports_of_interest.push_back(getId(argv[++i]));
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
void parse_flight(Flight& flight, string& line){
	vector<string> splittedLine;
	split_string(splittedLine, line, ';');
	if(splittedLine.size() == 7){
		flight.id = getId(splittedLine[0]);
		flight.from = getId(splittedLine[1]);
		flight.take_off_time = convert_string_to_timestamp(splittedLine[2]);
		flight.to = getId(splittedLine[3]);
		flight.land_time = convert_string_to_timestamp(splittedLine[4]);
		flight.cost = atof(splittedLine[5].c_str());
		flight.company = getId(splittedLine[6]);
	}
}

/**
 * \fn void parse_flights(vector<Flight>& flights, string filename)
 * \brief This function parses the flights from a file.
 * \param flights The vector of flights.
 * \param filename The name of the file containing the flights.
 */
float parse_flights(Flights& flights, string filename){
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
		parse_flight(flight, line);
		flights[flight.from].insert({{flight.take_off_time, flight}});
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
void parse_alliance(Alliances &alliances, string line){
	vector<string> splittedLine;
	split_string(splittedLine, line, ';');
	for(unsigned int i=0; i<splittedLine.size(); i++)
	for(unsigned int j=i+1; j<splittedLine.size(); j++){
		Id c1 = getId(splittedLine[i]);
		Id c2 = getId(splittedLine[j]);
		if (c1 == c2)
			continue;

		alliances.areAllied.insert(pair<Id,Id>(c1,c2));
		alliances.areAllied.insert(pair<Id,Id>(c2,c1));
	}
}

/**
 * \fn void parse_alliances(Alliances &alliances, string filename)
 * \brief This function parses a line containing alliances between companies.
 * \param alliances A 2D vector representing the alliances. Companies on the same line are in the same alliance.
 * \param filename The name of the file containing the alliances description.
 */
void parse_alliances(Alliances &alliances, string filename){
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
		parse_alliance(alliances, line);
	}
}

/**
 * \fn bool company_are_in_a_common_alliance(Id c1, Id c2, const Alliances& alliances)
 * \brief Check if 2 companies are in the same alliance.
 * \param c1 The first company's name.
 * \param c2 The second company's name.
 * \param alliances A 2D vector representing the alliances. Companies on the same line are in the same alliance.
 */
bool company_are_in_a_common_alliance(Id c1, Id c2, const Alliances& alliances){
	return alliances.areAllied.count(pair<Id,Id>(c1,c2)) != 0;
}

/**
 * \fn bool has_just_traveled_with_company(vector<Flight>& flights_before, Flight& current_flight)
 * \brief The 2 last flights are with the same company.
 * \param flight_before The first flight.
 * \param current_flight The second flight.
 * \return The 2 flights are with the same company
 */
bool has_just_traveled_with_company(const Flight& flight_before, const Flight& current_flight){
	return flight_before.company == current_flight.company;
}

/**
 * \fn bool has_just_traveled_with_alliance(Flight& flight_before, Flight& current_flight, vector<vector<string> >& alliances)
 * \brief The 2 last flights are with the same alliance.
 * \param flight_before The first flight.
 * \param current_flight The second flight.
 * \param alliances The alliances.
 * \return The 2 flights are with the same alliance.
 */
bool has_just_traveled_with_alliance(const Flight& flight_before, const Flight& current_flight, const Alliances& alliances){
	return company_are_in_a_common_alliance(current_flight.company, flight_before.company, alliances);
}

/**
 * \fn void print_alliances(vector<vector<string> > &alliances)
 * \brief Display the alliances on the standard output.
 * \param alliances The alliances.
 */
void print_alliances(const Alliances &alliances){
	cout<<"Alliance: NOT implemented"<<endl;
}

/**
 * \fn void print_flights(vector<Flight>& flights, ostream& output)
 * \brief Display the flights on the standard output.
 * \param flights The flights.
 */
void print_flights(const vector<const Flight *>& flights, const vector<float> &discounts, ostream& output){
	for(unsigned int i=0; i<flights.size(); i++)
		print_flight(*flights[i], discounts[i], output);
}

/**
 * \fn bool never_traveled_to(Travel travel, Id city)
 * \brief Indicates if the city has already been visited in the travel. This function is used to avoid stupid loops.
 * \param travel The travels.
 * \apram city The city.
 * \return The current travel has never visited the given city.
 */
bool never_traveled_to(Travel travel, Id city){
	for(unsigned int i=0; i<travel.flights.size(); i++)
		if(travel.flights[i]->from == city || travel.flights[i]->to == city)
			return false;
	return true;
}

/**
 * \fn void print_travel(Travel& travel, const Alliances&alliances)
 * \brief Display a travel on the standard output.
 * \param travel The travel.
 * \param alliances The alliances (used to compute the price).
 */
void print_travel(const Travel& travel, const Alliances&alliances, ostream& output){
	vector<float> discounts = apply_discount(travel, alliances);
	output<<"Price : "<<compute_cost(travel, alliances)<<endl;
	print_flights(travel.flights, discounts, output);
	output<<endl;
}

/**
 * \fn void output_play_hard(vector<Flight>& flights, Parameters& parameters, const Alliances& alliances)
 * \brief Display the solution of the "Play Hard" problem by solving it first.
 * \param flights The list of available flights.
 * \param parameters The parameters.
 * \param alliances The alliances between companies.
 */
void output_play_hard(Flights& flights, Parameters& parameters, const Alliances& alliances){
	ofstream output;
	output.open(parameters.play_hard_file.c_str());
	vector<Travel> travels = play_hard(flights, parameters, alliances);
	list<Id> cities = parameters.airports_of_interest;
	for(unsigned int i=0; i<travels.size(); i++){
		output<<"“Play Hard” Proposition "<<(i+1)<<" : "<<getName(cities.front())<<endl;
		print_travel(travels[i], alliances, output);
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
void output_work_hard(Flights& flights, Parameters& parameters, const Alliances& alliances){
	ofstream output;
	output.open(parameters.work_hard_file.c_str());
	Travel travel = work_hard(flights, parameters, alliances);
	output<<"“Work Hard” Proposition :"<<endl;
	print_travel(travel, alliances, output);
	output.close();
}

int main(int argc, char **argv) {
	TIMEME;
	//Declare variables and read the args
	Parameters parameters;
	Alliances alliances;
	read_parameters(parameters, argc, argv);
//	cout<<"Printing parameters..."<<endl;
//	print_params(parameters);
	Flights flights;
	TIMEME;
	parameters.highestCost = parse_flights(flights, parameters.flights_file);
	TIMEME;
//	cout<<"Printing flights..."<<endl;
//	print_flights(flights);
//	cout<<"flights printed "<<endl;
	parse_alliances(alliances, parameters.alliances_file);
//	cout<<"Printing alliances..."<<endl;
//	print_alliances(alliances);
	TIMEME;
	output_play_hard(flights, parameters, alliances);
	TIMEME;
	output_work_hard(flights, parameters, alliances);
	TIMEME;
}

//./run -from Paris -to Los\ Angeles -departure_time_min 11152012000000 -departure_time_max 11172012000000 -arrival_time_min 11222012000000 -arrival_time_max 11252012000000 -max_layover 100000 -vacation_time_min 432000 -vacation_time_max 604800 -vacation_airports Rio London Chicago -flights flights.txt -alliances alliances.txt
