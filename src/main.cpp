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

#include "ComputePathDirect.hh"
#include "ComputePathReverse.hh"
#include "PlannerMergers.hh"

namespace Planner {

vector<Travel> computePath(
	const Alliances& alliances,
	const Flights& flights,
	Id from, Id to,
	Time tMin, Time tMax,
	Time maxLayover,
	Cost maxDiscount)
{
	set<Id> destinationSet;
	destinationSet.insert(to);
	return computePath(alliances, flights,
		from, destinationSet,
		tMin, tMax,
		maxLayover,
		maxDiscount)[to];
}


Travel workHard(const Alliances& alliances, const Flights& flights, const Parameters& parameters, vector<Travel> &homeToConf, vector<Travel> &confToHome)
{
	/* Get most expensive flights from/to conference */
	float maxToConf   = priciestFlight(flights.landings(parameters.to, parameters.dep_time_min, parameters.dep_time_max));
	float maxFromConf = priciestFlight(flights.takeoffs(parameters.to, parameters.ar_time_min, parameters.ar_time_max));

	timeMe("start workHard");
	homeToConf = computePath(
		alliances, flights, /* description about the world */
		parameters.from, parameters.to, /* source, destination airport */
		parameters.dep_time_min, parameters.dep_time_max, /* interval of time during which to fly */
		parameters.max_layover_time, /* other trip parameters */
		(maxToConf + maxFromConf) * 0.3 /* pruning parameter */);
	timeMe("homeToConf");
	confToHome = computePath(
		alliances, flights,
		parameters.to, parameters.from,
		parameters.ar_time_min, parameters.ar_time_max,
		parameters.max_layover_time,
		(maxToConf + maxFromConf) * 0.3 /* pruning parameter */);
	timeMe("confToHome");

	return findCheapestAndMerge(alliances, homeToConf, confToHome);
}

map<Id, Travel> playHard(const Alliances &alliances, const Flights& flights, Parameters& parameters, const vector<Travel> &homeToConf, const vector<Travel> &confToHome)
{
	set<Id> vacations = set<Id>(parameters.airports_of_interest.begin(),
		parameters.airports_of_interest.end());

	unordered_map<Id, vector<Travel>>
		homeToVacations,
		vacationsToConf,
		confToVacations,
		vacationsToHome;

	/*
	 * The first part compute a travel from home -> vacation -> conference -> home
	 * We'll use the terminology A -> B -> C -> D and AB BC CD for travels
	 */
	{
		/* Cities */
		Id a = parameters.from;
		set<Id> bs = vacations;
		Id c = parameters.to;
		// Id d = parameters.from; /* unused, already computed in workHard */

		/* Compute valid travel times */
		Time tminAB = parameters.dep_time_min - parameters.vacation_time_max;
		Time tmaxAB = parameters.dep_time_min - parameters.vacation_time_min;
		Time tminBC = parameters.dep_time_min;
		Time tmaxBC = parameters.dep_time_max;
		Time tminCD = parameters.ar_time_min;
		Time tmaxCD = parameters.ar_time_max;
		Time maxLayover = parameters.max_layover_time;
	
		/* Get most expensive flights from/to vacation and conference(for pruning) */
		float maxToB   = 0;
		float maxFromB = 0;
		for (Id b : bs) {
			maxToB   = max(maxToB  , priciestFlight(flights.landings(b, tminAB, tmaxAB)));
			maxFromB = max(maxFromB, priciestFlight(flights.takeoffs(b, tminBC, tmaxBC)));
		}

		float maxToC   = priciestFlight(flights.landings(c, tminBC, tmaxBC));
		float maxFromC = priciestFlight(flights.takeoffs(c, tminCD, tmaxCD));

#pragma omp task untied shared(alliances, flights, homeToVacations)
		homeToVacations = std::move(computePath(
			alliances, flights, /* description about the world */
			a, bs, /* source, destination airport */
			tminAB, tmaxAB, /* interval of time during which to fly */
			maxLayover, /* other trip parameters */
			(maxToB + maxFromB) * 0.3 /* pruning parameter */));

#pragma omp task untied shared(alliances, flights, vacationsToConf)
		vacationsToConf = std::move(computePath(
			alliances, flights, /* description about the world */
			bs, c, /* source, destination airport */
			tminBC, tmaxBC, /* interval of time during which to fly */
			maxLayover, /* other trip parameters */
			(maxToB + maxFromB + maxToC + maxFromC) * 0.3 /* pruning parameter */));
	}
	
	/*
	 * The second part compute a travel from home -> conference -> vacation -> home
	 * We'll use the terminology A -> B -> C -> D and AB BC CD for travels
	 */
	{
		/* Cities */
		// Id a = parameters.from; /* unused, already computed in workHard */
		Id b = parameters.to;
		set<Id> cs = vacations;
		Id d = parameters.from;

		/* Compute valid travel times */
		Time tminAB = parameters.dep_time_min;
		Time tmaxAB = parameters.dep_time_max;
		Time tminBC = parameters.ar_time_min;
		Time tmaxBC = parameters.ar_time_max;
		Time tminCD = parameters.ar_time_max + parameters.vacation_time_min;
		Time tmaxCD = parameters.ar_time_max + parameters.vacation_time_max;
		Time maxLayover = parameters.max_layover_time;
	
		/* Get most expensive flights from/to vacation and conference(for pruning) */
		float maxToB   = priciestFlight(flights.landings(b, tminAB, tmaxAB));
		float maxFromB = priciestFlight(flights.takeoffs(b, tminBC, tmaxBC));

		float maxToC   = 0;
		float maxFromC = 0;
		for (Id c : cs) {
			maxToC   = max(maxToC, priciestFlight(flights.landings(c, tminBC, tmaxBC)));
			maxFromC = max(maxFromC, priciestFlight(flights.takeoffs(c, tminCD, tmaxCD)));
		}

#pragma omp task untied shared(alliances, flights, confToVacations)
		confToVacations = std::move(computePath(
			alliances, flights, /* description about the world */
			b, cs, /* source, destination airport */
			tminBC, tmaxBC, /* interval of time during which to fly */
			maxLayover, /* other trip parameters */
			(maxToB + maxFromB + maxToC + maxFromC) * 0.3 /* pruning parameter */));

#pragma omp task untied shared(alliances, flights, vacationsToHome)
		vacationsToHome = std::move(computePath(
			alliances, flights, /* description about the world */
			cs, d, /* source, destination airport */
			tminCD, tmaxCD, /* interval of time during which to fly */
			maxLayover, /* other trip parameters */
			(maxToC + maxFromC) * 0.3 /* pruning parameter */));

	}

#pragma omp taskwait
	timeMe("computePath");

	/*
	 * Run time-consuming findCheapest
	 */
	map<Id, Travel> results;
	for (Id vacation : vacations) {
		/* Make a copies, since findCheapestAndMerge is destructive */
		auto copyOfHomeToConf = homeToConf;
		auto copyOfConfToHome = confToHome;

		auto best1 = findCheapestAndMerge(alliances,
			homeToVacations[vacation],
			vacationsToConf[vacation],
			copyOfConfToHome);
		auto best2 = findCheapestAndMerge(alliances,
			copyOfHomeToConf,
			confToVacations[vacation],
			vacationsToHome[vacation]);

		if (best1.totalCost > best2.totalCost)
			results[vacation] = best2;
		else
			results[vacation] = best1;

		timeMe("findCheapest");
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
	parameters.nb_threads = 0;

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
	if (parameters.nb_threads > 0)
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
