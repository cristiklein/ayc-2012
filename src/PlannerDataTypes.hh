#ifndef PLANNERDATATYPES_HH
#define PLANNERDATATYPES_HH

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

} /* namespace Planner */

#endif /* PLANNERDATATYPES_HH */
