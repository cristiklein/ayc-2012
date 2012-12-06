#ifndef FLIGHTS_HH
#define FLIGHTS_HH

#include <map>
#include <unordered_map>

/*! \brief Stores information about a flight */
template<
	class _Id = int,
	class _Time = unsigned long,
	class _Airport = _Id,
	class _Company = _Id,
	class _Cost = float
>
struct Flight {
	typedef _Id Id; //!< Type of flight identifier
	typedef _Time Time; //!< Type of takeoff and landing time
	typedef _Airport Airport; //!< Type of airport identifier
	typedef _Company Company; //!< Type of company identifier
	typedef _Cost Cost; //!< Type of cost

	Id id; //!< Unique id of the flight
	Airport from; //!< Airport where the flight takes off
	Airport to; //!< Airport where the flight lands
	Time takeoffTime; //!< Landing time
	Time landTime; //!< Takeoff time
	Company company; //!< Company of the flight
	Cost cost; //!< Cost (without any discounts) of the flight
};

/*! \brief Stores flight schedules for the whole world
 *
 * It allows fast retrieval of flights from / to an airport
 * taking off or landing in a given time interval.
*/
template<class _Flight = Flight<>>
class Flights {
public:
	typedef _Flight Flight; //!< Type of flight
	typedef typename Flight::Airport Airport; //!< Type of airport
	typedef typename Flight::Time Time; //!< Type of time

	typedef std::multimap<Time, Flight> AirportSchedule; //!< Type of index for fast retrieval of flights given a time interval
	typedef std::unordered_map<Airport, AirportSchedule> Index; //!< Type of index for fast retrievel of flights given an airport

	//! Type of iterator pointing to one flight
	class Iterator {
	public:
		//! Indirection operator
		const Flight &operator*() const { return it->second; }
		//! Forward operator
		void operator++() { it++; }
		//! Non-equality operator
		bool operator!=(const Iterator &other) const { return this->it != other->id; }
		//! Equality operator
		bool operator==(const Iterator &other) const { return this->it == other->id; }
	private:
		typename AirportSchedule::iterator it; //!< Underlying iterator
	};
	typedef std::pair<Iterator, Iterator> Range; //!< Type of range of flights

	/*!
	 * \brief Add a flight to the world
	 */
	void addFlight(const Flight &&flight)
	{
		takeoffsSchedule[flight.from][flight.takeoffTime] = flight;
		landingsSchedule[flight.to][flight.landingTime  ] = flight;
	}

	/*!
	 * \brief Get all flight from an airport taking off between tMin and tMax
	 */
	Range takeoffs(Airport airport, Time tMin, Time tMax) const
	{
		const auto &takeoffsFromAirport = takeoffsSchedule[airport];
		auto lower = takeoffsFromAirport.lower_bound(tMin);
		auto upper = takeoffsFromAirport.upper_bound(tMax);
		return { lower, upper };
	}

	/*!
	 * \brief Get all flight from an airport taking off between tMin and tMax
	 */
	Range landings(Airport airport, Time tMin, Time tMax) const
	{
		const auto &landingsOnAirport = landingsSchedule[airport];
		auto lower = landingsOnAirport.lower_bound(tMin);
		auto upper = landingsOnAirport.upper_bound(tMax);
		return { lower, upper };
	}

private:
	Index takeoffsSchedule; //!< Flights index by airport, then by takeoff time
	Index landingsSchedule; //!< Flights index by airport, then by landing time
};

namespace std {

Flights<>::Iterator begin(const Flights<>::Range &range) { return range.first; }
Flights<>::Iterator end(const Flights<>::Range &range) { return range.second; }

};

#endif /* FLIGHTS_HH */
