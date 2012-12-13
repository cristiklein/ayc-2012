#ifndef COMPUTEPATHREVERSE_HH
#define COMPUTEPATHREVERSE_HH

#include "PlannerHelpers.hh"

namespace Planner {

using namespace std;

unordered_map<Id, vector<Travel>> computePath(
	const Alliances& alliances,
	const Flights& flights,
	set<Id> from, Id to,
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
	auto range = flights.landings(to,
		tMin,
		tMax);
	for(const Flight &newFlight : range) {
		if (newFlight.takeoffTime < tMin)
			continue;

		Segment *seg = new Segment();
		seg->prev = NULL;
		seg->flight = &newFlight;
		seg->prevTotalCost = 0;
		seg->totalCost = newFlight.cost;
		seg->discount = 1;
		queue.push(pair<float, const Segment *>(-seg->totalCost, seg));
	}

	unordered_map<Id, vector<const Segment *>> finalSegments;
	float lastCheapestTravel = INFINITY;
	while (!queue.empty()) {
		const Segment *seg = queue.top().second;
		queue.pop();

		/* Faster retrieval */
		Id airport = seg->flight->from;
		Time takeoffTime = seg->flight->takeoffTime;

		/* Pruning */
		if (seg->totalCost - maxDiscount > lastCheapestTravel)
			break;

		/* Check if we reached any destination */
		if (from.count(airport)) {
			finalSegments[airport].push_back(seg);
			if (finalSegments.size() == from.size())
				lastCheapestTravel = seg->totalCost;
		}

		/* Where can we go to from here? */
		auto range = flights.landings(airport,
			takeoffTime - maxLayover,
			takeoffTime);
		for(const Flight &newFlight : range) {
			/* Out of time? */
			if (newFlight.takeoffTime < tMin)
				continue;

			/* Avoid cycles */
			/* Implementation note: we tried accelerating this with an unordered_set,
			 * but we would actually observe a 2x slowdown */
			if (newFlight.from == to)
				continue;

			bool wasHereBefore = false;
			for (const Segment *s = seg; s != NULL && !wasHereBefore; s = s->prev)
				if (s->flight->from == newFlight.from)
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
	unordered_map<Id, vector<Travel>> travels;
	for (Id airport : from) {
		for (const Segment *seg : finalSegments[airport]) {
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
			travels[airport].push_back(travel);
		}
	}

	for (const Segment *s : allSegments)
		delete s;
	return travels;
}

} /* namespace Planner */

#endif /* COMPUTEPATHREVERSE_HH */
