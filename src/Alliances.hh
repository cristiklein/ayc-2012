#ifndef ALLIANCES_HH
#define ALLIANCES_HH

#include <unordered_set>

/*! \brief Determine fast if two companies are allied */
template<class _Id = int, class _IdPair = long>
class Alliances {
public:
	typedef _Id Id; //!< Type of identifiers
	typedef _IdPair IdPair; //!< Type of identifier-pair

	/*! \brief Add a new alliance between id1 and id2 */
	void add(Id id1, Id id2)
	{
		IdsToIdPair idp;
		idp.id1 = id1;
		idp.id2 = id2;
		alliances.insert(idp.id);
		
		idp.id1 = id2;
		idp.id2 = id1;
		alliances.insert(idp.id);
	}

	/*! \brief Return whether id1 and id2 are allied */
	bool areAllied(Id id1, Id id2) const
	{
		IdsToIdPair idp;
		idp.id1 = id1;
		idp.id2 = id2;
		return alliances.count(idp.id);
	}

private:
	//! "Magic" type to do conversion from a pair of identifier to an identifier-pair
	typedef union {
		struct {
			Id id1, id2;
		}; //!< identifiers
		IdPair id; //!< identifier-pair
	} IdsToIdPair;

	std::unordered_set<IdPair> alliances; //!< Stores the set of identifier-pairs which are allied
};

#endif /* ALLIANCES_HH */
