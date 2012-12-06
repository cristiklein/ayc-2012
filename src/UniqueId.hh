#ifndef UNIQUEID_HH
#define UNIQUEID_HH

#include <string>
#include <unordered_map>
#include <vector>

/*! \brief Associate bulky names (e.g., strings)
 * with light, unique identifiers (e.g., integers)
 */
template<class _Id = int, class _Name = std::string>
class UniqueId
{
public:
	typedef _Id Id; //!< Type of identifiers
	typedef _Name Name; //!< Type of names

	UniqueId()
	{
		/* The identifier zero is special, meaning identifier not yet registered */
		idToName.emplace_back("");
	}

	/*!
	 * \brief Retrieve the identifier or generate a new identifier for the given name
	*/
	Id getId(const Name &s)
	{
		/* Retrieve identifier from map */
		Id &id = nameToId[s];

		/* If name was not found */
		if (id == 0)
		{
			/* Generate a new identifier for it */
			Id newId = idToName.size();
			idToName.emplace_back(s);
			id = newId;
		}

		return id;
	}

	/*!
	 * \brief Retrieve the identifier or generate a new identifier for the given name
	*/
	Name getName(Id id) const
	{
		return idToName.at(id);
	}

private:
	std::unordered_map<Name, Id> nameToId; //!< Hashmap to store name to identifier association
	std::vector<Name> idToName; //!< Vector to store identifier to name association
};

#endif /* UNIQUE_ID */
