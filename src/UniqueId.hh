#ifndef UNIQUEID_HH
#define UNIQUEID_HH

#include <string>
#include <unordered_map>
#include <vector>

template<class Id = int, class Name = std::string>
class UniqueId
{
	std::unordered_map<Name, Id> nameToId;
	std::vector<Name> idToName;

	Id getId(Name s)
	{
		/* Check if Id is already registered */
		auto it = stringToId.find(s);
		if (it != std::end(stringToId))
			return it->second;

		/* Not found, add */
		Id newId = idToString.size();
		idToString.push_back(s);
		stringToId[s] = newId;
		return newId;
	}

	Name getString(Id id)
	{
		return idToString[id];
	}
};

#endif /* UNIQUE_ID */
