#ifndef UNIQUEID_HH
#define UNIQUEID_HH

#include <string>
#include <unordered_map>
#include <vector>

template<class _Id = int, class _Name = std::string>
class UniqueId
{
public:
	typedef _Id Id;
	typedef _Name Name;

	Id getId(Name s)
	{
		/* Check if Id is already registered */
		auto it = nameToId.find(s);
		if (it != std::end(nameToId))
			return it->second;

		/* Not found, add */
		Id newId = idToName.size();
		idToName.push_back(s);
		nameToId[s] = newId;
		return newId;
	}

	Name getName(Id id)
	{
		return idToName[id];
	}

private:
	std::unordered_map<Name, Id> nameToId;
	std::vector<Name> idToName;
};

#endif /* UNIQUE_ID */
