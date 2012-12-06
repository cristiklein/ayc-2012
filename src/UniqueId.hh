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

	UniqueId()
	{
		/* Id zero is special */
		idToName.emplace_back("");
	}

	Id getId(const Name &s)
	{
		Id &id = nameToId[s];

		/* Not found, add */
		if (id == 0)
		{
			Id newId = idToName.size();
			idToName.emplace_back(s);
			id = newId;
		}
		return id;
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
