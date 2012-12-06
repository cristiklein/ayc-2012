#include "UniqueId.hh"

#include <cassert>
#include <cstdio>

#include "TimeMe.hh"

int main()
{
	UniqueId<> ids;

	timeMe("start");
	for (int i = 0; i < 1000000; i++) {
		assert(ids.getId("Hello") == ids.getId("Hello"));
		assert(ids.getId("Hello") != ids.getId("Hello2"));
	}
	timeMe("end");

	return 0;
}
