#include "Alliances.hh"

#include <cstdio>
#include <cstdlib>

#include "TimeMe.hh"

int main()
{
	timeMe("start");

	Alliances<> alliances;
	timeMe("initialization");

	for (int i = 0; i < 1000000; i++)
	{
		int id1 = rand();
		int id2 = rand();
		alliances.add(id1, id2);
	}
	timeMe("insertions");

	for (int i = 0; i < 1000000; i++)
	{
		int id1 = rand();
		int id2 = rand();
		alliances.areAllied(id1, id2);
	}

	timeMe("retrievals");
}
