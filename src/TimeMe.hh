#ifndef TIMEME_HH
#define TIMEME_HH

#include <sys/time.h>
#include <time.h>

void inline timeMe(const char *desc)
{
	static struct timeval last = { 0, 0 };
	struct timeval now, diff;

	gettimeofday(&now, NULL);
	if (last.tv_sec != 0)
		timersub(&now, &last, &diff);
	else
		diff = { 0, 0 };
	last = now;

	fprintf(stderr, "+%2ld.%06ld %s\n", diff.tv_sec, diff.tv_usec, desc);
}

#endif /* TIMEME_HH */
