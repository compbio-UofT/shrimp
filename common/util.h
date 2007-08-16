/*	$Id$	*/

#define MAX(_a, _b) ((_a) > (_b) ? (_a) : (_b))
#define MIN(_a, _b) ((_a) < (_b) ? (_a) : (_b))

uint64_t rdtsc(void);
double   cpuhz(void);
u_int	 strchrcnt(const char *, const char);
