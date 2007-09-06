/*	$Id$	*/

#define MAX(_a, _b) ((_a) > (_b) ? (_a) : (_b))
#define MIN(_a, _b) ((_a) < (_b) ? (_a) : (_b))

uint64_t	rdtsc(void);
double		cpuhz(void);
u_int		strchrcnt(const char *, const char);
void	       *xmalloc(size_t);
char	       *xstrdup(const char *);
unsigned int	hash_string(const char *);
double		factorial(u_int);
double		ls_choose(int64_t, int64_t);
