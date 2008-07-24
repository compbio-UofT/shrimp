/*	$Id$	*/

#define DEF_PCHANCE_CUTOFF	0.05
#define DEF_PGENOME_CUTOFF	0.0
#define DEF_NORMODDS_CUTOFF	0.0
#define DEF_TOP_MATCHES		10


/* Stats stuff */

double maxCount(int ins, int dels, int len); /* max indel Z */
double minCount(int ins, int dels, int len); /* min indel Z */
double subCount(int subs, int len); /* subs Z */
double fastchoose(int n, int m); /* fast choose function: uses fastlchoose */
double fastlchoose(int n, int m); /* fast log choose function: uses lgamma */
void initStats(int maxlen); /* initiate statistics - build the lookup tables */
