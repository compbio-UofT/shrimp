/*	$Id$	*/

#define DEF_PCHANCE_CUTOFF	0.05
#define DEF_PGENOME_CUTOFF	0.0
#define DEF_NORMODDS_CUTOFF	0.0
#define DEF_TOP_MATCHES		10


/* Stats stuff */

double maxCount(int ins, int dels, int len);
double minCount(int ins, int dels, int len);
double subCount(int subs, int len);
double fact(double n);
double choose(int n, int m);
void initfastfact();
double fastfact(int n);
double fastchoose(int n, int m);
void initlookupchoose();
double lookupchoose(int n, int m) ;
