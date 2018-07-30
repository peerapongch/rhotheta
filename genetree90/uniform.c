#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
/* Long period random number generator from
   Marsagalia, G. and Zaman, A.
	 "A new class of random number generators",
	 Annals of Applied Probability, 1, 462-480.
	 Possibilities for base b=2^24 generators are
	 r     s
	 39    25
	 28     8
	 25    11
	 24    10
	 These generators have VERY long cycles.
*/

#define b 16777216
#define xb 16777216.0
#define r 39
#define s 25
#define WARMUP 1000
long table[r];
long c__=0L;

double uni01(long *idum) {
	long z=1L,k;
	static int position=0;
	int seed,i,j,l;
	/**************  Initialize generator ************************/

	if(*idum < 0) {
	/* Fill up table with system random numbers, quick and dirty */
		seed= -(int)(*idum);
		srand(seed);
		for(i=0;i<r;i++) {
			table[i]=((long)rand())%b;
			if(table[i]<=0) table[i] *= (-1);
		}
		/* Make a random permutation of initial table values */
		for(i=0;i<r;i++) {
			l=i+rand()%(r-i);
			/* Swap table[l] and table [i] */
			k=table[l];table[l]=table[i];table[i]=k;
		}
		/* Warmup main generator */
		for(i=0;i<WARMUP;i++) uni01(&z);
	}

	/************  Main code for random numbers ****************/

	/* Position is where the rth most recently generated number is in
		the array. Its a "wrap around array".
	*/
	if(position-s<0) z=table[position+r-s]-table[position]-c__ ;
	else z=table[position-s]-table[position]-c__ ;
	if(z>=0) c__=0L;
	else { z += b;c__ = 1L;}
	table[position]=*idum=z;
	/* Next oldest is to the right */
	position++;
	if(position==r) position=0;
	/* Really return from 0.5/2^24 to 1.0 - 0.5/2^24 */
	return (z+0.5)/xb;
}

#undef b
#undef xb
#undef r
#undef s
#undef WARMUP
