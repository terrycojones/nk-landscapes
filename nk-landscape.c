#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/time.h>
#ifndef NeXT
#include <malloc.h>
#endif

#include "nk_landscape.h"

extern const char *const sys_errlist[];
extern int sys_nerr;
extern int errno;

#define nk_uniform(n) ((int) (nk_random() * (double)(n)))
#define WHY errno >= sys_nerr ? "Reason unknown" : sys_errlist[errno]
#define JUMP_MULTIPLIER 128

static double nk_random();
static long nk_seed();
static void error();
static char *Malloc();

static void
error(u, v, w, x, y, z)
char *u, *v, *w, *x, *y, *z;
{
    extern void perror();

    fprintf(stderr, "nk_landscape: ");
    fprintf(stderr, u, v, w, x, y, z);
    fprintf(stderr, "\n");
    if (fflush(stderr) == EOF){
        perror("fflush");
    }

    exit(1);
}

static char *
Malloc(n)
int n;
{
    char *s;

    if (!(s = malloc(n))){
	error("could not malloc %d chars. (%s)", n, WHY);
    }
    
    return s;
}

NK_LANDSCAPE *
nk_create(n, k, a, base_allele, neighborhood_type, seed, show_epistasis, fp)
int n;
int k;
int a;
int neighborhood_type;
long *seed;
int show_epistasis;
FILE *fp;
{
    static int **get_influencers();
    register NK_LANDSCAPE *new;
    register int locus; 
    register int i; 
    int njumps;

    if (n <= 0){
	error("nk_create() called with n <= 0 (n=%d)", n);
    }
    
    if (k >= n){
	error("nk_create() called with k >= n (k=%d, n=%d)", k, n);
    }
    
    if (k < 0){
	error("nk_create() called with k < 0 (k=%d)", k);
    }
    
    if (a < 2){
	error("nk_create() called with a < 2 (a=%d)", a);
    }
    
    if (base_allele < 0){
	error("nk_create() called with base_allele < 0 (base_allele=%d)", base_allele);
    }
    
    new = (NK_LANDSCAPE *) Malloc(sizeof(NK_LANDSCAPE));

    njumps = k * JUMP_MULTIPLIER;
    
    if (njumps < n){
	njumps = n;
    }

    new->n           = n;
    new->k           = k;
    new->a           = a;
    new->base_allele = base_allele;
    new->njumps      = njumps;
    new->jump_table  = (NK_JUMP *) Malloc(njumps * sizeof(NK_JUMP));
    
    *seed = nk_seed(*seed);

    for (i = 0; i < njumps; i++){

	int allele;
	
#ifdef NK_STATS
	/* Create and initialize space to store hit counts per allele. */
	new->jump_table[i].hits = (double *) Malloc(a * sizeof(double));
	for (allele = 0; allele < a; allele++){
	    new->jump_table[i].hits[allele] = 0;
	}
#endif
	/* Create and initialize space to store jump destination per allele. */
	new->jump_table[i].location = (int *) Malloc(a * sizeof(int));
	for (allele = 0; allele < a; allele++){
	    new->jump_table[i].location[allele] = nk_uniform(njumps);
	}

	/* Create and initialize space to store fitness contribution per allele. */
	new->jump_table[i].fitness = (double *) Malloc(a * sizeof(double));
	for (allele = 0; allele < a; allele++){
	    new->jump_table[i].fitness[allele] = nk_random();
	}
    }
    
    if (k){
	new->influencers = get_influencers(n, k, neighborhood_type);
    }
    
    if (show_epistasis){

	if (k){
	    fprintf(fp, "The landscape's epistatic interactions are as follows:\n");
	
	    for (locus = 0; locus < n; locus++){
		fprintf(fp, "Locus %d: ", locus);
	    
		for (i = 0; i < k; i++){
		    fprintf(fp, "%d ", new->influencers[locus][i]);
		} 
	    
		putc('\n', fp);
	    } 
	}
	else {
	    fprintf(fp, "There are no epistatic interactions (k = 0).\n");
	}
    }

    new->neighbors = Malloc(k + 2);
    new->neighbors[k + 1] = '\0';
    
    return new;
}

void
nk_free(landscape)
NK_LANDSCAPE *landscape;
{
    register int njumps = landscape->njumps;
    register int i;
    
    for (i = 0; i < njumps; i++){

#ifdef NK_STATS
	free(landscape->jump_table[i].hits);
#endif
	free(landscape->jump_table[i].location);
	free(landscape->jump_table[i].fitness);
    }

    free(landscape->jump_table);
    free(landscape->neighbors);
    free(landscape);
}


double
nk_fitness(point, landscape)
char *point;
NK_LANDSCAPE *landscape;
{
    double fitness = 0.0;
    register int locus;
    register int n = landscape->n;
    register int k = landscape->k;
    register char *neighbors = landscape->neighbors;
    register int *influencers_this_locus;
    
    for (locus = 0; locus < n; locus++){

	static double nk_locus_fitness();
	
	neighbors[0] = point[locus];

	if (k){
	    register int i;
	    influencers_this_locus = landscape->influencers[locus];

	    for (i = 0; i < k; i++){
		neighbors[i + 1] = point[influencers_this_locus[i]];
	    }
	}
    
	fitness += nk_locus_fitness(landscape, locus, neighbors);
    } 
    
    return fitness / (double) n;
}

static double
nk_locus_fitness(landscape, locus, neighbors)
NK_LANDSCAPE *landscape;
int locus;
char *neighbors;
{
    register int k = landscape->k;
    register int a = landscape->a;
    register int base_allele = landscape->base_allele;
    double fitness = 0.0;
    register NK_JUMP *jumps = landscape->jump_table;
    register int i;

    for (i = 0; i <= k; i++){
	register int index = neighbors[i] - base_allele;
	fitness += jumps[locus].fitness[index];
	locus = jumps[locus].location[index];
#ifdef NK_STATS
	jumps[locus].hits[index]++;
#endif
    } 
 
    return fitness - floor(fitness);
}


#define MBIG                 1000000000
#define FAC                  (1.0 / MBIG)

static int inext;
static int inextp;
static long ma[56];

static double
nk_random()
{
    long mj;
    
    if (++inext == 56){
	inext = 1;
    }
    
    if (++inextp == 56){
	inextp = 1;
    }

    mj = ma[inext] - ma[inextp];
    
    if (mj < 0L){
	mj += MBIG;
    }

    ma[inext] = mj;
    return mj * FAC;
}

static long
nk_seed(seed)
long seed;
{
    long mj;
    long mk;
    register int i;
    register int k;

    if (seed < 0){
	extern int gettimeofday();

	struct timeval tp;
	if (gettimeofday(&tp, (struct timezone *)0) == -1){
	    error("Could not gettimeofday in nk_seed().");
	}
	
	seed = tp.tv_sec;
    }
    
    if (seed >= MBIG){
	error("Seed value too big (> %d) in nk_seed().", MBIG);
    }

    ma[55] = mj = seed;
    mk = 1;
    
    for (i = 1; i <= 54; i++){
	register int ii = (21 * i) % 55;
	ma[ii] = mk;
	mk = mj - mk;
	if (mk < 0){
	    mk += MBIG;
	}
	mj = ma[ii];
    }
    
    for (k = 0; k < 4; k++){
	for (i = 1; i <= 55; i++){
	    ma[i] -= ma[1 + (i + 30) % 55];
	    if (ma[i] < 0){
		ma[i] += MBIG;
	    }
	}
    }
    
    inext = 0;
    inextp = 31;
    
    return seed;
}

static int **
get_influencers(n, k, neighborhood_type)
register int n;
register int k;
int neighborhood_type;
{
    register int locus;
    static void pick_influencers();
    int **influencers = (int **) Malloc(n * sizeof(int *));

    for (locus = 0; locus < n; locus++){
	/* Get space for k others. */
	influencers[locus] = (int *)Malloc(k * sizeof(int)); 
	    
	/* Choose the k influencing, put into influencers[locus]. */
	pick_influencers(locus, n, k, neighborhood_type, influencers[locus]);
    }

    return influencers;
}

static void
pick_influencers(locus, n, k, neighborhood_type, influencers)
register int locus;
register int n;
register int k;
int neighborhood_type;
int *influencers;
{
    static char *picked = (char *)0;
    register int i;
    register int picked_index = 0;
    
    /* Small optimization if we know that all loci except locus must be influencers. */
    if (k == n - 1){
	for (i = 0; i < n; i++){
	    if (i != locus){
		influencers[picked_index] = i;
		picked_index++;
	    }
	}
	
	return;
    }

    switch (neighborhood_type){
	case RANDOM_NEIGHBORS:{
	    /*
	     * Knuth gives a better algorithm for this.
	     * it's simple, linear, and doesn't need malloc...
	     * I didn't want to put it in now, since I don't
	     * have time to check it properly for you at the
	     * moment. This works, but is slower than need be,
	     * especially as k approaches n - 1.
	     */
	     
	    int npicked = 0;
	    
	    if (picked == (char *)0){
		/* This is never freed, but that's no big deal. */
		picked = Malloc(n);
	    }
	    
	    for (i = 0; i < n; i++){
		picked[i] = 0;
	    }
	    
	    while (npicked < k) {
		int new = nk_uniform(n);
		
		if (new != locus && !picked[new]){
		    picked[new] = 1;
		    npicked++;
		    influencers[picked_index] = new;
		    picked_index++;
		}
	    }
	    
	    break;
	}

	case NEXTDOOR_NEIGHBORS:{
	    register int nleft;
	    register int nright;

	    nleft = nright = k >> 1;
	    
	    /* Put the extra one left or right, at random. */
	    if (k % 2){
		if (!nk_uniform(2)){
		    nright++;
		}
		else {
		    nleft++;
		}
	    }

	    /* Get the loci to the left, watch for the edge of the world. */
	    if (locus >= nleft){
		register int i; 
		for (i = locus - nleft; i < locus; i++){
		    influencers[picked_index] = i;
		    picked_index++;
		} 
	    }
	    else {
		for (i = 0; i < locus; i++){
		    influencers[picked_index] = i;
		    picked_index++;
		} 
		
		for (i = n - (nleft - locus); i < n; i++){
		    influencers[picked_index] = i;
		    picked_index++;
		}
	    }
	    
	    /* Get the loci to the right, watch for the edge of the world. */
	    if (locus + nright < n){
		for (i = locus + 1; i <= locus + nright; i++){
		    influencers[picked_index] = i;
		    picked_index++;
		}
	    }
	    else {
		for (i = locus + 1; i < n; i++){
		    influencers[picked_index] = i;
		    picked_index++;
		}
		
		for (i = 0; i <= locus + nright - n; i++){
		    influencers[picked_index] = i;
		    picked_index++;
		}
	    }

	    break;
	}

	default :{
	    error("unknown neighborhood_type found in pick_influencers.");
	}
    }
    
    return;
}


#ifdef NK_STATS
void
nk_stats(landscape, fp)
NK_LANDSCAPE *landscape;
FILE *fp;
{
    register int i;
    register int base_allele = landscape->base_allele;
    register int njumps = landscape->njumps;
    register int a = landscape->a;
    int width = 1 + (int) log10((double) njumps);
    
    printf("Jump table statistics:\n\
\n\
Format is one line per jump table entry. The line shows a [...] area for\n\
each possible allele. The [...] contains:\n\
[allele, destination, fitness, hits].\n\
\n\
The hits count is the number of times this particular location\n\
in the jump table was passed through in the course of evaluating\n\
fitnesses.\n\
\n\
I think this code exists because I was doing testing to see how usage\n\
of the jump table varied with size, and to see that things appeared to\n\
be happening fairly (with randomly created individuals). If you create\n\
enough random indidividuals and evaluate their fitnesses (once each),\n\
you should see an even distribution in the hits in the jump table. If\n\
you don't, make sure you did create your individuals randomly.\n\
\n\
Terry Jones (terry@cliffs.ucsd.edu).\n\
Sept 7, 1998.\n\
\n");

for (i = 0; i < njumps; i++){

	register int allele;
	
	fprintf(fp, "Entry %*d", width, i);

	for (allele = 0; allele < a; allele++){
	    fprintf(fp, " [%c, %*d, %f, %3d]",
		    base_allele + allele,
		    width, landscape->jump_table[i].location[allele],
		    landscape->jump_table[i].fitness[allele],
		    landscape->jump_table[i].hits[allele]);
	}

	fprintf(fp, "\n");
    } 
    
    return;
}
#endif
