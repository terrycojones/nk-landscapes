#define RANDOM_NEIGHBORS     0
#define NEXTDOOR_NEIGHBORS   1
#define DEF_NEIGHBORHOOD     RANDOM_NEIGHBORS

typedef struct {
    int *location;
    double *fitness;
#ifdef NK_STATS
    int *hits;
#endif
} NK_JUMP;

typedef struct {
    int n;
    int k;
    int a;
    int base_allele;
    int njumps;
    NK_JUMP *jump_table;
    int **influencers;
    char *neighbors;
} NK_LANDSCAPE;

extern NK_LANDSCAPE *nk_create();
extern double nk_fitness();

#ifdef NK_STATS
extern void nk_stats();
#endif
