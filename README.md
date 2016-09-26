# Efficient computation of NK landscapes

This uses something I call a *jump table* to allow the `k` parameter to be
very large without consuming much memory. I thought of this approach while
at the [Santa Fe Institute](http://www.santafe.edu/) in around 1994. I had
written several other implementations that were not efficient enough.

## Usage

In your code put

```c
    #include "nk-landscape.h"
```

and then you can do things like

```c
    int n = 10;
    int k = 7;
    int a = 2;
    long seed = -1L;

    /* Make a landscape. */
    NK_LANDSCAPE *landscape = nk_create(n, k, a, '0', RANDOM_NEIGHBORS, &seed, 1, stdout);

    /* Evaluate an individual. */
    double fitness = nk_fitness("1000101101", landscape);

    /* Free up the space associated with the landscape. */
    nk_free(landscape);
```

Which is about as simple as could be. Arguments are:

* The first two arguments to `nk_create` are (should be!) obvious.
* The next (`a`) indicates how many alleles you want (must be `>=2`). There
  is no default, so pass `2` if you want the normal NK landscape.
* I added this functionality on Sept 7, 1998 at the suggestion of
  Wim Hordijk (wim@santafe.edu).
* The next (`0` above) gives the base value you'll be putting in
  your genome. Thus if you pass `0` with `a = 2`, my code will use the
  ASCII characters '0' and '1' when trying to figure out what's in your
  genome. If you pass `0`, I'll use `0`, `1`, `2`, ..., `a`.
* The next (`RANDOM_NEIGHBORS` above) indicates how you want neighbors
  chosen (you have a choice between the immediate neighbors of each
  locus and neighbors chosen at random).
* The next (`&seed` above) is a pointer to a `long`. If you give a
  negative value, this will be filled in for you from the system
  clock. You can then print it etc.
* The next (`1` above) indicates whether you want to see the table of
  epistatic interactions (that I call influencers in the code)
  printed out. This will normally be `0` but if it is not, the table
  will be printed to the `FILE *` that is passed as the final
  argument.
* Finally, pass an already opened file descriptor. if you're not
  choosing to have the table of epistatic interactions printed, this
  can be `NULL`.

Errors are printed to `stderr`.

## Thanks

To Wim Hordijk (wim@santafe.edu) for suggestions and improvements, to Emily
Dickinson (sadly deceased), and to Bill Macready (wgm@santafe.edu) with
whom I spent time pondering how to efficiently implement NK landscapes with
large `n` and `k`.

## Notes

If you bother to work out how the code works, I'd be interested in your
reactions. I've thought about how to implement NK landscapes quite a lot. I
went through a few iterations before hitting on this 'jump table'
idea. This idea is the basis of the code used by the people who work for
Stu Kauffman now (1994). I did about a day's worth of fairly intensive
testing to see that it was random enough. If you want to do more, I'd like
to hear what you find. The random number generator in there is taken from
Knuth, and should be fine. It's fast and has been looked at a lot. I can
give you a reference to it if you like. It's in numerical recipes too, but
i think they have a small bug (and their C code is clearly written by a
fortran programmer :-)).

If you do figure out what I'm doing, you might like to

```c
#define NK_STATS
```

and call `nk_stats()` at the end of your code to see some statistics on how
the jump table was used during the run.


There is an inefficient part of my code in the building of the table of
influencers when you use random neighbors. I put a comment in there for
you. I was going to change it, there's a simple and elegant algorithm to do
what I clumsily do, but I didn't want to go through checking it to make
sure you could be confident. I'll put it in sometime. Again, I have a
reference to that part of Knuth if you're interested.

That was much more than I had intended to say...

Let me know if you have any problems, and I'll try to help. Hopefully
you'll be able to use large values of `N` and `K` in this code and not run
out of space. I'd be interested in knowing how large you manage to go. The
jump table in the `NK_LANSCAPE struct` is the largest thing malloc'd -
currently it's `128 * k * 2 * (sizeof(int) + sizeof(double))`, or something
like that. `128 = JUMP_MULTIPLIER`, which you can increase to make the
probability that the landscape is sufficiently random larger.


    Terry Jones (terry@santafe.edu)
    November 16, 1994.

    Terry Jones (terry@cliffs.ucsd.edu).
    September 7, 1998
    October 12, 1998

    README updated and Github repo created in September 2016.
