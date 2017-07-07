#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include "../libraries/gnuplot_i/src/gnuplot_i.h"

/*  TODO: Ints to char array
    TODO: Calculate bits needed
    TODO: Implement Error Tolerance (EPSILON)
    TODO: Show plot for standard deviation, avg fitness, max fitness for each GENERATIONS
    TODO: Pass max fitness to next gen without applying operations to it
    TODO: Evaluate multivariable equation
    TODO: Ask for length and limits for each variable
    TODO: Research De Jong Functions
    TODO: HowTo search minimum instead of maximum
    TODO: Read a lil bit David Goldberg book Genetic Algorithms in Search */

#define POPULATION_SIZE 10      /* SIZE OF THE` POPULATION                   */
#define MUTATION_RATE   0.01    /* PROBABILITY OF A MUTATION OCURRING       */
#define CROSSOVER_RATE  1       /* PROBABILITY OF A CROSSOVER HAPPENING     */
#define SLEEP_LGTH      1       /* GNUPLOT SLEEP TIMER THROUGH GENERATIONS  */

/* Chromosome structure */
struct chromosome{
    uint8_t x;              /* Phenotype / Genotype                     */
    int64_t fitness;        /* Result of f(x)                           */
    int64_t true_fitness;   /* True fitness value (No Normalization)    */
    double  probability;    /* Probability of selection                 */
    double  cdf;            /* Cumulative probability distribution      */
    double  expected_p;     /* Expected population                      */
};

/* Functions */
void                crossover           (uint8_t x1, uint8_t x2, struct chromosome *child1, struct chromosome *child2);
void                show_bits           (uint8_t x);
int64_t             f                   (uint8_t x);
void                first_gen           (struct chromosome p[POPULATION_SIZE]);
void                mutate              (uint8_t *x);
int64_t             get_fitness         (struct chromosome p[POPULATION_SIZE]);
void                get_probabilites    (struct chromosome p[POPULATION_SIZE], int64_t sof);
struct chromosome   *select_parent      (struct chromosome p[POPULATION_SIZE], int64_t *sof);
void                procreate           (struct chromosome p[POPULATION_SIZE], struct chromosome pnext[POPULATION_SIZE], int64_t *sof);
void                run                 (struct chromosome p[POPULATION_SIZE], uint16_t iters);
double              rng                 ();
double              cud                 (double x, double a, double b);
double              c_unif              (double x, double a, double b);
uint16_t            d_unif              (double x, double a, double b);


int main() {

    struct chromosome population[POPULATION_SIZE];
    uint16_t iterations;
    /* |||||||||||||||||||||||||||||||||||||||||||||||||| */
    /* |||||||||||||||||| TESTING AREA |||||||||||||||||| */
    /* |||||||||||||||||||||||||||||||||||||||||||||||||| */
    //printf("rng: %f\n", rng());
    //uint8_t cutpoint = d_unif(rng(), 1, sizeof(uint8_t)*8);
    //printf("cutpoint: %d\n", cutpoint);
    //uint8_t mask1 = pow(2, cutpoint) - 1;
    //printf("mask1: ");
    //show_bits(mask1);

    //printf("sizeof(uint8_t)*8: %lu\n", sizeof(uint8_t)*8);
    //uint8_t maxvalue = pow(2, sizeof(uint8_t)*8)-1;
    //printf("maxvalue: %d\n", maxvalue);

    //printf("discrete rng (1-5)[1,6]:\t%d\n", d_unif(rng(), 1, 6));
    //printf("continuous rng (1-5)[1,6]:\t%f\n", c_unif(rng(), 1, 6));
    /* |||||||||||||||||||||||||||||||||||||||||||||||||| */
    /* |||||||||||||||||||||||||||||||||||||||||||||||||| */
    /* |||||||||||||||||||||||||||||||||||||||||||||||||| */

    /*  UI */
    printf("\nWelcome to this Simple Genetic Algorithm Software written by Joe Vázquez-Mellado\n");
    printf("\nPlease enter number of iterations: ");
    scanf("%hd", &iterations);

    /* Create starting population */
    first_gen(population);
    /* Run algorithm */
    run(population, iterations);
    return 0;
}

/*  Run the genetic algorithm through a given population a given amount of times */
void run(struct chromosome p[POPULATION_SIZE], uint16_t iters) {
    int64_t sum_of_fitness;
    struct chromosome new_gen[POPULATION_SIZE];
    gnuplot_ctrl *h = gnuplot_init();

    /* Run generations */
    for (uint16_t i = 0; i < iters; i++) {
        printf("\n\n----------------- GENERATION %d -----------------\n", i);
        gnuplot_setstyle(h, "lines");
        gnuplot_cmd(h, "set xrange [0:255]");
        gnuplot_cmd(h, "set yrange [-50000:2000]");
        gnuplot_plot_equation(h, "60*x - x**2", "f(x) = 60x - x^2");
        gnuplot_setstyle(h, "points");
        /* Get generation fitness */
        sum_of_fitness = get_fitness(p);
        /* Get generation probabilities */
        get_probabilites(p, sum_of_fitness);
        /* Get new generation */
        procreate(p, new_gen, &sum_of_fitness);
        for (uint16_t j = 0; j < POPULATION_SIZE; j++) {
            double  px  = p[j].x,
                    ptf = p[j].true_fitness;

            gnuplot_plot_xy(h, &px, &ptf, 1, "");
            p[j].x          = new_gen[j].x;
            new_gen[j].x    = 0;
        }
        sleep(SLEEP_LGTH);
        gnuplot_resetplot(h);
    }
    printf("\n\n----------------- LAST GENERATION -----------------\n");
    gnuplot_setstyle(h, "lines");
    gnuplot_cmd(h, "set xrange [0:255]");
    gnuplot_cmd(h, "set yrange [-50000:2000]");
    gnuplot_plot_equation(h, "60*x - x**2", "f(x) = 60x - x^2");
    gnuplot_setstyle(h, "points");
    /* Get last generation fitness at this point, probabilities are useless */
    get_fitness(p);
    for (uint16_t i = 0; i < POPULATION_SIZE; i++){
        double  px  = p[i].x,
                ptf = p[i].true_fitness;
        gnuplot_plot_xy(h, &px, &ptf, 1, "");
        printf("p[%d]: x: %d\ttrue_fitness: %lld\tfitness: %llu\t - %p\n", i, p[i].x, p[i].true_fitness, p[i].fitness, &p[i]);
    }
    gnuplot_close(h);
}

/*  Cryptographically Secure Pseudornadom Number Generator
    Read bytes from /dev/random device file. Each byte from the file is
    a cryptographically random value from 0-255. This function concatenates
    those bytes to generate a random number of an arbitrary size, only to then
    get a point value out of it. */
double rng() {
    /* Establish connection with /dev/random */
    FILE *fp = fopen("/dev/random", "r");

    /* Abort if connection is lost */
    if (!fp) {
        perror("RNG");
        exit(-1);
    }
    uint16_t    true_value = 0;
    double      point_value;

    for (uint8_t i = 0; i < sizeof(true_value); i++) {
        true_value <<= 8;
        true_value |= fgetc(fp);
    }
    /* Get a point value from 0 to 1 */
    point_value = (double)true_value / (pow(2,(sizeof(true_value)*8)) - 1);
    /* Close connection with /dev/random */
    fclose(fp);
    //printf("\nRNG true_value: %u\tpoint_value: %f\n", true_value, point_value);
    return point_value;
}

/* Discrete Uniform Value.
    With a point value of x when 0 ≤ x ≤ 1, get an uniformly
    equivalent int value y as a ≤ y < b*/
uint16_t d_unif(double x, double a, double b) {
    return (x == 1) ? (uint16_t)b-1 : (uint16_t)floor(x*(b-a) + a );
}
/*  Continuous Uniform Value.
    With a point value of x when 0 ≤ x ≤ 1, get an uniformly
    equivalent point value y as a ≤ y < b */
double c_unif(double x, double a, double b) {
    return (x == 1) ? b-0.000001 : x*(b-a) + a;
}
/*  Continuous Uniform Distribution */
double cud(double x, double a, double b) {
    return (x-a) / (b-a) * (b-a);
}


/*  Create new generation */
void procreate(struct chromosome p[POPULATION_SIZE], struct chromosome pnext[POPULATION_SIZE], int64_t *sof) {
    struct chromosome *parent1, *parent2, child1, child2;
    printf("inside procreate p: %p\t new_gen: %p\n", &p[0], &pnext[0]);
    for (uint16_t i = 0; i < POPULATION_SIZE; i += 2) {
        parent1 = select_parent(p, sof);
        parent2 = select_parent(p, sof);

        //printf("\nparent1: x: %d\ttrue fitness: %lld\tfitness: %llu\tprobability: %f\tcdf: %f\t - %p\n", parent1->x, parent1->true_fitness,parent1->fitness, parent1->probability, parent1->cdf, &parent1);
        //printf("parent2: x: %d\ttrue fitness: %lld\tfitness: %llu\tprobability: %f\tcdf: %f\t - %p\n", parent2->x, parent2->true_fitness, parent2->fitness, parent2->probability, parent2->cdf, &parent2);
        //show_bits(parent1->x);
        //show_bits(parent2->x);

        /* Do crossover if CROSSOVER_RATE is met */
        if (rng() <= CROSSOVER_RATE)
            crossover(parent1->x, parent2->x, &child1, &child2);

        //printf("crossover after:\n");
        //show_bits(child1.x);
        //show_bits(child2.x);
        //printf("mutation after:\n");

        /* Do mutation if MUTATION_RATE is met */
        if (rng() <= MUTATION_RATE)
            mutate(&child1.x);
        if (rng() <= MUTATION_RATE )
            mutate(&child2.x);

        /* Assign children to new generation */
        pnext[i]    = child1;
        pnext[i+1]  = child2;

        //printf("\nchilds bits:\n");
        //printf("child1: x: %d - %p\n", child1.x, &child1);
        //printf("child2: x: %d - %p\n", child2.x, &child2);
        //show_bits(child1.x);
        //show_bits(child2.x);
    }
}


/*  Select chromosomes for next generation through Roulette Wheel selection */
struct chromosome *select_parent(struct chromosome p[POPULATION_SIZE], int64_t *sof) {
    double roulette_shot;
    uint16_t i      = 0;
    roulette_shot   = rng();

    //printf("\nroulette_shot: %f\n", roulette_shot);

    while (roulette_shot > p[i].cdf)
        i++;

    //printf("p[%d]: x: %d\ttrue_fitness: %lld\tfitness: %llu\tprobability: %f\tcdf: %f\t - %p\n", i, p[i].x, p[i].true_fitness, p[i].fitness, p[i].probability, p[i].cdf, &p[i]);

    return &p[i];
}


/*  Get probability of selection and Cumulative probability distribution of each chromosome */
void get_probabilites(struct chromosome p[POPULATION_SIZE], int64_t sof) {
    double sum_pr = 0;
    for (uint16_t i = 0; i < POPULATION_SIZE; i++) {
        p[i].probability    = (double)p[i].fitness / (double)sof;
        sum_pr              += p[i].probability;
        p[i].cdf            = sum_pr;

        printf("p[%d]: x: %d\ttrue_fitness: %lld\tfitness: %llu\tprobability: %f\tcdf: %f\t %p\n", i, p[i].x, p[i].true_fitness, p[i].fitness, p[i].probability, p[i].cdf, &p[i]);
    }
}


/*  Evaluate fitness of each chromosome, normalize to handle negative values, and return sum of fitness   */
int64_t get_fitness(struct chromosome p[POPULATION_SIZE]) {
    int64_t     sum = 0,
                min = 0;

    /* Get fitness of each chromosome and find minimum value for normalization */
    for (uint16_t i = 0; i < POPULATION_SIZE; i++) {
        p[i].fitness        = f(p[i].x);
        p[i].true_fitness   = p[i].fitness;
        if (p[i].fitness < min)
            min = p[i].fitness;

        //printf("p[%d] x: %d\tfitness: %lld\tmin: %lld\n", i, p[i].x, p[i].fitness, min);
    }

    /* Normalize fitness to handle negative values and obtain sum of all */
    for (uint16_t i = 0; i < POPULATION_SIZE; i++) {
        p[i].fitness    += labs(min);
        sum             += p[i].fitness;

        //printf("p[%d] x: %d\tfitness: %lld\tsum: %lld\n", i, p[i].x, p[i].fitness, sum);
    }

    //printf("\n\n SUM OF FITNESS: %lld\n", sum);

    return sum;
}


/*  Mutation Operator.
    Create a random gene of a given chromosome */
void mutate(uint8_t *x) {
    uint8_t bit_to_mutate   = d_unif(rng(), 0, sizeof(uint8_t)*8);
    uint8_t  mask           = 1 << bit_to_mutate;
    *x                      ^= mask;
}


/* Create starting population */
void first_gen(struct chromosome p[POPULATION_SIZE]) {
    /* Set maxvalue of each chromosome based off the chromosome's size */
    uint8_t maxvalue = pow(2, sizeof(uint8_t)*8 - 1);
    //printf("maxvalue: %d\n", maxvalue);
    for (int i = 0; i < POPULATION_SIZE; i++){
        uint8_t rnd = d_unif(rng(), 0, maxvalue+1);
    //    printf("rnd: %d", rnd);
        p[i].x = rnd;
        show_bits(p[i].x);
    }
}


/*  Crossover Operator.
    Get two new offsprings by crossing two given parents at a single point */
void crossover(uint8_t x1, uint8_t x2, struct chromosome *child1, struct chromosome *child2) {
    /*  Get random cutpoint's locus
        cutpoint >= 0 && cutpoint < chromosome's size */
    uint8_t cutpoint = d_unif(rng(), 1, sizeof(uint8_t)*8);
    //printf("cutpoint: %d", cutpoint);
    /*  Get mask with all 0's left to cutpoint's locus and all 1's to its right
        Used for & bit operation */
    uint8_t mask1 = pow(2, cutpoint) - 1;
    //printf("mask1: ");
    show_bits(mask1);
    /*  Get mask with all 1's left to cutpoint's locus and all 0s to its right
        Used for | bit operation */
    uint8_t mask2 = pow(2, sizeof(uint8_t)*8) - 1 - mask1;
    /* Crossover bit operation */
    child1->x = (x1 | (x2 & mask1)) & (x2 | mask2);
    child2->x = (x2 | (x1 & mask1)) & (x1 | mask2);

    //printf("cutpoint: %d\n", cutpoint);
}


/*  Function to evaluate: f(x) = 60x - x^2 */
int64_t f(uint8_t x) {
    return 60 * x - pow(x,2);
}


/*  Print bit array of a given number */
void show_bits(uint8_t x) {
    for (int i = (sizeof(uint8_t)*8) - 1; i >= 0; i--)
        (x & (1u << i)) ? putchar('1') : putchar('0');
	printf("\n");
}
