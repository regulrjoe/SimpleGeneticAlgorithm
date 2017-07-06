#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//#include "../libraries/gnuplot_i/src/gnuplot_i.h"

#define POPULATION_SIZE 10      /* SIZE OF THE POPULATION               */
#define MUTATION_RATE   0.01    /* PROBABILITY OF A MUTATION OCURRING   */
#define CROSSOVER_RATE  1       /* PROBABILITY OF A CROSSOVER HAPPENING */

/* Chromosome structure */
struct chromosome{
    uint8_t x;              /* Phenotype / Genotype                     */
    int64_t fitness;        /* Result of f(x)                           */
    int64_t true_fitness;   /* True fitness value (No Normalization)    */
    double  probability;    /* Probability of selection                 */
    double  cdf;            /* Cumulative probability distribution      */
    double  expected_p;     /* Expected population                      */
};


void                crossover           (uint8_t x1, uint8_t x2, struct chromosome *child1, struct chromosome *child2);
void                show_bits           (uint8_t x);
int64_t             f                   (uint8_t x);
void                first_gen           (struct chromosome p[POPULATION_SIZE]);
void                mutate              (uint8_t *x);
int64_t             get_fitness         (struct chromosome p[POPULATION_SIZE]);
void                get_probabilites    (struct chromosome p[POPULATION_SIZE], int64_t sof);
struct chromosome   *select             (struct chromosome p[POPULATION_SIZE], int64_t *sof);
void                procreate           (struct chromosome p[POPULATION_SIZE], struct chromosome pnext[POPULATION_SIZE], int64_t *sof);
void                run                 (struct chromosome p[POPULATION_SIZE], uint16_t iters);

int main() {
    struct chromosome population[POPULATION_SIZE];
    uint16_t iterations;

    /* Set seed for rand() */
    srand(time(NULL));

    /*  UI */
    printf("\nWelcome to this Simple Genetic Algorithm Software written by Joe Vázquez-Mellado\n");
    printf("\nPlease enter number of iterations: ");
    scanf("%hd", &iterations);

    /* Create starting population */

    first_gen(population);
    printf("before run p: %p\n", &population);
    run(population, iterations);
    return 0;
}

/*  Run the genetic algorithm through a given population a given amount of times */
void run(struct chromosome p[POPULATION_SIZE], uint16_t iters) {
    int64_t sum_of_fitness;
    struct chromosome new_gen[POPULATION_SIZE];

    for (uint16_t i = 0; i < iters; i++) {
        printf("\n\n----------------- GENERATION %d -----------------\n", i);
        sum_of_fitness = get_fitness(p);
        get_probabilites(p, sum_of_fitness);
        //printf("before procreate p: %p\t new_gen: %p\n", &p[0], &new_gen[0]);
        procreate(p, new_gen, &sum_of_fitness);
        for (uint16_t j = 0; j < POPULATION_SIZE; j++) {
            p[j].x = new_gen[j].x;
            new_gen[j].x = 0;
        }
    }
    printf("\n\n----------------- LAST GENERATION -----------------\n");

    get_fitness(p);
    for (uint16_t i = 0; i < POPULATION_SIZE; i++)
        printf("p[%d]: x: %d\ttrue_fitness: %lld\tfitness: %llu\t - %p\n", i, p[i].x, p[i].true_fitness, p[i].fitness, &p[i]);
}


/*  Create new generation */
void procreate(struct chromosome p[POPULATION_SIZE], struct chromosome pnext[POPULATION_SIZE], int64_t *sof) {
    struct chromosome *parent1, *parent2, child1, child2;
    printf("inside procreate p: %p\t new_gen: %p\n", &p[0], &pnext[0]);
    for (uint16_t i = 0; i < POPULATION_SIZE; i += 2) {
        parent1 = select(p, sof);
        parent2 = select(p, sof);
        printf("\nparent1: x: %d\ttrue fitness: %lld\tfitness: %llu\tprobability: %f\tcdf: %f\t - %p\n", parent1->x, parent1->true_fitness,parent1->fitness, parent1->probability, parent1->cdf, &parent1);
        printf("parent2: x: %d\ttrue fitness: %lld\tfitness: %llu\tprobability: %f\tcdf: %f\t - %p\n", parent2->x, parent2->true_fitness, parent2->fitness, parent2->probability, parent2->cdf, &parent2);
        show_bits(parent1->x);
        show_bits(parent2->x);

        if ( (double)(rand() % 1000) / 1000 <= CROSSOVER_RATE )
            crossover(parent1->x, parent2->x, &child1, &child2);
        //printf("crossover after:\n");
        //show_bits(child1.x);
        //show_bits(child2.x);

        //printf("mutation after:\n");
        if ( (double)(rand() % 1000) / 1000 <= MUTATION_RATE )
            mutate(&child1.x);
        if ( (double)(rand() % 1000) / 1000 <= MUTATION_RATE )
            mutate(&child2.x);
        printf("\nchilds bits:\n");
        pnext[i]    = child1;
        pnext[i+1]  = child2;
        printf("child1: x: %d - %p\n", child1.x, &child1);
        printf("child2: x: %d - %p\n", child2.x, &child2);
        show_bits(child1.x);
        show_bits(child2.x);
    }
}


/*  Select chromosomes for next generation through Roulette Wheel selection */
struct chromosome *select(struct chromosome p[POPULATION_SIZE], int64_t *sof) {
    double roulette_shot;
    uint16_t i      = 0;
    roulette_shot   = (double)(rand() % *sof) / (double)*sof;
    printf("\nroulette_shot: %f\n", roulette_shot);
    while (roulette_shot > p[i].cdf) {
       //printf("p[%d]: x: %d\ttrue_fitness: %lld\tfitness: %llu\tprobability: %f\tcdf: %f\t - %p\n", i, p[i].x, p[i].true_fitness, p[i].fitness, p[i].probability, p[i].cdf, &p[i]);
        i++;
    }
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
        p[i].fitness = f(p[i].x);
        p[i].true_fitness = p[i].fitness;
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
    uint8_t bit_to_mutate   = rand() % (sizeof(uint8_t)*8 -1);
    uint8_t  mask           = 1 << bit_to_mutate;
    *x                      ^= mask;
}


/* Create starting population */
void first_gen(struct chromosome p[POPULATION_SIZE]) {
    /* Set maxvalue of each chromosome based off the chromosome's size */
    uint8_t maxvalue = pow(2, sizeof(uint8_t)*8) - 1;
    for (int i = 0; i < POPULATION_SIZE; i++){
        p[i].x = rand() % maxvalue;
        show_bits(p[i].x);
    }
}


/*  Crossover Operator.
    Get two new offsprings by crossing two given parents at a single point */
void crossover(uint8_t x1, uint8_t x2, struct chromosome *child1, struct chromosome *child2) {
    /*  Get random cutpoint's locus
        cutpoint >= 0 && cutpoint < chromosome's size */
    uint8_t cutpoint = rand() % ( sizeof(uint8_t)*8 - 1 );
    /*  Get mask with all 0's left to cutpoint's locus and all 1's to its right
        Used for & bit operation */
    uint8_t mask1 = pow(2, cutpoint) - 1;
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
