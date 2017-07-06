#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//#include "../libraries/gnuplot_i/src/gnuplot_i.h"

#define POPULATION_SIZE 5 /* Population Size */

/* Chromosome structure */
struct chromosome{
    uint8_t    x;              /* Phenotype / Genotype */
    int64_t    fitness;        /* Result of f(x) */
    double      probability;    /* Probability of selection */
    double      cdf;            /* Cumulative probability distribution */
    double      expected_p;     /* Expected population */
};


void                crossover       (uint8_t chr1, uint8_t chr2, uint8_t *child1, uint8_t *child2);
void                show_bits       (uint8_t x);
int64_t             f               (uint8_t x);
void                new_population  (struct chromosome p[POPULATION_SIZE]);
void                mutation        (uint8_t *x);
int64_t             get_fitness     (struct chromosome p[POPULATION_SIZE]);
void                probabilities   (struct chromosome p[POPULATION_SIZE], int64_t sof);
struct chromosome   *selection      (struct chromosome p[POPULATION_SIZE], int64_t *sof);
void                next_gen        (struct chromosome p[POPULATION_SIZE], struct chromosome pnext[POPULATION_SIZE], int64_t *sof);
void                run             (struct chromosome p[POPULATION_SIZE]);

int main() {
    struct chromosome population[POPULATION_SIZE];
    uint16_t iters;

    /* Set seed for rand() */
    srand(time(NULL));

    /*  UI */
    printf("\nWelcome to this Simple Genetic Algorithm Software written by Joe Vázquez-Mellado\n");
    printf("\nPlease enter number of iterations: ");
    scanf("%hd", &iters);

    /* Create starting population */

    new_population(population);
    run(population);
    return 0;
}

void run(struct chromosome p[POPULATION_SIZE]) {
    int64_t sum_of_fitness;
    struct chromosome new_population[POPULATION_SIZE];


        sum_of_fitness = get_fitness(p);
        probabilities(p, sum_of_fitness);

        next_gen(p, new_population, &sum_of_fitness);
        //for (uint16_t i = 0; i < POPULATION_SIZE; i++)
        //    printf("p[%d]: fitness: %llu\tprobability: %f\tcdf: %f\n", i, p[i].fitness, p[i].probability, p[i].cdf);

}

void next_gen(struct chromosome p[POPULATION_SIZE], struct chromosome pnext[POPULATION_SIZE], int64_t *sof) {
    struct chromosome *parent1, *parent2, child1, child2;

    parent1 = selection(p, sof);
    parent2 = selection(p, sof);
    printf("\nparent1: x: %d\tfitness: %llu\tprobability: %f\tcdf: %f\n", parent1->x, parent1->fitness, parent1->probability, parent1->cdf);
    printf("parent2: x: %d\tfitness: %llu\tprobability: %f\tcdf: %f\n", parent2->x, parent2->fitness, parent2->probability, parent2->cdf);

}

/*  Selection of chromosomes for next generation */
struct chromosome *selection(struct chromosome p[POPULATION_SIZE], int64_t *sof) {
    double roulette_shot;
    uint16_t i      = 0;
    roulette_shot   = (double)(rand() % *sof) / (double)*sof;
    printf("\nroulette_shot: %f\n", roulette_shot);
    while (roulette_shot > p[i].cdf) {
        printf("p[%d]: x: %d\tfitness: %llu\tprobability: %f\tcdf: %f\n", i, p[i].x, p[i].fitness, p[i].probability, p[i].cdf);
        i++;
    }

    return &p[i];
}

/*  Get probability of selection and Cumulative probability distribution of each chromosome */
void probabilities(struct chromosome p[POPULATION_SIZE], int64_t sof) {
    double sum_pr = 0;
    for (uint16_t i = 0; i < POPULATION_SIZE; i++) {
        p[i].probability    = (double)p[i].fitness / (double)sof;
        sum_pr              += p[i].probability;
        p[i].cdf            = sum_pr;
        printf("p[%d]: x: %d\tfitness: %llu\tprobability: %f\tcdf: %f\n", i, p[i].x, p[i].fitness, p[i].probability, p[i].cdf);
    }
}


/*  Evaluate fitness of each chromosome, normalize to handle negative values, and return sum of fitness   */
int64_t get_fitness(struct chromosome p[POPULATION_SIZE]) {
    int64_t     sum = 0,
                min = 0;

    /* Get fitness of each chromosome and find minimum value for normalization */
    for (uint16_t i = 0; i < POPULATION_SIZE; i++) {
        p[i].fitness = f(p[i].x);
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
void mutation(uint8_t *x) {
    uint8_t bit_to_mutate   = rand() % (sizeof(uint8_t)*8 -1);
    uint8_t  mask           = 1 << bit_to_mutate;
    *x                      ^= mask;
}


/* Create starting population */
void new_population(struct chromosome p[POPULATION_SIZE]) {
    /* Set maxvalue of each chromosome based off the chromosome's size */
    uint8_t maxvalue = pow(2, sizeof(uint8_t)*8) - 1;
    for (int i = 0; i < POPULATION_SIZE; i++){
        p[i].x = rand() % maxvalue;
        show_bits(p[i].x);
    }
}


/*  Crossover Operator.
    Get two new offsprings by crossing two given parents at a single point */
void crossover(uint8_t chr1, uint8_t chr2, uint8_t *child1, uint8_t *child2) {
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
    *child1 = (chr1 | (chr2 & mask1)) & (chr2 | mask2);
    *child2 = (chr2 | (chr1 & mask1)) & (chr1 | mask2);
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
