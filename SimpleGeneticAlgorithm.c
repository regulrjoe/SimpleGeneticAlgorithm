#include <math.h>
#include <stdint.h>
#include <stdio.h>
//#include "../libraries/gnuplot_i/src/gnuplot_i.h"

#define POPULATION_SIZE 20 /* Population Size */

void        crossover       (uint16_t a, uint16_t b, uint16_t *child_a, uint16_t *child_b);
void        showbits        (uint16_t x);
uint32_t    f               (uint16_t x);
void        newpopulation   (uint16_t *p);
void        mutation        (uint16_t *x);

int main() {
    uint16_t population[POPULATION_SIZE];
    uint16_t iters;
    /* Set seed for rand() */
    srand(time(NULL));

    printf("\nWelcome to the Simple Genetic Algorithm Software written by Joe Vázquez-Mellado\n");
    printf("\nPlease enter number of iterations: ");
    scanf("%hd", &iters);

    newpopulation(&population);
    //mutation(&);
}

/*  Mutation Operator.
    Create a random gene of a given chromosome */
void mutation(uint16_t *x) {
    uint8_t bit_to_mutate = rand() % (sizeof(uint16_t)*8 -1);
    uint16_t mask = 1 << bit_to_mutate;
    *x ^= mask;
}

/* Create starting population */
void newpopulation(uint16_t *p) {
    /* Set maxvalue of each chromosome based off the chromosome's size */
    uint16_t maxvalue = pow(2, sizeof(uint16_t)*8) - 1;

    for (int i = 0; i < POPULATION_SIZE; i++)
        p[i] = rand() % maxvalue;
}

/*  Crossover Operator.
    Get two new offsprings by crossing two given parents at a single point */
void crossover(uint16_t a, uint16_t b, uint16_t *child_a, uint16_t *child_b) {
    /*  Get random cutpoint's locus
        cutpoint >= 0 && cutpoint < chromosome's size */
    uint8_t cutpoint = rand() % ( sizeof(uint16_t)*8 - 1 );
    /*  Get mask with all 0's left to cutpoint's locus and all 1's to its right
        Used for & bit operation */
    uint16_t mask1 = pow(2, cutpoint) - 1;
    /*  Get mask with all 1's left to cutpoint's locus and all 0s to its right
        Used for | bit operation */
    uint16_t mask2 = pow(2, sizeof(uint16_t)*8) - 1 - mask1;
    /* Crossover bit operation */
    *child_a = (a | (b & mask1)) & (b | mask2);
    *child_b = (b | (a & mask1)) & (a | mask2);
}

/*  Print bit array of a given number */
void showbits(uint16_t x) {
        for (int i = (sizeof(uint16_t)*8) - 1; i >= 0; i--)
            (x & (1u << i)) ? putchar('1') : putchar('0');
    	printf("\n");
}

/*  Function to evaluate: f(x) = 4x^3 - 48x^2 + 144x */
uint32_t f(uint16_t x) {
    return (144 * x) - (48 * pow(x,2)) + (4 * pow(x,3));
}
