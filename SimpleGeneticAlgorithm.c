#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include "../libraries/gnuplot_i/src/gnuplot_i.h"

/*
    TODO: Codificación Gray
    TODO: Implement Error Tolerance (EPSILON)
    TODO: Show plot for standard deviation, avg fitness, max fitness for each GENERATIONS
    TODO: Research De Jong Functions
    TODO: HowTo search minimum instead of maximum
    TODO: Read a lil bit David Goldberg book Genetic Algorithms in Search
    TODO: Obtener variancia de la población
    TODO: Use Tournament Selection
*/

#define POPULATION_SIZE 10      /* SIZE OF THE POPULATION                   */
#define MUTATION_RATE   0.05    /* PROBABILITY OF A MUTATION OCURRING       */
#define CROSSOVER_RATE  1       /* PROBABILITY OF A CROSSOVER HAPPENING     */
#define SLEEP_LGTH      1       /* GNUPLOT SLEEP TIMER THROUGH GENERATIONS  */

/* Variable structure */
struct variable{
    uint8_t     position_in_genotype;   /* Position of given variable in genes array    */
    int32_t     lower_lim;              /* Lower limit of variable                      */
    int32_t     upper_lim;              /* Upper limit of variable                      */
    double      precision;              /* Precision of variable value                  */
    uint8_t     bits_len;               /* Length of bits used for variable             */
    double      value;                  /* Current value of the variable                */
};

/* Chromosome structure */
struct chromosome{
    uint8_t         *genotype;      /* Chromosome genes array                   */
    int64_t         phenotype;      /* Decimal value of Genotype                */
    int64_t         fitness;        /* Result of f(phenotype)                   */
    int64_t         true_fitness;   /* True fitness value (No Normalization)    */
    double          probability;    /* Probability of selection                 */
    double          cdf;            /* Cumulative probability distribution      */
    double          expected_p;     /* Expected population                      */
    struct variable *vars;          /* Array of variables in genomea            */
};

/* Functions */
uint8_t             get_bits_len        (int32_t xl, int32_t xu, double precision);
uint8_t             define_vars         (struct variable *vars, uint8_t n);
void                assign_vars         (struct chromosome *chrome, struct variable *vars, uint8_t n);
void                eval_all_vars       (struct chromosome p[POPULATION_SIZE], uint8_t n);
void                eval_var            (struct chromosome chrome, uint8_t i);
void                init_population     (struct chromosome p[POPULATION_SIZE], struct variable *vars, uint8_t chrome_len, uint8_t n_vars);
void                eval_population     (struct chromosome p[POPULATION_SIZE], uint8_t n_vars);
double              eval_bin_array      (uint8_t *bin_array, uint8_t len);
int64_t             get_fitness         (struct chromosome p[POPULATION_SIZE]);
void                get_probabilites    (struct chromosome p[POPULATION_SIZE], int64_t sof);
void                crossover           (uint8_t *geno1, uint8_t *geno2, struct chromosome *child1, struct chromosome *child2, uint8_t chrome_len);
struct chromosome   *select_parent      (struct chromosome p[POPULATION_SIZE]);
void                mutate              (uint8_t *geno, uint8_t chrome_len);
void                run                 (struct chromosome p[POPULATION_SIZE], uint16_t iters, uint8_t chrome_len, struct variable *vars, uint8_t n_vars);
void                procreate           (struct chromosome p[POPULATION_SIZE], struct chromosome pnext[POPULATION_SIZE], uint8_t chrome_len, uint8_t n_vars);
int64_t             f                   (int64_t x);
double              rng                 ();
double              cud                 (double x, double a, double b);
double              c_unif              (double x, double a, double b);
uint16_t            d_unif              (double x, double a, double b);


int main() {
    struct chromosome   *population = (struct chromosome*) calloc(POPULATION_SIZE, sizeof(struct chromosome));
    uint16_t            iterations;
    uint8_t             num_vars;

    /*  UI */
    printf("\nWelcome to this Simple Genetic Algorithm Software written by Joe Vázquez-Mellado\n");
    printf("\nHow many variables do you wish to work with?: ");
    scanf("%d", &num_vars);

    struct variable *genome_vars = (struct variable*) calloc(num_vars, sizeof(struct variable));

    uint8_t chromosome_len = define_vars(genome_vars, num_vars);

    printf("\nPlease enter number of iterations: ");
    scanf("%hd", &iterations);

    /* Create starting population */
    init_population(population, genome_vars, chromosome_len, num_vars);
    /* Run algorithm */
    run(population, iterations, chromosome_len, genome_vars, num_vars);
    return 0;
}

/*  Define the variables metainformation to use in chromosome */
uint8_t define_vars(struct variable *vars, uint8_t n) {
    uint8_t len = 0;
    for (uint8_t i = 0; i < n; i++) {
        printf("\nVAR %d\nlower limit: ", i);
        scanf("%u", &vars[i].lower_lim);
        printf("upper limit: ");
        scanf("%u", &vars[i].upper_lim);
        printf("precision: ");
        scanf("%lf", &vars[i].precision);

        vars[i].bits_len = get_bits_len(vars[i].lower_lim, vars[i].upper_lim, vars[i].precision);
        vars[i].position_in_genotype = len;
        len += vars[i].bits_len;
    }
    return len;
}

/* Assign given variables to given chromosomes  */
void assign_vars(struct chromosome *chrome, struct variable *vars, uint8_t n) {
    chrome->vars = (struct variable*) calloc(n, sizeof(struct variable));
    for (uint8_t i = 0; i < n; i++) {
        chrome->vars[i].position_in_genotype    = vars[i].position_in_genotype;
        chrome->vars[i].lower_lim               = vars[i].lower_lim;
        chrome->vars[i].upper_lim               = vars[i].upper_lim;
        chrome->vars[i].precision               = vars[i].precision;
        chrome->vars[i].bits_len                = vars[i].bits_len;
    }
}

/* Evalute all variables of all chromosomes */
void eval_all_vars(struct chromosome p[POPULATION_SIZE], uint8_t n) {
    for (uint8_t i = 0; i < POPULATION_SIZE; i++)
        for (uint8_t j = 0; j < n; j++)
            eval_var(p[i], j);
}

/* Evaluate a single variable */
void eval_var(struct chromosome chrome, uint8_t i) {
    chrome.vars[i].value = eval_bin_array(&chrome.genotype[chrome.vars[i].position_in_genotype], chrome.vars[i].bits_len);
}

/* Create starting population */
void init_population(struct chromosome p[POPULATION_SIZE], struct variable *vars, uint8_t chrome_len, uint8_t n_vars) {
    for (int i = 0; i < POPULATION_SIZE; i++){
        p[i].genotype = (uint8_t*) calloc(chrome_len, sizeof(uint8_t));
        p[i].vars = (struct variable*) calloc(n_vars, sizeof(struct variable));

        for (uint8_t j = 0; j < chrome_len; j++)
            p[i].genotype[j] = d_unif(rng(), 0, 2);

        p[i].phenotype = (int64_t)eval_bin_array(p[i].genotype, chrome_len);

        for (uint8_t k = 0; k < n_vars; k++) {
            p[i].vars[k].position_in_genotype   = vars[k].position_in_genotype;
            p[i].vars[k].lower_lim              = vars[k].lower_lim;
            p[i].vars[k].upper_lim              = vars[k].upper_lim;
            p[i].vars[k].precision              = vars[k].precision;
            p[i].vars[k].bits_len               = vars[k].bits_len;
            eval_var(p[i], k);
        }
    }
}

/*  Get amount of bits needed to cover search space with given precision
    ceil( log10(xu-xl) / log10(precision) )*/
uint8_t get_bits_len(int32_t xl, int32_t xu, double precision) {
    return (uint8_t)ceil( log2( ((double)(xu-xl)/precision) - 1.0) );
}


/* Evaluate binary int array */
double eval_bin_array(uint8_t *bin_array, uint8_t len){
    double value = 0;
    for (int8_t i = len-1; i >= 0; i--)
        value += bin_array[i] * pow(2, len-i-1);

    return value;
}

/* Evaluate population */
void eval_population(struct chromosome p[POPULATION_SIZE], uint8_t n_vars) {
    /* Get generation fitness */
    int64_t sum_of_fitness = get_fitness(p);
    /* Get generation probabilities */
    get_probabilites(p, sum_of_fitness);
    /* Get generation vars values */
    eval_all_vars(p, n_vars);
}

/*  Run the genetic algorithm through a given population a given amount of times */
void run(struct chromosome p[POPULATION_SIZE], uint16_t iters, uint8_t chrome_len, struct variable *vars, uint8_t n_vars) {
    int64_t sum_of_fitness;
    struct chromosome new_gen[POPULATION_SIZE+1];
    // gnuplot_ctrl *h = gnuplot_init();

    for (uint16_t i = 0; i < POPULATION_SIZE+1; i++) {
        new_gen[i].genotype = (uint8_t*) calloc(chrome_len, sizeof(uint8_t));
        assign_vars(&new_gen[i], vars, n_vars);
    }

    /* Run generations */
    for (uint16_t i = 0; i < iters; i++) {
        printf("\n\n----------------- GENERATION %d -----------------\n", i);
        gnuplot_setstyle(h, "lines");
        gnuplot_cmd(h, "set xrange [0:255]");
        gnuplot_cmd(h, "set yrange [-50000:2000]");
        gnuplot_plot_equation(h, "60*x - x**2", "f(x) = 60x - x^2");
        gnuplot_setstyle(h, "points");


        /* Evaluate fitness, probabilities and var values */
        eval_population(p, n_vars);

        /* Get new generation */
        procreate(p, new_gen, chrome_len, n_vars);

        /* Assign new generation's properties to current population */
        for (uint16_t j = 0; j < POPULATION_SIZE; j++) {
            double  px  = p[j].phenotype,
                    ptf = p[j].true_fitness;
            gnuplot_plot_xy(h, &px, &ptf, 1, "");

            /* assign new_gen genotype to current_gen */
            for (uint8_t k = 0; k < chrome_len; k++)
                p[j].genotype[k] = new_gen[j].genotype[k];
            /* assign new_gen phenotype to current_gen */
            p[j].phenotype = (int64_t)eval_bin_array(p[j].genotype, chrome_len);
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

    /* Get last generation fitness, at this point probabilities are useless */
    get_fitness(p);
    for (uint16_t i = 0; i < POPULATION_SIZE; i++){
        double  px  = p[i].phenotype,
                ptf = p[i].true_fitness;
        gnuplot_plot_xy(h, &px, &ptf, 1, "");
    }
    gnuplot_close(h);
}

/*  Evaluate fitness of each chromosome, normalize to handle negative values, and return sum of fitness   */
int64_t get_fitness(struct chromosome p[POPULATION_SIZE]) {
    int64_t     sum = 0,
                min = 0;

    /* Get fitness of each chromosome and find minimum value for normalization */
    for (uint16_t i = 0; i < POPULATION_SIZE; i++) {
        p[i].fitness        = f(p[i].phenotype);
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

/*  Get probability of selection and Cumulative probability distribution of each chromosome */
void get_probabilites(struct chromosome p[POPULATION_SIZE], int64_t sof) {
    double sum_pr = 0;
    for (uint16_t i = 0; i < POPULATION_SIZE; i++) {
        p[i].probability    = (double)p[i].fitness / (double)sof;
        sum_pr              += p[i].probability;
        p[i].cdf            = sum_pr;
        printf("p[%d]: phenotype: %llu\ttrue_fitness: %lld\tfitness: %llu\tprobability: %f\tcdf: %f\t %p\n", i, p[i].phenotype, p[i].true_fitness, p[i].fitness, p[i].probability, p[i].cdf, &p[i]);
    }
}


/*  Create new generation */
void procreate(struct chromosome p[POPULATION_SIZE], struct chromosome pnext[POPULATION_SIZE], uint8_t chrome_len, uint8_t n_vars) {
    struct chromosome *parent1, *parent2, child1, child2;

    child1.genotype     = (uint8_t*) calloc(chrome_len, sizeof(uint8_t));
    child2.genotype     = (uint8_t*) calloc(chrome_len, sizeof(uint8_t));
    child1.vars         = (struct variable*) calloc(n_vars, sizeof(struct variable));
    child2.vars         = (struct variable*) calloc(n_vars, sizeof(struct variable));
    uint16_t i_fittest  = 0;
    int64_t fittest     = 0;

    /* Pass fittest chromosome straight to next gen */
    for (uint16_t i = 0; i < POPULATION_SIZE; i++)
        if (p[i].fitness > fittest) {
            i_fittest = i;
            fittest = p[i].fitness;
        }
    for (uint8_t j = 0; j < chrome_len; j++)
        pnext[0].genotype[j] = p[i_fittest].genotype[j];

    for (uint16_t i = 1; i < POPULATION_SIZE; i += 2) {
        parent1 = select_parent(p);
        parent2 = select_parent(p);
        for (uint8_t j = 0; j < chrome_len; j++) {
            child1.genotype[j] = parent1->genotype[j];
            child2.genotype[j] = parent2->genotype[j];
        }
        /* Do crossover if CROSSOVER_RATE is met */
        if (rng() <= CROSSOVER_RATE)
            crossover(parent1->genotype, parent2->genotype, &child1, &child2, chrome_len);
        /* Do mutation if MUTATION_RATE is met */
        if (rng() <= MUTATION_RATE)
            mutate(child1.genotype, chrome_len);
        if (rng() <= MUTATION_RATE )
            mutate(child2.genotype, chrome_len);

        /* Assign children to new generation */
        for (uint8_t j = 0; j < chrome_len; j++) {
            pnext[i].genotype[j]    = child1.genotype[j];
            pnext[i+1].genotype[j]  = child2.genotype[j];
        }
    }
}

/*  Crossover Operator.
    Get two new offsprings by crossing two given parents at a single point */
void crossover(uint8_t *geno1, uint8_t *geno2, struct chromosome *child1, struct chromosome *child2, uint8_t chrome_len) {
    /*  Get random cutpoint's locus
        cutpoint >= 1 && cutpoint < chromosome's size */
    uint8_t cutpoint = d_unif(rng(), 1, chrome_len);

    for (uint8_t i = cutpoint; i < chrome_len; i++) {
        child1->genotype[i] = geno2[i];
        child2->genotype[i] = geno1[i];
    }
}


/*  Mutation Operator.
    Create a random gene of a given chromosome */
void mutate(uint8_t *geno, uint8_t chrome_len) {
    uint8_t locus   = d_unif(rng(), 0, chrome_len);
    geno[locus] = (geno[locus] == 0) ? 1 : 0;
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

/*  Select chromosomes for next generation through Roulette Wheel selection */
struct chromosome *select_parent(struct chromosome p[POPULATION_SIZE]) {
    double roulette_shot;
    uint16_t i      = 0;
    roulette_shot   = rng();
    //printf("\nroulette_shot: %f\n", roulette_shot);
    while (roulette_shot > p[i].cdf)
        i++;
    //printf("p[%d]: x: %d\ttrue_fitness: %lld\tfitness: %llu\tprobability: %f\tcdf: %f\t - %p\n", i, p[i].x, p[i].true_fitness, p[i].fitness, p[i].probability, p[i].cdf, &p[i]);
    return &p[i];
}

/*  Function to evaluate: f(x) = 60x - x^2 */
int64_t f(int64_t x) {
    return 60 * x - pow(x,2);
}
