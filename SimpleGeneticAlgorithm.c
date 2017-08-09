#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include "../libraries/gnuplot_i/src/gnuplot_i.h"

/*
    BUGS

    PENDING
    TODO: Codificación Gray
    TODO: Implement Error Tolerance (EPSILON)
    TODO: Show plot for standard deviation, avg fitness, max fitness for each GENERATIONS
    TODO: Read a lil bit David Goldberg book Genetic Algorithms in Search
    TODO: Obtener variancia de la población
    TODO: Two-point crossover
    TODO: define_vars creates and returns pointer to vars instead of using pointer as argument
    TODO: add chrome_len attribute to chromosome struct
    TODO: Restructure Code into modal implementation.
    TODO: Producir bitácora de configuración y resultados por corrida.
    TODO: Función de selección Boltzmann
*/

/* Constants */
#define POPULATION_SIZE 30      /* SIZE OF THE POPULATION                       */
#define MUTATION_RATE   0.2     /* PROBABILITY OF A MUTATION OCURRING           */
#define CROSSOVER_RATE  0.8     /* PROBABILITY OF A CROSSOVER HAPPENING         */
#define SLEEP_LENGTH    1       /* GNUPLOT SLEEP TIMER THROUGH GENERATIONS      */
#define PASSTHROUGH_PCT 0.1     /* % OF FITTEST CHROMES WHO GO STRAIGHT THROUGH */
#define TOURNEY_RATE    0.8     /* PROBABILITY OF SELECTING FITTEST OF ROUND    */


/* Variable structure */
struct variable{
    uint8_t pos;        /* Position of given variable in genes array    */
    float   xl;         /* Lower limit of variable                      */
    float   xu;         /* Upper limit of variable                      */
    float   pr;         /* Precision of variable value                  */
    uint8_t bits_len;   /* Length of bits used for variable             */
};

/* Chromosome structure */
struct chromosome{
    uint8_t *genotype;      /* Chromosome bit array                     */
    int64_t phenotype;      /* Decimal value of Genotype                */
    double  fitness;        /* Result of f(phenotype)                   */
    double  probability;    /* Probability of selection                 */
    double  cdf;            /* Cumulative probability distribution      */
    double  expected_p;     /* Expected population                      */
    double  *values;        /* Decoded value of genes                   */
};

/* Population structure */
struct population{
    struct chromosome   chromosomes[POPULATION_SIZE];
    double              sum_of_fitness;
    double              avg_fitness;
    double              fittest;
    uint8_t             chromosome_len;
    uint8_t             num_vars;
    struct variable     *vars;
}

/* Functions */
uint8_t             define_vars         (struct variable *vars, uint8_t nv);
uint8_t             get_genotype_len    (float xl, float xu, float pr);
void                init_population     (struct chromosome p[POPULATION_SIZE], struct variable *vars, uint8_t chrome_len, uint8_t nv);
void                run                 (struct chromosome p[POPULATION_SIZE], uint16_t iters, uint8_t chrome_len, struct variable *vars, uint8_t nv);
void                eval_population     (struct chromosome p[POPULATION_SIZE], struct variable *vars, uint8_t nv);
void                eval_vars           (struct chromosome chrome, struct variable *vars, uint8_t nv);
double              eval_bits           (uint8_t *bits, uint8_t len);
double              get_fitness         (struct chromosome p[POPULATION_SIZE]);
void                get_probabilites    (struct chromosome p[POPULATION_SIZE], double sof);
void                procreate           (struct chromosome p[POPULATION_SIZE], struct chromosome new_gen[POPULATION_SIZE], uint8_t chrome_len, uint8_t nv);
void                crossover           (uint8_t *geno1, uint8_t *geno2, struct chromosome *child1, struct chromosome *child2, uint8_t chrome_len);
struct chromosome   *select_parent      (struct chromosome p[POPULATION_SIZE]);
void                mutate              (uint8_t *geno, uint8_t chrome_len);
struct chromosome   *tournament         (struct chromosome p[POPULATION_SIZE]);
struct chromosome   *roulette           (struct chromosome p[POPULATION_SIZE]);
void                print_chromosome    (struct chromosome chrome, uint8_t nv);
void                print_genotype      (uint8_t *genotype, uint8_t chrome_len);
double              f                   (double *values);
double              rng                 ();
double              cud                 (double x, double a, double b);
double              c_unif              (double x, double a, double b);
uint16_t            d_unif              (double x, double a, double b);
void                plot_population     (gnuplot_ctrl *h, struct chromosome p[POPULATION_SIZE], struct variable *vars, uint8_t nv);
void                save_datapoints     (struct chromosome p[POPULATION_SIZE], uint8_t nv);

int main() {
    struct chromosome   *population = (struct chromosome*) calloc(POPULATION_SIZE, sizeof(struct chromosome));
    uint16_t            iterations;
    uint8_t             nv; /*  Number of variables */

    /*  UI */
    printf("\nWelcome to this Simple Genetic Algorithm Software written by Joe Vázquez-Mellado\n");
    printf("\nHow many variables do you wish to work with?: ");
    scanf("%d", &nv);

    struct variable *genome_vars = (struct variable*) calloc(nv, sizeof(struct variable));
    uint8_t chrome_len = define_vars(genome_vars, nv);

    printf("\nEnter number of generations: ");
    scanf("%hd", &iterations);

    /* Create starting population */
    init_population(population, genome_vars, chrome_len, nv);
    /* Run algorithm */
    run(population, iterations, chrome_len, genome_vars, nv);
    return 0;
}

/*  Define the variables metainformation to use in chromosome */
uint8_t define_vars(struct variable *vars, uint8_t nv) {
    uint8_t len = 0;
    for (uint8_t i = 0; i < nv; i++) {
        printf("\nVAR %d\nlower limit: ", i);
        scanf("%f", &vars[i].xl);
        printf("upper limit: ");
        scanf("%f", &vars[i].xu);
        printf("precision: ");
        scanf("%f", &vars[i].pr);

        vars[i].bits_len = get_genotype_len(vars[i].xl, vars[i].xu, vars[i].pr);
        vars[i].pos = len;
        len += vars[i].bits_len;
    }
    return len;
}

/*  Get amount of bits needed to cover search space with given precision.
    Equivalent to: ceil( log10(xu-xl) / log10(precision) )*/
uint8_t get_genotype_len(float xl, float xu, float pr) {
    return (uint8_t)ceil( log2((xu-xl)/pr - 1.0) );
}

/* Create starting population */
void init_population(struct chromosome p[POPULATION_SIZE], struct variable *vars, uint8_t chrome_len, uint8_t nv) {
    for (int i = 0; i < POPULATION_SIZE; i++){
        /* Memory allocation */
        p[i].genotype   = (uint8_t*) calloc(chrome_len, sizeof(uint8_t));
        p[i].values     = (double*) calloc(nv, sizeof(double));

        /* Genotype intialization with random values 0-1 */
        for (uint8_t j = 0; j < chrome_len; j++)
            p[i].genotype[j] = d_unif(rng(), 0, 2);

        /* Get genotype decimal value */
        p[i].phenotype = (int64_t)eval_bits(p[i].genotype, chrome_len);
    }
}

/*  Run the genetic algorithm through a given population a given amount of times */
void run(struct chromosome p[POPULATION_SIZE], uint16_t iters, uint8_t chrome_len, struct variable *vars, uint8_t nv) {
    /* Initialize gnuplot handler */
    gnuplot_ctrl *h = gnuplot_init();

    /* Memory allocation */
    struct chromosome *new_gen = (struct chromosome*) calloc(POPULATION_SIZE+1, sizeof(struct chromosome));
    for (uint16_t i = 0; i < POPULATION_SIZE+1; i++) {
        new_gen[i].genotype = (uint8_t*) calloc(chrome_len, sizeof(uint8_t));
        new_gen[i].values   = (double*) calloc(nv, sizeof(double));
    }

    /* Run generations */
    for (uint16_t i = 0; i < iters; i++) {
        printf("\n\n----------------- GENERATION %d -----------------\n", i);

        /* Evaluate fitness, probabilities and var values */
        eval_population(p, vars, nv);

        /* Print chromosomes */
        for (uint8_t i = 0; i < POPULATION_SIZE; i++) {
            print_chromosome(p[i], nv);
        }

        /* Get new generation */
        procreate(p, new_gen, chrome_len, nv);

        /* Assign new generation's properties to current population */
        for (uint16_t j=0; j < POPULATION_SIZE; j++) {
            // double  px  = p[j].phenotype,
            //         ptf = p[j].fitness;
            //  gnuplot_plot_xy(h, &px, &ptf, 1, "");

            /* assign new_gen genotype to current_gen */
            for (uint8_t k=0; k < chrome_len; k++)
                p[j].genotype[k] = new_gen[j].genotype[k];
            /* Get new phenotype */
            p[j].phenotype = (int64_t)eval_bits(p[j].genotype, chrome_len);
        }
        save_datapoints(p, nv);
        plot_population(h, p, vars, nv);
    }
    // printf("\n\n----------------- LAST GENERATION -----------------\n");
    // // gnuplot_setstyle(h, "lines");
    // // gnuplot_cmd(h, "set xrange [0:255]");
    // // gnuplot_cmd(h, "set yrange [-50000:2000]");
    // // gnuplot_plot_equation(h, "60*x - x**2", "f(x) = 60x - x^2");
    // // gnuplot_setstyle(h, "points");
    //
    // /* Get last generation fitness, at this point probabilities are useless */
    // eval_population(p, vars, nv);
    // /* Print chromosomes */
    // for (uint8_t i = 0; i < POPULATION_SIZE; i++) {
    //     print_chromosome(p[i], nv);
    // }
    //
    // // for (uint16_t i = 0; i < POPULATION_SIZE; i++){
    // //     double  px  = p[i].phenotype,
    // //             ptf = p[i].fitness;
    // //      gnuplot_plot_xy(h, &px, &ptf, 1, "");
    // // }
    gnuplot_close(h);
}

/* Plot equation surface if 3D */
void plot_population(gnuplot_ctrl *h, struct chromosome p[POPULATION_SIZE], struct variable *vars, uint8_t nv) {

    char *xrange = (char *)malloc(30 * sizeof(char));
    char *yrange = (char *)malloc(30 * sizeof(char));
    sprintf(xrange, "set xrange [%f:%f]", vars[0].xl, vars[0].xu);
    sprintf(yrange, "set yrange [%f:%f]", vars[1].xl, vars[1].xu);

    gnuplot_setstyle(h, "lines");
    gnuplot_cmd(h, xrange);
    gnuplot_cmd(h, yrange);

    gnuplot_cmd(h, "set xlabel 'var 0'");
    gnuplot_cmd(h, "set ylabel 'var 1'");
    gnuplot_cmd(h, "set zlabel 'fitness'");
    gnuplot_cmd(h, "set ticslevel 0");
    gnuplot_cmd(h, "set hidden3d");
    gnuplot_cmd(h, "set isosample 70");
    gnuplot_cmd(h, "splot (1+cos(12*sqrt(x**2+y**2)))/(0.5*(x**2+y**2)+2), 'data/datapoints.csv' every::1" );


    sleep(SLEEP_LENGTH);
    gnuplot_resetplot(h);
}

/* Save datapoints to data/datapoints.data file */
void save_datapoints(struct chromosome p[POPULATION_SIZE], uint8_t nv) {
    /* Open or create data/datapoints.data for writing */
    FILE *fp = fopen("data/datapoints.csv", "w");
    /* Abort if connection is lost */
    if (!fp) {
        perror("Save Datapoints");
        exit(-1);
    }

    uint8_t i, j;
    for (i = 0; i < nv; i++)
        fprintf(fp, "var[%d], ", i);
    fprintf(fp, "fitness");
    fprintf(fp, "\n");

    for (i = 0; i < POPULATION_SIZE; i++) {
        /* save value[0], value[1], fitness to file.data */
        for (j = 0; j < nv; j++)
            fprintf(fp, "%lf,", p[i].values[j]);
        fprintf(fp, "%lf", p[i].fitness);
        fprintf(fp, "\n");
    }

    fclose(fp);
}

/* Evaluate population */
void eval_population(struct chromosome p[POPULATION_SIZE], struct variable *vars, uint8_t nv) {
    /* Get generation vars values */
    for (uint8_t i = 0; i < POPULATION_SIZE; i++)
        eval_vars(p[i], vars, nv);
    /* Get generation fitness */
    double sum_of_fitness = get_fitness(p);
    /* Get and print population's average fitness */
    double avg_fitness = sum_of_fitness /  POPULATION_SIZE;
    printf("Avg fitness: %lf\n", avg_fitness);
    /* Get generation probabilities */
    get_probabilites(p, sum_of_fitness);
}

/* Evaluate a single variable */
void eval_vars(struct chromosome chrome, struct variable *vars, uint8_t nv) {
    double decimal;
    for (uint8_t i = 0; i < nv; i++) {
        decimal = eval_bits(&chrome.genotype[vars[i].pos], vars[i].bits_len);

        /* xl + [(xu-xl)/(2^bl - 1)] * decimal */
        chrome.values[i] = vars[i].xl + ( (vars[i].xu - vars[i].xl) / (pow(2,vars[i].bits_len)-1) ) * decimal;
    }
}

/* Evaluate binary int array */
double eval_bits(uint8_t *bits, uint8_t len){
    double value = 0;
    for (int8_t i = len-1; i >= 0; i--)
        value += bits[i] * pow(2, len-i-1);

    return value;
}

/*  Evaluate fitness of each chromosome and return sum of fitness   */
double get_fitness(struct chromosome p[POPULATION_SIZE]) {
    double  sum = 0;
    /* Get fitness of each chromosome */
    for (uint16_t i = 0; i < POPULATION_SIZE; i++) {
        p[i].fitness = f(p[i].values);
        sum += p[i].fitness;
    }
    return sum;
}

/*  Get probability of selection and Cumulative probability distribution of each chromosome */
void get_probabilites(struct chromosome p[POPULATION_SIZE], double sof) {
    double sum_pr = 0;
    for (uint16_t i = 0; i < POPULATION_SIZE; i++) {
        p[i].probability    = p[i].fitness / sof;
        p[i].expected_p     = p[i].probability * POPULATION_SIZE;
        sum_pr              += p[i].probability;
        p[i].cdf            = sum_pr;
    }
}

/*  Create new generation */
void procreate(struct chromosome p[POPULATION_SIZE], struct chromosome new_gen[POPULATION_SIZE], uint8_t chrome_len, uint8_t nv) {

    struct chromosome *parent1, *parent2, child1, child2;
    /* Children memory allocation */
    child1.genotype = (uint8_t*)    calloc(chrome_len, sizeof(uint8_t));
    child1.values   = (double*)     calloc(nv, sizeof(double));
    child2.genotype = (uint8_t*)    calloc(chrome_len, sizeof(uint8_t));
    child2.values   = (double*)     calloc(nv, sizeof(double));

    uint16_t    i_fittest       = 0;        /* Fittest Chromosome's position    */
    double      fittest         = -DBL_MAX; /* Fittest Chromosome's fitness     */
    double      last_fittest    = DBL_MAX;  /* Last check's fittest Chromosome  */
    uint16_t    ngi             = 0;        /* Next Gen's index                 */

    for (uint16_t i = 0; i < ceil(POPULATION_SIZE * PASSTHROUGH_PCT); i++) {
        /* Select fittest chromosome */
        for (uint16_t k = 0; k < POPULATION_SIZE; k++)
            if (p[k].fitness > fittest && p[k].fitness < last_fittest) {
                i_fittest   = k;
                fittest     = p[k].fitness;
            }

        last_fittest = fittest;
        fittest = -DBL_MAX;

        /* Pass fittest chromosome straight to next gen */
        for (uint8_t j = 0; j < chrome_len; j++)
            new_gen[ngi].genotype[j] = p[i_fittest].genotype[j];

        ngi++;
    }

    while (ngi < POPULATION_SIZE) {
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
            new_gen[ngi].genotype[j]    = child1.genotype[j];
            new_gen[ngi+1].genotype[j]  = child2.genotype[j];
        }
        ngi = ngi+2;
    }
    /* Garbage collection */
    free(child1.genotype);
    free(child1.values);
    free(child2.genotype);
    free(child2.values);
}

/*  Select chromosomes for next generation through Roulette Wheel selection */
struct chromosome *select_parent(struct chromosome p[POPULATION_SIZE]) {
    /* Tournament Selection */
    return tournament(p);
    /* Roulette Shot Selection */
    //return roulette(p);
}

/*  Tournament Selection
    Get two random chromosomes
    Select the fittest if t_pct is met,
    Select weaker otherwise. */
struct chromosome *tournament(struct chromosome p[POPULATION_SIZE]) {
    /* Randomly select Chromosomes to fight */
    uint16_t    champ1  = floor(rng() * POPULATION_SIZE),
                champ2  = floor(rng() * POPULATION_SIZE),
                strong, weak;

    if (p[champ1].fitness > p[champ2].fitness) {
        strong  = champ1;
        weak    = champ2;
    } else {
        weak    = champ1;
        strong  = champ2;
    }

    double shot = rng();
    if (shot < TOURNEY_RATE)
        return &p[strong];
    else
        return &p[weak];
}

/*  Roulette Shot Selection
    Generate a random point value,
    get nearest chromosome with cdf > the shot */
struct chromosome *roulette(struct chromosome p[POPULATION_SIZE]) {
    double shot =   rng();
    uint16_t i  =   0;

    //printf("\nroulette: %f\n", roulette);
    while (shot > p[i].cdf)
        i++;
    //printf("p[%d]: x: %d\tfitness: %llu\tprobability: %f\tcdf: %f\t - %p\n", i, p[i].x, p[i].fitness, p[i].probability, p[i].cdf, &p[i]);
    return &p[i];
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

/* Print chromosome's genotype */
void print_genotype(uint8_t *genotype, uint8_t chrome_len) {
    for (uint8_t i = 0; i < chrome_len; i++) {
        putchar(genotype[i]);
    }
    printf("\n");
}

/* Print chromosome's relevant information */
void print_chromosome(struct chromosome chrome, uint8_t nv) {

    for (uint8_t i = 0; i < nv; i++) {
        printf("val%d: %lf\t", i, chrome.values[i]);
    }
    printf("fitness: %lf\tprob_selec: %lf\texpected_pop: %lf\n", chrome.fitness, chrome.probability, chrome.expected_p);
}

/*  De-Jong Drop-Wave Function
    (1+cos(12*sqrt(x^2+y^2))) / (0.5*(x^2+y^2)+2)*/
double f(double *values) {
    double x = values[0];
    double y = values[1];
    return ( 1 + cos( 12 * sqrt(pow(x,2) + pow(y,2)) ) ) / ( 0.5 * (pow(x,2) + pow(y,2)) + 2 );
}

/* Discrete Uniform Value.
    With a point value of x when 0 ≤ x ≤ 1, get an uniformly
    equivalent int value y as a ≤ y < b*/
uint16_t d_unif(double x, double a, double b) {
    return (x == 1) ? (uint16_t)b-1 : (uint16_t)floor(x*(b-a) + a );
}

/* ####### NOT USED ####### */
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
/* ######################### */

/*  Cryptographically Secure Pseudornadom Number Generator
    Read bytes from /dev/random device file. Each byte from the file is
    a cryptographically random value from 0-255. This function concatenates
    those bytes to generate a random number of an arbitrary size, only to then
    get a 0..1 point value out of it. */
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
