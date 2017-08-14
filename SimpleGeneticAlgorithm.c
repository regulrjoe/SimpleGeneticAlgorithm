#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include "../libraries/gnuplot_i/src/gnuplot_i.h"

/*
    PENDING
    TODO: Codificación Gray
    TODO: Implement Error Tolerance (EPSILON)
    TODO: Show plot for standard deviation, avg fitness, max fitness for each GENERATIONS
    TODO: Read a lil bit David Goldberg book Genetic Algorithms in Search
    TODO: Obtener variancia de la población
    TODO: Buscar mínimo o máximo
    TODO: Two-point crossover
    TODO: define_vars creates and returns pointer to vars instead of using pointer as argument
    TODO: add chrome_len attribute to chromosome struct
    TODO: Restructure Code into modal implementation.
    TODO: Producir bitácora de configuración y resultados por corrida.
    TODO: Función de selección Boltzmann
    TODO: Mutar por bit

    TODO: Estudiar programación lineal para entender restricciones
        "Programación lineal restricciones planteamiento de problemas"
*/

/* Constants */
#define POPULATION_S        50      /* SIZE OF THE POPULATION                       */
#define MUTATION_THRESH     0.2     /* PROBABILITY OF A MUTATION OCURRING           */
#define CROSSOVER_THRESH    0.8     /* PROBABILITY OF A CROSSOVER HAPPENING         */
#define PASSTHROUGH_THRESH  0.1     /* % OF FITTEST CHROMES WHO GO STRAIGHT THROUGH */
#define TOURNEY_THRESH      0.8     /* PROBABILITY OF SELECTING FITTEST OF ROUND    */


/* Enum to specify if it's looking for a min or a max optimum */
typedef enum {MIN, MAX} optimum_enum;



/* Variable structure */
typedef struct {
    uint8_t pos;        /* Position of given variable in genes array    */
    float   xl;         /* Lower limit of variable                      */
    float   xu;         /* Upper limit of variable                      */
    float   pr;         /* Precision of variable value                  */
    uint8_t bits_len;   /* Length of bits used for variable             */
} variable_t;

/* Chromosome structure */
typedef struct {
    uint8_t *genotype;      /* Chromosome bit array                     */
    int64_t phenotype;      /* Decimal value of Genotype                */
    double  fitness;        /* Result of f(phenotype)                   */
    double  probability;    /* Probability of selection                 */
    double  cdf;            /* Cumulative probability distribution      */
    double  expected_p;     /* Expected population                      */
    double  *values;        /* Decoded value of genes                   */
} chromosome_t;

/* Population structure */
typedef struct {
    chromosome_t    chromosomes[POPULATION_S];
    double          sum_of_fitness;
    double          avg_fitness;
    chromosome_t    *fittest_chromosome;
    uint8_t         chromosome_len;
    uint8_t         num_of_vars;
    variable_t      *vars;
    optimum_enum    optimum;
    double          std_deviation;
} population_t;



/* Functions */
uint8_t         define_vars         (variable_t *vars, uint8_t nv);
uint8_t         get_genotype_len    (float xl, float xu, float pr);
void            init_population     (population_t *p);
void            run                 (population_t *p, uint16_t iters);
void            eval_population     (population_t *p);
void            eval_vars           (population_t *p);
double          eval_bits           (uint8_t *bits, uint8_t len);
void            get_fitness         (population_t *p);
void            get_probabilites    (population_t *p);
void            procreate           (population_t *p, chromosome_t new_gen[POPULATION_S+1]);
chromosome_t    *select_parent      (population_t *p);
chromosome_t    *tournament         (population_t *p);
chromosome_t    *roulette           (population_t *p);
void            crossover           (uint8_t *parent1_geno, uint8_t *parent2_geno, chromosome_t *child1, chromosome_t *child2, uint8_t chrome_len);
void            mutate              (uint8_t *geno, uint8_t chrome_len);
void            print_chromosome    (chromosome_t chrome, uint8_t nv);
void            print_genotype      (uint8_t *genotype, uint8_t chrome_len);
void            save_datapoints     (population_t *p);
void            save_statistics     (chromosome_t fittest, uint16_t gen);
void            plot_avg_fitness    (gnuplot_ctrl *h, uint16_t iters);
void            plot_datapoints     (gnuplot_ctrl *h, variable_t *vars);
void            plot_fittest        (gnuplot_ctrl *h, uint16_t iters);
void            plot_std_deviation  (gnuplot_ctrl *h, uint16_t iters);
double          std_deviation       (population_t *p);
double          f                   (double *values);
double          rng                 ();
double          cud                 (double x, double a, double b);
double          c_unif              (double x, double a, double b);
uint16_t        d_unif              (double x, double a, double b);




int main() {
    population_t    p;          /* Population */
    uint16_t        iters; /* Generations */

    /*  UI */
    printf("\nWelcome to this Simple Genetic Algorithm Software written by Joe Vázquez-Mellado\n");
    printf("\nHow many variables do you wish to work with?: ");
    scanf("%d", &p.num_of_vars);

    p.vars = (variable_t*) calloc(p.num_of_vars, sizeof(variable_t));
    p.chromosome_len = define_vars(p.vars, p.num_of_vars);

    printf("\nEnter number of generations: ");
    scanf("%hd", &iters);

    /* Create starting population */
    init_population(&p);
    /* Run algorithm */
    run(&p, iters);

    return 0;
}



/*  Define the variables metainformation to use in chromosome */
uint8_t define_vars(variable_t *vars, uint8_t nv) {
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
void init_population(population_t *p) {
    for (int i = 0; i < POPULATION_S; i++){
        /* Memory allocation */
        p->chromosomes[i].genotype   = (uint8_t*) calloc(p->chromosome_len, sizeof(uint8_t));
        p->chromosomes[i].values     = (double*) calloc(p->num_of_vars, sizeof(double));

        /* Genotype intialization with random values 0-1 */
        for (uint8_t j = 0; j < p->chromosome_len; j++)
            p->chromosomes[i].genotype[j] = d_unif(rng(), 0, 2);

        /* Get genotype decimal value */
        p->chromosomes[i].phenotype = (int64_t)eval_bits(p->chromosomes[i].genotype, p->chromosome_len);
    }
}



/*  Run the genetic algorithm through a given population a given amount of times */
void run(population_t *p, uint16_t iters) {
    /* Initialize gnuplot handler */
    gnuplot_ctrl *h1 = gnuplot_init();
    gnuplot_ctrl *h2 = gnuplot_init();

    /* Memory allocation */
    chromosome_t new_gen[POPULATION_S+1];
    for (uint16_t i = 0; i < POPULATION_S+1; i++) {
        new_gen[i].genotype = (uint8_t*) calloc(p->chromosome_len, sizeof(uint8_t));
        new_gen[i].values   = (double*) calloc(p->num_of_vars, sizeof(double));
    }

    /* Run generations */
    for (uint16_t g = 0; g < iters; g++) {
        printf("\n\n----------------- GENERATION %d -----------------\n", g);
        /* Evaluate fitness, probabilities and var values */
        eval_population(p);
        /* Save to file and plot: datapoints, fittest chromosome, avg fitness & standard deviation */
        save_fittest(*p->fittest_chromosome, g, p->num_of_vars);
        save_datapoints(p);
        plot_fittest(h2, iters);
        plot_datapoints(h1, p->vars);

        /* Print chromosomes */
        for (uint8_t i = 0; i < POPULATION_S; i++) {
            print_chromosome(p->chromosomes[i], p->num_of_vars);
        }

        /* Get new generation */
        procreate(p, new_gen);

        /* Assign new generation's properties to current population */
        for (uint16_t j=0; j < POPULATION_S; j++) {
            /* assign new_gen genotype to current_gen */
            for (uint8_t k=0; k < p->chromosome_len; k++)
                p->chromosomes[j].genotype[k] = new_gen[j].genotype[k];
            /* Get new phenotype */
            p->chromosomes[j].phenotype = (int64_t)eval_bits(p->chromosomes[j].genotype, p->chromosome_len);
        }
    }
    gnuplot_close(h1);
    gnuplot_close(h2);
}



/* Evaluate population */
void eval_population(population_t *p) {
    /* Get generation vars values */
    eval_vars(p);
    /* Get generation fitness, avg fitness & sum of fitness */
    get_fitness(p);

    printf("Avg fitness: %lf\n", p->avg_fitness);
    /* Get generation probabilities */
    get_probabilites(p);
}



/* Evaluate a single variable */
void eval_vars(population_t *p) {
    double decimal;
    for (uint8_t i = 0; i < POPULATION_S; i++)
        for (uint8_t j = 0; j < p->num_of_vars; j++) {
            decimal = eval_bits(&p->chromosomes[i].genotype[p->vars[j].pos], p->vars[j].bits_len);

            /* xl + [(xu-xl)/(2^bl - 1)] * decimal */
            p->chromosomes[i].values[j] = p->vars[j].xl + ( (p->vars[j].xu - p->vars[j].xl) / (pow(2,p->vars[j].bits_len)-1) ) * decimal;
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
void get_fitness(population_t *p) {
    double      sum             = 0,        /* Sum of fitness                   */
                fittest         = -DBL_MAX; /* Fittest Chromosome's fitness     */
    uint16_t    i_fittest       = 0;        /* Fittest Chromosome's position    */

    /* Get fitness of each chromosome */
    for (uint16_t i = 0; i < POPULATION_S; i++) {
        p->chromosomes[i].fitness = f(p->chromosomes[i].values);
        sum += p->chromosomes[i].fitness;

        if (p->chromosomes[i].fitness > fittest) {
            i_fittest = i;
            fittest = p->chromosomes[i].fitness;
        }
    }
    p->fittest_chromosome = &p->chromosomes[i_fittest];
    p->sum_of_fitness = sum;
    p->avg_fitness = sum/POPULATION_S;
}



/*  Get probability of selection and Cumulative probability distribution of each chromosome */
void get_probabilites(population_t *p) {
    double sum_pr = 0;
    for (uint16_t i = 0; i < POPULATION_S; i++) {
        p->chromosomes[i].probability    = p->chromosomes[i].fitness / p->sum_of_fitness;
        p->chromosomes[i].expected_p     = p->chromosomes[i].probability * POPULATION_S;
        sum_pr                          += p->chromosomes[i].probability;
        p->chromosomes[i].cdf            = sum_pr;
    }
}



/*  Create new generation */
void procreate(population_t *p, chromosome_t new_gen[POPULATION_S+1]) {

    chromosome_t *parent1, *parent2, child1, child2;
    /* Children memory allocation */
    child1.genotype = (uint8_t*)    calloc(p->chromosome_len, sizeof(uint8_t));
    child1.values   = (double*)     calloc(p->num_of_vars, sizeof(double));
    child2.genotype = (uint8_t*)    calloc(p->chromosome_len, sizeof(uint8_t));
    child2.values   = (double*)     calloc(p->num_of_vars, sizeof(double));

    uint16_t    i_fittest       = 0,        /* Fittest Chromosome's position    */
                ngi             = 0;        /* Next Gen's index                 */
    double      fittest         = -DBL_MAX, /* Fittest Chromosome's fitness     */
                last_fittest    = DBL_MAX;  /* Last check's fittest Chromosome  */

    for (uint16_t i = 0; i < ceil(POPULATION_S * PASSTHROUGH_THRESH); i++) {
        /* Select fittest chromosome */
        for (uint16_t k = 0; k < POPULATION_S; k++) {
            if (p->chromosomes[k].fitness > fittest && p->chromosomes[k].fitness < last_fittest) {
                i_fittest = k;
                fittest = p->chromosomes[k].fitness;
            }
        }

        last_fittest = fittest;
        fittest = -DBL_MAX;

        /* Pass fittest chromosome straight to next gen */
        for (uint8_t j = 0; j < p->chromosome_len; j++)
            new_gen[ngi].genotype[j] = p->chromosomes[i_fittest].genotype[j];
        ngi++;
    }

    while (ngi < POPULATION_S) {
        parent1 = select_parent(p);
        parent2 = select_parent(p);
        for (uint8_t j = 0; j < p->chromosome_len; j++) {
            child1.genotype[j] = parent1->genotype[j];
            child2.genotype[j] = parent2->genotype[j];
        }
        /* Do crossover if CROSSOVER_THRESH is met */
        if (rng() <= CROSSOVER_THRESH)
            crossover(parent1->genotype, parent2->genotype, &child1, &child2, p->chromosome_len);
        /* Do mutation if MUTATION_THRESH is met */
        if (rng() <= MUTATION_THRESH)
            mutate(child1.genotype, p->chromosome_len);
        if (rng() <= MUTATION_THRESH )
            mutate(child2.genotype, p->chromosome_len);

        /* Assign children to new generation */
        for (uint8_t j = 0; j < p->chromosome_len; j++) {
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
chromosome_t *select_parent(population_t *p) {
    /* Tournament Selection */
    return tournament(p);
    /* Roulette Shot Selection */
    //return roulette(p);
}



/*  Tournament Selection
    Get two random chromosomes
    Select the fittest if t_pct is met,
    Select weaker otherwise. */
chromosome_t *tournament(population_t *p) {
    /* Randomly select Chromosomes to fight */
    uint16_t    champ1  = floor(rng() * POPULATION_S),
                champ2  = floor(rng() * POPULATION_S),
                strong, weak;

    if (p->chromosomes[champ1].fitness > p->chromosomes[champ2].fitness) {
        strong  = champ1;
        weak    = champ2;
    } else {
        weak    = champ1;
        strong  = champ2;
    }

    double shot = rng();
    if (shot < TOURNEY_THRESH)
        return &p->chromosomes[strong];
    else
        return &p->chromosomes[weak];
}



/*  Roulette Shot Selection
    Generate a random point value,
    get nearest chromosome with cdf > the shot */
chromosome_t *roulette(population_t *p) {
    double shot =   rng();
    uint16_t i  =   0;

    //printf("\nroulette: %f\n", roulette);
    while (shot > p->chromosomes[i].cdf)
        i++;
    //printf("p[%d]: x: %d\tfitness: %llu\tprobability: %f\tcdf: %f\t - %p\n", i, p[i].x, p[i].fitness, p[i].probability, p[i].cdf, &p[i]);
    return &p->chromosomes[i];
}



/*  Crossover Operator.
    Get two new offsprings by crossing two given parents at a single point */
void crossover(uint8_t *parent1_geno, uint8_t *parent2_geno, chromosome_t *child1, chromosome_t *child2, uint8_t chrome_len) {
    /*  Get random cutpoint's locus cutpoint >= 1 && cutpoint < chromosome's size */
    uint8_t cutpoint = d_unif(rng(), 1, chrome_len);

    for (uint8_t i = cutpoint; i < chrome_len; i++) {
        child1->genotype[i] = parent2_geno[i];
        child2->genotype[i] = parent1_geno[i];
    }
}



/* Mutation Operator.
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
void print_chromosome(chromosome_t chrome, uint8_t nv) {

    for (uint8_t i = 0; i < nv; i++) {
        printf("val%d: %lf\t", i, chrome.values[i]);
    }
    printf("fitness: %lf\tprob_selec: %lf\texpected_pop: %lf\n", chrome.fitness, chrome.probability, chrome.expected_p);
}



/* Save datapoints to data/datapoints.csv file */
void save_datapoints(population_t *p) {
    /* Open or create data/datapoints.csv for writing */
    FILE *fp = fopen("data/datapoints.csv", "w");
    /* Abort if connection is lost */
    if (!fp) {
        perror("Save Datapoints");
        exit(-1);
    }

    uint8_t i, j;
    for (i = 0; i < p->num_of_vars; i++)
        fprintf(fp, "var[%d], ", i);
    fprintf(fp, "fitness");
    fprintf(fp, "\n");

    for (i = 0; i < POPULATION_S; i++) {
        /* save value[0], value[1], fitness to file.data */
        for (j = 0; j < p->num_of_vars; j++)
            fprintf(fp, "%lf, ", p->chromosomes[i].values[j]);
        fprintf(fp, "%lf", p->chromosomes[i].fitness);
        fprintf(fp, "\n");
    }
    fclose(fp);
}



/* Save per-generation statistics to data/stats.csv file
    • fittest chromosome's fitness and values
    • average fitness
    • standard deviation */
void save_statistics(population_t *p, uint16_t gen) {
    FILE *fp;
    uint8_t i;

    /*  If first generation open file for writing and add headers
        Otherwise open file for appending */
    if (gen == 0) {
        /* Open or create data/stats.csv for writing
            Abort if connection is lost */
        if (!(fp = fopen("data/stats.csv", "w")) ) {
            perror("Save Generation [%d]'s statistics", gen);
            exit(-1);
        }

        fprintf(fp, "generation, ");
        fprintf(fp, "fittest fitness, ");
        for (i = 0; i < p->num_of_vars; i++) {
            fprintf(fp, "fittest val[%d], ",i);
        }

        fprintf(fp, "avg_fitness, ");
        fprintf(fp, "std_deviation");
        fprintf(fp, "\n");
    } else {
        /* Open data/stats.csv for appending */
        fp = fopen("data/generations_stats.csv", "a");
        /* Abort if connection is lost */
        if (!fp) {
            perror("Save Generation [%d]'s statistics", gen);
            exit(-1);
        }
    }

    fprintf(fp, "%d, ", gen);   /* save generation */
    fprintf(fp, "%lf, ", p->fittest_chromosome.fitness); /* save fittest chromosome's fitness */
    for (i = 0; i < p->num_of_vars; i++) {
        fprintf(fp, "%lf, ", p->fittest_chromosome.values[i]);  /* save fittest_chromosome's value[i] */
    }
    fprintf(fp, "%lf, ", p->avg_fitness);
    fprintf(fp, "%lf", p->std_deviation);
    fprintf(fp, "\n");
    fclose(fp);
}



/* Plot average fitness of population */
void plot_avgfitness(gnuplot_ctrl *h, uint16_t iters) {}



/* Plot equation surface and its datapoints on top of it */
void plot_datapoints(gnuplot_ctrl *h, variable_t *vars) {

    char *xrange = (char *)malloc(30 * sizeof(char));
    char *yrange = (char *)malloc(30 * sizeof(char));
    sprintf(xrange, "set xrange [%f:%f]", vars[0].xl, vars[0].xu);
    sprintf(yrange, "set yrange [%f:%f]", vars[1].xl, vars[1].xu);

    gnuplot_cmd(h, "set title 'Population Fitness'");
    gnuplot_setstyle(h, "lines");
    gnuplot_cmd(h, xrange);
    gnuplot_cmd(h, yrange);
    gnuplot_cmd(h, "set xlabel 'var 0'");
    gnuplot_cmd(h, "set ylabel 'var 1'");
    gnuplot_cmd(h, "set zlabel 'fitness'");
    gnuplot_cmd(h, "set ticslevel 0");
    gnuplot_cmd(h, "set hidden3d");
    gnuplot_cmd(h, "set isosample 70");
    /* Plot fitness function surface, and datapoints from data/datapoints.cvs */
    gnuplot_cmd(h, "splot (1+cos(12*sqrt(x**2+y**2)))/(0.5*(x**2+y**2)+2), 'data/datapoints.csv' every::1" );

    gnuplot_resetplot(h);
}



/* Plot fittest chromosome of population */
void plot_fittest(gnuplot_ctrl *h, uint16_t iters) {

    char *xrange = (char *)malloc(30 * sizeof(char));
    sprintf(xrange, "set xrange[0:%d]", iters);

    gnuplot_cmd(h, "set title 'Fittest Chromosome'");
    gnuplot_setstyle(h, "lines");
    gnuplot_cmd(h, xrange);
    gnuplot_cmd(h, "set yrange [0:]");
    gnuplot_cmd(h, "set xlabel 'generation'");
    gnuplot_cmd(h, "set ylabel 'highest fitness'");
    gnuplot_cmd(h, "plot 'data/fittest.csv' using 1:2 with lines");

    gnuplot_resetplot(h);

}



/* Plot standard deviation of population */
void plot_stddev(gnuplot_ctrl *h, uint16_t iters) {}



/* De-Jong Drop-Wave Function
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
