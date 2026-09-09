/*
 Trait-based coevolutionary dynamics on the empirical antagonistic networks
 published in the /networks directory of this repository (network_001 ...
 network_122, see networks/README.md).

 adjacency_matrix.csv convention (see networks/README.md):
   - rows    = exploiter species (M_Exploiters)
   - columns = victim species (M_Victims)
   - f_ij = 1 if exploiter i interacts with victim j, 0 otherwise

 Compile:
   gcc -o coevolution_antagonistic_networks coevolution_antagonistic_networks.c -O3 -lm -Wall -Wextra

 Run (from the repository root, so that ./networks/ and ./Files/ resolve):
   ./coevolution_antagonistic_networks net_id flag_temporal

   net_id        -> network identifier, 1..122 (matches networks/network_<net_id padded to 3 digits>)
   flag_temporal -> 0: write only the final (t -> TMAX) state of each simulation
                     1: write the full time series of every simulation

 Output directories `files/networks/` and `files/networks/temporal/` must
 exist before running (this program does not create them).
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

double PI = 3.14159265;

/*######################################################

RAN3 ALGORITHM: RANDOM NUMBER GENERATOR BETWEEN (0,1)

########################################################*/

#define MBIG 1000000000
#define MZ 0
#define FAC (1.0 / MBIG)
int MSEED;
long idum = -1;
// Any large MBIG, and any smaller (but still large) MSEED can be substituted for the above values.

void ini_ran()
{
    // MSEED = (int)(time(NULL));
    MSEED = 1234567890;
}

double Random()
{
    // Set idum to any negative value to initialize or reinitialize the sequence.
    static int inext, inextp;
    static long ma[56]; // The value 56 (range ma[1..55]) is special and should not be modified; see Knuth.
    static int iff = 0;
    long mj, mk;
    int i, ii, k;

    if (idum < 0 || iff == 0)
    { // Initialization
        iff = 1;
        mj = labs(MSEED - labs(idum)); // Initialize ma[55] using the seed idum and the large number MSEED.
        mj %= MBIG;
        ma[55] = mj;
        mk = 1;

        for (i = 1; i <= 54; i++)
        { // Now initialize the rest of the table in a slightly random order with numbers that are not especially random.
            ii = (21 * i) % 55;
            ma[ii] = mk;
            mk = mj - mk;
            if (mk < MZ)
                mk += MBIG;
            mj = ma[ii];
        }
        for (k = 1; k <= 4; k++) // We randomize them by "warming up the generator"
            for (i = 1; i <= 55; i++)
            {
                ma[i] -= ma[1 + (i + 30) % 55];
                if (ma[i] < MZ)
                    ma[i] += MBIG;
            }

        inext = 0;   // Prepare indices for our first generated number.
        inextp = 31; // The constant 31 is special; see Knuth.
        idum = 1;
    }

    // Here is where we start, except on initialization.
    if (++inext == 56)
        inext = 1; // Increment inext and inextp, wrapping around 56 to 1.
    if (++inextp == 56)
        inextp = 1;
    mj = ma[inext] - ma[inextp]; // Generate a new random number subtractively.
    if (mj < MZ)
        mj += MBIG; // Be sure that it is in range.
    ma[inext] = mj; // Store it and output the derived uniform deviate.
    return mj * FAC;
}

/*######################################################

BOX-MULLER ALGORITHM: RANDOM NUMBER GENERATOR WITH NORMAL DISTRIBUTION

########################################################*/

void Box_Muller(double *g1, double *g2)
{
    double d1, d2;
    d1 = sqrt(-2.0 * log(Random()));
    d2 = 2.0 * PI * Random();
    *g1 = -d1 * cos(d2);
    *g2 = -d1 * sin(d2);
}
/*################################################

SPECIES STRUCTS DEFINITIONS, WITH ALL VARIABLES USED

##################################################*/
//

typedef struct
{
    int *w_E;             // indices of the neighbouring exploiters, size = k_victim
    int k_victim;         // number of connections in the network
    double z_trait;       // value of trait
    double theta;         // stabilizing selection, theta = z_trait(0)
    double xi_s;          // intensity of environmental selection
    double xi_d;          // intensity of interaction patterns
    double S_i;           // partial selection by environment
    double *M_ij;         // partial selection by interaction patterns, size = k_victim
    double fitness_V;     // fitness of the species
    double fitness_amb_V; // ambiental term of fitness
    double fitness_int_V; // interaction term of fitness
} Victim;                 //

typedef struct
{
    int *w_V;             // indices of the neighbouring victims, size = k_exploiter
    int k_exploiter;      // number of connections in the network
    double z_trait;       // value of trait
    double theta;         // stabilizing selection, theta = z_trait(0)
    double xi_s;          // intensity of environmental selection
    double xi_d;          // intensity of interaction patterns
    double S_i;           // environmental selection
    double *M_ij;         // partial selection by interaction patterns, size = k_exploiter
    double fitness_E;     // fitness of the species
    double fitness_amb_E; // ambiental term of fitness
    double fitness_int_E; // interaction term of fitness
} Exploiter;              //

/**####################

FUNCTION DECLARATION

#######################*/

int ini_species();              // reads networks/network_<net_id>/adjacency_matrix.csv and builds the species lists
int ini_species_distribution(); // initial distribution of species' trait each simulation
void Runge_Kutta(double time);
double interaction_function_exploiters(Exploiter *exploiters_net, Victim *victims_net, int id_exploiter, double time); // interaction of the agents
double interaction_function_victims(Exploiter *exploiters_net, Victim *victims_net, int id_victim, double time);       // interaction of the agents
void measure_fitness();

void write_header_global_temporal(char name_1[200], char name_2[200]); // writes a header with important data in the output files
void write_header_global_final(char name_1[200]);                      // writes a header with important data in the output files
void write_data_time(int simul, double t_iter);                        // writes data about temporal information
void write_data_final(int simul);                                      // writes data about temporal information

/**###########################

GLOBAL VARIABLES DEFINITION

##############################*/

int net_id, flag_temporal; // input variables
double xi_d_V, xi_d_E;     // input variables
double xi_d_V_min, xi_d_V_max, xi_d_V_delta;
double xi_d_E_min, xi_d_E_max, xi_d_E_delta;

int N_Sim;      // number of simulation
int N_Max_Iter; // max iteration of each simulation

int M_Victims;    // number of victims in the network (columns of adjacency_matrix.csv)
int M_Exploiters; // number of exploiters in the network (rows of adjacency_matrix.csv)

double h, TMAX; // time step - TMAX
int write_time_window; // writing time window
double epsilon; // evolutionary response on victim species
double alpha;   // exploiter preference
double delta;   // integration constant to avoid dividing by 0

FILE *f_out_victims, *f_out_exploiters, *f_out_final;

Victim *Victim_Species_List;       // List of species with all the information
Exploiter *Exploiter_Species_List; // List of species with all the information

/**############

MAIN PROGRAM

###############*/

int main(int argc, char **argv)
{

    clock_t begin = clock();

    ini_ran();
    Random(); // initialization of the RNG

    /**----------------------------------------------------------

                    CHANGE WHEN NECESSARY

    ----------------------------------------------------------*/

    N_Sim = 100;

    alpha = 0.1;
    TMAX = 300;
    delta = 0.000001;
    h = 0.01;
    write_time_window = 50;
    N_Max_Iter = (int)TMAX / h;

    xi_d_V_min = 0.1;
    xi_d_V_max = 0.9;
    xi_d_V_delta = 0.1;

    xi_d_E_min = 0.1;
    xi_d_E_max = 0.9;
    xi_d_E_delta = 0.1;

    /**----------------------------------------------------------

    Read the input arguments:
        - net_id: network identifier (1..122), matches networks/network_<net_id padded to 3 digits>
        - flag_temporal: 0 (final state only) or 1 (full time series)

    ----------------------------------------------------------*/

    if (argc < 3)
    {
        fprintf(stderr, "Usage: %s net_id flag_temporal\n", argv[0]);
        exit(1);
    }

    double variables[2];
    int i;
    for (i = 1; i < argc; i++)
    {
        sscanf(argv[i], " %lf", &variables[i - 1]);
    }

    net_id = (int)variables[0];
    flag_temporal = (int)variables[1];

    /*----------------------------------------------------------

    SPECIES NETWORK INITIALIZATION
    (reads networks/network_<net_id>/adjacency_matrix.csv and sets
     M_Exploiters, M_Victims, Exploiter_Species_List, Victim_Species_List)

    ----------------------------------------------------------*/

    ini_species();

    /*----------------------------------------------------------

    SIMULATION

    ----------------------------------------------------------*/

    if (flag_temporal == 1)
    {
        char filename_victims[200], filename_exploiters[200];
        sprintf(filename_victims, "files/networks/temporal/network_%03d_temporal_V_alpha=%.1lf.csv", net_id, alpha);
        sprintf(filename_exploiters, "files/networks/temporal/network_%03d_temporal_E_alpha=%.1lf.csv", net_id, alpha);
        write_header_global_temporal(filename_victims, filename_exploiters);
    }
    else
    {
        char filename_final[200];

        sprintf(filename_final, "files/networks/network_%03d_final_alpha=%.1lf.csv", net_id, alpha);

        write_header_global_final(filename_final);
    }

    for (xi_d_V = xi_d_V_min; xi_d_V < xi_d_V_max + xi_d_V_delta / 10.; xi_d_V += xi_d_V_delta)
    {
        for (xi_d_E = xi_d_E_min; xi_d_E < xi_d_E_max + xi_d_E_delta / 10.; xi_d_E += xi_d_E_delta)
        {
            int sim, iter;

            for (sim = 0; sim < N_Sim; sim++)
            {

                ini_species_distribution();
                measure_fitness();

                iter = 0;
                if (flag_temporal == 1)
                {
                    write_data_time(sim, iter * h);
                }

                do
                {
                    Runge_Kutta(iter * h);
                    measure_fitness();

                    iter++;
                    if (flag_temporal == 1 && iter % write_time_window == 0)
                    {
                        write_data_time(sim, iter * h);
                    }

                } while (iter < N_Max_Iter);

                if (flag_temporal == 0)
                {
                    write_data_final(sim);
                }

                if (flag_temporal == 1)
                {
                    fprintf(f_out_victims, "\n");
                    fprintf(f_out_exploiters, "\n");
                }
            }
        }
    }

    if (flag_temporal == 1)
    {
        fclose(f_out_victims);
        fclose(f_out_exploiters);
    }

    if (flag_temporal == 0)
    {
        fclose(f_out_final);
    }

    /* ------------------------------------- */

    int i_victims, i_exploiters;

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        free(Victim_Species_List[i_victims].w_E);
        free(Victim_Species_List[i_victims].M_ij);
    }
    free(Victim_Species_List);

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        free(Exploiter_Species_List[i_exploiters].w_V);
        free(Exploiter_Species_List[i_exploiters].M_ij);
    }
    free(Exploiter_Species_List);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("Time of simulation: %.2lf\n", time_spent);

    return 0;
}

/**########################

FUNCTION DEFINITIONS

#########################*/

int ini_species()
{
    /*
     * Loads networks/network_<net_id>/adjacency_matrix.csv (rows = exploiters,
     * columns = victims, see networks/README.md) and builds Exploiter_Species_List
     * and Victim_Species_List from it, mirroring the neighbour-list format
     * expected by the rest of the simulation (w_V / w_E and k_exploiter / k_victim).
     */

    char net_filename[200];
    sprintf(net_filename, "networks/network_%03d/adjacency_matrix.csv", net_id);

    FILE *matrix_file = fopen(net_filename, "r");
    if (matrix_file == NULL)
    {
        fprintf(stderr, "Error: could not open network file %s\n", net_filename);
        exit(1);
    }

    /* First pass: number of rows (M_Exploiters) and columns (M_Victims) */
    char line[8192];
    M_Exploiters = 0;
    M_Victims = 0;

    while (fgets(line, sizeof(line), matrix_file) != NULL)
    {
        if (line[0] == '\n' || line[0] == '\0')
            continue;

        if (M_Exploiters == 0)
        {
            char line_copy[8192];
            strcpy(line_copy, line);
            char *tok = strtok(line_copy, ",\n");
            while (tok != NULL)
            {
                M_Victims++;
                tok = strtok(NULL, ",\n");
            }
        }
        M_Exploiters++;
    }

    /* Second pass: fill in the dense adjacency matrix */
    int **adjacency = (int **)malloc(M_Exploiters * sizeof(int *));
    int i_exploiters, i_victims;
    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
        adjacency[i_exploiters] = (int *)malloc(M_Victims * sizeof(int));

    rewind(matrix_file);
    i_exploiters = 0;
    while (fgets(line, sizeof(line), matrix_file) != NULL)
    {
        if (line[0] == '\n' || line[0] == '\0')
            continue;

        char *tok = strtok(line, ",\n");
        i_victims = 0;
        while (tok != NULL && i_victims < M_Victims)
        {
            adjacency[i_exploiters][i_victims] = atoi(tok);
            i_victims++;
            tok = strtok(NULL, ",\n");
        }
        i_exploiters++;
    }
    fclose(matrix_file);

    /* Build the victim species list (columns of the adjacency matrix) */
    Victim_Species_List = (Victim *)malloc(M_Victims * sizeof(Victim));

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        int k = 0;
        for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
            if (adjacency[i_exploiters][i_victims])
                k++;

        Victim_Species_List[i_victims].k_victim = k;
        Victim_Species_List[i_victims].w_E = (int *)calloc(k, sizeof(int));
        Victim_Species_List[i_victims].M_ij = (double *)calloc(k, sizeof(double));

        int idx = 0;
        for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
            if (adjacency[i_exploiters][i_victims])
                Victim_Species_List[i_victims].w_E[idx++] = i_exploiters;
    }

    /* Build the exploiter species list (rows of the adjacency matrix) */
    Exploiter_Species_List = (Exploiter *)malloc(M_Exploiters * sizeof(Exploiter));

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        int k = 0;
        for (i_victims = 0; i_victims < M_Victims; i_victims++)
            if (adjacency[i_exploiters][i_victims])
                k++;

        Exploiter_Species_List[i_exploiters].k_exploiter = k;
        Exploiter_Species_List[i_exploiters].w_V = (int *)calloc(k, sizeof(int));
        Exploiter_Species_List[i_exploiters].M_ij = (double *)calloc(k, sizeof(double));

        int idx = 0;
        for (i_victims = 0; i_victims < M_Victims; i_victims++)
            if (adjacency[i_exploiters][i_victims])
                Exploiter_Species_List[i_exploiters].w_V[idx++] = i_victims;
    }

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
        free(adjacency[i_exploiters]);
    free(adjacency);

    return 0;
}

int ini_species_distribution()
{
    double gauss1, gauss2;
    gauss1 = 0.;
    gauss2 = 0.;

    int i_victims, i_exploiters;
    int idx_victims, idx_exploiters;

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        if (i_victims % 2 != 0)
        {
            Victim_Species_List[i_victims].z_trait = 0.5 * gauss1;
        }
        else
        {
            Box_Muller(&gauss1, &gauss2);

            Victim_Species_List[i_victims].z_trait = 0.5 * gauss2;
        }

        Victim_Species_List[i_victims].theta = Victim_Species_List[i_victims].z_trait;
        Victim_Species_List[i_victims].S_i = 0.;
        Victim_Species_List[i_victims].xi_d = xi_d_V;
        Victim_Species_List[i_victims].xi_s = 1. - Victim_Species_List[i_victims].xi_d;

        for (idx_exploiters = 0; idx_exploiters < Victim_Species_List[i_victims].k_victim; idx_exploiters++)
        {
            Victim_Species_List[i_victims].M_ij[idx_exploiters] = 0.;
        }
    }

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        if (i_exploiters % 2 != 0)
        {
            Exploiter_Species_List[i_exploiters].z_trait = 0.5 * gauss1;
        }
        else
        {
            Box_Muller(&gauss1, &gauss2);

            Exploiter_Species_List[i_exploiters].z_trait = 0.5 * gauss2;
        }

        Exploiter_Species_List[i_exploiters].theta = Exploiter_Species_List[i_exploiters].z_trait;
        Exploiter_Species_List[i_exploiters].S_i = 0.;
        Exploiter_Species_List[i_exploiters].xi_d = xi_d_E;
        Exploiter_Species_List[i_exploiters].xi_s = 1. - Exploiter_Species_List[i_exploiters].xi_d;

        for (idx_victims = 0; idx_victims < Exploiter_Species_List[i_exploiters].k_exploiter; idx_victims++)
        {
            Exploiter_Species_List[i_exploiters].M_ij[idx_victims] = 0.;
        }
    }

    return 0;
}

void Runge_Kutta(double time)
{

    /* Runge Kutta Algorithm */
    Exploiter *aux_exploiters_net;
    Victim *aux_victims_net;

    aux_exploiters_net = (Exploiter *)malloc(M_Exploiters * sizeof(Exploiter));
    aux_victims_net = (Victim *)malloc(M_Victims * sizeof(Victim));

    int i_exploiters, i_victims;

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        aux_exploiters_net[i_exploiters].z_trait = Exploiter_Species_List[i_exploiters].z_trait;
        aux_exploiters_net[i_exploiters].k_exploiter = Exploiter_Species_List[i_exploiters].k_exploiter;
        aux_exploiters_net[i_exploiters].theta = Exploiter_Species_List[i_exploiters].theta;
        aux_exploiters_net[i_exploiters].xi_d = Exploiter_Species_List[i_exploiters].xi_d;
        aux_exploiters_net[i_exploiters].xi_s = Exploiter_Species_List[i_exploiters].xi_s;
        aux_exploiters_net[i_exploiters].w_V = (int *)malloc(aux_exploiters_net[i_exploiters].k_exploiter * sizeof(int));
        aux_exploiters_net[i_exploiters].M_ij = (double *)malloc(aux_exploiters_net[i_exploiters].k_exploiter * sizeof(double));
        int idx_victims;
        for (idx_victims = 0; idx_victims < aux_exploiters_net[i_exploiters].k_exploiter; idx_victims++)
        {
            aux_exploiters_net[i_exploiters].w_V[idx_victims] = Exploiter_Species_List[i_exploiters].w_V[idx_victims];
        }
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        aux_victims_net[i_victims].z_trait = Victim_Species_List[i_victims].z_trait;
        aux_victims_net[i_victims].k_victim = Victim_Species_List[i_victims].k_victim;
        aux_victims_net[i_victims].theta = Victim_Species_List[i_victims].theta;
        aux_victims_net[i_victims].xi_d = Victim_Species_List[i_victims].xi_d;
        aux_victims_net[i_victims].xi_s = Victim_Species_List[i_victims].xi_s;
        aux_victims_net[i_victims].w_E = (int *)malloc(aux_victims_net[i_victims].k_victim * sizeof(int));
        aux_victims_net[i_victims].M_ij = (double *)malloc(aux_victims_net[i_victims].k_victim * sizeof(double));
        int idx_exploiters;
        for (idx_exploiters = 0; idx_exploiters < aux_victims_net[i_victims].k_victim; idx_exploiters++)
        {
            aux_victims_net[i_victims].w_E[idx_exploiters] = Victim_Species_List[i_victims].w_E[idx_exploiters];
        }
    }

    double *k1_exploiters, *k1_victims;
    double *k2_exploiters, *k2_victims;
    double *k3_exploiters, *k3_victims;
    double *k4_exploiters, *k4_victims;

    k1_exploiters = (double *)malloc(M_Exploiters * sizeof(double));
    k2_exploiters = (double *)malloc(M_Exploiters * sizeof(double));
    k3_exploiters = (double *)malloc(M_Exploiters * sizeof(double));
    k4_exploiters = (double *)malloc(M_Exploiters * sizeof(double));
    k1_victims = (double *)malloc(M_Victims * sizeof(double));
    k2_victims = (double *)malloc(M_Victims * sizeof(double));
    k3_victims = (double *)malloc(M_Victims * sizeof(double));
    k4_victims = (double *)malloc(M_Victims * sizeof(double));

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        k1_exploiters[i_exploiters] = h * (interaction_function_exploiters(aux_exploiters_net, aux_victims_net, i_exploiters, time));
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        k1_victims[i_victims] = h * (interaction_function_victims(aux_exploiters_net, aux_victims_net, i_victims, time));
    }

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        aux_exploiters_net[i_exploiters].z_trait = Exploiter_Species_List[i_exploiters].z_trait + 0.5 * k1_exploiters[i_exploiters];
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        aux_victims_net[i_victims].z_trait = Victim_Species_List[i_victims].z_trait + 0.5 * k1_victims[i_victims];
    }

    // obtain k2 for the new networks

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        k2_exploiters[i_exploiters] = h * (interaction_function_exploiters(aux_exploiters_net, aux_victims_net, i_exploiters, time));
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        k2_victims[i_victims] = h * (interaction_function_victims(aux_exploiters_net, aux_victims_net, i_victims, time));
    }

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        aux_exploiters_net[i_exploiters].z_trait = Exploiter_Species_List[i_exploiters].z_trait + 0.5 * k2_exploiters[i_exploiters];
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        aux_victims_net[i_victims].z_trait = Victim_Species_List[i_victims].z_trait + 0.5 * k2_victims[i_victims];
    }

    // obtain k3 for the new networks

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        k3_exploiters[i_exploiters] = h * (interaction_function_exploiters(aux_exploiters_net, aux_victims_net, i_exploiters, time));
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        k3_victims[i_victims] = h * (interaction_function_victims(aux_exploiters_net, aux_victims_net, i_victims, time));
    }

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        aux_exploiters_net[i_exploiters].z_trait = Exploiter_Species_List[i_exploiters].z_trait + k3_exploiters[i_exploiters];
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        aux_victims_net[i_victims].z_trait = Victim_Species_List[i_victims].z_trait + k3_victims[i_victims];
    }

    // obtain k4 for the new networks

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        k4_exploiters[i_exploiters] = h * (interaction_function_exploiters(aux_exploiters_net, aux_victims_net, i_exploiters, time));
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        k4_victims[i_victims] = h * (interaction_function_victims(aux_exploiters_net, aux_victims_net, i_victims, time));
    }

    // update trait values for h + 1

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        Exploiter_Species_List[i_exploiters].z_trait = Exploiter_Species_List[i_exploiters].z_trait + 1. / 6 * (k1_exploiters[i_exploiters] + 2. * k2_exploiters[i_exploiters] + 2. * k3_exploiters[i_exploiters] + k4_exploiters[i_exploiters]);
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        Victim_Species_List[i_victims].z_trait = Victim_Species_List[i_victims].z_trait + 1. / 6 * (k1_victims[i_victims] + 2. * k2_victims[i_victims] + 2. * k3_victims[i_victims] + k4_victims[i_victims]);
    }

    free(k1_exploiters);
    free(k2_exploiters);
    free(k3_exploiters);
    free(k4_exploiters);
    free(k1_victims);
    free(k2_victims);
    free(k3_victims);
    free(k4_victims);

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        free(aux_exploiters_net[i_exploiters].w_V);
        free(aux_exploiters_net[i_exploiters].M_ij);
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        free(aux_victims_net[i_victims].w_E);
        free(aux_victims_net[i_victims].M_ij);
    }

    free(aux_exploiters_net);
    free(aux_victims_net);
}

double interaction_function_exploiters(Exploiter *exploiters_net, Victim *victims_net, int id_exploiter, double time)
{
    /*
    Function which describes the function of the exploiters
    */

    int idx_victims;

    exploiters_net[id_exploiter].S_i = exploiters_net[id_exploiter].xi_s * (exploiters_net[id_exploiter].theta - exploiters_net[id_exploiter].z_trait);
    double aux_denominator_Ei = 0.;

    for (idx_victims = 0; idx_victims < exploiters_net[id_exploiter].k_exploiter; idx_victims++) // loop for the denominator
    {
        int id_victim;
        id_victim = exploiters_net[id_exploiter].w_V[idx_victims];

        double difference = victims_net[id_victim].z_trait - exploiters_net[id_exploiter].z_trait;

        aux_denominator_Ei += exp(-alpha * difference * difference);
    }

    double aux_selection_sum = 0.;
    for (idx_victims = 0; idx_victims < exploiters_net[id_exploiter].k_exploiter; idx_victims++)
    {
        double aux_p_ij = 0.;
        int id_victim;
        id_victim = exploiters_net[id_exploiter].w_V[idx_victims];

        double difference = victims_net[id_victim].z_trait - exploiters_net[id_exploiter].z_trait + delta;
        double abs_difference = fabs(difference);

        aux_p_ij = exp(-alpha * difference * difference);
        exploiters_net[id_exploiter].M_ij[idx_victims] = aux_p_ij;

        aux_selection_sum += exploiters_net[id_exploiter].M_ij[idx_victims] * difference / abs_difference / (aux_denominator_Ei + delta);
    }

    return exploiters_net[id_exploiter].S_i + exploiters_net[id_exploiter].xi_d * aux_selection_sum; // equation for exploiters
}

double interaction_function_victims(Exploiter *exploiters_net, Victim *victims_net, int id_victim, double time)
{

    int idx_exploiters;

    victims_net[id_victim].S_i = victims_net[id_victim].xi_s * (victims_net[id_victim].theta - victims_net[id_victim].z_trait);

    double aux_denominator_Vi = 0.;

    for (idx_exploiters = 0; idx_exploiters < victims_net[id_victim].k_victim; idx_exploiters++) // loop for the denominator
    {
        int id_exploiter;
        id_exploiter = victims_net[id_victim].w_E[idx_exploiters];

        double difference = victims_net[id_victim].z_trait - exploiters_net[id_exploiter].z_trait;

        aux_denominator_Vi += exp(-alpha * difference * difference);
    }

    double aux_selection_sum = 0.;

    for (idx_exploiters = 0; idx_exploiters < victims_net[id_victim].k_victim; idx_exploiters++)
    {
        int id_exploiter;
        id_exploiter = victims_net[id_victim].w_E[idx_exploiters];
        double difference = victims_net[id_victim].z_trait - exploiters_net[id_exploiter].z_trait + delta;
        double abs_difference = fabs(difference);

        victims_net[id_victim].M_ij[idx_exploiters] = exp(-alpha * difference * difference);
        aux_selection_sum += victims_net[id_victim].M_ij[idx_exploiters] * difference / abs_difference / (aux_denominator_Vi + delta);
    }

    return victims_net[id_victim].S_i + victims_net[id_victim].xi_d * aux_selection_sum;
}

void measure_fitness()
{
    int i_exploiters, i_victims, idx_victims, idx_exploiters;
    double fitness_amb_E, fitness_int_E, diff_amb_E, diff_sel_E, diff_amb_EV;
    double fitness_amb_V, fitness_int_V, diff_amb_V, diff_sel_V;
    double max_xi_E = xi_d_E / (1. - xi_d_E);
    double max_xi_V = xi_d_V / (1. - xi_d_V);

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        diff_amb_E = fabs(Exploiter_Species_List[i_exploiters].z_trait - Exploiter_Species_List[i_exploiters].theta);

        fitness_amb_E = (1. - diff_amb_E / max_xi_E);
        int k_E = Exploiter_Species_List[i_exploiters].k_exploiter;
        fitness_int_E = 0.;
        for (idx_victims = 0; idx_victims < k_E; idx_victims++)
        {
            int id_victim;
            id_victim = Exploiter_Species_List[i_exploiters].w_V[idx_victims];
            diff_amb_EV = fabs(Victim_Species_List[id_victim].theta - Exploiter_Species_List[i_exploiters].theta);
            diff_sel_E = fabs(Victim_Species_List[id_victim].z_trait - Exploiter_Species_List[i_exploiters].z_trait);
            fitness_int_E += 1. / k_E * (1. - diff_sel_E / (max_xi_E + max_xi_V + diff_amb_EV));
        }
        Exploiter_Species_List[i_exploiters].fitness_amb_E = fitness_amb_E;
        Exploiter_Species_List[i_exploiters].fitness_int_E = fitness_int_E;
        Exploiter_Species_List[i_exploiters].fitness_E = (1. - xi_d_E) * fitness_amb_E + xi_d_E * fitness_int_E;
    }

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        diff_amb_V = fabs(Victim_Species_List[i_victims].z_trait - Victim_Species_List[i_victims].theta);

        fitness_amb_V = (1. - diff_amb_V / max_xi_V);
        int k_V = Victim_Species_List[i_victims].k_victim;
        fitness_int_V = 0.;
        for (idx_exploiters = 0; idx_exploiters < k_V; idx_exploiters++)
        {
            int id_exploiter;
            id_exploiter = Victim_Species_List[i_victims].w_E[idx_exploiters];
            diff_amb_EV = fabs(Victim_Species_List[i_victims].theta - Exploiter_Species_List[id_exploiter].theta);
            diff_sel_V = fabs(Victim_Species_List[i_victims].z_trait - Exploiter_Species_List[id_exploiter].z_trait);
            fitness_int_V += 1. / k_V * (diff_sel_V / (max_xi_E + max_xi_V + diff_amb_EV));
        }
        Victim_Species_List[i_victims].fitness_amb_V = fitness_amb_V;
        Victim_Species_List[i_victims].fitness_int_V = fitness_int_V;
        Victim_Species_List[i_victims].fitness_V = (1. - xi_d_V) * fitness_amb_V + xi_d_V * fitness_int_V;
    }
}

/**########################

TEXTFILE FUNCTIONS

#########################*/

void write_header_global_temporal(char name_1[200], char name_2[200])
{
    /*
        File with the full time series for each simulation
    */

    int i_exploiters, i_victims;

    time_t curtime;
    time(&curtime);

    f_out_victims = fopen(name_1, "w");

    fprintf(f_out_victims, "#Date of simulation: %s", ctime(&curtime));
    fprintf(f_out_victims, "#Seed: %d\n#\n", MSEED);
    fprintf(f_out_victims, "#Number of simulations: %d\n#\n", N_Sim);
    fprintf(f_out_victims, "#----------------------------------------------------------\n");
    fprintf(f_out_victims, "#net_id,xi_E,xi_V,sim,time");
    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        fprintf(f_out_victims, ",theta_V_%d,z_V_%d,fitness_V_%d", i_victims, i_victims, i_victims);
    }
    fprintf(f_out_victims, "\n");

    f_out_exploiters = fopen(name_2, "w");

    fprintf(f_out_exploiters, "#Date of simulation: %s", ctime(&curtime));
    fprintf(f_out_exploiters, "#Seed: %d\n#\n", MSEED);
    fprintf(f_out_exploiters, "#Number of simulations: %d\n#\n", N_Sim);
    fprintf(f_out_exploiters, "#----------------------------------------------------------\n");
    fprintf(f_out_exploiters, "#net_id,xi_E,xi_V,sim,time");
    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        fprintf(f_out_exploiters, ",theta_E_%d,z_E_%d,fitness_E_%d", i_exploiters, i_exploiters, i_exploiters);
    }
    fprintf(f_out_exploiters, "\n");
}

void write_data_time(int simul, double t_iter)
{
    int i_exploiters, i_victims;

    fprintf(f_out_victims, "%03d,%.3lf,%.3lf,%d,%.5lf", net_id, xi_d_E, xi_d_V, simul, t_iter);
    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        fprintf(f_out_victims, ",%.5lf,%.5lf,%.5lf", Victim_Species_List[i_victims].theta, Victim_Species_List[i_victims].z_trait, Victim_Species_List[i_victims].fitness_V);
    }
    fprintf(f_out_victims, "\n");

    fprintf(f_out_exploiters, "%03d,%.3lf,%.3lf,%d,%.5lf", net_id, xi_d_E, xi_d_V, simul, t_iter);
    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        fprintf(f_out_exploiters, ",%.5lf,%.5lf,%.5lf", Exploiter_Species_List[i_exploiters].theta, Exploiter_Species_List[i_exploiters].z_trait, Exploiter_Species_List[i_exploiters].fitness_E);
    }
    fprintf(f_out_exploiters, "\n");
}

void write_header_global_final(char name_1[200])
{
    /*
        File with the final (t -> TMAX) state of each simulation
    */

    int i_exploiters, i_victims;

    time_t curtime;
    time(&curtime);

    f_out_final = fopen(name_1, "w");

    fprintf(f_out_final, "#Date of simulation: %s", ctime(&curtime));
    fprintf(f_out_final, "#Seed: %d\n#\n", MSEED);
    fprintf(f_out_final, "#Number of simulations: %d\n#\n", N_Sim);
    fprintf(f_out_final, "#----------------------------------------------------------\n");
    fprintf(f_out_final, "#net_id,xi_E,xi_V,sim");

    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        fprintf(f_out_final, ",theta_V_%d,z_V_%d_tinf,fitness_V_%d,fitness_amb_V_%d,fitness_int_V_%d", i_victims, i_victims, i_victims, i_victims, i_victims);
    }

    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        fprintf(f_out_final, ",theta_E_%d,z_E_%d_tinf,fitness_E_%d,fitness_amb_E_%d,fitness_int_E_%d", i_exploiters, i_exploiters, i_exploiters, i_exploiters, i_exploiters);
    }

    fprintf(f_out_final, "\n");
}

void write_data_final(int simul)
{
    int i_exploiters, i_victims;
    fprintf(f_out_final, "%03d,%.3lf,%.3lf,%d", net_id, xi_d_E, xi_d_V, simul);
    for (i_victims = 0; i_victims < M_Victims; i_victims++)
    {
        fprintf(f_out_final, ",%.5lf,%.5lf,%.5lf,%.5lf,%.5lf", Victim_Species_List[i_victims].theta, Victim_Species_List[i_victims].z_trait,
                Victim_Species_List[i_victims].fitness_V, Victim_Species_List[i_victims].fitness_amb_V, Victim_Species_List[i_victims].fitness_int_V);
    }
    for (i_exploiters = 0; i_exploiters < M_Exploiters; i_exploiters++)
    {
        fprintf(f_out_final, ",%.5lf,%.5lf,%.5lf,%.5lf,%.5lf", Exploiter_Species_List[i_exploiters].theta, Exploiter_Species_List[i_exploiters].z_trait,
                Exploiter_Species_List[i_exploiters].fitness_E, Exploiter_Species_List[i_exploiters].fitness_amb_E, Exploiter_Species_List[i_exploiters].fitness_int_E);
    }
    fprintf(f_out_final, "\n");
}
