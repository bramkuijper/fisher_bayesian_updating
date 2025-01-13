//      Fisherian sexual selection 
//      with Bayesian updating / learning during mate choice
//

// load some libraries that we need
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <cmath>
#include <cassert>
#include <random>



// random number generators
//
// we need something to 'seed' our random
// number generators (i.e., which sequence 
// of random numbers will we pick)
std::random_device rd; 
unsigned seed = rd();

// the actual random number generator
std::mt19937 rng_r(seed);

// allocate a uniform [0,1] random number distribution
std::uniform_real_distribution <double> uniform(0.0,1.0);

// standard normal distribution with mean 0 and stddev = 1.0
std::normal_distribution <double> standard_normal;

// allocate a bernoulli distribution for allele 
// segregration during inheritance
std::bernoulli_distribution segregator(0.5);



// parameters (some of these values are changed upon initialization
// through the command line, see the argv and argc arrays in the 
// init function down below)
const int N = 5000; // population size
const int N_mate_sample = 10; // number of mates sampled
const int clutch_size = 10; // number of offspring produced

double init_t = 0.0; // initial value for ornament
double init_p = 0.0; // initial value for preference
double init_prior_var = 1.0; // initial value for preference
double a = 1.0; // efficacy of sexual selection
double b = 0.5; // cost of preference 
double c = 0.5; // cost of trait
double biast = 0.5; // mutation bias 

double mu_p 	  = 0.05;            // mutation rate preference
double mu_t 	  = 0.05;            // mutation rate ornament
double mu_prior_mean_t 	  = 0.05;            // mutation rate ornament
double mu_prior_sigma_t 	  = 0.05;            // mutation rate ornament
double sdmu         = 0.4;			 // standard deviation mutation stepsize
double sdmu_prior         = 0.05;			 // mutational probability of the prior
                                         //
const double NumGen = 150000; // number of generations
const int skip = 10; // n generations interval before data is printed

// stats for mate assessment
int n_mates_assessed_mean = 0;
int n_mates_assessed_ss = 0;

double diff_mean_prior_posterior = 0;
double diff_var_prior_posterior = 0;

// default name of the output file
std::string file_name = "output.csv";

// preference function:
// 0: open-ended preferences
// 1: absolute
// 2: relative
int pref = 0;

// statistics variables
double meanornsurv = 0;

int popsize = N; // population size between 
bool do_stats = 0;

int generation = 0;
int Nfemales = N / 2, Nmales = N / 2;
int msurvivors = 0;
int fsurvivors = 0;

int father_eggs[N];
int mother_eggs[N];


// the components of an actual individual
struct Individual
{
	double t[2]; // diploid, additive loci for t,p
	double p[2]; // this is now the choice threshold against which the mean of the
                 // prior or posterior ornament distribution is compared.
    double prior_mean_t[2]; // prior for the mean male ornament, expressed by females
    double prior_sigma_t[2]; // prior for the variance in male ornaments, expressed by females
    double t_expr; // and store their expressed values
    double p_expr;

    double posterior_mean_t; //posterior mean
    double posterior_sigma_t; // posterior variance

    double cost; // cost of mate assessment (when female)
    int mate;

};

// generate the population
typedef Individual Population[N];

// declare various arrays of N individuals 
Population Females, Males, FemaleSurvivors, MaleSurvivors;

// make an array to store the index values of the parents
// to generate offspring (rather than to make an array with lots of
// offspring an array with lots of indices is cheaper).
int Parents[N*clutch_size][2]; 


// declaring parameters and arrays done, now let's move
// on to declaring functions

// function which obtains arguments from the command line
// for definitions of the various parameters see top of the file
void initArguments(int argc, char *argv[])
{
	a = std::stod(argv[1]);
	b = std::stod(argv[2]);
	c = std::stod(argv[3]);
	biast = std::stod(argv[4]);
	mu_p = std::stod(argv[5]);
	mu_t = std::stod(argv[6]);
	sdmu = std::stod(argv[7]);
	sdmu_prior = std::stod(argv[8]);
    file_name = argv[9];
} // end initArguments

// mutation function:
// mutate an allele G given a particular mutation rate mu
// a standard deviation of the mutational effect distribution
// of sdmu and a certain mutational bias
void mutate(double &G, double mu, double sdmu, double bias=0.5)
{
    if (uniform(rng_r) < mu)
    {
        double effect = sdmu * uniform(rng_r);
        G+= uniform(rng_r) < bias ? -effect : effect;
    }
}

void mutate_mean(double &G, double mu, double sdmu)
{
    if (uniform(rng_r) < mu)
    {
        G+= sdmu * standard_normal(rng_r);
    }
}

void mutate_var(double &G, double mu, double sdmu)
{
    if (uniform(rng_r) < mu)
    {
        G+= sdmu * standard_normal(rng_r);

        if (G < 0.0)
        {
            G = 0.0;
        }
    }
}

// write the parameters to the DataFile
void WriteParameters(std::ofstream &DataFile)
{
	DataFile << std::endl
		<< std::endl
		<< "type;" << "gonochorist_fisherian" << ";" << std::endl
		<< "popsize_init;" << N << ";" << std::endl
		<< "n_mate_sample;" << N_mate_sample << ";"<< std::endl
		<< "init_t;" << init_t << ";"<< std::endl
		<< "init_p;" << init_p << ";"<< std::endl
		<< "a;" <<  a << ";"<< std::endl
		<< "b;" <<  b << ";"<< std::endl
		<< "c;" <<  c << ";"<< std::endl
		<< "pref;" <<  pref << ";"<< std::endl
		<< "mu_p;" <<  mu_p << ";"<< std::endl
		<< "mu_t;" <<  mu_t << ";"<< std::endl
		<< "mu_prior_mean_t;" <<  mu_prior_mean_t << ";"<< std::endl
		<< "mu_prior_sigma_t;" <<  mu_prior_sigma_t << ";"<< std::endl
		<< "mu_std;" <<  sdmu << ";"<< std::endl
		<< "mu_std_prior;" <<  sdmu_prior << ";"<< std::endl
		<< "biast;" <<  biast << ";"<< std::endl
		<< "seed;" << seed << ";"<< std::endl;
} // end WriteParameters()

// initialize all the phenotypes
void Init()
{
	// initialize the whole populatin
	for (int i = 0; i < Nfemales; ++i)
	{
        // initialize both diploid loci
		for (int j = 0; j < 2; ++j)
		{
			Females[i].t[j] = init_t;
			Females[i].p[j] = init_p;
			Females[i].prior_mean_t[j] = init_t;
			Females[i].prior_sigma_t[j] = init_prior_var;
		}
        
        // and the expressed values
        Females[i].t_expr = init_t;
        Females[i].p_expr = init_p;
			
	}

    // initialize the male part of the population
	for (int i = 0; i < Nmales; ++i)
	{
		for (int j = 0; j < 2; ++j)
		{
			Males[i].t[j] = init_t;
			Males[i].p[j] = init_p;
			Males[i].prior_mean_t[j] = init_t;
			Males[i].prior_sigma_t[j] = init_prior_var;
		}
			
        Males[i].t_expr = init_t;
        Males[i].p_expr = init_p;
	}
} // end Init

// create an offspring 
void Create_Kid(int const mother, int const father, Individual &kid)
{
	assert(mother >= 0);
    assert(mother < fsurvivors);
	assert(father >= 0); 
    assert(father < msurvivors);

    // inherit ornament
	kid.t[0] = FemaleSurvivors[mother].t[segregator(rng_r)];
    mutate(kid.t[0], mu_t, sdmu, biast);

	kid.t[1] = MaleSurvivors[father].t[segregator(rng_r)];
    mutate(kid.t[1], mu_t, sdmu, biast);

    // inherit preference
	kid.p[0] = FemaleSurvivors[mother].p[segregator(rng_r)];
    mutate(kid.p[0], mu_p, sdmu);

	kid.p[1] = MaleSurvivors[father].p[segregator(rng_r)];
    mutate(kid.p[1], mu_p, sdmu);

    // inherit prior mean
    kid.prior_mean_t[0] = FemaleSurvivors[mother].prior_mean_t[segregator(rng_r)];
    mutate_mean(kid.prior_mean_t[0], mu_prior_mean_t, sdmu_prior);

    kid.prior_mean_t[1] = MaleSurvivors[mother].prior_mean_t[segregator(rng_r)];
    mutate_mean(kid.prior_mean_t[1], mu_prior_mean_t, sdmu_prior);
    
    // inherit prior
    kid.prior_sigma_t[0] = FemaleSurvivors[mother].prior_sigma_t[segregator(rng_r)];
    mutate_var(kid.prior_sigma_t[0], mu_prior_sigma_t, sdmu_prior);

    kid.prior_sigma_t[1] = MaleSurvivors[mother].prior_sigma_t[segregator(rng_r)];
    mutate_var(kid.prior_sigma_t[1], mu_prior_sigma_t, sdmu_prior);
} // end Create_Kid

// survival stage
void Survive(std::ofstream &DataFile)
{
    // store individual fitness values
    double w = 0; 

    // keep track of the 
    // mean ornament size of surviving males
    // necessary for absolute/relative preference
    // functions
    meanornsurv = 0;

    // all females simply survive,
    // as survival costs will be incurred prior to reproduction
    // because females now pay a cost per assessment

    fsurvivors = 0;
    msurvivors = 0;
    
    // allow females to survive
	for (int i = 0; i < Nfemales; ++i)
	{
        FemaleSurvivors[fsurvivors++] = Females[i];
	}

    // male survival
	for (int i = 0; i < Nmales; ++i)
	{
		double t_expr = Males[i].t_expr;

		w = exp(-c*t_expr*t_expr);
        
        if (uniform(rng_r) < w)
        {
            // in case of relative preferences get the mean ornament
            meanornsurv += t_expr;
            MaleSurvivors[msurvivors++] = Males[i];
        }
	}

    // extinction?
    if (msurvivors == 0)
    {
        WriteParameters(DataFile);

        exit(1);
    }

    // take the average of the surviving male trait value
    meanornsurv /= msurvivors;

    assert(msurvivors > 0);
    assert(msurvivors < popsize);
} // end survival stage

// mate choice

// have a female choose among males
void Choose(Individual &female) 
{
    female.mate = -1;

	// check if there are enough other males to choose from
	// otherwise restrict the sample to the amount of individuals present
    // note that N_mate_sample is the max here,
    // once threshold of a female hits
	int current_mate_sample = N_mate_sample > msurvivors ? msurvivors : N_mate_sample;

    std::uniform_int_distribution <int> msurvivor_sampler(0, msurvivors - 1);

    // develop the prior mean
    double prior_mean = female.prior_mean_t[0] + female.prior_mean_t[1];
    double prior_sigma = female.prior_sigma_t[0] + female.prior_sigma_t[1];
    double prior_sigma_sq_inv = 1.0 / (prior_sigma * prior_sigma);

    double posterior_mean = 0.0;
    double posterior_sigma = 0.0;
    double posterior_variance = 0.0;
    double posterior_var_inv = 0.0;

    // auxiliary variable that contains a sample from the prior distribution 
    // (before mate assessment)
    // and otherwise a sample from the posterior distribution
    double estimate;

    double preference = female.p_expr;

    // aux variables for updating things
    double mean_ornaments = 0.0;
    double ss_ornaments = 0.0;

    unsigned n_mates_assessed = 0;

    // mate choice among the sample of males
	for (int j = 0; j < current_mate_sample; ++j)
	{
		// get a random surviving male
		int random_mate = msurvivor_sampler(rng_r);

        assert(random_mate >= 0);
        assert(random_mate < msurvivors);

        // obtain a male's ornament
		double trait = MaleSurvivors[random_mate].t_expr;

        if (j == 0)
        {
            estimate = prior_mean + prior_sigma * standard_normal(rng_r);
        } 
        else
        {
            estimate = posterior_mean + 
                sqrt(posterior_variance) * standard_normal(rng_r);
        }

        // check whether current estimate exceeds 
        // the female preference threshold
        if (estimate >= preference)
        {
            // ok we have to think about this, because if the 
            // estimate >= preference clause below then this will 
            // also be automatically false, but it is what it is.
            female.mate = random_mate;
            break;
        } 
        else // ok estimate does not exceed preference, let's assess this male
        {
            ++n_mates_assessed;

            posterior_var_inv = 1.0 / (posterior_sigma * posterior_sigma);

            // update estimate based on mean and variance
            // this is eq. (11) of Luttbeg (1996)
            posterior_mean = (
                    trait * prior_sigma_sq_inv + 
                        posterior_mean * 1.0 / posterior_var_inv
                    )/ (posterior_var_inv + prior_sigma_sq_inv);
            
            posterior_variance = 1.0 / prior_sigma_sq_inv;

            mean_ornaments += trait;
            ss_ornaments += trait * trait;

            if (n_mates_assessed > 1)
            {
                // this is rho in eq. (12)
                double ornament_variance_sampled_males = ss_ornaments / n_mates_assessed - 
                    mean_ornaments / n_mates_assessed;

                if (ornament_variance_sampled_males <= 0)
                {
                    ornament_variance_sampled_males = 1.0 / prior_sigma_sq_inv;
                }

                // yielding rho_new in eq. (12)
                posterior_variance = 1.0 / (prior_sigma_sq_inv + 
                        1.0/ornament_variance_sampled_males);

            }

            // now sample estimate after assessing male
            estimate = posterior_mean + 
                sqrt(posterior_variance) * standard_normal(rng_r);

            if (estimate >= preference)
            {
                female.mate = random_mate;
                break;
            }
        }
	} // end for (int j = 0 

    // female exhausted all options, mate randomly
    if (female.mate < 0 && n_mates_assessed == current_mate_sample)
    {
        female.mate = msurvivor_sampler(rng_r);
    }

    assert(n_mates_assessed <= N_mate_sample);

    n_mates_assessed_mean += n_mates_assessed;
    n_mates_assessed_ss += n_mates_assessed * n_mates_assessed;

    // calculate cost of female choice
    female.cost = n_mates_assessed * b;
} // end ChooseMates


// produce the next generation
void NextGen()
{
    // do some stats of the adult reproductive
    // population
    if (do_stats)
    {
        for (int i = 0; i < msurvivors; ++i)
        {
            father_eggs[i] = 0;
        }

        for (int i = 0; i < fsurvivors; ++i)
        {
            mother_eggs[i] = 0;
        }
    } // end if(do_stats)
    
    std::vector<double> costs_of_choice;

    n_mates_assessed_mean = 0;
    n_mates_assessed_ss = 0;

    // let the females choose a mate
	for (int i = 0; i < Nfemales; ++i)
	{
		Choose(FemaleSurvivors[i]);

		assert(FemaleSurvivors[i].mate >= 0 && 
               FemaleSurvivors[i].mate < msurvivors);

        // for a fitness weighting, we need to select
        // females that have the highest fitness value
        // (i.e., the lowest cost weighting). Obvious
        // choice would be to make costs negative, but weightings
        // of a discrete_distribution always have to be positive
        // Hence what we do is:
        // max possible cost is N_mate_sample * b. 
        // from this we subtract the actual cost, so that individuals
        // who minimize cost have the largest payoff
        costs_of_choice.push_back(N_mate_sample * b - FemaleSurvivors[i].cost);
    }

    n_mates_assessed_mean /= Nfemales;
    n_mates_assessed_ss /= Nfemales;

    std::discrete_distribution <int> fitness_distribution(
        costs_of_choice.begin(),
            costs_of_choice.end());

    int sons = 0;
    int daughters = 0;

    int mother_idx;

    // replace the next generation
    for (int i = 0; i < popsize; ++i)
    {
        // create an offspring
        Individual Kid;

        mother_idx = fitness_distribution(rng_r);

        assert(FemaleSurvivors[mother_idx].mate >= 0);
        assert(FemaleSurvivors[mother_idx].mate < msurvivors);

        // randomly sample an offspring to replace the population
        Create_Kid(mother_idx, 
                    FemaleSurvivors[mother_idx].mate, 
                    Kid);
            
        double t = 0.5 * ( Kid.t[0] + Kid.t[1]);
        double p = 0.5 * ( Kid.p[0] + Kid.p[1]);

        Kid.t_expr = t; 
        Kid.p_expr = p; 

        // it's a boy
        if (uniform(rng_r) < 0.5)
        {
            Males[sons++] = Kid;
        }
        else
        {
            Females[daughters++] = Kid;
        }
    }

    Nmales = sons;
    Nfemales = daughters;
}



// write the data
void WriteData(std::ofstream &DataFile)
{
    // in case population is extinct,
    // quit
	if (Nmales == 0 || Nfemales == 0)
	{
		WriteParameters(DataFile);
		exit(1);
	}

    double meanp = 0.0;
    double meant = 0.0; 
    double ssp = 0.0;
    double sst = 0.0;
    double spt = 0.0;

    double mean_mean_prior = 0.0;
    double ss_mean_prior = 0.0;
    double mean_sigma_prior = 0.0;
    double ss_sigma_prior = 0.0;

    double p,t,meanmrs,meanfrs,varfrs,varmrs,prior_mean,prior_sigma;
	double ssmrs = 0, ssfrs = 0, summrs=0, sumfrs=0;

    // calculate means and variances for the males
	for (int i = 0; i < Nmales; ++i)
	{
		p = 0.5 * ( Males[i].p[0] + Males[i].p[1]);
		t = 0.5 * ( Males[i].t[0] + Males[i].t[1]);

        meanp += p;
        meant += t;

        ssp += p * p;
        sst += t * t;
        spt += t * p;

        prior_mean = 0.5 * (Males[i].prior_mean_t[0] + Males[i].prior_mean_t[1]);
        prior_sigma = 0.5 * (Males[i].prior_sigma_t[0] + Males[i].prior_sigma_t[1]);

        mean_mean_prior += prior_mean;
        mean_sigma_prior += prior_sigma;
        ss_mean_prior += prior_mean * prior_mean;
        ss_sigma_prior += prior_sigma * prior_sigma;

        if (i < msurvivors)
        {
            summrs += father_eggs[i];
            ssmrs += father_eggs[i] * father_eggs[i];
        }
	}

    // calculate means and variances for the females
	for (int i = 0; i < Nfemales; ++i)
	{
		p = 0.5 * ( Females[i].p[0] + Females[i].p[1]);
		t = 0.5 * ( Females[i].t[0] + Females[i].t[1]);

        meanp += p;
        meant += t;

        ssp += p * p;
        sst += t * t;
        spt += t * p;
        
        prior_mean = 0.5 * (Females[i].prior_mean_t[0] + Females[i].prior_mean_t[1]);
        prior_sigma = 0.5 * (Females[i].prior_sigma_t[0] + Females[i].prior_sigma_t[1]);

        mean_mean_prior += prior_mean;
        mean_sigma_prior += prior_sigma;
        ss_mean_prior += prior_mean * prior_mean;
        ss_sigma_prior += prior_sigma * prior_sigma;

        
        if (i < fsurvivors)
        {
            sumfrs += mother_eggs[i];
            ssfrs += mother_eggs[i] * mother_eggs[i];
        }
	} 

    double varp = 0; 
    double vart = 0; 
    double covpt = 0;

    double sum_sexes = Nmales + Nfemales;

    meant /= sum_sexes;
    meanp /= sum_sexes;
    mean_mean_prior /= sum_sexes;
    mean_sigma_prior /= sum_sexes;
    ss_mean_prior /= sum_sexes;
    ss_sigma_prior /= sum_sexes;

    vart = sst / sum_sexes - meant * meant;
    varp = ssp / sum_sexes - meanp * meanp;
    covpt = spt / sum_sexes - meanp * meant;

    double var_mean_prior = ss_mean_prior - mean_mean_prior * mean_mean_prior;
    double var_sigma_prior = ss_sigma_prior - mean_sigma_prior * mean_sigma_prior;

    meanfrs = sumfrs / Nfemales;
    varfrs = ssfrs / Nfemales - meanfrs * meanfrs;

    meanmrs = summrs / Nmales;
    varmrs = ssmrs / Nmales - meanmrs * meanmrs;

    // calculate variance in the numbers of males assessed
    double n_mates_assessed_var = 
        n_mates_assessed_ss - 
            n_mates_assessed_mean * n_mates_assessed_mean;

    // output of all the statistics
	DataFile << generation 
		<< ";" << meanp
		<< ";" << meant
		
        << ";" << meanfrs 
        << ";" << meanmrs 
		
        << ";" << varp
		<< ";" << vart
		<< ";" << covpt 
		
        << ";" << varfrs 
        << ";" << varmrs 

        << ";" << n_mates_assessed_mean
        << ";" << n_mates_assessed_var

        << ";" << mean_mean_prior
        << ";" << mean_sigma_prior
        
        << ";" << var_mean_prior
        << ";" << var_sigma_prior

        << ";" << msurvivors
        << ";" << fsurvivors
		<< ";" << sum_sexes
        << ";" << std::endl;

} // end WriteData

// headers of the datafile
void WriteDataHeaders(std::ofstream &DataFile)
{
	DataFile << "generation" 
		<< ";mean_p" 
		<< ";mean_t" 
        << ";meanfrs"
        << ";meanmrs"
        << ";varp"
        << ";vart"
        << ";covpt"
        << ";varfrs"
        << ";varmrs"
        << ";nmates_assessed"
        << ";nmates_assessed_var"
        << ";mean_mean_prior"
        << ";mean_sigma_prior"
        << ";var_mean_prior"
        << ";var_sigma_prior"
        << ";surviving_males"
        << ";surviving_females"
		<< ";N;"
		<< std::endl;
} // end WriteDataHeaders

// the core part of the code
int main(int argc, char ** argv)
{
    // initialize parameters based on command line arguments
	initArguments(argc, argv);

    // initialize output file
    std::ofstream output_file(file_name.c_str());
	

    // write the headers to the data file
	WriteDataHeaders(output_file);

    // initialize the population
	Init();

    // loop through each generation
	for (generation = 0; generation <= NumGen; ++generation)
	{
        // decide whether we want to write out stats or not
		do_stats = generation % skip == 0;

        // individuals survive (or not)
		Survive(output_file);
        
        // individuals reproduce (or not)
        // female choice, so all females reproduce (by definition)
        // but not all males.
        NextGen();
        
        if (do_stats)
		{
			WriteData(output_file);
		}
	}

    WriteParameters(output_file);
} // end main
