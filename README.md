# GENETIC_ALGORITHM
# DESCRIPTION:
The genetic algorithm (GA) here is understood as an extension to weighted gene co-expression
network analysis (WGCNA) by Peter Langfelder and Steve Horvath (BMC Bioinformatics, 2008, 9:559).
whereby the genes expression data corresponding to the genes identified by WGCNA as gene modules
is provided to GA as a parameter as well as the trait of interest.
The GA optimizes the gene module to trait relationship by gradually increasing the correlation
between the trait and a subset of genes of the gene module
1.) It does so by creating an inital population of chromosomes, which are binary vectors of length gene
    module. The values within a chromosome are set randomly to '1' or '0' goverened by the parameter number.of.genes.
    For each chromosome an expression matrix is created, where the values within the vector are set to '1'.
2.) As the next step the fitness value for each expression matrix is determined. For that the 1st principal
    component of each epxression matrix is computed and then correlated (Pearson correlation) to the trait
    of interest. 
3.) Chromosomes with a higher fitness value have greater chances to contribute to the next generation. The next
    generation is generated via a recombination event, where chromosome_father x chromosome_mother. The recombination
    event is performed at a random location/locus L of the chromosome, so that chromosome_father[1:L] + 
   chromosome_mother[L+1:length(chromosome_mother)] result in chromosome_offspring.
4.) As the last step mutation events can occur randomly at each position of the chromosome_offspring flipping 
    the value from '0' to '1' or from '1' to '0', respectively
    Then the GA is reiterated. The GA is repeated for n generations.
   For a full description of the GA, we refer the user to the accompanying publication.


# PARAMETERS:
   data: a dataset containing gene expression data
   population.size: the number of individuals (here chromosomes) the population is composed of
   number.of.genes: the number of genes that should be on the chromosome in the initial population
   mutation.rate: the chance for a mutation to occur on either of the cells/genes in the chromosome vector
   the mutation changes the value of 0 to 1 and 1 to 0
   crossover.events: during the reproduction phase the individual chromosomes of the father and mother
   will crossover up to the times specified by crossover.events
   generations: how many generations should be simulated to improve fitness
   traits: contains the other traits the genes will be correlated to
   stats: stats is an optional parameter that will create a file containing statisics about the GA
   run in real time. If a file is defined, the statistics will be written to it
   iteration: iteration is associated with parameter stats. It is adivsable to run the GA for several
   iterations. If the GA is run with multiple iterations only one stats file will be created

