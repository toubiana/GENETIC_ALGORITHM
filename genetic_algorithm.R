#/
## GENETIC ALGORITHM
## AUTHOR: David Toubiana
## email: dtoubiana@ucdavis.edu
## SHORT DESCRIPTION: Using the genetic algorithm to test a set of genes with the strongest correlation
## to a trait of interest
#/

#/
## ----------------------------------------------DISCLAIMER----------------------------------------------------##
## This SOFTWARE PRODUCT is provided by THE PROVIDER "as is" and "with all faults."
## THE PROVIDER makes no representations or warranties of any kind concerning the safety,
## suitability, lack of viruses, inaccuracies, typographical errors, or other harmful
## components of this SOFTWARE PRODUCT. There are inherent dangers in the use of any software,
## and you are solely responsible for determining whether this SOFTWARE PRODUCT is compatible
## with your equipment and other software installed on your equipment.
## You are also solely responsible for the protection of your equipment and backup of your data,
## and THE PROVIDER will not be liable for any damages you may suffer in connection with using,
## modifying, or distributing this SOFTWARE PRODUCT.
## ---------------------------------------------DISCLAIMER END--------------------------------------------------##
#/

#/
##----------------------------------------------DESCRIPTION-----------------------------------------------------##
## the genetic algorithm (GA) here is understood as an extension to weighted gene co-expression
## network analysis (WGCNA) by Peter Langfelder and Steve Horvath (BMC Bioinformatics, 2008, 9:559).
## whereby the genes expression data corresponding to the genes identified by WGCNA as gene modules
## is provided to GA as a parameter as well as the trait of interest.
## The GA optimizes the gene module to trait relationship by gradually increasing the correlation
## between the trait and a subset of genes of the gene module
## 1.) It does so by creating an inital population of chromosomes, which are binary vectors of length gene
## module. The values within a chromosome are set randomly to '1' or '0' goverened by the parameter number.of.genes.
## For each chromosome an expression matrix is created, where the values within the vector are set to '1'.
## 2.) As the next step the fitness value for each expression matrix is determined. For that the 1st principal
## component of each epxression matrix is computed and then correlated (Pearson correlation) to the trait
## of interest. 
## 3.) Chromosomes with a higher fitness value have greater chances to contribute to the next generation. The next
## generation is generated via a recombination event, where chromosome_father x chromosome_mother. The recombination
## event is performed at a random location/locus L of the chromosome, so that chromosome_father[1:L] + 
## chromosome_mother[L+1:length(chromosome_mother)] result in chromosome_offspring.
## 4.) As the last step mutation events can occur randomly at each position of the chromosome_offspring flipping 
## the value from '0' to '1' or from '1' to '0', respectively
## Then the GA is reiterated. The GA is repeated for n generations.
## For a full description of the GA, we refer the user to the accompanying publication.
##----------------------------------------------DESCRIPTION END-------------------------------------------------##
#/ 

#/
## ------------------------------------------------PARAMETERS---------------------------------------------------##
## data: a dataset containing gene expression data
## population.size: the number of individuals (here chromosomes) the population is composed of
## number.of.genes: the number of genes that should be on the chromosome in the initial population
## mutation.rate: the chance for a mutation to occur on either of the cells/genes in the chromosome vector
## the mutation changes the value of 0 to 1 and 1 to 0
## crossover.events: during the reproduction phase the individual chromosomes of the father and mother
## will crossover up to the times specified by crossover.events
## generations: how many generations should be simulated to improve fitness
## traits: contains the other traits the genes will be correlated to
## stats: stats is an optional parameter that will create a file containing statisics about the GA
## run in real time. If a file is defined, the statistics will be written to it
## iteration: iteration is associated with parameter stats. It is adivsable to run the GA for several
## iterations. If the GA is run with multiple iterations only one stats file will be created
## ----------------------------------------------PARAMETERS END--------------------------------------------------##
#/

# License
# 
# The following license governs the use of the source code for GENETIC_ALGORIHM in academic and educational environments. 
# Commercial use requires a commercial license from the author directly: 
#     David Toubiana
# email: dtoubiana@ucdavis.edu or david.toubiana@gmail.com
# 
# ACADEMIC PUBLIC LICENSE
# 
# Original license written by Andras Varga and modified by David Toubiana (license text is in public domain)
# 
# Preamble
# 
# This license contains the terms and conditions of using the source code for GENETIC_ALGORIHM in noncommercial settings: at academic
# institutions for teaching and research use, and at non-profit research organizations. You will find that this license provides 
# noncommercial users of the source code for the GENETIC_ALGORIHM with rights that are similar to the well-known GNU General Public License, 
# yet it retains the possibility for the source code for the GENETIC_ALGORIHM authors to financially support the development by selling 
#commercial licenses. In fact, if you intend to use the source code for the GENETIC_ALGORIHM in a “for-profit” environment, where research
# is conducted to develop or enhance a product, is used in a commercial service offering, or when a commercial company uses the source 
# code for the GENETIC_ALGORIHM to participate in a research project (for example government-funded or EU-funded research projects), 
# then you need to obtain a commercial license for the source code for the GENETIC_ALGORIHM. In that case, please contact the Author to
# inquire about commercial licenses.
# 
# What are the rights given to noncommercial users? Similarly to GPL, you have the right to use the software, to distribute copies, to 
# receive source code, to change the software and distribute your modifications or the modified software. Also similarly to the GPL, 
# if you distribute verbatim or modified copies of this software, they must be distributed under this license.
# 
# By modeling the GPL, this license guarantees that you’re safe when using the source code for the GENETIC_ALGORIHM in your work, for 
# teaching or research. This license guarantees that the source code for the GENETIC_ALGORIHM will remain available free of charge for 
# nonprofit use. You can modify the source code for the GENETIC_ALGORIHM to your purposes, and you can also share your modifications. 
# Even in the unlikely case of the authors abandoning the source code for the GENETIC_ALGORIHM entirely, this license permits anyone to 
# continue developing it from the last release, and to create further releases under this license.
# 
# We believe that the combination of noncommercial open-source and commercial licensing will be beneficial for the whole user community, 
# because income from commercial licenses will enable faster development and a higher level of software quality, while further enjoying 
# the informal, open communication and collaboration channels of open source development.
# 
# The precise terms and conditions for using, copying, distribution and modification follow.
# 
# TERMS AND CONDITIONS FOR USE, COPYING, DISTRIBUTION AND MODIFICATION
# Definitions
# 
# “Program” means a copy of the source code for the GENETIC_ALGORIHM, which is said to be distributed under this Academic Public License.
# 
# “Work based on the Program” means either the Program or any derivative work under copyright law: that is to say, a work containing the 
# Program or a portion of it, either verbatim or with modifications and/or translated into another language. (Hereinafter, translation is 
# included without limitation in the term “modification”.)
# 
# “Using the Program” means any act of creating executables that contain or directly use libraries that are part of the Program, running 
# any of the tools that are part of the Program, or creating works based on the Program.
# 
# Each licensee is addressed as “you”.
# 
# §1. Permission is hereby granted to use the Program free of charge for any noncommercial purpose, including teaching and research at 
# universities, colleges and other educational institutions, research at non-profit research institutions, and personal non-profit 
# purposes. For using the Program for commercial purposes, including but not restricted to consulting activities, design of commercial
# hardware or software networking products, and a commercial entity participating in research projects, you have to contact the Author 
# for an appropriate license. Permission is also granted to use the Program for a reasonably limited period of time for the purpose of
# evaluating its usefulness for a particular purpose.
# 
# §2. You may copy and distribute verbatim copies of the Program’s source code as you receive it, in any medium, provided that you 
# conspicuously and appropriately publish on each copy an appropriate copyright notice and disclaimer of warranty; keep intact all the
# notices that refer to this License and to the absence of any warranty; and give any other recipients of the Program a copy of this 
# License along with the Program.
# 
# §3. You may modify your copy or copies of the Program or any portion of it, thus forming a work based on the Program, and copy and 
# distribute such modifications or work under the terms of Section 2 above, provided that you also meet all of these conditions:
#     
# a) You must cause the modified files to carry prominent notices stating that you changed the files and the date of any change.
# 
# b) You must cause any work that you distribute or publish, that in whole or in part contains or is derived from the Program or any part 
#    thereof, to be licensed as a whole at no charge to all third parties under the terms of this License.
# 
# These requirements apply to the modified work as a whole. If identifiable sections of that work are not derived from the Program, and
# can be reasonably considered independent and separate works in themselves, then this License, and its terms, do not apply to those 
# sections when you distribute them as separate works. But when you distribute the same sections as part of a whole which is a work based 
# on the Program, the distribution of the whole must be on the terms of this License, whose regulations for other licensees extend to the 
# entire whole, and thus to each and every part regardless of who wrote it. (If the same, independent sections are distributed as part of 
# a package that is otherwise reliant on, or is based on the Program, then the distribution of the whole package, including but not 
# restricted to the independent section, must be on the unmodified terms of this License, regadless of who the author of the included 
# sections was.)
# 
# Thus, it is not the intent of this section to claim rights or contest your rights to work written entirely by you; rather, the intent
# is to exercise the right to control the distribution of derivative or collective works based or reliant on the Program.
# 
# In addition, mere aggregation of another work not based on the Program with the Program (or with a work based on the Program) on a 
# volume of storage or distribution medium does not bring the other work under the scope of this License.
# 
# §4. You may copy and distribute the Program (or a work based on it, under Section 3) in object code or executable form under the terms
#     of Sections 2 and 3 above provided that you also do one of the following:
#     
# a) Accompany it with the complete corresponding machine-readable source code, which must be distributed under the terms of Sections 2
#    and 3 above on a medium customarily used for software interchange; or,
# 
# b) Accompany it with a written offer, valid for at least three years, to give any third party, for a charge no more than your cost of
#    physically performing source distribution, a complete machine-readable copy of the corresponding source code, to be distributed 
#    under the terms of Sections 2 and 3 above on a medium customarily used for software interchange; or,
# 
# c) Accompany it with the information you received as to the offer to distribute corresponding source code. (This alternative is allowed 
#    only for noncommercial distribution and only if you received the program in object code or executable form with such an offer, in
#    accord with Subsection b) above.)
# 
# The source code for a work means the preferred form of the work for making modifications to it. For an executable work, complete source
# code means all the source code for all modules it contains, plus any associated interface definition files, plus the scripts used to
# control compilation and installation of the executable. However, as a special exception, the source code distributed need not include 
# anything that is normally distributed (in either source or binary form) with the major components (compiler, kernel, and so on) of the
# operating system on which the executable runs, unless that component itself accompanies the executable.
# 
# If distribution of executable or object code is made by offering access to copy from a designated place, then offering equivalent access
# to copy the source code from the same place counts as distribution of the source code, even though third parties are not compelled to 
# copy the source along with the object code.
# 
# §5. You may not copy, modify, sublicense, or distribute the Program except as expressly provided under this License. Any attempt 
# otherwise to copy, modify, sublicense or distribute the Program is void, and will automatically terminate your rights under this License. 
# However, parties who have received copies, or rights, from you under this License will not have their licenses terminated so long as 
# such parties remain in full compliance.
# 
# §6. You are not required to accept this License, since you have not signed it. Nothing else grants you permission to modify or 
# distribute the Program or its derivative works; law prohibits these actions if you do not accept this License. Therefore, by modifying
# or distributing the Program (or any work based on the Program), you indicate your acceptance of this License and all its terms and 
# conditions for copying, distributing or modifying the Program or works based on it, to do so.
# 
# §7. Each time you redistribute the Program (or any work based on the Program), the recipient automatically receives a license from the
# original licensor to copy, distribute or modify the Program subject to these terms and conditions. You may not impose any further 
# restrictions on the recipients’ exercise of the rights granted herein. You are not responsible for enforcing compliance by third parties
# to this License.
# 
# §8. If, as a consequence of a court judgment or allegation of patent infringement or for any other reason (not limited to patent 
# issues), conditions are imposed on you (whether by court order, agreement or otherwise) that contradict the conditions of this License,
# they do not excuse you from the conditions of this License. If you cannot distribute so as to satisfy simultaneously your obligations 
# under this License and any other pertinent obligations, then as a consequence you may not distribute the Program at all. For example, 
# if a patent license would not permit royalty-free redistribution of the Program by all those who receive copies directly or indirectly 
# through you, then the only way you could satisfy both it and this License would be to refrain entirely from distribution of the Program.
# 
# If any portion of this section is held invalid or unenforceable under any particular circumstance, the balance of the section is 
# intended to apply and the section as a whole is intended to apply in other circumstances.
# 
# §9. If the distribution and/or use of the Program are restricted in certain countries either by patents or by copyrighted interfaces,
# the original copyright holder who places the Program under this License may add an explicit geographical distribution limitation 
# excluding those countries, so that distribution is permitted only in or among countries not thus excluded. In such case, this License
# incorporates the limitation as if written in the body of this License.
# 
# NO WARRANTY
# 
# §10. BECAUSE THE PROGRAM IS LICENSED FREE OF CHARGE, THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY APPLICABLE LAW. 
# EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM “AS IS” WITHOUT WARRANTY OF ANY 
# KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
# PURPOSE. THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM IS WITH YOU. SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME 
# THE COST OF ALL NECESSARY SERVICING, REPAIR OR CORRECTION.
# 
# §11. IN NO EVENT UNLESS REQUIRED BY APPLICABLE LAW OR AGREED ON IN WRITING WILL ANY COPYRIGHT HOLDER, OR ANY OTHER PARTY WHO MAY MODIFY 
# AND/OR REDISTRIBUTE THE PROGRAM AS PERMITTED ABOVE, BE LIABLE TO YOU FOR DAMAGES, INCLUDING ANY GENERAL, SPECIAL, INCIDENTAL OR 
# CONSEQUENTIAL DAMAGES ARISING OUT OF THE USE OR INABILITY TO USE THE PROGRAM INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA BEING
# RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH ANY OTHER PROGRAMS), EVEN 
# IF SUCH HOLDER OR OTHER PARTY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
# 
# END OF TERMS AND CONDITIONS
# 
# In case the above text differs from the license file in the source distribution, the latter is the valid one.

## -------------------------------------------------LICENSE END-------------------------------------------------------------------------##

## --------------------------------------------------FUNCTIONS--------------------------------------------------------------------------##
## this is the container function, which is called externally
genetic_algorithm <- function(data,population.size = 1000, chromosome.size = ncol(data), number.of.genes = 10,
         mutation.rate = 0.0001,crossover.events = 1,generations = 1000,traits,
         stats = '',iteration =  NULL){
  
  
  # initial population
  parent.population <- init_population(population.size, chromosome.size, number.of.genes,crossover.events,
                                       stats,iteration)
  child.population <- vector(mode="list",length=length(parent.population))
  
  ## the GA can be modified to intake more traits and change the fitness computations
  ## modifications: February 18th 2019
  ## if the variable traits contains multiple vectors, then calculate its first principal component
  ## if not keep it the way it is
  t <- NULL
  if(ncol(traits)>1){
    t <- prcomp(t(traits))$rotation[,1]
  } else{
    t <- traits
  }
  fitness <- sapply(parent.population,FUN = estimate_fitness,data=data,trait.1=t) #,
  #trait.2=traits[,2],trait.3=traits[,3])
  avg.fitness.pergeneration <- mean(fitness)
  avg.fitness <- mean(avg.fitness.pergeneration)
  if(exists('global.stats')){
    if(is.null(iteration)){
      global.stats <<- paste('Parent_generation,',mean(fitness),',',
                             median(fitness),',',max(fitness),',',min(fitness),',',
                             sd(fitness),sep="")  
    }else{
      global.stats <<- paste(paste('Iteration_',iteration,sep=""),
                             ',Parent_generation,',mean(fitness),',',
                             median(fitness),',',max(fitness),',',min(fitness),',',
                             sd(fitness),sep="") 
    }
    
  }
  
  for(i in 1:generations){
    print(paste("Generation: ",i,sep=""))
    ## uncomment to following section to create plot during run time
    # if(i==1){
    #     plot(y=avg.fitness.pergeneration[i],x=i,xlim=c(0,generations),ylim=c(0,1),xlab="Generations",
    #          ylab="Average Fitness")
    # }
    # else{
    #     points(y=avg.fitness.pergeneration[i],x=i)
    # }
    ## end section to be commented/uncommented
    
    ## only considering positive correlations = elite parents
    ## negative correlations are considered lethal factor
    # determine order
    elite <- sort(fitness,decreasing = T)
    elite <- elite[elite>0]
    elite.ordered <- order(fitness,decreasing = T)
    elite.ordered <- elite.ordered[1:length(elite)]
    
    # now reproduce with fitness
    if(exists('mutation.counter')){
      ## set it to zero
      mutation.counter <<- 0
    }
    for(j in 1:length(parent.population)){
      
      ## only elite parents can contribute to the next generation
      father <- sample(elite.ordered,size=1,prob=elite)
      mother <- sample(elite.ordered,size=1,prob=elite)
      # parents cannot be the same
      repeat{
        if(mother != father)
          break
        else
          mother <- sample(elite.ordered,size=1,prob=elite)
      } # end repeat
      # reproduce
      child.population[[j]] <- reproduce(parent.population[[father]],
                                         parent.population[[mother]],crossover.events)
      # mutate
      child.population[[j]] <- mutate(child.population[[j]],mutation.rate)
    } # end nested for
    ## prepare output for stats
    if(exists('global.stats')){
      if(i>1){
        global.stats <<- paste(global.stats,length(elite),mutation.counter,sep=",")
        ## making the actual stats output
        cat(global.stats,file=stats,sep='\n',append=T)
      }
    }
    parent.population <- lapply(child.population,unlist)
    
    # determine fitness of each individual in population
    fitness <- sapply(parent.population,FUN = estimate_fitness,data=data,trait.1=t)#,
    # trait.2=traits[,2],trait.3=traits[,3])
    
    avg.fitness.pergeneration <- c(avg.fitness.pergeneration,mean(fitness))
    avg.fitness <- mean(avg.fitness.pergeneration)
    ## write to stats output if parameter exists
    if(exists('global.stats')){
      if(i==1){
        global.stats <<- paste(global.stats,length(elite),mutation.counter,sep=",")
        ## making the actual stats output
        cat(global.stats,file=stats,sep='\n',append=T)
      } 
      if(is.null(iteration)){
        global.stats <<- paste(paste('Generation_',i,sep=""),',',mean(fitness),',',
                               median(fitness),',',max(fitness),',',min(fitness),',',
                               sd(fitness),sep="")  
      }else{
        global.stats <<- paste(paste('Iteration_',iteration,',',sep=""),
                               paste('Generation_',i,sep=""),',',mean(fitness),',',
                               median(fitness),',',max(fitness),',',min(fitness),',',
                               sd(fitness),sep="") 
      }
    }
    gc()
  } # end for
  
  # make a list with all the relavent information, which contains 
  # all the individuals/chromosomes of the last population
  last.generation <- vector(mode='list',length=4)
  last.generation$lastpopulation <- parent.population
  last.generation$fitness <- fitness
  last.generation$average_fitness_per_generation <- avg.fitness.pergeneration
  last.generation$average_fitness <- avg.fitness
  
  return(last.generation)
} ## end function

init_population <- function(population.size, chromosome.size, number.of.genes,crossover.events,stats,iteration){
  
  ## Firts, make sure there are not more genes than the size of the chromosomes
  if(number.of.genes > chromosome.size){
    stop("There cannot be more genes than the size of chromosomes! \n Ending program",call.=F)
  } # end if
  else if(population.size < 10){
    stop("Population size should be at least 10! \n Ending program",call.=F)
  }
  else if(number.of.genes < 1){
    stop("Number of genes should be at least 1! \n Ending program",call.=F)
  }
  else if(crossover.events < 1){
    stop("Crossover rate should be at least 1! \n Ending program",call.=F)
  }
  else if(crossover.events > chromosome.size){
    stop("Crossover rate should not exceed chromosome size! \n Ending program",call.=F)
  }
  else{
    ## create the initial population
    if(stats!=''){
      warning('A file name for parameter stats was entered! A (large) stats file will be generated.')
      ## create a global variable that will be available in all subroutines of the GA
      global.stats <<- NULL
      mutation.counter <<- NULL
      ## print header line to file defined in parameter 'stats'
      ## stats file is created as a comma separated value file
      if(is.null(iteration)){
        cat(paste('Generation,Average correlation,Median correlation,Max correlation,',
                  'Min correlation,SD correlation,Number of elite individuals,',
                  'Number of mutations',sep=""),file=stats,sep='\n',append=F)
      } else{
        if(iteration==1){
          cat(paste('Iteration,Generation,Average correlation,Median correlation,Max correlation,',
                    'Min correlation,SD correlation,Number of elite individuals,',
                    'Number of mutations',sep=""),file=stats,sep='\n',append=F)
        }
      }
      
    }
    return(lapply(1:population.size,function(x) 
      sample(c(rep(1,number.of.genes),rep(0,chromosome.size-number.of.genes)))))
  }
  
} ## end function

reproduce <- function(father,mother,crossover.events){
  ## parameter crossover.events determines the maximum number of cross over events
  ## the actual crossover may lie anywhere between 1 and number in the parameter
  cr <- sample(c(1:crossover.events),1)
  ## determine where the crossover/s should occur
  ## mother and father have the same length of chromosome
  loci <- sort(sample(c(1:length(father)-1),cr))
  # construct child chromosome
  child <- father
  for(i in (seq(1,length(loci),2))){
    if(i<length(loci)){
      child[loci[i]:loci[i+1]] <- mother[loci[i]:loci[i+1]] 
    }
    else{
      child[loci[i]:length(child)] <- mother[loci[i]:length(mother)]
    }
  } # end for
  return(child)
} ## end function

mutate <- function(chromosome,mutation.rate){
  ## create a binary vector of size chromosome determining if and where
  ## mutations should occur
  mutation.loci <- sample(0:1, length(chromosome), 
                          replace=T,prob=c(1-mutation.rate,mutation.rate))
  number.of.mutations <- sum(mutation.loci)
  
  
  if(number.of.mutations > 0){
    chromosome[mutation.loci] <- abs(chromosome[mutation.loci]-1)
    ## if mutation occurred increase mutation.counter
    if(exists('mutation.counter')){
      mutation.counter <<- mutation.counter+number.of.mutations
    }
  } # end if
  return(list(chromosome))
} ## end function

estimate_fitness <- function(data,chromosome,trait.1,trait.2=NULL,trait.3=NULL){
  # finding indeces/loci where genes are located
  chromosome.ind <- which(chromosome==1)
  ## no genes were in the chromsome
  fitness <- NULL
  if(length(chromosome.ind)==0){
    fitness <- NA
  }
  else{
    # computing eigengene
    eigengene <- prcomp(t(data[,chromosome.ind]))$x[,1]
    fitness <- cor(trait.1,eigengene)
  }
  ## change here if you want to modify how the fitness is extimated
  ## e.g.
  #cor.trait.1 <- cor(trait.1,eigengene)
  #cor.trait.2 <- cor(trait.2,eigengene)
  #cor.trait.3 <- cor(trait.3,eigengene)
  #fitness <- cor.trait.1 * cor.trait.2 * -(cor.trait.3)
  
  ## in these cases we have a positive outcome although we don't want to
  ## multiply fitness by -1 to negate the results
  # if(cor.trait.1 < 0 & cor.trait.2 < 0 & cor.trait.3 < 0)
  #     fitness <- fitness * -1
  # else if(cor.trait.1 > 0 & cor.trait.2 < 0 & cor.trait.3 > 0)
  #     fitness <- fitness * -1
  # else if(cor.trait.1 < 0 & cor.trait.2 > 0 & cor.trait.3 > 0)
  #     fitness <- fitness * -1
  
  return(fitness)
} ## end function
## --------------------------------------------------FUNCTIONS END-----------------------------------------------------------------------##
