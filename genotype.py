class Genotype:
    """A genotype class

    A genotype is composed of one or more alleles, depending on ploidy. Alleles
    must be represented by integers. Typically the number of alleles is 2, for
    diploid organisms.
    In addition to the constituent alleles, each genotype has a number of
    characteristics describing its fitness. These are the relative fitness of
    the genotype in the selected environment (e.g. exposure to a PIP-expressing
    crop), its maximumn growth rate and maximum population density. These
    latter two can be used to model fitness costs of an otherwise favoured
    genotype.    
    """

    def __init__(self,
                 alleles,
                 fitness=1.0,
                 growthRate=1.0,
                 maxDensity=1.0):
        assert all(isinstance(x, int) for x in alleles),\
            "Alleles must be a list of integers"
        self.alleles = alleles
        self.alleles.sort()
        self.fitness = float(fitness)
        assert (self.fitness >= 0) and (self.fitness <= 1.0),\
            "Fitness must be between 0 and 1"
        self.growthRate = float(growthRate)
        self.maxDensity = float(maxDensity)
        assert (self.maxDensity > 0.0) and (self.maxDensity <= 1.0),\
            "Maximum density must be greater than 0 and no greater than 1"

    def toString(self):
        """Gives string representation of the genoptype.
        
        Genotypes are represented by their alleles separated by slashes
        """
        return "/".join([str(x) for x in self.alleles])

    def alleleCount(self, allele):
        """Count the copies of the specified allele"""
        return sum([x == allele for x in self.alleles])

    def alleleFreq(self, allele):
        """Get the requency of the specified allele

        For example, in the case of the the genotype 0/1, the frequency of
        allele 1 is 0.5
        """
        return float(self.alleleCount(allele)) / len(self.alleles)

    def logisticGrowth(self, density, totalDensity):
        """Calculate the change in density of the genotype according to a
        logistic growth model.
        
        Density of the genotype is increased by
        
        D_t+1 = lambda D_t
        
        where

        lambda = 1+r - ((r/K)*N)

        and

        r, K and N are the genotype-specific growth rate, the genotype-specific
        maximum density and the total population density, respectively.
        """
        assert (float(totalDensity) > 0.0) and (float(totalDensity) <= 1.0),\
            "Total population density out of range"
        lam = 1.0 + self.growthRate
        lam = lam - ((self.growthRate / self.maxDensity) * float(totalDensity))
        return density * lam

    
