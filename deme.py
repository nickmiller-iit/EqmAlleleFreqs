"""Represents a deme"""

from genotype import Genotype


class Neighbour:
    """Contains information about a deme's neighbour

    Demes can send migrants to their neighbours at a specified rate.
    To avoid side effects, Neighbour instance should not directly modify
    genotype densities either in the parent population that contains them
    nor in the neighbouring population to which they point.
    """
    def __init__(self,
                 neighbour,
                 migrationRate):
        assert isinstance(neighbour, Deme),\
            "Neighbours must be an instance of Deme"
        assert isinstance(migrationRate, float),\
            "Migration rate should be a float between 0.0 and 1.0"
        assert (migrationRate >= 0.0) and (migrationRate <= 1.0),\
            "Migration rate should be a float between 0.0 and 1.0"
        self.neighbour = neighbour
        self.migrationRate = migrationRate
        self.migrantPool = []

    def fillPool(self,
                 sourceDensities):
        """Fills the migrant pool

        sourceDensities is the current list of genotype densities in the
        source population.
        """
        self.migrantPool = [x * self.migrationRate for x in sourceDensities]

    def sendMigrants(self):
        self.neighbour.receiveMigrants(self.migrantPool)

    def emptyPool(self):
        self.migrantPool = []


class Deme:
    """A deme contains a population of individuals of different genotypes
    
    """

    def __init__(self,
                 genotypeList,
                 genotypeDensities,
                 refugeProportion=0.0):
        """Creates a new deme
    
        Creating a new deme requires a reference to a list of Genotype objects
        and a list of initial genotype densities. The two lists should be in
        the same order, ie genotypeDensities[0] represents the initial density
        of genotypeList[0]
        """
        assert all(isinstance(g, Genotype) for g in genotypeList),\
            "Elements of genotypeList must be instances of Genotype"
        self.genotypes = genotypeList
        assert all(isinstance(f,
                              float)
                   for f in genotypeDensities),\
            "Elements of genotypeDensities must be type float"
        assert sum(genotypeDensities) <= 1.0,\
            "Total genotype densities must not exceed 1.0"
        self.genotypeDensities = genotypeDensities
        assert isinstance(refugeProportion, float),\
            "Refuge proprotion must be a float"
        assert (refugeProportion >= 0.0) and (refugeProportion <= 1.0),\
            "Refuge proportion must be between 0.0 and 1.0"
        self.refugeProportion = refugeProportion
        # set up list of alleles and initial frequencies
        # Note that because allele states are integers we simply
        # use the alleles integer value to index its frequency
        # in the frequencies list
        self.alleles = []
        for g in self.genotypes:
            self.alleles += g.alleles
        self.alleles = sorted(set(self.alleles))
        self.alleleFreqs = [0.0 for x in range(max(self.alleles) + 1)]
        self.updateAlleleFreqs()
        # Stuff that needs to be set up later
        self.neighbours = []

    def currentDensity(self):
        """Get the total population density"""
        return sum(self.genotypeDensities)
        
    def genotypeFreqs(self):
        """Get the relative frequencies of each genotype"""
        tot = self.currentDensity()
        return [x / tot for x in self.genotypeDensities]

    def updateAlleleFreqs(self):
        """Updates the allele frequencies from the current genotype frequencies
        """
        gFreqs = self.genotypeFreqs()
        for a in self.alleles:
            freq = 0.0
            for g in range(len(self.genotypes)):
                freq += self.genotypes[g].alleleFreq(a) * gFreqs[g]
            self.alleleFreqs[a] = freq

    def randomMating(self):
        """Carry out a round of random mating
        
        IMPORTANT: This assumes diploid genotypes.
        Genotype densities are updated based on allele frequencies
        and H-W. The total density is not changed as a result of
        mating.
        """
        self.updateAlleleFreqs()
        totDensity = self.currentDensity()
        gIdx = [x for x in range(len(self.genotypes))]
        newDensities = [0.0 for g in gIdx]
        for g in gIdx:
            a = self.genotypes[g].alleles
            if a[0] == a[1]:
                newDensities[g] = self.alleleFreqs[a[0]] ** 2
            else:
                newDensities[g] = self.alleleFreqs[a[0]] *\
                                  self.alleleFreqs[a[1]] * 2
        newDensities = [x * totDensity for x in newDensities]
        for g in gIdx:
            self.genotypeDensities[g] = newDensities[g]

    def logisticGrowth(self):
        """Increase density of each genotype according to a logistic growth \
        model

        New density of each genotype is
        
        D_t+1 = lambda D_t
        
        where

        lambda = 1+r - ((r/K)*N)

        and

        r, K and N are the genotype-specific growth rate, the genotype-specific
        maximum density and the total population density, respectively.
        """
        totDensity = self.currentDensity()
        for g in range(len(self.genotypes)):
            lam = 1.0 + self.genotypes[g].growthRate
            lam = lam - ((self.genotypes[g].growthRate /
                          self.genotypes[g].maxDensity) * totDensity)
            self.genotypeDensities[g] = self.genotypeDensities[g] * lam

    def selection(self):
        """Applies selection to genotypes

        A simple viability selection model is applied. The proportion of
        survivors is simply the genotypes fitness. If a fraction of the deme
        is refuge then that proportion of each genotype escapes selection.
        """
        selected = [x * (1.0 - self.refugeProportion)
                    for x in self.genotypeDensities]
        unselected = [x * self.refugeProportion
                      for x in self.genotypeDensities]
        gIdx = [x for x in range(len(self.genotypes))]
        for g in gIdx:
            selected[g] = selected[g] * self.genotypes[g].fitness
        newDensities = [sum(x) for x in zip(selected, unselected)]
        for g in gIdx:
            self.genotypeDensities[g] = newDensities[g]

    def addNeighbour(self,
                     neighbour,
                     migrationRate):
        """Add a new neighbour"""
        self.neighbours.append(Neighbour(neighbour,
                                         migrationRate))
        
    def receiveMigrants(self,
                        migrants):
        """Adds incoming migrants to the population

        The migrants argument should be a list of genotype densities
        in the same order as the deme's genotypeDensities list.
        """
        assert len(self.genotypeDensities) == len(migrants),\
            "Mismatch between migrants and destination"
        for g in range(len(self.genotypes)):
            self.genotypeDensities[g] += migrants[g]

    def fillMigrantPools(self):
        """Fills migrant pools for all neighbours

        The migrant pools of any neighbours are filled and the migrants
        are removed from the population.
        """
        for n in range(len(self.neighbours)):
            self.neighbours[n].fillPool(self.genotypeDensities)
        for n in range(len(self.neighbours)):
            for g in range(len(self.genotypeDensities)):
                self.genotypeDensities[g] -= self.neighbours[n].migrantPool[g]

    def sendMigrants(self):
        """Sends migrants to all neighbours

        Migrants are sent to each neighbour. Migrant pools are then emptied.
        """
        for n in range(len(self.neighbours)):
            self.neighbours[n].sendMigrants()
            self.neighbours[n].emptyPool()
            

