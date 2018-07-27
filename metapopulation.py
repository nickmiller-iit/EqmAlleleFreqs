"""Classes to represent metapopulations of connected demes"""

from genotype import Genotype
from deme import Deme


class TwoDTorus:
    """A 2-dimensional stepping stone model, arranged on a torus

    At present, migration rates must be uniform in all directions and
    between all demes.
    """

    def __init__(self,
                 xlen,
                 ylen,
                 genotypes,
                 genotypeDensities,
                 migrationRate,
                 refugeProportion=0.0):
        """Create a new 2D torus steping stone model

        At a minimum, the following arguments must be provided:
        1. xlen, the number of demes aling the x dimension
        2. ylen, the number of demes along the y dimension
        3. genotypes, a list of genotype objects representing all
        possible genotypes in the simulation
        4. genotypDensities, a list of lists of genotype frequencies of length
        xlen * ylen, giving the initial genotype frequencies for each deme.
        5. migrationRate, currently a uniform migration rate between all demes
        in all directions must be specified

        Optionally, refugeProportion can specify the proportion of each deme
        that occupies untreated refuge. A single float will generate uniform
        refuge size for all demes, a list of floats of length xlen * ylen
        can be used to set variable refuge sizes
        """
        # Check parameter types
        assert isinstance(xlen, int),\
            "xlen should be type int"
        assert isinstance(ylen, int),\
            "ylen should be type int"
        assert all(isinstance(g, Genotype) for g in genotypes),\
            "Elements of genotypes must be instances of Genotype"
#        assert all(isinstance(f, float) for f in genotypeDensities),\
#            "Genotype densities should all be floats"
        assert isinstance(migrationRate, float),\
            "migrationRate should be type float"
        assert isinstance(refugeProportion, float) or\
            all(isinstance(x, float) for x in refugeProportion),\
            "refugeProportion should be type float or list of floats"
        self.xlen = xlen
        self.ylen = ylen
        # Set up refuge proportions
        if isinstance(refugeProportion, float):
            self.refugeProps = [refugeProportion for i in range(xlen * ylen)]
        else:
            assert len(refugeProportion) == xlen * ylen,\
                "list of refuge proportions does not match size of metapop"
            self.refugeProps = refugeProportion
        # Create demes
        d = []
        for x in range(xlen * ylen):
            d.append(Deme(genotypes,
                          genotypeDensities[x],
                          self.refugeProps[x]))
        # Place demes in a matrix
        self.demes = []
        for x in range(0, xlen * ylen, xlen):
            self.demes.append(d[x: x + xlen])
        # Set up neighbours and migration rates
        #
        # Start with rows
        for r in range(ylen):
            for c in range(xlen):
                # left neighbour
                if c == 0:
                    self.demes[r][c].addNeighbour(self.demes[r][xlen - 1],
                                                  migrationRate)
                else:
                    self.demes[r][c].addNeighbour(self.demes[r][c - 1],
                                                  migrationRate)
                # right neighbour
                if c == xlen - 1:
                    self.demes[r][c].addNeighbour(self.demes[r][0],
                                                  migrationRate)
                else:
                    self.demes[r][c].addNeighbour(self.demes[r][c + 1],
                                                  migrationRate)
        # Then columns
        for c in range(xlen):
            for r in range(ylen):
                # upper neighbour
                if r == 0:
                    self.demes[r][c].addNeighbour(self.demes[ylen - 1][c],
                                                  migrationRate)
                else:
                    self.demes[r][c].addNeighbour(self.demes[r - 1][c],
                                                  migrationRate)
                # lower neighbour
                if r == ylen - 1:
                    self.demes[r][c].addNeighbour(self.demes[0][c],
                                                  migrationRate)
                else:
                    self.demes[r][c].addNeighbour(self.demes[r + 1][c],
                                                  migrationRate)

    def randomMating(self):
        """Carry ouy random mating in each deme"""
        for r in range(self.ylen):
            for c in range(self.xlen):
                self.demes[r][c].randomMating()

    def logisticGrowth(self):
        """Carry out logistic population growth in each deme"""
        for r in range(self.ylen):
            for c in range(self.xlen):
                self.demes[r][c].logisticGrowth()

    def selection(self):
        """Carry out selection in each deme"""
        for r in range(self.ylen):
            for c in range(self.xlen):
                self.demes[r][c].selectio()

    def migration(self):
        """Carry out migration between demes"""
        for r in range(self.ylen):
            for c in range(self.xlen):
                self.demes[r][c].fillMigrantPools()
        for r in range(self.ylen):
            for c in range(self.xlen):
                self.demes[r][c].sendMigrants()

    def getGenotypeDensities(self):
        """Get genotype densities from each deme

        Returns a list of lists. Each element contains the deme's x position
        the deme's y position and the genotype frequencies
        """
        # TO DO
        pass




