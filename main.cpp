#include "pearson.h"

#include <array>
#include <chrono>
#include <iostream>
#include <random>

using namespace std;

class Allele {
public:
    uint type;

    /// Creates an allele with a random type, the possible number of types
    /// depend on the number of types parameter.
    explicit Allele(uint numberOfTypes) :type( rand() % numberOfTypes ) {}

    /// Define equal operator on Allele, returns true if type is the same
    bool operator==(const Allele & other) const { return other.type == type; }
};


struct Locus {
public:
    array<Allele, 2> alleli;

    /// Creates a locus where the alleles have the specified number of types
    explicit Locus( int numberOfTypes ) : alleli({
            Allele(numberOfTypes),
            Allele(numberOfTypes)
        }){}

    /// Returns wethere the 2 alleli are equal in this locus
    bool isEterozigote() const {
        return !(alleli[0] == alleli[1]);
    }

    /// Computes normalized (0..1) similarity between 2 loci
    double similarity( const Locus & l ) const {
        int c = 0;
        c += alleli[0].type == l.alleli[0].type;
        c += alleli[0].type == l.alleli[1].type;
        c += alleli[1].type == l.alleli[0].type;
        c += alleli[1].type == l.alleli[1].type;
        return c / 4.0;
    }
};


struct Entity {
public:
    /// To change loci definition change loci count and modify the array
    /// accordingly, Locus constructor get the number of types of Allele
    /// for that Locus.
    static constexpr uint lociCount = 30;
    std::array<Locus, lociCount> loci = {
        Locus(3), Locus(3), Locus(3), Locus(3), Locus(3), Locus(3),
        Locus(3), Locus(3), Locus(3), Locus(3), Locus(3), Locus(3),
        Locus(3), Locus(3), Locus(3), Locus(3), Locus(3), Locus(3),
        Locus(3), Locus(3), Locus(3), Locus(3), Locus(3), Locus(3),
        Locus(3), Locus(3), Locus(3), Locus(3), Locus(3), Locus(3)
    };

    /// Returns a normalized (0..1) value representing how many loci are
    /// eterozygotes compared to the total number of loci in this entity
    double eterozigosity() const {
        double count = 0;
        for (const auto & locus : loci ) {
            count += locus.isEterozigote();
        }
        return count/(double)lociCount;
    }

    /// Return a normalized (0..1) value representing the average locus
    /// similarity between 2 entities
    double similarity( const Entity & e2 ) const {
        auto e1 = *this;
        double count = 0;
        for (uint i = 0; i < lociCount; i++) {
            auto locus1 = e1.loci[i];
            auto locus2 = e2.loci[i];
            count += locus1.similarity(locus2);
        }
        return count / (double)lociCount;
    }
};


struct Population {
public:
    uint size;
    std::vector<Entity> entities;

    /// Constructor: Creates a population of the given size
    explicit Population(uint _size) : size(_size) {
        entities.resize(size);
    }

    /// Performs a tests where the given number of Entities pairs are selected
    /// from the population and the correlation between eterozigosity and
    /// similarity is returned.
    double eteroSimilarityCorrelation( uint selectCount ) const {
        vector<double> eterozigosity;
        vector<double> similarity;
        eterozigosity.reserve(selectCount);
        similarity.reserve(selectCount);

        for (uint i = 0; i < selectCount; ) {
            const Entity & e1 = entities[rand() % size];
            const Entity & e2 = entities[rand() % size];
            if (&e1 == &e2) continue; // Skip comparison between the same entity
            eterozigosity.push_back(e1.eterozigosity());
            similarity.push_back(e1.similarity(e2));
            i++;
        }
        return pearsoncoeff(eterozigosity, similarity);
    }

    /// Perform the specified number of correlation tests where the specified
    /// number of entities are randomly selected in the population and returns
    /// a vector of correlations, each value represents the correlation in each
    /// performed test.
    vector<double> performEteroSimilarityTests(uint numberOfTests, uint selectCount) const {
        vector<double> testCorrelations;
        testCorrelations.reserve(numberOfTests);
        for (uint j = 0; j < numberOfTests; j++ ) {
            testCorrelations.push_back( eteroSimilarityCorrelation(selectCount) );
        }
        return testCorrelations;
    }
};


int main( int argc, char *argv[] )
{
    /// Seed random number generator with time so it changes at every run
    //srand(time(NULL));

    uint populationSize =       argc < 2 ? 1000 : stoi(argv[1]);
    uint numberOfTests =        argc < 3 ? 1000 : stoi(argv[2]);
    uint numberOfSelections =   argc < 4 ? 1000 : stoi(argv[3]);

    Population population(populationSize);

    vector<double> testsCorrelation = population.performEteroSimilarityTests(numberOfTests, numberOfSelections);
    cout << "Correlation mean  : " << mean(testsCorrelation) << endl;
    cout << "Correlation stdev : " << stdev(testsCorrelation) << endl;
}

