#include "pearson.h"
#include <array>
#include <chrono>
#include <iostream>
#include <random>

using namespace std;


class Allele {
public:
    explicit Allele(int types) {
        type = rand() % types;
    }

    int type;
};

struct Locus {
public:
    array<Allele, 2> alleli;

    explicit Locus( int types ) : alleli({Allele(types), Allele(types)}) {}

    bool isEterozigote() const {
        return alleli[0].type != alleli[1].type;
    }

    double similarity( const Locus & l ) const {
        // AA AA -> Returns 4 -> OK?
        // AB AA -> Returns 2 -> OK?
        // AB AB -> Returns 2 -> OK?
        // BA AB -> Returns 2 -> OK?
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
    static constexpr uint lociCount = 30;

    std::array<Locus, lociCount> loci = {
        Locus(3), Locus(3), Locus(3), Locus(3), Locus(3), Locus(3),
        Locus(3), Locus(3), Locus(3), Locus(3), Locus(3), Locus(3),
        Locus(3), Locus(3), Locus(3), Locus(3), Locus(3), Locus(3),
        Locus(3), Locus(3), Locus(3), Locus(3), Locus(3), Locus(3),
        Locus(3), Locus(3), Locus(3), Locus(3), Locus(3), Locus(3)
    };

    double eterozigosity() const {
        double count = 0;
        for (const auto & locus : loci ) {
            count += locus.isEterozigote();
        }
        return count/(double)lociCount;
    }

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
    static constexpr uint size = 10000;
    std::array<Entity, size> entities;

    void geneticSimilarity() const {
        /*for (const Entity & e1 : entities) {
            for (const Entity & e2 : entities) {
                if (&e1 == &e2) continue;
                cout << e1.eterozigosity() << "\t" << e2.eterozigosity() << "\t" << e1.similarity(e2) << "\n";
            }
        }*/
        vector<double> pearson;
        pearson.reserve(1000);
        for (int j = 0; j < 1000; j++ ) {
            vector<double> eterozigosity;
            vector<double> similarity;
            eterozigosity.reserve(10000);
            similarity.reserve(10000);

            for (uint i = 0; i < 10000; i++) {
                const Entity & e1 = entities[rand() % size];
                const Entity & e2 = entities[rand() % size];
                if (&e1 == &e2) continue;
                eterozigosity.push_back(e1.eterozigosity());
                similarity.push_back(e1.similarity(e2));

                // Print data of sigle comparison:
                //cout << e1.eterozigosity() << "\t" << e1.similarity(e2) << "\n";
            }
            cout << "PEARSON: " << pearsoncoeff(eterozigosity, similarity) << "\n";
            pearson.push_back( pearsoncoeff(eterozigosity, similarity));
        }

        cout << "AVG: " << mean(pearson) << endl;
        cout << "STD: " << stdev(pearson) << endl;
    }
};

int main()
{
    srand(time(NULL));
    Population e;
    e.geneticSimilarity();
}

