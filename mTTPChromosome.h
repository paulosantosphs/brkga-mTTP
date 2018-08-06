//
// Created by paulo on 28/02/18.
//

#ifndef BRKGA_MTTP_MTTPCHROMOSOME_H
#define BRKGA_MTTP_MTTPCHROMOSOME_H


#include <vector>

class mTTPChromosome {
public:

    mTTPChromosome();

    mTTPChromosome(const std::vector<double> &chromosome, const std::vector<int> &permutation,
                   const std::vector<std::vector<int>> &haps, double fitness);

    void setChromosome(const std::vector<double> &chromosome);

    void setPermutation(const std::vector<int> &permutation);

    void setHaps(const std::vector<std::vector<int>> &haps);

    void setFitness(double fitness);

    const std::vector<double> &getChromosome() const;

    const std::vector<int> &getPermutation() const;

    const std::vector<std::vector<int>> &getHaps() const;

    long getBeforeVND() const;

    void setBeforeVND(long beforeVND);

    long getAfterVND() const;

    void setAfterVND(long afterVND);

    int getCount() const;

    void setCount(int count);

    double getFitness() const;

private:
    std::vector<double> chromosome;
    std::vector<int> permutation;
    std::vector<std::vector<int>> haps;
    double fitness;
    long int beforeVND;
    long int afterVND;
    int count;
};


#endif //BRKGA_MTTP_MTTPCHROMOSOME_H
