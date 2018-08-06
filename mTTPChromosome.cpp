//
// Created by paulo on 28/02/18.
//

#include "mTTPChromosome.h"

mTTPChromosome::mTTPChromosome() {
    count = 0;
    beforeVND = 0;
    afterVND = 0;
}

const std::vector<double> &mTTPChromosome::getChromosome() const {
    return chromosome;
}

const std::vector<int> &mTTPChromosome::getPermutation() const {
    return permutation;
}

const std::vector<std::vector<int>> &mTTPChromosome::getHaps() const {
    return haps;
}

double mTTPChromosome::getFitness() const {
    return fitness;
}

void mTTPChromosome::setChromosome(const std::vector<double> &chromosome) {
    mTTPChromosome::chromosome = chromosome;
}

void mTTPChromosome::setPermutation(const std::vector<int> &permutation) {
    mTTPChromosome::permutation = permutation;
}

void mTTPChromosome::setHaps(const std::vector<std::vector<int>> &haps) {
    mTTPChromosome::haps = haps;
}

void mTTPChromosome::setFitness(double fitness) {
    mTTPChromosome::fitness = fitness;
}

mTTPChromosome::mTTPChromosome(const std::vector<double> &chromosome, const std::vector<int> &permutation,
                               const std::vector<std::vector<int>> &haps, double fitness) : chromosome(chromosome),
                                                                                            permutation(permutation),
                                                                                            haps(haps),
                                                                                            fitness(fitness) {
    count = 0;
    beforeVND = 0;
    afterVND = 0;
}

long mTTPChromosome::getBeforeVND() const {
    return beforeVND;
}

void mTTPChromosome::setBeforeVND(long beforeVND) {
    mTTPChromosome::beforeVND = beforeVND;
}

long mTTPChromosome::getAfterVND() const {
    return afterVND;
}

void mTTPChromosome::setAfterVND(long afterVND) {
    mTTPChromosome::afterVND = afterVND;
}

int mTTPChromosome::getCount() const {
    return count;
}

void mTTPChromosome::setCount(int count) {
    mTTPChromosome::count = count;
}


