//
// Created by paulo on 17/02/18.
//

#ifndef BRKGA_MTTP_MTTPDECODER_H
#define BRKGA_MTTP_MTTPDECODER_H

#include <list>
#include <vector>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include "brkgaAPI/MTRand.h"
#include "mTTPChromosome.h"

class mTTPDecoder {
public:

    mTTPDecoder(const std::vector<std::vector<int>> _distances, mTTPChromosome &hapChromosome, MTRand &rng,
                unsigned lsPattern);

    ~mTTPDecoder();

    double decode(std::vector<double> &chromosome) const;

    double applyLocalSearch(std::vector<double> &chromosome, int fitness) const;

    int totalDistanceTravelled(std::vector<std::vector<int>> table, std::vector<std::vector<int>> distances) const;

private:
    std::vector<std::vector<int>> distances;

    mTTPChromosome &hapChromosome;

    MTRand &rng;

    unsigned lsPattern;

    std::vector<std::vector<int>> polygon(std::vector<int> permutation) const;

    std::vector<std::vector<int>> tableSRR(std::vector<std::vector<int>> rounds) const;

    std::vector<std::vector<int>> hapsAssignment(std::vector<std::vector<int>> srr) const;


    bool verifyFeasibility(std::vector<std::vector<int>> mdrr) const;

    std::vector<std::vector<int>> hapsMatrix(std::vector<std::vector<int>> mdrr) const;

    std::vector<std::vector<int>>
    haSwap(std::vector<std::vector<int>> mdrr1, int &fitness, int team, int round, int r) const;

    std::vector<std::vector<int>>
    roundSwap(std::vector<std::vector<int>> mdrr1, int &fitness, int r1, int r2, int r) const;

    std::vector<std::vector<int>>
    teamSwap(std::vector<std::vector<int>> mdrr1, int &fitness, int team1, int team2, int t, int r) const;

    std::vector<std::vector<int>> partialRoundSwap(std::vector<std::vector<int>> mdrr1, int &fitness, int r1, int r2, int team, int t, int r) const;

    std::vector<std::vector<int>>
    partialTeamSwap(std::vector<std::vector<int>> mdrr1, int &fitness, int t1, int t2, int round,int t, int r) const;

    std::vector<std::vector<int>>
    haSwapMove(std::vector<std::vector<int>> mdrr1, int team, int round, int r) const;

    std::vector<std::vector<int>>
    roundSwapMove(std::vector<std::vector<int>> mdrr1, int r1, int r2, int r) const;

    std::vector<std::vector<int>>
    teamSwapMove(std::vector<std::vector<int>> mdrr1, int team1, int team2, int t, int r) const;

    std::vector<std::vector<int>> partialRoundSwapMove(std::vector<std::vector<int>> mdrr1, int &fitness, int r1, int r2, int team, int t, int r) const;

    std::vector<std::vector<int>>
    partialTeamSwapMove(std::vector<std::vector<int>> mdrr1, int &fitness, int t1, int t2, int round,int t, int r) const;

    std::vector<int> polygonInverse(std::vector<std::vector<int>> mdrr) const;
};


#endif //BRKGA_MTTP_MTTPDECODER_H
