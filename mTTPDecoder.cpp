//
// Created by paulo on 17/02/18.
//

#include <iostream>
#include <tuple>
#include "mTTPDecoder.h"


mTTPDecoder::mTTPDecoder(const std::vector<std::vector<int>> _distances, mTTPChromosome &hapChromosome, MTRand &rng,
                         unsigned lsPattern)
        : distances(_distances), hapChromosome(hapChromosome), rng(rng), lsPattern(lsPattern) {
}

mTTPDecoder::~mTTPDecoder() {}

double mTTPDecoder::decode(std::vector<double> &chromosome) const {
    int myFitness = 0.0;
    std::vector<double> caux = chromosome;
    typedef std::pair<double, unsigned> ValueKeyPair;
    std::vector<ValueKeyPair> rank(chromosome.size());


    for (unsigned i = 0; i < chromosome.size(); ++i) {
        rank[i] = ValueKeyPair(chromosome[i], i);
    }


    // Here we sort 'permutation', which will then produce a permutation of [n]
    // stored in ValueKeyPair::second:
    std::sort(rank.begin(), rank.end());
    //std::cout << "1:" << std::endl;
    // permutation[i].second is in {0, ..., n - 1}; a permutation can be obtained as follows
    std::vector<int> permutation;
    for (std::vector<ValueKeyPair>::const_iterator i = rank.begin(); i != rank.end(); ++i) {
        permutation.push_back(i->second + 1);
        //std::cout<<i->second + 1<<" ";
    }
    //std::cout << std::endl;

    //std::cout << "Poligono:" << std::endl;
    std::vector<std::vector<int>> rounds = polygon(permutation);
    //std::cout << "Tabela SRR: " << std::endl;
    std::vector<std::vector<int>> table = tableSRR(rounds);

    std::vector<std::vector<int>> mdrr = hapsAssignment(table);

    //Repete em média apenas 1 ou 2 vezes quando não encontra um mdrr viável de primeira
    while (!verifyFeasibility(mdrr)) {
        mdrr = hapsAssignment(table);
    }

    myFitness = totalDistanceTravelled(mdrr, distances);

    if (myFitness < hapChromosome.getFitness()) {
        hapChromosome.setFitness(myFitness);
        hapChromosome.setHaps(mdrr);
        hapChromosome.setPermutation(permutation);
        hapChromosome.setChromosome(chromosome);
    }

    return myFitness;
}

std::vector<std::vector<int>> mTTPDecoder::polygon(std::vector<int> round) const {
    std::vector<std::vector<int>> rounds;
    rounds.push_back(round);
    for (unsigned i = 1; i < round.size() - 1; i++) {
        int aux = round[1];
        for (unsigned j = 2; j < round.size(); j++) {
            round[j - 1] = round[j];
        }
        round[round.size() - 1] = aux;
        rounds.push_back(round);
    }
    return rounds;
}

std::vector<int> mTTPDecoder::polygonInverse(std::vector<std::vector<int>> mdrr) const {
    int k;
    bool var;
    int t1;
    int t2;
    int t3;
    int t4;
    int t5;
    int t6;
    std::vector<std::vector<int>> chrs;
    for (int i = 0; i < mdrr.size(); i++) {
        k = mdrr.size() - 2;
        var = true;
        for (int j = 1; j <= (mdrr.size() - 2) / 2; j++) {
            t1 = mdrr[i][j];
            if (t1 < 0) {
                t1 = t1 * -1;
            }
            t2 = mdrr[i][k];
            if (t2 < 0) {
                t2 = t2 * -1;
            }
            t3 = mdrr[t1 - 1][0];
            if (t3 < 0) {
                t3 = t3 * -1;
            }
            t4 = mdrr[t2 - 1][0];
            if (t4 < 0) {
                t4 = t4 * -1;
            }

            if (!((t3 == t2) && (t4 == t1))) {
                var = false;
            }
            k--;
        }
        if (var == true) {
            std::vector<int> chr;
            chr.push_back(i + 1);
            for (int c = 0; c < mdrr.size() - 1; c++) {
                t5 = mdrr[i][c];
                if (t5 < 0) {
                    t5 = t5 * -1;
                }
                chr.push_back(t5);
            }
            chrs.push_back(chr);
        }
    }
    std::vector<int> chr;
    std::vector<std::vector<int>> srr;
    if (chrs.size() == 0) {
        std::vector<int> p;
        p.push_back(0);
        return p;
    }
    for (int i = 0; i < chrs.size(); i++) {
        int c = 0;
        chr = chrs[i];
        srr = tableSRR(polygon(chr));
        for (int r = 0; r < srr.size() - 1; r++) {
            bool v = true;
            for (int t = 0; t < srr.size(); t++) {
                t6 = mdrr[t][r];
                if (t6 < 0) {
                    t6 = t6 * -1;
                }
                if (t6 != srr[t][r]) {
                    v = false;
                }
            }
            if (v == true) {
                c++;
            }
        }
        if (chrs.size() == 1) {
            return chr;
        } else if (chrs.size() == 2 && c == mdrr.size() - 1) {
            return chr;
        } else return chr;
    }
}

std::vector<std::vector<int>> mTTPDecoder::tableSRR(std::vector<std::vector<int>> rounds) const {
    std::vector<std::vector<int>> table;
    int n = rounds.size() + 1;
    //times
    for (int i = 0; i < n; i++) {
        std::vector<int> games;
        //rounds
        for (unsigned j = 0; j < rounds.size(); j++) {
            std::vector<int> round = rounds[j];
            //Verificar em cada round posição do time i

            int positionI = 0;
            int k = 0;
            for (auto &r: round) {
                if (r == i + 1) {
                    positionI = k;
                }
                k++;
            }

            //pegar o elemento da posição certa que joga com i

            if (positionI == 0) {
                games.push_back(round[1]);
            } else if (positionI == 1) {
                games.push_back(round[0]);
            } else if (positionI > 1 && positionI <= ((n / 2))) {
                games.push_back(round[n - positionI + 1]);
            } else if (positionI > ((n / 2))) {
                games.push_back(round[n - positionI + 1]);
            }
        }
        table.push_back(games);
    }
    return table;
}

std::vector<std::vector<int>> mTTPDecoder::hapsAssignment(std::vector<std::vector<int>> mdrr) const {
    std::vector<int> n;
    for (unsigned r = 0; r < mdrr.size() - 1; r++) {
        for (unsigned t = 0; t < mdrr.size(); t++) {
            bool v = false;
            for (int k = t; k >= 0; k--) {
                unsigned aux;
                if (mdrr[k][r] < 0) {
                    aux = mdrr[k][r] * -1;
                } else {
                    aux = mdrr[k][r];
                }
                if ((aux - 1) == t) {
                    v = true;
                }
            }

            int aux;
            if (mdrr[t][r] < 0) {
                aux = mdrr[t][r] * -1;
            } else {
                aux = mdrr[t][r];
            }

            std::vector<int> ht = mdrr[t];
            std::vector<int> hta = mdrr[aux - 1];
            if (r == 1) {
                for (unsigned i = 0; i < mdrr.size(); i++) {
                    n.push_back(1);
                }
            } else if (r == 2) {
                if ((ht[r - 1] < 0 && ht[r - 2] < 0) || (ht[r - 1] > 0 && ht[r - 2] > 0)) {
                    n[t] = 2;
                } else {
                    n[t] = 1;
                }

                if ((hta[r - 1] < 0 && hta[r - 2] < 0) || (hta[r - 1] > 0 && hta[r - 2] > 0)) {
                    n[aux - 1] = 2;
                } else {
                    n[aux - 1] = 1;
                }

            } else if (r > 2) {
                if ((ht[r - 1] < 0 && ht[r - 2] < 0 && ht[r - 3] < 0) ||
                    (ht[r - 1] > 0 && ht[r - 2] > 0 && ht[r - 3] > 0)) {
                    n[t] = 3;
                } else if ((ht[r - 1] < 0 && ht[r - 2] < 0) || (ht[r - 1] > 0 && ht[r - 2] > 0)) {
                    n[t] = 2;
                } else {
                    n[t] = 1;
                }

                if ((hta[r - 1] < 0 && hta[r - 2] < 0 && hta[r - 3] < 0) ||
                    (hta[r - 1] > 0 && hta[r - 2] > 0 && hta[r - 3] > 0)) {
                    n[aux - 1] = 3;
                } else if ((hta[r - 1] < 0 && hta[r - 2] < 0) || (hta[r - 1] > 0 && hta[r - 2] > 0)) {
                    n[aux - 1] = 2;
                } else {
                    n[aux - 1] = 1;
                }
            }

            if (!v && (mdrr[t][r] > 0)) {
                if (r == 0) {
                    if (rng.randInt(10) < 5) {
                        mdrr[t][r] = mdrr[t][r] * -1;
                    } else {
                        mdrr[aux - 1][r] = mdrr[aux - 1][r] * -1;
                    }
                } else if (r == (mdrr.size() - 2)) {
                    if (n[t] > n[aux - 1]) {
                        if (mdrr[t][r - 1] > 0) {
                            int ax = mdrr[t][r];
                            mdrr[t][r] = mdrr[t][r] * -1;
                            int r1 = mdrr[t][0] * -1;
                            int r2 = mdrr[t][1] * -1;
                            int r3 = mdrr[t][2] * -1;
                            int rbbl = mdrr[t][r - 2];
                            int rbl = mdrr[t][r - 1];
                            int rl = mdrr[t][r];

                            if ((r1 > 0 && r2 > 0 && rbl > 0 && rl > 0) || (r1 < 0 && r2 < 0 && rbl < 0 && rl < 0) ||
                                (r1 > 0 && r2 > 0 && r3 > 0 && rl > 0) || (r1 < 0 && r2 < 0 && r3 < 0 && rl < 0) ||
                                (r1 > 0 && rbbl > 0 && rbl > 0 && rl > 0) ||
                                (r1 < 0 && rbbl < 0 && rbl < 0 && rl < 0)) {
                                mdrr[t][r] = ax;
                                mdrr[aux - 1][r] = mdrr[aux - 1][r] * -1;
                            }

                        } else {
                            int ax = mdrr[aux - 1][r];
                            mdrr[aux - 1][r] = mdrr[aux - 1][r] * -1;
                            int r1 = mdrr[t][0] * -1;
                            int r2 = mdrr[t][1] * -1;
                            int r3 = mdrr[t][2] * -1;
                            int rbbl = mdrr[t][r - 2];
                            int rbl = mdrr[t][r - 1];
                            int rl = mdrr[t][r];

                            if ((r1 > 0 && r2 > 0 && rbl > 0 && rl > 0) || (r1 < 0 && r2 < 0 && rbl < 0 && rl < 0) ||
                                (r1 > 0 && r2 > 0 && r3 > 0 && rl > 0) || (r1 < 0 && r2 < 0 && r3 < 0 && rl < 0) ||
                                (r1 > 0 && rbbl > 0 && rbl > 0 && rl > 0) ||
                                (r1 < 0 && rbbl < 0 && rbl < 0 && rl < 0)) {
                                mdrr[aux - 1][r] = ax;
                                mdrr[t][r] = mdrr[t][r] * -1;
                            }
                        }
                    } else if (n[t] < n[aux - 1]) {
                        if (mdrr[aux - 1][r - 1] > 0) {
                            int ax = mdrr[aux - 1][r];
                            mdrr[aux - 1][r] = mdrr[aux - 1][r] * -1;
                            int r1 = mdrr[t][0] * -1;
                            int r2 = mdrr[t][1] * -1;
                            int r3 = mdrr[t][2] * -1;
                            int rbbl = mdrr[t][r - 2];
                            int rbl = mdrr[t][r - 1];
                            int rl = mdrr[t][r];

                            if ((r1 > 0 && r2 > 0 && rbl > 0 && rl > 0) || (r1 < 0 && r2 < 0 && rbl < 0 && rl < 0) ||
                                (r1 > 0 && r2 > 0 && r3 > 0 && rl > 0) || (r1 < 0 && r2 < 0 && r3 < 0 && rl < 0) ||
                                (r1 > 0 && rbbl > 0 && rbl > 0 && rl > 0) ||
                                (r1 < 0 && rbbl < 0 && rbl < 0 && rl < 0)) {
                                mdrr[aux - 1][r] = ax;
                                mdrr[t][r] = mdrr[t][r] * -1;
                            }

                        } else {
                            int ax = mdrr[t][r];
                            mdrr[t][r] = mdrr[t][r] * -1;
                            int r1 = mdrr[t][0] * -1;
                            int r2 = mdrr[t][1] * -1;
                            int r3 = mdrr[t][2] * -1;
                            int rbbl = mdrr[t][r - 2];
                            int rbl = mdrr[t][r - 1];
                            int rl = mdrr[t][r];

                            if ((r1 > 0 && r2 > 0 && rbl > 0 && rl > 0) || (r1 < 0 && r2 < 0 && rbl < 0 && rl < 0) ||
                                (r1 > 0 && r2 > 0 && r3 > 0 && rl > 0) || (r1 < 0 && r2 < 0 && r3 < 0 && rl < 0) ||
                                (r1 > 0 && rbbl > 0 && rbl > 0 && rl > 0) ||
                                (r1 < 0 && rbbl < 0 && rbl < 0 && rl < 0)) {
                                mdrr[t][r] = ax;
                                mdrr[aux - 1][r] = mdrr[aux - 1][r] * -1;
                            }
                        }

                    } else {
                        if ((mdrr[t][r - 1] < 0) && (mdrr[aux - 1][r - 1] > 0)) {
                            int ax = mdrr[aux - 1][r];
                            mdrr[aux - 1][r] = mdrr[aux - 1][r] * -1;
                            int r1 = mdrr[t][0] * -1;
                            int r2 = mdrr[t][1] * -1;
                            int r3 = mdrr[t][2] * -1;
                            int rbbl = mdrr[t][r - 2];
                            int rbl = mdrr[t][r - 1];
                            int rl = mdrr[t][r];

                            if ((r1 > 0 && r2 > 0 && rbl > 0 && rl > 0) || (r1 < 0 && r2 < 0 && rbl < 0 && rl < 0) ||
                                (r1 > 0 && r2 > 0 && r3 > 0 && rl > 0) || (r1 < 0 && r2 < 0 && r3 < 0 && rl < 0) ||
                                (r1 > 0 && rbbl > 0 && rbl > 0 && rl > 0) ||
                                (r1 < 0 && rbbl < 0 && rbl < 0 && rl < 0)) {
                                mdrr[aux - 1][r] = ax;
                                mdrr[t][r] = mdrr[t][r] * -1;
                            }

                        } else if ((mdrr[t][r - 1] > 0) && (mdrr[aux - 1][r - 1] < 0)) {
                            int ax = mdrr[t][r];
                            mdrr[t][r] = mdrr[t][r] * -1;
                            int r1 = mdrr[t][0] * -1;
                            int r2 = mdrr[t][1] * -1;
                            int r3 = mdrr[t][2] * -1;
                            int rbbl = mdrr[t][r - 2];
                            int rbl = mdrr[t][r - 1];
                            int rl = mdrr[t][r];

                            if ((r1 > 0 && r2 > 0 && rbl > 0 && rl > 0) || (r1 < 0 && r2 < 0 && rbl < 0 && rl < 0) ||
                                (r1 > 0 && r2 > 0 && r3 > 0 && rl > 0) || (r1 < 0 && r2 < 0 && r3 < 0 && rl < 0) ||
                                (r1 > 0 && rbbl > 0 && rbl > 0 && rl > 0) ||
                                (r1 < 0 && rbbl < 0 && rbl < 0 && rl < 0)) {
                                mdrr[t][r] = ax;
                                mdrr[aux - 1][r] = mdrr[aux - 1][r] * -1;
                            }

                        } else {
                            if (rng.randInt(10) < 5) {
                                int ax = mdrr[t][r];
                                mdrr[t][r] = mdrr[t][r] * -1;
                                int r1 = mdrr[t][0] * -1;
                                int r2 = mdrr[t][1] * -1;
                                int r3 = mdrr[t][2] * -1;
                                int rbbl = mdrr[t][r - 2];
                                int rbl = mdrr[t][r - 1];
                                int rl = mdrr[t][r];

                                if ((r1 > 0 && r2 > 0 && rbl > 0 && rl > 0) ||
                                    (r1 < 0 && r2 < 0 && rbl < 0 && rl < 0) || (r1 > 0 && r2 > 0 && r3 > 0 && rl > 0) ||
                                    (r1 < 0 && r2 < 0 && r3 < 0 && rl < 0) ||
                                    (r1 > 0 && rbbl > 0 && rbl > 0 && rl > 0) ||
                                    (r1 < 0 && rbbl < 0 && rbl < 0 && rl < 0)) {
                                    mdrr[t][r] = ax;
                                    mdrr[aux - 1][r] = mdrr[aux - 1][r] * -1;
                                }
                            } else {
                                int ax = mdrr[aux - 1][r];
                                mdrr[aux - 1][r] = mdrr[aux - 1][r] * -1;
                                int r1 = mdrr[t][0] * -1;
                                int r2 = mdrr[t][1] * -1;
                                int r3 = mdrr[t][2] * -1;
                                int rbbl = mdrr[t][r - 2];
                                int rbl = mdrr[t][r - 1];
                                int rl = mdrr[t][r];

                                if ((r1 > 0 && r2 > 0 && rbl > 0 && rl > 0) ||
                                    (r1 < 0 && r2 < 0 && rbl < 0 && rl < 0) || (r1 > 0 && r2 > 0 && r3 > 0 && rl > 0) ||
                                    (r1 < 0 && r2 < 0 && r3 < 0 && rl < 0) ||
                                    (r1 > 0 && rbbl > 0 && rbl > 0 && rl > 0) ||
                                    (r1 < 0 && rbbl < 0 && rbl < 0 && rl < 0)) {
                                    mdrr[aux - 1][r] = ax;
                                    mdrr[t][r] = mdrr[t][r] * -1;
                                }
                            }
                        }
                    }

                } else {
                    if (n[t] > n[aux - 1]) {
                        if (mdrr[t][r - 1] > 0) {
                            mdrr[t][r] = mdrr[t][r] * -1;

                        } else {
                            mdrr[aux - 1][r] = mdrr[aux - 1][r] * -1;

                        }
                    } else if (n[t] < n[aux - 1]) {
                        if (mdrr[aux - 1][r - 1] > 0) {
                            mdrr[aux - 1][r] = mdrr[aux - 1][r] * -1;

                        } else {
                            mdrr[t][r] = mdrr[t][r] * -1;
                        }

                    } else {
                        if ((mdrr[t][r - 1] < 0) && (mdrr[aux - 1][r - 1] > 0)) {
                            mdrr[aux - 1][r] = mdrr[aux - 1][r] * -1;

                        } else if ((mdrr[t][r - 1] > 0) && (mdrr[aux - 1][r - 1] < 0)) {
                            mdrr[t][r] = mdrr[t][r] * -1;

                        } else {
                            if (rng.randInt(10) < 5) {
                                mdrr[t][r] = mdrr[t][r] * -1;
                            } else {
                                mdrr[aux - 1][r] = mdrr[aux - 1][r] * -1;
                            }
                        }
                    }
                }
            }
        }
    }
    for (unsigned i = 0; i < mdrr.size(); i++) {
        std::vector<int> t = mdrr[i];
        std::vector<int> aux = t;
        for (unsigned j = 0; j < t.size(); j++) {
            aux.push_back(t[j] * -1);
        }
        mdrr[i] = aux;
    }
    return mdrr;
}

int mTTPDecoder::totalDistanceTravelled(std::vector<std::vector<int>> table,
                                        std::vector<std::vector<int>> distances) const {
    int totalDistance = 0;
    int i = 0;
    for (auto &team: table) {
        int distance = 0;
        unsigned j = 0;
        for (j = 0; j < team.size(); j++) {
            if ((j == 0)) {
                if (team[0] < 0) {
                    distance = distances[i][-1 * team[j] - 1];
                }
            } else if (team[j - 1] < 0) {
                if (team[j] < 0) {
                    distance = distance + distances[(-1 * team[j - 1]) - 1][(-1 * team[j]) - 1];
                } else {
                    distance = distance + distances[(-1 * team[j - 1]) - 1][i];
                }
            } else if (team[j - 1] > 0) {
                if (team[j] < 0) {
                    distance = distance + (distances[(-1 * team[j]) - 1])[i];
                }
            }

            if (j == (team.size() - 1)) {
                if (team[j] < 0) {
                    distance = distance + distances[(-1 * team[j]) - 1][i];
                }
            }
        }
        totalDistance = totalDistance + distance;
        i++;
    }
    return totalDistance;
}

bool mTTPDecoder::verifyFeasibility(std::vector<std::vector<int>> mdrr) const {
    bool var = true;
    for (auto &t: mdrr) {
        int c = 1;
        for (unsigned i = 1; i < t.size(); i++) {
            if ((t[i] < 0 && t[i - 1] < 0) || (t[i] > 0 && t[i - 1] > 0)) {
                c++;
                if (c > 3) {
                    var = false;
                }
                //std::cout<<"R: "<<i+1<<" c: "<<c<<std::endl;
            } else {
                c = 1;
            }
        }
    }
    return var;
}

std::vector<std::vector<int>> mTTPDecoder::hapsMatrix(std::vector<std::vector<int>> mdrr) const {
    std::vector<std::vector<int>> haps;
    for (auto &t: mdrr) {
        std::vector<int> hap;
        for (auto &r: t) {
            if (r > 0) {
                hap.push_back(1);
            } else {
                hap.push_back(-1);
            }
        }
        haps.push_back(hap);
    }
    return haps;
}

std::vector<std::vector<int>>
mTTPDecoder::haSwap(std::vector<std::vector<int>> mdrr1, int &fitness, int team, int round, int r) const {
    std::vector<std::vector<int>> mdrr = mdrr1;
    std::vector<std::vector<int>> mdrra = mdrr1;
    int hasFitness = fitness;
    int team2;
    int fit;

    team2 = mdrr[team][round];
    if (team2 < 0) {
        team2 = team2 * -1;
    }
    team2--;
    mdrr[team][round] = mdrr[team][round] * -1;
    mdrr[team2][round] = mdrr[team2][round] * -1;

    mdrr[team][round + r] = mdrr[team][round + r] * -1;
    mdrr[team2][round + r] = mdrr[team2][round + r] * -1;

    fit = totalDistanceTravelled(mdrr, distances);
    if (fit < hasFitness && verifyFeasibility(mdrr)) {
        hasFitness = fit;
        mdrra = mdrr;
    } else {
        mdrr = mdrra;
    }

    fitness = hasFitness;
    return mdrr;
}

std::vector<std::vector<int>>
mTTPDecoder::roundSwap(std::vector<std::vector<int>> mdrr1, int &fitness, int r1, int r2, int r) const {
    std::vector<std::vector<int>> mdrr = mdrr1;
    std::vector<std::vector<int>> mdrra = mdrr1;
    int hasFitness = fitness;
    int aux;
    int fit;


    for (unsigned i = 0; i < mdrr.size(); i++) {
        aux = mdrr[i][r1];
        mdrr[i][r1] = mdrr[i][r2];
        mdrr[i][r2] = aux;

    }

    for (unsigned i = 0; i < mdrr.size(); i++) {
        for (int j = 0; j < r; j++) {
            mdrr[i][j + r] = mdrr[i][j] * -1;
        }
    }

    fit = totalDistanceTravelled(mdrr, distances);
    if (fit < hasFitness && verifyFeasibility(mdrr)) {
        hasFitness = fit;
        mdrra = mdrr;
    } else {
        mdrr = mdrra;
    }

    fitness = hasFitness;
    return mdrr;
}

std::vector<std::vector<int>>
mTTPDecoder::teamSwap(std::vector<std::vector<int>> mdrr1, int &fitness, int team1, int team2, int t, int r) const {
    std::vector<std::vector<int>> mdrr = mdrr1;
    std::vector<std::vector<int>> mdrra = mdrr1;
    int hasFitness = fitness;
    int fit;
    int team3;
    std::vector<int> aux;

    for (int i = 0; i < t; i++) {
        for (int j = 0; j < r; j++) {
            team3 = mdrr[i][j];
            if (team1 == i + 1) {
                if (team3 > 0) {
                    if (team3 == team2) {
                        mdrr[i][j] = team1;
                        mdrr[i][j + r] = mdrr[i][j] * -1;
                    }
                } else {
                    team3 = team3 * -1;
                    if (team3 == team2) {
                        mdrr[i][j] = team1 * -1;
                        mdrr[i][j + r] = mdrr[i][j] * -1;
                    }

                }
            } else if (team2 == i + 1) {
                if (team3 > 0) {
                    if (team3 == team1) {
                        mdrr[i][j] = team2;
                        mdrr[i][j + r] = mdrr[i][j] * -1;
                    }
                } else {
                    team3 = team3 * -1;
                    if (team3 == team1) {
                        mdrr[i][j] = team2 * -1;
                        mdrr[i][j + r] = mdrr[i][j] * -1;
                    }

                }
            } else {
                if (team3 > 0) {
                    if (team3 == team1) {
                        mdrr[i][j] = team2;
                        mdrr[i][j + r] = mdrr[i][j] * -1;
                    } else if (team3 == team2) {
                        mdrr[i][j] = team1;
                        mdrr[i][j + r] = mdrr[i][j] * -1;
                    }
                } else {
                    team3 = team3 * -1;
                    if (team3 == team1) {
                        mdrr[i][j] = team2 * -1;
                        mdrr[i][j + r] = mdrr[i][j] * -1;
                    } else if (team3 == team2) {
                        mdrr[i][j] = team1 * -1;
                        mdrr[i][j + r] = mdrr[i][j] * -1;
                    }
                }

            }
        }
    }
    aux = mdrr[team1 - 1];
    mdrr[team1 - 1] = mdrr[team2 - 1];
    mdrr[team2 - 1] = aux;

    fit = totalDistanceTravelled(mdrr, distances);
    if (fit < hasFitness && verifyFeasibility(mdrr)) {
        hasFitness = fit;
        mdrra = mdrr;
    } else {
        mdrr = mdrra;
    }

    fitness = hasFitness;
    return mdrr;

}

std::vector<std::vector<int>> mTTPDecoder::partialRoundSwap(std::vector<std::vector<int>> mdrr1, int &fitness, int r1, int r2, int team, int t, int r) const {
    std::vector<std::vector<int>> mdrr = mdrr1;
    std::vector<std::vector<int>> mdrra = mdrr1;
    int hasFitness = fitness;
    int fit;
    int team1;
    int team2;
    int taux;
    std::vector<int> aux;



    team1 = mdrr[team][r1];
    team2 = mdrr[team][r2];
    if (team1 < 0) {
        team1 = team1 * -1;
    }
    if (team2 < 0) {
        team2 = team2 * -1;
    }

    for (int i = 0; i < t; i++) {
        for (int j = 0; j < r; j++) {
            taux = mdrr[i][j];
            if (taux < 0) {
                taux = taux * -1;
            }
            if (taux == team1) {
                if (mdrr[i][j] < 0) {
                    mdrr[i][j] = team2 * -1;
                    mdrr[i][j + r] = team2;
                } else {
                    mdrr[i][j] = team2;
                    mdrr[i][j + r] = team2 * -1;
                }
            }
            if (taux == team2) {
                if (mdrr[i][j] < 0) {
                    mdrr[i][j] = team1 * -1;
                    mdrr[i][j + r] = team1;
                } else {
                    mdrr[i][j] = team1;
                    mdrr[i][j + r] = team1 * -1;
                }
            }
        }

    }

    aux = mdrr[team1 - 1];
    mdrr[team1 - 1] = mdrr[team2 - 1];
    mdrr[team2 - 1] = aux;

    fit = totalDistanceTravelled(mdrr, distances);
    if (fit < hasFitness && verifyFeasibility(mdrr)) {
        hasFitness = fit;
        mdrra = mdrr;
    } else {
        mdrr = mdrra;
    }

    //std::cout<<"Partial Round: "<<count<<std::endl;
    fitness = hasFitness;
    return mdrr;

}

std::vector<std::vector<int>> mTTPDecoder::partialTeamSwap(std::vector<std::vector<int>> mdrr1, int &fitness, int t1, int t2, int round,int t, int r) const {
    std::vector<std::vector<int>> mdrr = mdrr1;
    std::vector<std::vector<int>> mdrra = mdrr1;
    int hasFitness = fitness;
    int fit;
    int team;
    int team1;
    int team2;
    int team3;
    int team4;
    int taux;
    std::vector<int> aux;
    std::vector<int> aux2;



    team1 = mdrr[t1][round];
    team2 = mdrr[t2][round];
    team3 = t1 + 1;
    team4 = t2 + 1;

    if (team1 < 0) {
        team1 = team1 * -1;
    }
    if (team2 < 0) {
        team2 = team2 * -1;
    }
    if (team1 != team3 && team1 != team4 && team2 != team3 && team2 != team4) {
        for (int i = 0; i < t; i++) {
            for (int j = 0; j < r; j++) {
                taux = mdrr[i][j];
                if (taux < 0) {
                    taux = taux * -1;
                }
                if (taux == team1) {
                    if (mdrr[i][j] < 0) {
                        mdrr[i][j] = team2 * -1;
                        mdrr[i][j + r] = team2;
                    } else {
                        mdrr[i][j] = team2;
                        mdrr[i][j + r] = team2 * -1;
                    }
                }
                if (taux == team2) {
                    if (mdrr[i][j] < 0) {
                        mdrr[i][j] = team1 * -1;
                        mdrr[i][j + r] = team1;
                    } else {
                        mdrr[i][j] = team1;
                        mdrr[i][j + r] = team1 * -1;
                    }
                }
                if (taux == team3) {
                    if (mdrr[i][j] < 0) {
                        mdrr[i][j] = team4 * -1;
                        mdrr[i][j + r] = team4;
                    } else {
                        mdrr[i][j] = team4;
                        mdrr[i][j + r] = team4 * -1;
                    }
                }
                if (taux == team4) {
                    if (mdrr[i][j] < 0) {
                        mdrr[i][j] = team3 * -1;
                        mdrr[i][j + r] = team3;
                    } else {
                        mdrr[i][j] = team3;
                        mdrr[i][j + r] = team3 * -1;
                    }
                }
            }

        }

        aux = mdrr[team1 - 1];
        mdrr[team1 - 1] = mdrr[team2 - 1];
        mdrr[team2 - 1] = aux;

        aux2 = mdrr[t1];
        mdrr[t1] = mdrr[t2];
        mdrr[t2] = aux2;
    }

    fit = totalDistanceTravelled(mdrr, distances);
    if (fit < hasFitness && verifyFeasibility(mdrr)) {
        hasFitness = fit;
        mdrra = mdrr;
    } else {
        mdrr = mdrra;
    }

    //std::cout<<"Partial Team: "<<count<<std::endl;
    fitness = hasFitness;
    return mdrr;

}

std::vector<std::vector<int>>
mTTPDecoder::haSwapMove(std::vector<std::vector<int>> mdrr1, int team, int round, int r) const {
    std::vector<std::vector<int>> mdrr = mdrr1;
    int team2;

    team2 = mdrr[team][round];
    if (team2 < 0) {
        team2 = team2 * -1;
    }
    team2--;
    mdrr[team][round] = mdrr[team][round] * -1;
    mdrr[team2][round] = mdrr[team2][round] * -1;

    mdrr[team][round + r] = mdrr[team][round + r] * -1;
    mdrr[team2][round + r] = mdrr[team2][round + r] * -1;

    if(verifyFeasibility(mdrr)){
        return mdrr;
    }else{
        return mdrr1;
    }
}

std::vector<std::vector<int>>
mTTPDecoder::roundSwapMove(std::vector<std::vector<int>> mdrr1, int r1, int r2, int r) const {
    std::vector<std::vector<int>> mdrr = mdrr1;
    int aux;

    for (unsigned i = 0; i < mdrr.size(); i++) {
        aux = mdrr[i][r1];
        mdrr[i][r1] = mdrr[i][r2];
        mdrr[i][r2] = aux;

    }

    for (unsigned i = 0; i < mdrr.size(); i++) {
        for (int j = 0; j < r; j++) {
            mdrr[i][j + r] = mdrr[i][j] * -1;
        }
    }

    if(verifyFeasibility(mdrr)){
        return mdrr;
    }else{
        return mdrr1;
    }
}

std::vector<std::vector<int>>
mTTPDecoder::teamSwapMove(std::vector<std::vector<int>> mdrr1, int team1, int team2, int t, int r) const {
    std::vector<std::vector<int>> mdrr = mdrr1;
    int team3;
    std::vector<int> aux;

    for (int i = 0; i < t; i++) {
        for (int j = 0; j < r; j++) {
            team3 = mdrr[i][j];
            if (team1 == i + 1) {
                if (team3 > 0) {
                    if (team3 == team2) {
                        mdrr[i][j] = team1;
                        mdrr[i][j + r] = mdrr[i][j] * -1;
                    }
                } else {
                    team3 = team3 * -1;
                    if (team3 == team2) {
                        mdrr[i][j] = team1 * -1;
                        mdrr[i][j + r] = mdrr[i][j] * -1;
                    }

                }
            } else if (team2 == i + 1) {
                if (team3 > 0) {
                    if (team3 == team1) {
                        mdrr[i][j] = team2;
                        mdrr[i][j + r] = mdrr[i][j] * -1;
                    }
                } else {
                    team3 = team3 * -1;
                    if (team3 == team1) {
                        mdrr[i][j] = team2 * -1;
                        mdrr[i][j + r] = mdrr[i][j] * -1;
                    }

                }
            } else {
                if (team3 > 0) {
                    if (team3 == team1) {
                        mdrr[i][j] = team2;
                        mdrr[i][j + r] = mdrr[i][j] * -1;
                    } else if (team3 == team2) {
                        mdrr[i][j] = team1;
                        mdrr[i][j + r] = mdrr[i][j] * -1;
                    }
                } else {
                    team3 = team3 * -1;
                    if (team3 == team1) {
                        mdrr[i][j] = team2 * -1;
                        mdrr[i][j + r] = mdrr[i][j] * -1;
                    } else if (team3 == team2) {
                        mdrr[i][j] = team1 * -1;
                        mdrr[i][j + r] = mdrr[i][j] * -1;
                    }
                }

            }
        }
    }
    aux = mdrr[team1 - 1];
    mdrr[team1 - 1] = mdrr[team2 - 1];
    mdrr[team2 - 1] = aux;

    if(verifyFeasibility(mdrr)){
        return mdrr;
    }else{
        return mdrr1;
    }
}

std::vector<std::vector<int>> mTTPDecoder::partialRoundSwapMove(std::vector<std::vector<int>> mdrr1, int &fitness, int r1, int r2, int team, int t, int r) const {
    std::vector<std::vector<int>> mdrr = mdrr1;
    std::vector<std::vector<int>> mdrra = mdrr1;
    int team1;
    int team2;
    int taux;
    std::vector<int> aux;



    team1 = mdrr[team][r1];
    team2 = mdrr[team][r2];
    if (team1 < 0) {
        team1 = team1 * -1;
    }
    if (team2 < 0) {
        team2 = team2 * -1;
    }

    for (int i = 0; i < t; i++) {
        for (int j = 0; j < r; j++) {
            taux = mdrr[i][j];
            if (taux < 0) {
                taux = taux * -1;
            }
            if (taux == team1) {
                if (mdrr[i][j] < 0) {
                    mdrr[i][j] = team2 * -1;
                    mdrr[i][j + r] = team2;
                } else {
                    mdrr[i][j] = team2;
                    mdrr[i][j + r] = team2 * -1;
                }
            }
            if (taux == team2) {
                if (mdrr[i][j] < 0) {
                    mdrr[i][j] = team1 * -1;
                    mdrr[i][j + r] = team1;
                } else {
                    mdrr[i][j] = team1;
                    mdrr[i][j + r] = team1 * -1;
                }
            }
        }

    }

    aux = mdrr[team1 - 1];
    mdrr[team1 - 1] = mdrr[team2 - 1];
    mdrr[team2 - 1] = aux;

    verifyFeasibility(mdrr);
    fitness = totalDistanceTravelled(mdrr, distances);
    return mdrr;

}

std::vector<std::vector<int>> mTTPDecoder::partialTeamSwapMove(std::vector<std::vector<int>> mdrr1, int &fitness, int t1, int t2, int round,int t, int r) const {
    std::vector<std::vector<int>> mdrr = mdrr1;
    std::vector<std::vector<int>> mdrra = mdrr1;
    int hasFitness = fitness;
    int fit;
    int team;
    int team1;
    int team2;
    int team3;
    int team4;
    int taux;
    std::vector<int> aux;
    std::vector<int> aux2;



    team1 = mdrr[t1][round];
    team2 = mdrr[t2][round];
    team3 = t1 + 1;
    team4 = t2 + 1;

    if (team1 < 0) {
        team1 = team1 * -1;
    }
    if (team2 < 0) {
        team2 = team2 * -1;
    }
    if (team1 != team3 && team1 != team4 && team2 != team3 && team2 != team4) {
        for (int i = 0; i < t; i++) {
            for (int j = 0; j < r; j++) {
                taux = mdrr[i][j];
                if (taux < 0) {
                    taux = taux * -1;
                }
                if (taux == team1) {
                    if (mdrr[i][j] < 0) {
                        mdrr[i][j] = team2 * -1;
                        mdrr[i][j + r] = team2;
                    } else {
                        mdrr[i][j] = team2;
                        mdrr[i][j + r] = team2 * -1;
                    }
                }
                if (taux == team2) {
                    if (mdrr[i][j] < 0) {
                        mdrr[i][j] = team1 * -1;
                        mdrr[i][j + r] = team1;
                    } else {
                        mdrr[i][j] = team1;
                        mdrr[i][j + r] = team1 * -1;
                    }
                }
                if (taux == team3) {
                    if (mdrr[i][j] < 0) {
                        mdrr[i][j] = team4 * -1;
                        mdrr[i][j + r] = team4;
                    } else {
                        mdrr[i][j] = team4;
                        mdrr[i][j + r] = team4 * -1;
                    }
                }
                if (taux == team4) {
                    if (mdrr[i][j] < 0) {
                        mdrr[i][j] = team3 * -1;
                        mdrr[i][j + r] = team3;
                    } else {
                        mdrr[i][j] = team3;
                        mdrr[i][j + r] = team3 * -1;
                    }
                }
            }

        }

        aux = mdrr[team1 - 1];
        mdrr[team1 - 1] = mdrr[team2 - 1];
        mdrr[team2 - 1] = aux;

        aux2 = mdrr[t1];
        mdrr[t1] = mdrr[t2];
        mdrr[t2] = aux2;
    }

    verifyFeasibility(mdrr);
    fitness = totalDistanceTravelled(mdrr, distances);
    return mdrr;

}

double mTTPDecoder::applyLocalSearch(std::vector<double> &chromosome, int fitness) const {
    int myFitness = 0.0;
    std::vector<double> caux = chromosome;
    typedef std::pair<double, unsigned> ValueKeyPair;
    std::vector<ValueKeyPair> rank(chromosome.size());


    for (unsigned i = 0; i < chromosome.size(); ++i) {
        rank[i] = ValueKeyPair(chromosome[i], i);
    }


    // Here we sort 'permutation', which will then produce a permutation of [n]
    // stored in ValueKeyPair::second:
    std::sort(rank.begin(), rank.end());
    //std::cout << "1:" << std::endl;
    // permutation[i].second is in {0, ..., n - 1}; a permutation can be obtained as follows
    std::vector<int> permutation;
    for (std::vector<ValueKeyPair>::const_iterator i = rank.begin(); i != rank.end(); ++i) {
        permutation.push_back(i->second + 1);
        //std::cout<<i->second + 1<<" ";
    }
    //std::cout << std::endl;

    //std::cout << "Poligono:" << std::endl;

    std::vector<std::vector<int>> rounds = polygon(permutation);
    /*for (auto &r: rounds) {
        for (auto &rr: r) {
            std::cout << rr<<" ";
        }
        std::cout << std::endl;
    }
    std::cout << "Tabela SRR: " << std::endl;*/
    std::vector<std::vector<int>> table = tableSRR(rounds);
    /*int i = 1;
    for (auto &t: table) {
        std::cout << "Games de " << i << ": " << std::endl;
        for (auto &g: t) {
            std::cout << g<<" ";
        }
        i++;
        std::cout << std::endl;
    }*/
    std::vector<std::vector<int>> mdrr = hapsAssignment(table);

    //Repete em média apenas 1 ou 2 vezes quando não encontra um mdrr viável de primeira
    while (!verifyFeasibility(mdrr)) {
        mdrr = hapsAssignment(table);
    }


    // Return the fitness:
    int myFitnessOld = totalDistanceTravelled(mdrr, distances);

    int t = mdrr.size();
    int r = mdrr.size() - 1;


    std::vector<std::pair<int, int>> pairsHa;
    for (int i = 0; i < t; i++) {
        for (int j = 0; j < r; j++) {
            pairsHa.push_back(std::make_pair(i, j));
        }
    }
    auto rand = std::default_random_engine(time(NULL));

    std::shuffle(std::begin(pairsHa), std::end(pairsHa), rand);

    std::vector<std::pair<int, int>> pairsRs;
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < r; j++) {
            if (i > j) {
                pairsRs.push_back(std::make_pair(i, j));
            }
        }
    }
    std::shuffle(std::begin(pairsRs), std::end(pairsRs), rand);

    std::vector<std::pair<int, int>> pairsTs;
    for (int i = 1; i <= t; i++) {
        for (int j = 1; j <= t; j++) {
            if (i > j) {
                pairsTs.push_back(std::make_pair(i, j));
            }
        }
    }
    std::shuffle(std::begin(pairsTs), std::end(pairsTs), rand);

    std::vector<std::tuple<int, int, int>> tuplesPrs;
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < r; j++) {
            for (int k = 0; k < t; k++) {
                if (i > j) {
                    std::tuple<int, int, int> tri(i, j, k);
                    tuplesPrs.push_back(tri);
                }
            }
        }
    }
    std::shuffle(std::begin(tuplesPrs), std::end(tuplesPrs), rand);

    std::vector<std::tuple<int, int, int>> tuplesPts;
    for (int i = 0; i < t; i++) {
        for (int j = 0; j < t; j++) {
            for (int k = 0; k < r; k++) {
                if (i > j) {
                    std::tuple<int, int, int> tri(i, j, k);
                    tuplesPts.push_back(tri);
                }
            }
        }
    }
    std::shuffle(std::begin(tuplesPts), std::end(tuplesPts), rand);

    int bestFit = myFitnessOld;
    std::vector<std::vector<int>> bestMdrr = mdrr;

    int fitBL = myFitnessOld;
    std::vector<std::vector<int>> mdrrBL = mdrr;

    //Para piorias
    int fit = myFitnessOld;
    std::vector<std::vector<int>> mdrrP = mdrr;

    unsigned iterT = 0;
    unsigned k = 1;
    unsigned reheat = 0;

    while(reheat <= 3){

        double tstart = (-0.8*bestFit)/log(0.5);
        double tend = (-0.45*bestFit)/log(0.5);
        double temp = tstart;

        fit = bestFit;
        mdrrP = bestMdrr;

        while(temp > 0.01){
            while (iterT < 300){
                iterT++;

                fitBL = fit;
                mdrrBL = mdrrP;

                int n = rng.randInt(2) + 1;
                if(n==1){
                    int iham = rng.randInt(pairsHa.size() - 1);
                    int teamm = pairsHa[iham].first;
                    int roundm = pairsHa[iham].second;
                    mdrrBL = haSwapMove(mdrrBL,teamm, roundm,r);
                    fitBL = totalDistanceTravelled(mdrrBL, distances);
                    std::shuffle(std::begin(pairsHa), std::end(pairsHa), rand);
                }
                if(n==2){
                    int irsm = rng.randInt(pairsRs.size() - 1);
                    int r1m = pairsRs[irsm].first;
                    int r2m = pairsRs[irsm].second;
                    mdrrBL = roundSwapMove(mdrrBL, r1m, r2m, r);
                    fitBL = totalDistanceTravelled(mdrrBL, distances);
                    std::shuffle(std::begin(pairsRs), std::end(pairsRs), rand);
                }
                if(n==3){
                    int itsm = rng.randInt(pairsTs.size() - 1);
                    int team1m = pairsTs[itsm].first;
                    int team2m = pairsTs[itsm].second;
                    mdrrBL = teamSwapMove(mdrrBL,team1m,team2m,t,r);
                    fitBL = totalDistanceTravelled(mdrrBL, distances);
                    std::shuffle(std::begin(pairsTs), std::end(pairsTs), rand);
                }

                float delta = fitBL - fit;

                if(delta < 0){
                    fit = fitBL;
                    mdrrP = mdrrBL;
                    if(fitBL < bestFit){
                        bestFit = fitBL;
                        bestMdrr = mdrrBL;
                        tend = (-0.45*bestFit)/log(0.5);
                    }
                }else{
                    double x = rng.rand();
                    if (x < exp((-1*delta)/temp)){
                        fit = fitBL;
                        mdrrP = mdrrBL;
                    }
                }
            }
            iterT = 0;
            double kk = 1.0/k;
            double resfri = (pow((tend/tstart),kk));
            temp = temp * resfri;
            if(resfri > 0.999){
                temp = 0;
            }
            k++;
        }
        reheat++;

        /*VND:

        unsigned n = 1;
        unsigned iha = 0;
        unsigned irs = 0;
        unsigned its = 0;

        bool flag = true;

        while(flag) {
            int oldFit = fit;
            if (n == 1) {
                if (iha == pairsHa.size()) {
                    n = 2;
                } else {
                    int team = pairsHa[iha].first;
                    int round = pairsHa[iha].second;
                    mdrrP = haSwap(mdrrP, fit, team, round, r);
                    if (fit < oldFit) {
                        n = 1;
                        iha = 0;
                        std::shuffle(std::begin(pairsHa), std::end(pairsHa), rand);
                    } else if (iha <= pairsHa.size() - 1) {
                        n = 1;
                        iha++;
                    }
                }
            } else if (n == 2) {
                if (irs == pairsRs.size()) {
                    n = 3;
                } else {
                    int r1 = pairsRs[irs].first;
                    int r2 = pairsRs[irs].second;
                    mdrrP = roundSwap(mdrrP, fit, r1, r2, r);
                    if (fit < oldFit) {
                        n = 1;
                        iha = 0;
                        irs = 0;
                        std::shuffle(std::begin(pairsHa), std::end(pairsHa), rand);
                        std::shuffle(std::begin(pairsRs), std::end(pairsRs), rand);
                    } else if (irs <= pairsRs.size() - 1) {
                        n = 2;
                        irs++;
                    }
                }
            } else if (n == 3) {
                if (its == pairsTs.size()) {
                    flag = false;
                } else {
                    int team1 = pairsTs[its].first;
                    int team2 = pairsTs[its].second;
                    mdrrP = teamSwap(mdrrP, fit, team1, team2, t, r);
                    if (fit < oldFit) {
                        n = 1;
                        iha = 0;
                        irs = 0;
                        its = 0;
                        std::shuffle(std::begin(pairsHa), std::end(pairsHa), rand);
                        std::shuffle(std::begin(pairsRs), std::end(pairsRs), rand);
                        std::shuffle(std::begin(pairsTs), std::end(pairsTs), rand);
                    } else if (its <= pairsTs.size() - 1) {
                        n = 3;
                        its++;
                    }
                }
            }
        }*/

        if(fit <= bestFit){
            bestFit = fit;
            bestMdrr = mdrrP;
        }
    }

    myFitness = bestFit;
    mdrr = bestMdrr;

    if (!verifyFeasibility(mdrr)) {
        std::cout << "Inviavel" << std::endl;
    }

    if (myFitness <= fitness * 1.05) {
        if (myFitness < hapChromosome.getFitness()) {
            std::cout<<myFitness<<std::endl;
            hapChromosome.setFitness(myFitness);
            hapChromosome.setHaps(mdrr);
            hapChromosome.setPermutation(permutation);
            hapChromosome.setChromosome(chromosome);
            hapChromosome.setBeforeVND(hapChromosome.getBeforeVND() + myFitnessOld);
            hapChromosome.setAfterVND(hapChromosome.getAfterVND() + myFitness);
            hapChromosome.setCount(hapChromosome.getCount() + 1);
        }

        /*std::cout << "Tabela MDRR: " << std::endl;
        int i = 1;
        for (auto &t:mdrr) {
            std::cout << "Games de " << i << ": " << std::endl;
            for (auto &g: t) {
                std::cout << g<<" ";
            }
            i++;
            std::cout << std::endl;
        }*/


        std::vector<int> newPermutation = polygonInverse(mdrr);
        if (newPermutation.size() != 1) {
            int j = 0;
            for (std::vector<ValueKeyPair>::const_iterator i = rank.begin(); i != rank.end(); ++i) {
                chromosome[newPermutation[j] - 1] = caux[i->second];
                j++;
            }
            /*std::cout << "2:" << std::endl;
            for (auto &g:newPermutation ) {
                std::cout << g<<" ";
            }
            std::cout << std::endl;
            std::vector<ValueKeyPair> ranki(chromosome.size());


            for (unsigned i = 0; i < chromosome.size(); ++i) {
                ranki[i] = ValueKeyPair(chromosome[i], i);
            }

            std::sort(ranki.begin(), ranki.end());
            std::cout << "3:" << std::endl;
            for (std::vector<ValueKeyPair>::const_iterator i = ranki.begin(); i != ranki.end(); ++i) {
                std::cout<<(i->second + 1)<<" ";
            }
            std::cout << std::endl;*/
        }
        return myFitness;
    } else {
        return fitness;
    }

}




