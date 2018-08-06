#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <cstring>
#include <time.h>
#include "brkgaAPI/BRKGA.h"
#include "brkgaAPI/MTRand.h"

#include "mTTPInstance.h"
#include "mTTPDecoder.h"


using namespace std::chrono;

double standardDeviation(std::vector<double> results, double &mean);

double gap(int best, int myBest);

double mean(std::vector<double> elements);

int main() {
    std::vector<std::string> instancesNames;
    instancesNames.push_back("nl4.txt");
    instancesNames.push_back("nl6.txt");
    instancesNames.push_back("nl8.txt");
    instancesNames.push_back("nl10.txt");
    instancesNames.push_back("nl12.txt");
    instancesNames.push_back("nl14.txt");
    instancesNames.push_back("nl16.txt");
    instancesNames.push_back("circ4.txt");
    instancesNames.push_back("circ6.txt");
    instancesNames.push_back("circ8.txt");
    instancesNames.push_back("circ10.txt");
    instancesNames.push_back("circ12.txt");
    instancesNames.push_back("circ14.txt");
    instancesNames.push_back("circ16.txt");
    instancesNames.push_back("circ18.txt");
    instancesNames.push_back("circ20.txt");
    std::string instanceFile;

    std::vector<int> bests;
    //NL
    bests.push_back(8276);
    bests.push_back(26588);
    bests.push_back(41928);
    bests.push_back(63832);
    bests.push_back(119608);
    bests.push_back(199363);
    bests.push_back(278305);

    //Circ
    bests.push_back(20);
    bests.push_back(72);
    bests.push_back(140);
    bests.push_back(272);
    bests.push_back(432);
    bests.push_back(672);
    bests.push_back(968);
    bests.push_back(1306);
    bests.push_back(1852);

    int iterations = 5;

    for (int k = 0; k < instancesNames.size(); k++) {

        std::vector<double> solutions;
        std::vector<double> times;
        std::vector<double> convs;
        std::vector<double> rhoevec;
        std::vector<double> pevec;
        std::vector<double> pmvec;

        mTTPChromosome bestSolution;
        bestSolution.setFitness(999999999999);

        for (int y = 0; y < iterations; y++) {
            instanceFile = "../instances/" + instancesNames[k];
            std::cout << "Instance file: " << instanceFile << std::endl;

            mTTPInstance instance(instanceFile);
            mTTPChromosome hapChromosome;
            hapChromosome.setFitness(999999999999);

            std::cout << "Numero de times: " << instance.getDistances().size() << std::endl;

            const unsigned n = instance.getDistances().size();// size of chromosomes
            const unsigned p = n * 5;// size of population

            double pe = 0.15;        // fraction of population to be the elite-set
            double pm = 0.15;        // fraction of population to be replaced by mutants
            double rhoe = 0.95;// probability that offspring inherit an allele from elite parent


            const unsigned K = 4;        // number of independent populations
            const unsigned MAXT = 1;    // number of threads for parallel decoding

            unsigned lsPattern = 1;

            const long unsigned rngSeed = time(NULL);    // seed to the random number generator
            MTRand rng(rngSeed);                // initialize the random number generator

            mTTPDecoder decoder(instance.getDistances(), hapChromosome, rng,
                                lsPattern);                // initialize the decoder
            // initialize the BRKGA-based heuristic
            BRKGA<mTTPDecoder, MTRand> algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);

            unsigned generation = 0;        // current generation
            const unsigned X_INTVL = 5;    // exchange best individuals at every 100 generations
            const unsigned X_NUMBER = 2;    // exchange top 2 best
            const unsigned MAX_GENS = 50;    // run for 1000 gens

            int convergence = 0;
            int count = 0;
            double beforeFit = 9999999999999;
            float adaptationFactor = 0.03;

            double runningTime;

            high_resolution_clock::time_point t1 = high_resolution_clock::now();
            std::cout << "Running for " << MAX_GENS << " generations..." << std::endl;
            do {
                algorithm.evolve();    // evolve the population for one generation
                if (algorithm.getBestFitness() == hapChromosome.getFitness()) {
                    if (beforeFit != algorithm.getBestFitness()) {
                        convergence = generation;
                        beforeFit = algorithm.getBestFitness();
                    }
                }

                //Adaptação automática de parametros
                /*if(convergence !=0){
                    if(((generation - convergence) > (MAX_GENS*adaptationFactor)) && count > MAX_GENS*adaptationFactor){
                        count = 0;
                        if(rhoe>0){
                            rhoe = rhoe - rhoe*adaptationFactor;
                            algorithm.setRhoe(rhoe);
                        }
                        if(pe*p<=0){
                            std::cout<<"Elite 0, geração: "<<generation<<std::endl;
                        }

                        if((pm + pm*adaptationFactor)*p+(algorithm.getPe()) <= p){
                            pm = pm + pm*adaptationFactor;
                            algorithm.setPm(pm);
                        }
                        if((pe + pe*adaptationFactor)*p+(algorithm.getPm()) <= p){
                            pe = pe - pe*adaptationFactor;
                            algorithm.setPe(pe);
                        }
                    }
                }
                count++;*/

                /*if (generation == 10 || generation == 20 || generation == 30 || generation == 40 ||
                    generation == 500 || generation == 600) {
                    std::cout << "Fitness of the top individuals of each population. Generation: " << generation
                              << std::endl;
                    const unsigned bound = std::min(p, unsigned(p * pe));    // makes sure we have 10 individuals
                    for (unsigned i = 0; i < K; ++i) {
                        std::cout << "Population #" << i << ":" << std::endl;
                        for (unsigned j = 0; j < bound; ++j) {
                            std::cout << "\t" << j << ") "
                                      << algorithm.getPopulation(i).getFitness(j) << std::endl;
                        }
                    }

                }*/

                if ((++generation) % X_INTVL == 0) {
                    algorithm.exchangeElite(X_NUMBER);    // exchange top individuals
                }
            } while (generation < MAX_GENS);


            high_resolution_clock::time_point t2 = high_resolution_clock::now();                //time taking
            duration<double> time_span = duration_cast<duration<double> >(t2 - t1);
            runningTime = time_span.count();

            // print the fitness of the top 10 individuals of each population:
            std::cout << "Fitness of the top 10 individuals of each population:" << std::endl;
            const unsigned bound = std::min(p, unsigned(p));    // makes sure we have 10 individuals
            for (unsigned i = 0; i < K; ++i) {
                std::cout << "Population #" << i << ":" << std::endl;
                for (unsigned j = 0; j < bound; ++j) {
                    std::cout << "\t" << j << ") "
                              << algorithm.getPopulation(i).getFitness(j) << std::endl;
                }
            }

            std::cout << "Best solution found has objective value = "
                      << algorithm.getBestFitness() << std::endl;

            std::vector<std::vector<int>> mdrr = hapChromosome.getHaps();
            int i = 1;
            for (auto &t: mdrr) {
                std::cout << "Games de " << i << ": " << std::endl;
                for (auto &g: t) {
                    std::cout << g << " ";
                }
                i++;
                std::cout << std::endl;
            }
            std::cout << decoder.totalDistanceTravelled(mdrr,instance.getDistances()) << std::endl;
            solutions.push_back(hapChromosome.getFitness());

            if (hapChromosome.getFitness() < bestSolution.getFitness()) {
                bestSolution.setFitness(hapChromosome.getFitness());
                bestSolution.setChromosome(hapChromosome.getChromosome());
                bestSolution.setHaps(hapChromosome.getHaps());
                bestSolution.setPermutation(hapChromosome.getPermutation());
            }
            bestSolution.setBeforeVND(hapChromosome.getBeforeVND() + bestSolution.getBeforeVND());
            bestSolution.setAfterVND(hapChromosome.getAfterVND() + bestSolution.getAfterVND());
            bestSolution.setCount(hapChromosome.getCount() + bestSolution.getCount());

            rhoevec.push_back(algorithm.getRhoe());
            pmvec.push_back((double)algorithm.getPm()/(n*5));
            pevec.push_back((double)algorithm.getPe()/(n*5));

            std::cout << "Convergencia:" << convergence << std::endl;
            convs.push_back(convergence);
            std::cout << "Tempo :" << runningTime << std::endl;
            times.push_back(runningTime);
        }

        double meanSolutions;
        double sd = standardDeviation(solutions, meanSolutions);
        double meanTime = mean(times);
        double meanConv = mean(convs);
        double gapValue = gap(bests[k], bestSolution.getFitness());
        double meanBeforeVND = bestSolution.getBeforeVND()/bestSolution.getCount();
        double meanAfterVND = bestSolution.getAfterVND()/bestSolution.getCount();
        double sumRhoe = 0;
        double sumPe = 0;
        double sumPm = 0;
        double meanRhoe;
        double meanPe;
        double meanPm;

        for(auto &v: rhoevec){
            sumRhoe = sumRhoe + v;
        }

        for(auto &v: pevec){
            sumPe = sumPe + v;
        }

        for(auto &v: pmvec){
            sumPm = sumPm + v;
        }

        meanRhoe = sumRhoe/rhoevec.size();
        meanPe = sumPe/pevec.size();
        meanPm = sumPm/pmvec.size();



        std::vector<std::vector<int>> mdrr = bestSolution.getHaps();

        FILE *fpSolution;

        char outputFileName[256];

        strncpy(outputFileName, instancesNames[k].c_str(), sizeof(outputFileName));
        outputFileName[sizeof(outputFileName) - 1] = 0;

        fpSolution = fopen(outputFileName,
                           "w");                //file that contains the information about the solution of a problem instance

        fprintf(fpSolution,
                "Best Total Distance Travelled : %f\nGAP: %f\nMean Solution: %f\nStandardDeviation: %f\nMean Running time: %lf seconds\nMean Convergence: %f\nMean Before VND: %f\nMean After VND: %f\nMean Rhoe: %f\nMean Pe: %f\nMean Pm: %f\n",
                bestSolution.getFitness(), gapValue, meanSolutions, sd, meanTime, meanConv, meanBeforeVND, meanAfterVND, meanRhoe, meanPe, meanPm);

        fprintf(fpSolution, "\n\n");
        fprintf(fpSolution, "Table: \n");

        fprintf(fpSolution, "Teams:    ");
        for (int j = 0; j < 2 * (mdrr.size() - 1); j++) {
            if (j < 9) {
                fprintf(fpSolution, "R%d:  ", j + 1);
            } else {
                fprintf(fpSolution, "R%d: ", j + 1);
            }

        }
        fprintf(fpSolution, "\n");
        int i = 1;
        for (auto &t: mdrr) {
            fprintf(fpSolution, "Team %d:", i);
            for (auto &g: t) {
                if (g < 0 && g > 9) {
                    fprintf(fpSolution, " %d", g);
                } else if (g < 0 || g > 9) {
                    fprintf(fpSolution, "   %d", g);
                } else {
                    fprintf(fpSolution, "    %d", g);
                }

            }
            i++;
            fprintf(fpSolution, "\n");
        }

        fclose(fpSolution);
    }
    return 0;
}

double standardDeviation(std::vector<double> results, double &mean) {
    float sum = 0, standardDeviation = 0;

    for (int i = 0; i < results.size(); ++i) {
        sum += results[i];
    }

    mean = sum / (results.size());

    for (int i = 0; i < results.size(); ++i)
        standardDeviation += pow(results[i] - mean, 2);

    return sqrt(standardDeviation / results.size());
}

double gap(int best, int myBest) {
    double gap = myBest - best;
    gap = gap / best;
    return (gap * 100);
}

double mean(std::vector<double> elements) {
    double sum = 0;
    for (int i = 0; i < elements.size(); i++) {
        sum = sum + elements[i];
    }
    return (sum / elements.size());
}
