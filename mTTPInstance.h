//
// Created by paulo on 17/02/18.
//

#ifndef BRKGA_MTTP_MTTPINSTANCE_H
#define BRKGA_MTTP_MTTPINSTANCE_H

#include <cmath>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <iterator>


class mTTPInstance {
public:

    mTTPInstance(const std::string &instanceFile);

    std::vector<std::vector<int>> getDistances();

private:
    std::vector<std::vector<int>> distances;

    std::vector<std::vector<int>> readFile(const std::string &filename);

    void setDistances(std::vector<std::vector<int>> distances);
};


#endif //BRKGA_MTTP_MTTPINSTANCE_H
