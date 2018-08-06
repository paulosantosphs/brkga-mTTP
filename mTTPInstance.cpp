//
// Created by paulo on 17/02/18.
//

#include <iostream>
#include "mTTPInstance.h"

mTTPInstance::mTTPInstance(const std::string &instanceFile) {

    std::vector<std::vector<int> > distances = readFile(instanceFile);
    setDistances(distances);
}

std::vector<std::vector<int>> mTTPInstance::readFile(const std::string &instanceFile) {
    std::ifstream source;
    source.open(instanceFile);
    std::string line;

    std::vector<std::vector<int> > distances;

    while (std::getline(source, line)) {
        //std::cout<<line<<std::endl;
        std::istringstream buf(line);
        std::istream_iterator<std::string> beg(buf), end;

        std::vector<std::string> column(beg, end);
        std::vector<int> col;
        for (auto &c: column) {
            int dis = atoi(c.c_str());
            col.push_back(dis);
        }

        distances.push_back(col);
    }

    return distances;
}

void mTTPInstance::setDistances(std::vector<std::vector<int>> distances) {
    this->distances = distances;
}

std::vector<std::vector<int> > mTTPInstance::getDistances() {
    return distances;
}
