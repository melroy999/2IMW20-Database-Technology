//
// Created by Nikolay Yakovets on 2018-01-31.
//

#include "SimpleGraph.h"

SimpleGraph::SimpleGraph(uint32_t n)   {
    setNoVertices(n);
}

uint32_t SimpleGraph::getNoVertices() const {
    return V;
}

void SimpleGraph::setDataStructureSizes() {
    adj.resize(getNoLabels());
    reverse_adj.resize(getNoLabels());

    for (auto &l : adj)
        l.resize(V);

    for (auto &l : reverse_adj)
        l.resize(V);
}

void SimpleGraph::setNoVertices(uint32_t n) {
    V = n;
}

uint32_t SimpleGraph::getNoEdges() const {
    uint32_t sum = 0;
    for (const auto &l : adj)
        for(const auto &v : l)
            sum += v->size();
    return sum;
}

uint32_t SimpleGraph::getNoDistinctEdges() const {

    uint32_t sum = 0;

    if(adj_ptr) {
        for (auto sourceVec : *adj_ptr) {
            if(sourceVec) {
                std::sort(sourceVec->begin(), sourceVec->end());

                uint32_t prevTarget = 0;
                bool first = true;

                for (const auto &target : *sourceVec) {
                    if (first || prevTarget != target) {
                        first = false;
                        sum++;
                        prevTarget = target;
                    }
                }
            }
        }
    }

    for (const auto &labelVec : adj) {
        for (auto sourceVec : labelVec) {
            if(sourceVec) {
                std::sort(sourceVec->begin(), sourceVec->end());

                uint32_t prevTarget = 0;
                bool first = true;

                for (const auto &target : *sourceVec) {
                    if (first || prevTarget != target) {
                        first = false;
                        sum++;
                        prevTarget = target;
                    }
                }
            }
        }
    }

    return sum;
}

uint32_t SimpleGraph::getNoLabels() const {
    return L;
}

void SimpleGraph::setNoLabels(uint32_t noLabels) {
    L = noLabels;
}

void SimpleGraph::addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel) {
    if(from >= V || to >= V || edgeLabel >= L)
        throw std::runtime_error(std::string("Edge data out of bounds: ") +
                                         "(" + std::to_string(from) + "," + std::to_string(to) + "," +
                                         std::to_string(edgeLabel) + ")");

    if(!adj[edgeLabel][from]) adj[edgeLabel][from] = new std::vector<uint32_t>();
    if(!reverse_adj[edgeLabel][to]) reverse_adj[edgeLabel][to] = new std::vector<uint32_t>();

    adj[edgeLabel][from]->emplace_back(to);
    reverse_adj[edgeLabel][to]->emplace_back(from);
}

void SimpleGraph::addEdges(std::shared_ptr<SimpleGraph> &in, uint32_t projectLabel, bool isInverse) {

    auto _adj = &in->adj[projectLabel];
    auto _reverse_adj = &in->reverse_adj[projectLabel];

    if(!isInverse) {
        adj_ptr = _adj;
        reverse_adj_ptr = _reverse_adj;
    } else {
        adj_ptr = _reverse_adj;
        reverse_adj_ptr = _adj;
    }
}

void SimpleGraph::readFromContiguousFile(const std::string &fileName) {

    std::fstream file(fileName, std::ios_base::in);

    std::string graphStructure;
    file >> graphStructure;

    std::istringstream iss(graphStructure);
    std::string line;
    uint32_t noNodes = 0;
    uint32_t noEdges = 0;
    uint32_t noLabels = 0;
    int j = 0;
    while(std::getline(iss, line, ',')) {
        if(j == 0) {
            noNodes = (uint32_t) std::stoul(line);
        }  else if(j == 1) {
            noEdges = (uint32_t) std::stoul(line);
        } else if(j == 2) {
            noLabels = (uint32_t) std::stoul(line);
        }
        j++;
    }

    setNoLabels(noLabels);
    setNoVertices(noNodes);
    setDataStructureSizes();

    for(int i = 0; i < noEdges; i++) {
        uint32_t subject, predicate, object;

        file >> subject >> predicate >> object;
        file.ignore(INT32_MAX, '\n');

        addEdge(subject, object, predicate);
    }

    file.close();
}