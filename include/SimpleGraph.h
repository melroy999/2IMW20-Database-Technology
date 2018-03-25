//
// Created by Nikolay Yakovets on 2018-01-31.
//

#ifndef QS_SIMPLEGRAPH_H
#define QS_SIMPLEGRAPH_H

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <iostream>
#include <regex>
#include <fstream>
#include "Graph.h"

class SimpleGraph : public Graph {
public:
    std::vector<std::vector<std::vector<uint32_t>>> adj;
    std::vector<std::vector<std::vector<uint32_t>>> reverse_adj;

    // Pointers to vectors in another vector, used to speed up certain operations.
    std::vector<std::vector<uint32_t>> *adj_ptr = nullptr;
    std::vector<std::vector<uint32_t>> *reverse_adj_ptr = nullptr;

    int L_left = -1;
    int L_right = -1;

protected:
    uint32_t V;
    std::vector<uint32_t> numEdges;
    uint32_t E;
    uint32_t L;

public:
    SimpleGraph() : V(0), E(0), L(0) {};
    ~SimpleGraph() = default;
    explicit SimpleGraph(uint32_t n);

    uint32_t getNoVertices() const override ;
    uint32_t getNoEdges() const override ;
    uint32_t getNoDistinctEdges() const override ;
    uint32_t getNoLabels() const override ;

    void addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel) override ;
    void readFromContiguousFile(const std::string &fileName) override ;

    void setNoVertices(uint32_t n);
    void setNoLabels(uint32_t noLabels);

    void setDataStructureSizes();

    void addEdges(std::shared_ptr<SimpleGraph> &in, uint32_t projectLabel, bool isInverse);
};

#endif //QS_SIMPLEGRAPH_H
