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

#define CHECK_BIT(var,pos) ((var) & (1ULL<<(pos)))
#define SET_BIT(pos) (1ULL << ((pos) % 64))

static uint32_t countBitsSet(std::vector<uint64_t>* result) {
    uint32_t sum = 0;
    for(auto v : *result) {
        if(v == 0) continue;
        sum += __builtin_popcountll(v);
    }
    return sum;
}

static std::vector<uint64_t> doAnd(const std::vector<uint64_t> *t, const std::vector<uint64_t> *s) {

    // If the sizes do not correspond, return an empty join.
    std::vector<uint64_t> result(t->size());
    if(t->size() != s->size()) {
        return result;
    }

    for(unsigned long i = t->size() ; i -- > 0 ; ) {
        result[i] = (*t)[i] & (*s)[i];
    }
    return result;
}

class SimpleGraph : public Graph {
public:
    // The adjacency list. The reverse adjacency list is unused in joins.
    std::vector<std::vector<std::vector<uint32_t>>> adj;
    std::vector<std::vector<std::vector<uint32_t>>> reverse_adj;

    // The sources and targets of the graph, represented as a bitmap for each individual label.
    std::vector<std::vector<uint64_t>> sources;
    std::vector<std::vector<uint64_t>> targets;

    // Pointers to vectors in another vector, used to speed up certain operations.
    std::vector<std::vector<uint32_t>> *adj_ptr = nullptr;
    std::vector<uint64_t> *sources_ptr = nullptr;
    std::vector<uint64_t> *targets_ptr = nullptr;

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

    void setDataStructureSizes(bool isJoin);

    void addEdges(std::shared_ptr<SimpleGraph> &in, uint32_t projectLabel, bool isInverse);
};

#endif //QS_SIMPLEGRAPH_H
