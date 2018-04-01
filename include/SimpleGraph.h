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

static inline uint32_t countBitsSet(std::vector<uint64_t> &result) {
    uint32_t sum = 0;
    for(const auto &v : result) {
        sum += __builtin_popcountll(v);
    }
    return sum;
}

static inline std::vector<uint64_t> doAnd(const std::vector<uint64_t> &t, const std::vector<uint64_t> &s) {

    // If the sizes do not correspond, return an empty join.
    std::vector<uint64_t> result(t.size());
    if(t.size() != s.size()) {
        return result;
    }

    for(unsigned long i = t.size() ; i -- > 0 ; ) {
        result[i] = t[i] & s[i];
    }
    return result;
}

class SimpleGraphStorage {
public:
    SimpleGraphStorage(uint32_t n, uint32_t N);
    ~SimpleGraphStorage();

    // The adjacency list. The reverse adjacency list is unused in joins.
    std::vector<std::vector<uint32_t>*> adj;
    std::vector<std::vector<uint32_t>*> reverse_adj;

    // The sources and targets of the graph, represented as a bitmap for each individual label.
    std::vector<uint64_t> sources;
    std::vector<uint64_t> targets;

    // The number of source and target vertices, plus the number of edges.
    uint32_t noEdges{};
    uint32_t noSources{};
    uint32_t noTargets{};

    // Store the minimum and maximum source/target numbers encountered.
    uint32_t minSource;
    uint32_t maxSource;
    uint32_t minTarget;
    uint32_t maxTarget;

    // The sizes of the storage.
    uint32_t n;
    uint32_t N;

    void setNoVertices(uint32_t n);
    void setNoBuckets(uint32_t n);
    void finalize();

    void addEdge(uint32_t from, uint32_t to);
};

class SimpleJoinStorage {
public:
    SimpleJoinStorage(uint32_t n, uint32_t N);
    SimpleJoinStorage(SimpleGraphStorage &storage, bool isInverse);

    // The adjacency list.
    std::vector<std::vector<uint32_t>*> adj;

    // The sources and targets of the graph, represented as a bitmap for each individual label.
    std::vector<uint64_t> sources;
    std::vector<uint64_t> targets;

    // A pointer to the adjacency list, used when the storage is a leaf.
    std::vector<std::vector<uint32_t>*> *adj_ptr{};

    // Pointers to the source and targets vectors, used when the storage is a leaf.
    std::vector<uint64_t> *sources_ptr{};

    // The number of source and target vertices, plus the number of edges.
    uint32_t noEdges{};
    uint32_t noSources{};
    uint32_t noTargets{};

    // Store the minimum and maximum source/target numbers encountered.
    uint32_t minSource;
    uint32_t maxSource;
    uint32_t minTarget;
    uint32_t maxTarget;

    // The sizes of the storage.
    uint32_t n;
    uint32_t N;

    /**
     * Add a collection of edges in one sweep.
     *
     * @param from The node the edges originate from.
     * @param data The target nodes.
     * @param end An iterator pointing to the last object that should be added.
     */
    void addEdges(uint32_t from, std::vector<uint32_t> &data, std::vector<uint32_t>::iterator &end);

    void setNoVertices(uint32_t n);
    void setNoBuckets(uint32_t N);
    void finalize();
};

class SimpleGraph : public Graph {
public:
    // The simple graph storage areas for each of the labels.
    std::vector<SimpleGraphStorage> graphs;

protected:
    uint32_t V;
    uint32_t E;
    uint32_t L;

public:
    SimpleGraph() : V(0), E(0), L(0) {};
    explicit SimpleGraph(uint32_t n);

    uint32_t getNoVertices() const override ;
    uint32_t getNoEdges() const override ;
    uint32_t getNoDistinctEdges() const override ;
    uint32_t getNoLabels() const override ;

    void readFromContiguousFile(const std::string &fileName) override ;

    void addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel) override;

    void setNoVertices(uint32_t n);
    void setNoLabels(uint32_t noLabels);

    void setDataStructureSizes();
};

#endif //QS_SIMPLEGRAPH_H
