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

#define SET_BIT(pos) (1ULL << (pos))

static uint32_t countBitsSet(std::vector<uint64_t> &result) {
    uint32_t sum = 0;
    for(auto v : result) {
        if(v == 0) continue;
        sum += __builtin_popcountll(v);
    }
    return sum;
}

static std::vector<uint64_t> doAnd(const std::vector<uint64_t> &t, const std::vector<uint64_t> &s) {

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

struct Entry {
    uint32_t i;
    uint32_t masks;
    uint64_t v;
};

struct BlockGraph {
    // The adjacency list and reverse adjacency list for this graph.
    std::vector<std::vector<Entry>*> adj;
    std::vector<std::vector<Entry>*> reverse_adj;

    // The sources and targets of the graph, represented as a bitmap for each individual label.
    std::vector<uint64_t> sources;
    std::vector<uint64_t> targets;

    // The number of sources, targets and edges.
    uint32_t noSources{};
    uint32_t noTargets{};
    uint32_t noEdges{};

    void setDataStructureSizes(uint32_t V) {
        adj = std::vector<std::vector<Entry>*>((V + 7) >> 3);
        reverse_adj = std::vector<std::vector<Entry>*>((V + 7) >> 3);
        sources.resize((V + 63) >> 6);
        targets.resize((V + 63) >> 6);
    }

    void addEdge(uint32_t i, uint32_t j, Entry x, Entry y) {
        if(x.v != 0ULL) {
            if(!adj[i]) {
                adj[i] = new std::vector<Entry>();
            }
            adj[i]->emplace_back(x);

            if(!reverse_adj[j]) {
                reverse_adj[j] = new std::vector<Entry>();
            }
            reverse_adj[j]->emplace_back(y);

            // Enable the bits of the sources and targets.
            // Note that i >> 3 == (i << 3) >> 6.
            // Isolate the first 8 bits, and shift to the appropriate spot.
            sources[i >> 3] |= (uint64_t(x.masks) & 255) << (8 * (i & 7));
            targets[j >> 3] |= (uint64_t(y.masks) & 255) << (8 * (j & 7));

            // Increment the number of edges.
            noEdges += __builtin_popcountll(x.v);
        }
    }

    uint32_t getSizeInBytes() {
        // The last part of this equation is the minimal size of the lists.
        uint32_t size = sizeof(noSources) + sizeof(noTargets) + sizeof(noEdges)
                        + 2 * (sizeof(std::vector<uint64_t>) + sizeof(uint64_t) * sources.size())
                        + 2 * (sizeof(std::vector<std::vector<Entry>*>) + sizeof(std::vector<Entry>*) * adj.size());

        for(const auto &v : adj) {
            if(v) {
                size += sizeof(std::vector<Entry>) + sizeof(Entry) * v->size();
            }
        }

        for(const auto &v : reverse_adj) {
            if(v) {
                size += sizeof(std::vector<Entry>) + sizeof(Entry) * v->size();
            }
        }

        return size;
    }

    void finalize() {
        noSources = countBitsSet(sources);
        noTargets = countBitsSet(targets);
    }

    void reportData(uint32_t id) {
        std::cout << "Graph for label " << id << ":" << std::endl;
        std::cout << "\t - |S| = " << noSources << std::endl;
        std::cout << "\t - |T| = " << noTargets << std::endl;
        std::cout << "\t - |E| = " << noEdges << std::endl;
    }

    virtual ~BlockGraph() {
        for(const auto &v : adj) {
            delete(v);
        }

        for(const auto &v : reverse_adj) {
            delete(v);
        }
    }
};

class SimpleGraph : public Graph {
public:
    // The individual graphs for all of the labels.
    std::vector<BlockGraph> graphs;

protected:
    uint32_t V;
    uint32_t E;
    uint32_t L;

public:
    SimpleGraph() : V(0), E(0), L(0) {};
    ~SimpleGraph() = default;
    explicit SimpleGraph(uint32_t n);

    uint32_t getNoVertices() const override ;
    uint32_t getNoEdges() const override ;
    uint32_t getNoDistinctEdges() const override { return getNoEdges(); };
    uint32_t getNoLabels() const override ;

    void addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel) override {
        std::runtime_error(std::string("This function is unused in our implementation."));
    };

    void readFromContiguousFile(const std::string &fileName) override ;

    void setNoVertices(uint32_t n);
    void setNoLabels(uint32_t noLabels);
};

#endif //QS_SIMPLEGRAPH_H
