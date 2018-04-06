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

    bool operator<(const Entry& rhs) const {
        return i < rhs.i;
    }
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

    // The number of sources and targets having a degree of one.
    uint32_t noSourcesDeg1{};
    uint32_t noTargetsDeg1{};


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
        uint32_t size = sizeof(noSources) + sizeof(noTargets) + sizeof(noEdges) + sizeof(noSourcesDeg1) + sizeof(noTargetsDeg1)+ sizeof(noTargets) + sizeof(noEdges)
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
        std::cout << "\t - memory = " << getSizeInBytes() << std::endl;
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

/**
 * A join graph, which holds resulting joins/projections.
 */
struct JoinGraph {
    // The adjacency list for this graph.
    std::vector<std::vector<Entry>*> adj;

    // The sources and targets of the graph, represented as a bitmap for each individual label.
    std::vector<uint64_t> sources;
    std::vector<uint64_t> targets;

    // Pointers to an adjacency list and sources bitmap used in projections.
    std::vector<std::vector<Entry>*>* adj_ptr = nullptr;
    std::vector<uint64_t>* sources_ptr = nullptr;

    // The number of sources, targets and edges.
    uint32_t noSources{};
    uint32_t noTargets{};
    uint32_t noEdges{};

    /**
     * Constructor used to create a projection of the specified blockGraph.
     * @param g The graph to copy the structure from.
     * @param isInverse Whether we want to take the inverse.
     */
    explicit JoinGraph(BlockGraph &g, bool isInverse) {
        if(!isInverse) {
            adj_ptr = &g.adj;
            sources_ptr = &g.sources;
            noSources = g.noSources;
            noTargets = g.noTargets;
        } else {
            adj_ptr = &g.reverse_adj;
            sources_ptr = &g.targets;
            noSources = g.noTargets;
            noTargets = g.noSources;
        }

        noEdges = g.noEdges;
    }

    /**
     * Constructor used during join operations.
     * @param left The left side of the join.
     * @param right The right side of the join.
     */
    JoinGraph(JoinGraph &left, JoinGraph &right) {
        adj = std::vector<std::vector<Entry>*>(left.adj_ptr ? left.adj_ptr->size() : left.adj.size());
        sources.resize(left.sources_ptr ? left.sources_ptr->size() : left.sources.size());
        targets.resize(left.sources_ptr ? left.sources_ptr->size() : left.sources.size());

        // Instead of doing the join in the evaluator, we do it here locally, to avoid constant function calls.
        join(left, right);

        // Finalize the class.
        finalize();
    }

    /**
     * Join the left and the right graphs.
     * @param left The left side of the join, in which the sources originate.
     * @param right The right side of the join, in which the targets terminate.
     */
    void join(JoinGraph &left, JoinGraph &right) {
        // We do the join block by block. First, instantiate pointers to the actual adjacency matrices to use.
        std::vector<std::vector<Entry>*>* leftAdj = left.adj_ptr ? left.adj_ptr : &left.adj;
        std::vector<std::vector<Entry>*>* rightAdj = right.adj_ptr ? right.adj_ptr : &right.adj;

        // We want to make sure that the result of the join is sorted in block order, without duplicates.
        // By using this assumption, we know that our input is always in sorted order as well.
        std::vector<Entry> targets;

        // Start iterating over all blocks.
        for(uint32_t i = 0; i < leftAdj->size(); ++i) {
            // Do a null check, as the value might be uninitialized.
            if((*leftAdj)[i]) {

                // Find all block locations we can go to from source block i.
                for(const Entry &commonBlock : *(*leftAdj)[i]) {

                    // Where can we go to from the common block?
                    std::vector<Entry>* targetBlocks = (*rightAdj)[commonBlock.i];

                    // Count the number of non zero insertions.
                    uint32_t count = 0;

                    // Do a null check, as the value might not exist.
                    if(targetBlocks) {

                        // Create the masks for our common block.
                        uint64_t m1 = commonBlock.v & 0x0101010101010101; 
                        m1 |= m1 << 1; m1 |= m1 << 2; m1 |= m1 << 4;
                        uint64_t m2 = (commonBlock.v >> 1) & 0x0101010101010101;
                        m2 |= m2 << 1; m2 |= m2 << 2; m2 |= m2 << 4;
                        uint64_t m3 = (commonBlock.v >> 2) & 0x0101010101010101;
                        m3 |= m3 << 1; m3 |= m3 << 2; m3 |= m3 << 4;
                        uint64_t m4 = (commonBlock.v >> 3) & 0x0101010101010101;
                        m4 |= m4 << 1; m4 |= m4 << 2; m4 |= m4 << 4;
                        uint64_t m5 = (commonBlock.v >> 4) & 0x0101010101010101;
                        m5 |= m5 << 1; m5 |= m5 << 2; m5 |= m5 << 4;
                        uint64_t m6 = (commonBlock.v >> 5) & 0x0101010101010101;
                        m6 |= m6 << 1; m6 |= m6 << 2; m6 |= m6 << 4;
                        uint64_t m7 = (commonBlock.v >> 6) & 0x0101010101010101;
                        m7 |= m7 << 1; m7 |= m7 << 2; m7 |= m7 << 4;
                        uint64_t m8 = (commonBlock.v >> 7) & 0x0101010101010101;
                        m8 |= m8 << 1; m8 |= m8 << 2; m8 |= m8 << 4;
                        
                        /* Note that if bit k is set in row m of commonBlock.v, it indicates that vertex 8 * i + m
                         * has vertex 8 * commonBlock.i + k as a target. Thus, to do the join, we should join on
                         * the columns of the common block, instead of the rows.
                         */

                        // Iterate over all target blocks, and construct the new blocks.
                        for(const Entry &targetBlock : *targetBlocks) {

                            // What if, the masks indicate that we have no vertices in common?
                            if(((commonBlock.masks >> 8) & targetBlock.masks) == 0 && ((commonBlock.masks) & targetBlock.masks >> 8) == 0) {

                                // The result will be zero, so better proceed.
                                continue;
                            }

                            // Create the masks for this specific target block.
                            uint64_t t1 = targetBlock.v & 0xff;
                            t1 |= t1 << 8; t1 |= t1 << 16; t1 |= t1 << 32;
                            uint64_t t2 = (targetBlock.v >> 8) & 0xff;
                            t2 |= t2 << 8; t2 |= t2 << 16; t2 |= t2 << 32;
                            uint64_t t3 = (targetBlock.v >> 16) & 0xff;
                            t3 |= t3 << 8; t3 |= t3 << 16; t3 |= t3 << 32;
                            uint64_t t4 = (targetBlock.v >> 24) & 0xff;
                            t4 |= t4 << 8; t4 |= t4 << 16; t4 |= t4 << 32;
                            uint64_t t5 = (targetBlock.v >> 32) & 0xff;
                            t5 |= t5 << 8; t5 |= t5 << 16; t5 |= t5 << 32;
                            uint64_t t6 = (targetBlock.v >> 40) & 0xff;
                            t6 |= t6 << 8; t6 |= t6 << 16; t6 |= t6 << 32;
                            uint64_t t7 = (targetBlock.v >> 48) & 0xff;
                            t7 |= t7 << 8; t7 |= t7 << 16; t7 |= t7 << 32;
                            uint64_t t8 = (targetBlock.v >> 56) & 0xff;
                            t8 |= t8 << 8; t8 |= t8 << 16; t8 |= t8 << 32;

                            // Set the v value using the masks above.
                            uint64_t v = m1 & t1 | m2 & t2 | m3 & t3 | m4 & t4 | m5 & t5 | m6 & t6 | m7 & t7 | m8 & t8;

                            // Create the new entry.
                            if(v != 0ULL) {
                                count++;

                                // For now, the mask is zero. This will be altered when inserting the results.
                                targets.push_back({targetBlock.i, 0, v});
                            }
                        }

                        // Do an in place merge.
                        std::inplace_merge(targets.begin(), targets.end() - count, targets.end());
                    }
                }

                if(!targets.empty()) {
                    Entry* entry_ptr = &targets.front();
                    bool first = true;

                    for(auto &entry : targets) {
                        if(first) {
                            // Do nothing.
                            first = false;
                        } else if(entry_ptr->i != entry.i) {
                            // We know we are done with the current entry pointer, so set the mask.
                            // The first eight bits denote if a value in a row exist.
                            // The next eight bits denote if a value in a column exist.
                            entry_ptr->masks |= ((entry_ptr->v &  0xFF) != 0)
                                                | ((entry_ptr->v &  0xFF00) != 0) << 1
                                                | ((entry_ptr->v &  0xFF0000) != 0) << 2
                                                | ((entry_ptr->v &  0xFF000000) != 0) << 3
                                                | ((entry_ptr->v &  0xFF00000000) != 0) << 4
                                                | ((entry_ptr->v &  0xFF0000000000) != 0) << 5
                                                | ((entry_ptr->v &  0xFF000000000000) != 0) << 6
                                                | ((entry_ptr->v &  0xFF00000000000000) != 0) << 7;

                            entry_ptr->masks |= ((entry_ptr->v & 0x0101010101010101) != 0) << 8
                                                | ((entry_ptr->v & 0x0202020202020202) != 0) << 9
                                                | ((entry_ptr->v & 0x0404040404040404) != 0) << 10
                                                | ((entry_ptr->v & 0x0808080808080808) != 0) << 11
                                                | ((entry_ptr->v & 0x1010101010101010) != 0) << 12
                                                | ((entry_ptr->v & 0x2020202020202020) != 0) << 13
                                                | ((entry_ptr->v & 0x4040404040404040) != 0) << 14
                                                | ((entry_ptr->v & 0x8080808080808080) != 0) << 15;

                            // The i value has changed, so change the pointer.
                            entry_ptr = &entry;
                        } else {
                            // We found another occurrence of the same i value, so merge the flags.
                            entry_ptr->v |= entry.v;
                        }
                    }

                    // Set the mask for the last encountered entry.
                    entry_ptr->masks |= ((entry_ptr->v &  0xFF) != 0)
                                        | ((entry_ptr->v &  0xFF00) != 0) << 1
                                        | ((entry_ptr->v &  0xFF0000) != 0) << 2
                                        | ((entry_ptr->v &  0xFF000000) != 0) << 3
                                        | ((entry_ptr->v &  0xFF00000000) != 0) << 4
                                        | ((entry_ptr->v &  0xFF0000000000) != 0) << 5
                                        | ((entry_ptr->v &  0xFF000000000000) != 0) << 6
                                        | ((entry_ptr->v &  0xFF00000000000000) != 0) << 7;

                    entry_ptr->masks |= ((entry_ptr->v & 0x0101010101010101) != 0) << 8
                                        | ((entry_ptr->v & 0x0202020202020202) != 0) << 9
                                        | ((entry_ptr->v & 0x0404040404040404) != 0) << 10
                                        | ((entry_ptr->v & 0x0808080808080808) != 0) << 11
                                        | ((entry_ptr->v & 0x1010101010101010) != 0) << 12
                                        | ((entry_ptr->v & 0x2020202020202020) != 0) << 13
                                        | ((entry_ptr->v & 0x4040404040404040) != 0) << 14
                                        | ((entry_ptr->v & 0x8080808080808080) != 0) << 15;

                    // Insert all the edges.
                    auto it = std::unique(targets.begin(), targets.end(), [](const Entry& a, const Entry& n) {return a.i == n.i;});
                    adj[i] = new std::vector<Entry>(static_cast<unsigned long>(std::distance(targets.begin(), it)));
                    std::copy(targets.begin(), it, adj[i]->begin());

                    // Set the number of edges plus the source and target flags.
                    for(const auto &entry : *adj[i]) {
                        noEdges += __builtin_popcountll(entry.v);
                        JoinGraph::sources[i >> 3] |= (uint64_t(entry.masks) & 255) << (8 * (i & 7));
                        JoinGraph::targets[entry.i >> 3] |= (uint64_t(entry.masks >> 8) & 255) << (8 * (entry.i & 7));
                    }

                    // Clear the cache of targets.
                    targets.clear();
                }
            }
        }
    }

    uint32_t getSizeInBytes() {
        // The last part of this equation is the minimal size of the lists.
        uint32_t size = sizeof(noSources) + sizeof(noTargets) + sizeof(noEdges) + sizeof(adj_ptr) + sizeof(sources_ptr)
                        + 2 * (sizeof(std::vector<uint64_t>) + sizeof(uint64_t) * sources.size())
                        + sizeof(std::vector<std::vector<Entry>*>) + sizeof(std::vector<Entry>*) * adj.size();

        for(const auto &v : adj) {
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

    virtual ~JoinGraph() {
        for(const auto &v : adj) {
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
