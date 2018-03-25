//
// Created by Nikolay Yakovets on 2018-01-31.
//

#ifndef QS_SIMPLEGRAPH_H
#define QS_SIMPLEGRAPH_H

#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <memory>
#include "Graph.h"

#define CHECK_BIT(var,pos) ((var) & (1ULL<<(pos)))
#define SET_BIT(pos) (1ULL << ((pos) % 64))
#define CLEAR_BIT(pos) (~(1ULL << ((pos) % 64)))

static uint32_t popcountVector(std::vector<uint64_t>* result) {
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

static void doAndIP(std::vector<uint64_t> *s, const std::vector<uint64_t> *t) {

    for(unsigned long i = t->size() ; i -- > 0 ; ) {
        (*s)[i] &= (*t)[i];
    }
}

static void doOrIP(std::vector<uint64_t> *s, const std::vector<uint64_t> *t) {

    for(unsigned long i = t->size() ; i -- > 0 ; ) {
        (*s)[i] |= (*t)[i];
    }
}

// Represents nodes within the K2Tree.
template <uint32_t k>
struct Node {

    // We store bits as bytes, since it is the minimum size that is allocatable.
    uint64_t v{};

    // First of all, we only want a pointer to a child to exist if the corresponding bit n the array is 1.
    // Otherwise, we will waste a lot of space unintentionally on objects or pointers, of which the size >> bits.
    std::vector<std::pair<uint8_t, Node<k>*>>* children{};

    /**
     * Add an edge between the given source and the target.
     *
     * @param s The source of the edge.
     * @param t The target of the edge.
     * @param n The current width of the sub block.
     */
    void addEdge(uint32_t s, uint32_t t, uint32_t n) {
        // If we are at a leaf, we should just set the appropriate bit.
        if(n == k) {
            v |= SET_BIT(s * k + t);
        } else {
            // Find the appropriate quadrant, or expand it if it does not exist.
            uint32_t s_i = s / (n / k);
            uint32_t t_i = t / (n / k);

            // We want to mode to block s_i, t_i, which is at position s_i * k + t_i.
            auto i = static_cast<uint8_t>(s_i * k + t_i);
            if(!CHECK_BIT(v, i)) {
                // There does not exist a node yet for this subtree. Create it.
                v |= SET_BIT(i);
                std::pair<uint8_t, Node<k>*> pair = {i, new Node<k>()};

                if(!children) {
                    children = new std::vector<std::pair<uint8_t, Node<k>*>>();
                }
                children->push_back(pair);
            }

            // Find the note we want to add the edge to.
            // Note: keep the typename here, unless you want to cause the apocalypse.
            typename std::vector<std::pair<uint8_t, Node<k>*>>::iterator it = std::find_if(
                    children->begin(), children->end(), [i](std::pair<uint8_t, Node<k>*> const& obj){return obj.first == i;});

            // We found our target. Call recursively with the reduced s and t.
            it->second->addEdge(s % (n / k), t % (n / k), n / k);
        }
    }

    void getDirect(uint32_t s, uint32_t x, uint32_t n, std::vector<uint32_t> &results) {
        // If we are in a leaf, report the found vertices.
        if(n == k) {
            for(uint32_t i = 0; i < k; i++) {
                if(CHECK_BIT(v, k * (s % k) + i)) {
                    results.push_back(x + i);
                }
            }
        } else {
            if(children) {
                // Otherwise, keep exploring. Which blocks contain the row that we want?
                uint32_t s_i_b = s / (n / k);

                for(const std::pair<uint8_t, Node<k>*> &pair : *children) {
                    if(pair.first / k == s_i_b) {
                        pair.second->getDirect(s % (n / k), x + (pair.first % k) * (n / k), n / k, results);
                    }
                }
            }
        }
    }

    void getReverse(uint32_t t, uint32_t x, uint32_t n, std::vector<uint32_t> &results) {
        // If we are in a leaf, report the found vertices.
        if(n == k) {
            for(uint32_t i = 0; i < k; i++) {
                if(CHECK_BIT(v, t % k + i * k)) {
                    auto sum = x + i;
                    results.push_back(x + i);
                }
            }
        } else {
            if(children) {
                // Otherwise, keep exploring. Which blocks contain the column that we want?
                uint32_t t_i_b = t / (n / k);

                for(const std::pair<uint8_t, Node<k>*> &pair : *children) {
                    if(pair.first % k == t_i_b) {
                        pair.second->getReverse(t % (n / k), x + (pair.first / k) * (n / k), n / k, results);
                    }
                }
            }
        }
    }

    void getNoEdges(uint32_t& currentSize) {
        if(!children) {
            currentSize += __builtin_popcountll(v);
        } else {
            for(const std::pair<uint8_t, Node<k>*> &pair : *children) {
                pair.second->getNoEdges(currentSize);
            }
        }
    }

    ~Node() {
        if(children) {
            for(const std::pair<uint8_t, Node<k>*> &pair : *children) {
                delete pair.second;
            }
            delete children;
        }
    }

};

// Object that contains statistical data about the K2Tree, including an access point to the main node.
class K2Tree {
    // The root of the K2Tree.
    Node<8> root;

    // The logical size of the tree, which is a power of 2.
    uint32_t N;

    // The size of the tree, which corresponds to the number of vertices in the tree.
    uint32_t n;

    // The maximum depth of the tree.
    uint32_t h;

    // The number of unique source and target vertices.
    std::vector<uint64_t> sources;
    std::vector<uint64_t> targets;

    // The number of sources and targets.
    uint32_t noSources{};
    uint32_t noTargets{};

    // The sources and targets with degree one.
    std::vector<uint64_t> sourcesDeg1;
    std::vector<uint64_t> targetsDeg1;

    // The number of sources and targets with degree one.
    uint32_t noSourcesDeg1{};
    uint32_t noTargetsDeg1{};

    // The size of the tree, which we can fetch after construction.
    uint32_t size{};

public:
    explicit K2Tree(uint32_t n) : n(n) {
        h = static_cast<uint32_t>(std::ceil(std::log(n) / std::log(8)));
        N = static_cast<uint32_t>(std::pow(8, h));
        root = Node<8>();
        sources = std::vector<uint64_t>(static_cast<unsigned long>(std::ceil((double) n / 64)));
        targets = std::vector<uint64_t>(static_cast<unsigned long>(std::ceil((double) n / 64)));
        sourcesDeg1 = std::vector<uint64_t>(static_cast<unsigned long>(std::ceil((double) n / 64)), std::numeric_limits<uint64_t>::max());
        targetsDeg1 = std::vector<uint64_t>(static_cast<unsigned long>(std::ceil((double) n / 64)), std::numeric_limits<uint64_t>::max());
    }

    /**
     * Add an edge between the given source and the target.
     *
     * @param s The source of the edge.
     * @param t The target of the edge.
     */
    void addEdge(uint32_t s, uint32_t t) {
        root.addEdge(s, t, N);

        if(CHECK_BIT(sources[s / 64], s % 64)) {
            sourcesDeg1[s / 64] &= CLEAR_BIT(s);
        } else {
            sources[s / 64] |= SET_BIT(s);
        }

        if(CHECK_BIT(targets[t / 64],t % 64)) {
            targetsDeg1[t / 64] &= CLEAR_BIT(t);
        } else {
            targets[t / 64] |= SET_BIT(t);
        }
    }

    void finalizeTree() {
        root.getNoEdges(size);
        noSources = popcountVector(&sources);
        noTargets = popcountVector(&targets);
        doAndIP(&sourcesDeg1, &sources);
        doAndIP(&targetsDeg1, &targets);
        noSourcesDeg1 = popcountVector(&sourcesDeg1);
        noTargetsDeg1 = popcountVector(&targetsDeg1);
    }

    uint32_t getNoEdges() const {
        return size;
    }

    uint32_t getNoSources() const {
        return noSources;
    }

    uint32_t getNoTargets() const {
        return noTargets;
    }

    uint32_t getNoSourcesDeg1() const {
        return noSourcesDeg1;
    }

    uint32_t getNoTargetsDeg1() const {
        return noTargetsDeg1;
    }

    std::vector<uint32_t> getDirect(uint32_t s) {
        std::vector<uint32_t> results;
        root.getDirect(s, 0, N, results);
        return results;
    }

    std::vector<uint32_t> getReverse(uint32_t t) {
        std::vector<uint32_t> results;
        root.getReverse(t, 0, N, results);
        return results;
    }

    const std::vector<uint64_t> &getSources() const {
        return sources;
    }

    const std::vector<uint64_t> &getTargets() const {
        return targets;
    }
};

class SimpleGraph : public Graph {
public:
    std::vector<K2Tree> trees;

    // Pointer to a specific tree in the original graph.
    K2Tree *tree = nullptr;
    bool inverse{false};

    // The number of sources and targets.
    uint32_t noSources{};
    uint32_t noTargets{};

protected:
    uint32_t V;
    uint32_t E;
    uint32_t L;
    const uint32_t k = 8;

    // The number of unique source and target vertices, together with a bitmap representation of the vertices.
    std::vector<uint64_t> sources;
    std::vector<uint64_t> targets;

    uint32_t S;
    uint32_t T;

public:

    SimpleGraph() : V(0), L(0), E(0), S(0), T(0) {};
    ~SimpleGraph() = default;
    explicit SimpleGraph(uint32_t n, uint32_t l);
    explicit SimpleGraph(const std::shared_ptr<SimpleGraph> &g, uint32_t l);

    uint32_t getNoVertices() const override ;
    uint32_t getNoEdges() const override ;
    uint32_t getNoDistinctEdges() const override ;
    uint32_t getNoLabels() const override ;
    uint32_t getNoSources() const;
    uint32_t getNoTargets() const;

    const std::vector<uint64_t> &getSources() const;
    const std::vector<uint64_t> &getTargets() const;

    void addEdge(uint32_t from, uint32_t to, uint32_t edgeLabel) override ;
    void readFromContiguousFile(const std::string &fileName) override ;

    void setNoVertices(uint32_t n);
    void setNoLabels(uint32_t noLabels);

    void setDataStructureSizes();

    void finalizeTree();
};

#endif //QS_SIMPLEGRAPH_H
