//
// Created by Nikolay Yakovets on 2018-02-01.
//

#ifndef QS_SIMPLEESTIMATOR_H
#define QS_SIMPLEESTIMATOR_H

#include <set>
#include <cmath>
#include <iomanip>
#include "Estimator.h"
#include "SimpleGraph.h"

struct exCardStat {
    double noOut;
    double noPaths;
    double noIn;

    // The vertices that have been processed in this statistics collection.
    std::vector<uint32_t> vertices;

    void print() {
        std::cout << "(" << noOut << ", " << noPaths << ", " << noIn << ")" << std::endl;
    }

    explicit operator cardStat() {
        return cardStat {static_cast<uint32_t>(noOut), static_cast<uint32_t>(noPaths), static_cast<uint32_t>(noIn)};
    }
};

struct joinStat;
struct labelStat;

class SimpleEstimator : public Estimator {

    std::shared_ptr<SimpleGraph> graph;

    std::vector<labelStat> labelData = std::vector<labelStat>();
    std::vector<std::vector<joinStat>> joinData;

public:
    explicit SimpleEstimator(std::shared_ptr<SimpleGraph> &g);
    ~SimpleEstimator() = default;

    void prepare() override ;
    cardStat estimate(RPQTree *q) override ;
    exCardStat doEstimation(RPQTree *q);

    exCardStat estimateLeafNode(RPQTree *q);

    exCardStat estimateSimpleJoin(exCardStat *leftStat, exCardStat *rightStat);

    exCardStat estimateJoin(exCardStat *leftStat, exCardStat *rightStat);
};

struct labelStat {

    // Whether this stat is a twin.
    bool isInverse;

    // The twin of the current label.
    labelStat* twin{};

    // The id of the label.
    uint32_t uid;
    uint32_t id;

private:

    // A reference to the appropriate KTree.
    // It contains the sources, targets, number of sources and number of targets, and the number of edges.
    K2Tree* tree{};

public:

    // The transitions, encoded as bit collections.
    // Here, we limit the pair value to char, as it is a value of at most 64.
    labelStat(uint32_t uid, uint32_t id, bool isInverse, unsigned long noBins) {

        labelStat::uid = uid;
        labelStat::id = id;
        labelStat::isInverse = isInverse;
    }

    void setTwin(labelStat* twin) {
        labelStat::twin = twin;
        twin->twin = this;
    }

    const std::vector<uint64_t> &getSources() const { return isInverse ? twin->getTargets() : tree->getSources(); }
    const std::vector<uint64_t> &getTargets() const { return isInverse ? twin->getSources() : tree->getTargets(); }
    const uint32_t getNoSources() const { return isInverse ? twin->getNoTargets() : tree->getNoSources(); }
    const uint32_t getNoTargets() const { return isInverse ? twin->getNoSources() : tree->getNoTargets(); }
    const uint32_t getNumEdges() const { return isInverse ? twin->getNumEdges() : tree->getNoEdges(); }
    const uint32_t getSourceOut1() const { return isInverse ? twin->getTargetIn1() : tree->getNoSourcesDeg1(); }
    const uint32_t getTargetIn1() const { return isInverse ? twin->getSourceOut1() : tree->getNoTargetsDeg1(); }

    void setTree(K2Tree *tree) {
        labelStat::tree = tree;
    }

    const void print() const {
        std::cout << "Label id=" << id << (isInverse ? "-" : "+") << ", uid=" << uid << ": #sources=" << std::left << std::setw(8)
                  << getNoSources() << "#targets=" << std::setw(8) << getNoTargets() << "#edges="
                  << std::setw(8) << getNumEdges() << "#sourceOut1=" << std::setw(8) << getSourceOut1()
                  << "#targetIn1=" << std::setw(8) << getTargetIn1() << std::endl;
    }
};

struct joinStat {
    // The id of the label we join.
    uint32_t source;
    uint32_t target;

    // The nodes they have in common.
    std::vector<uint64_t> commonNodes;
    std::vector<uint64_t> sourceNodes;
    std::vector<uint64_t> targetNodes;
    uint32_t numCommonNodes;
    uint32_t numSourceNodes;
    uint32_t numTargetNodes;

    // The number of (non-distinct) paths between the source and the target label.
    uint32_t numPaths;
    uint32_t numSourceEdges;
    uint32_t numTargetEdges;
    uint32_t uniquePaths;

    void changeBucketSize(unsigned long noBins) {
        sourceNodes.resize(noBins);
        targetNodes.resize(noBins);
    }

    void join(const labelStat* source, const labelStat* target) {
        joinStat::source = source->uid;
        joinStat::target = target->uid;

        auto targets = &source->getTargets();
        auto sources = &target->getSources();
        commonNodes = doAnd(&source->getTargets(), &target->getSources());
        numCommonNodes = popcountVector(&commonNodes);

        // If the set of common nodes is empty, free the space of the vector.
        if(numCommonNodes == 0) {
            commonNodes.resize(0);
        }
    }

    const void print() const {
        std::cout << "Join " << source << "x" << target << ": #numCommonNodes=" << std::left << std::setw(8)
                  << numCommonNodes << "#numSourceNodes=" << std::setw(10) << numSourceNodes << "#numTargetNodes="
                  << std::setw(8) << numTargetNodes << "#numPaths="  << std::setw(8) << numPaths << "#numSourceEdges="
                  << std::setw(8) << numSourceEdges << "#numTargetEdges="  << std::setw(8) << numTargetEdges << "#uniquePaths="  << std::setw(8) << uniquePaths << std::endl;

    }
};


#endif //QS_SIMPLEESTIMATOR_H
