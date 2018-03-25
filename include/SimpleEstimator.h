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

public:
    explicit SimpleEstimator(std::shared_ptr<SimpleGraph> &g);
    ~SimpleEstimator() = default;

    void prepare() override ;
    cardStat estimate(RPQTree *q) override ;
    exCardStat doEstimation(RPQTree *q);

    static std::vector<uint64_t> doAnd(const std::vector<uint64_t> *t, const std::vector<uint64_t> *s);

    static uint32_t countBitsSet(std::vector<uint64_t> *result);

    exCardStat estimateLeafNode(RPQTree *q);

    exCardStat estimateSimpleJoin(exCardStat *leftStat, exCardStat *rightStat);

    exCardStat estimateJoin(exCardStat *leftStat, exCardStat *rightStat);

    std::vector<std::vector<joinStat>> joinData;
    std::vector<labelStat> labelData = std::vector<labelStat>();
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

    // The number of unique edges that use this label.
    uint32_t numEdges{};

    // The sources and targets of the label, encoded as a bit collection.
    std::vector<uint64_t> sources;
    std::vector<uint64_t> targets;
    uint32_t numSources{};
    uint32_t numTargets{};

    // The number of sources of the labelStat with source out degree 1 and target out degree 1.
    uint32_t sourceOut1{};
    uint32_t targetIn1{};

public:

    // The transitions, encoded as bit collections.
    // Here, we limit the pair value to char, as it is a value of at most 64.
    labelStat(uint32_t uid, uint32_t id, bool isInverse, unsigned long noBins) {

        labelStat::uid = uid;
        labelStat::id = id;
        labelStat::isInverse = isInverse;
        sources = std::vector<uint64_t>(noBins);
        targets = std::vector<uint64_t>(noBins);
    }

    void setTwin(labelStat* twin) {
        labelStat::twin = twin;
        twin->twin = this;
    }

    const std::vector<uint64_t> &getSources() const { return isInverse ? twin -> getTargets() : sources; }
    const std::vector<uint64_t> &getTargets() const { return isInverse ? twin -> getSources() : targets; }
    const uint32_t getNumSources() const { return isInverse ? twin -> getNumTargets() : numSources; }
    const uint32_t getNumTargets() const { return isInverse ? twin -> getNumSources() : numTargets; }
    const uint32_t getNumEdges() const { return isInverse ? twin -> getNumEdges() : numEdges; }

    void calculateSize() {
        numSources = SimpleEstimator::countBitsSet(&sources);
        numTargets = SimpleEstimator::countBitsSet(&targets);
    }

    void insertEdge(uint32_t s, uint32_t t) {
        sources[s / 64] |= 1ULL << (s % 64);
        targets[t / 64] |= 1ULL << (t % 64);
        numEdges++;
    }

    void incrementSourceOut1() {
        sourceOut1++;
    }

    void incrementTargetIn1() {
        targetIn1++;
    }

    uint32_t getSourceOut1() const {
        return isInverse ? twin -> targetIn1 : sourceOut1;
    }

    uint32_t getTargetIn1() const {
        return isInverse ? twin -> sourceOut1 : targetIn1;
    }

    const void print() const {
        std::cout << "Label id=" << id << (isInverse ? "-" : "+") << ", uid=" << uid << ": #sources=" << std::left << std::setw(8)
                  << getNumSources() << "#targets=" << std::setw(8) << getNumTargets() << "#edges="
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
        commonNodes = SimpleEstimator::doAnd(&source->getTargets(), &target->getSources());
        numCommonNodes = SimpleEstimator::countBitsSet(&commonNodes);

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
