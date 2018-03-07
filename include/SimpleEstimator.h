//
// Created by Nikolay Yakovets on 2018-02-01.
//

#ifndef QS_SIMPLEESTIMATOR_H
#define QS_SIMPLEESTIMATOR_H

#include <set>
#include "Estimator.h"
#include "SimpleGraph.h"

struct labelStat {

    // Whether this stat is a twin.
    bool isInverse;

    // The twin of the current label.
    labelStat* twin;

    // The id of the label.
    uint32_t id;

    // The number of unique edges that use this label.
    uint32_t numEdges;

    // The sources and targets of the label, encoded as a bit collection.
    std::vector<uint64_t> sources;
    std::vector<uint64_t> targets;
    uint32_t numSources;
    uint32_t numTargets;

    // The transitions, encoded as bit collections.
    // Here, we limit the pair value to char, as it is a value of at most 64.
    labelStat(uint32_t id, bool isInverse, unsigned long noBins) {

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
    uint32_t getNumSources() const { return isInverse ? twin -> getNumTargets() : numSources; }
    uint32_t getNumTargets() const { return isInverse ? twin -> getNumSources() : numTargets; }
    uint32_t getNumEdges() const { return isInverse ? twin -> getNumEdges() : numEdges; }
};

struct joinStat;

class SimpleEstimator : public Estimator {

    std::shared_ptr<SimpleGraph> graph;

    std::vector<labelStat> labelData = std::vector<labelStat>();
    std::vector<std::vector<joinStat>> joinData;

public:
    explicit SimpleEstimator(std::shared_ptr<SimpleGraph> &g);
    ~SimpleEstimator() = default;

    void prepare() override ;
    cardStat estimate(RPQTree *q) override ;

    static std::vector<uint64_t> doAnd(const std::vector<uint64_t> *t, const std::vector<uint64_t> *s);

    static uint32_t countBitsSet(std::vector<uint64_t> *result);
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

    void changeBucketSize(unsigned long noBins) {
        sourceNodes.resize(noBins);
        targetNodes.resize(noBins);
    }

    void join(const labelStat* source, const labelStat* target) {
        joinStat::source = source->id;
        joinStat::target = target->id;

        auto targets = &source->getTargets();
        auto sources = &target->getSources();
        commonNodes = SimpleEstimator::doAnd(&source->getTargets(), &target->getSources());
        numCommonNodes = SimpleEstimator::countBitsSet(&commonNodes);
    }
};


#endif //QS_SIMPLEESTIMATOR_H
