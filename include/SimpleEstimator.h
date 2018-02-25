//
// Created by Nikolay Yakovets on 2018-02-01.
//

#ifndef QS_SIMPLEESTIMATOR_H
#define QS_SIMPLEESTIMATOR_H

#include "Estimator.h"
#include "SimpleGraph.h"

// A data structure holding information about vertices in the graph.
struct vertexStat {
    // The in and out degree of the vertex.
    uint32_t outDegree;
    uint32_t inDegree;

    // Also keep the in and out degree of each individual label.
    std::map<uint32_t, uint32_t> labelOutDegrees;
    std::map<uint32_t, uint32_t> labelInDegrees;
};

// A data structure holding information about labels in the graph.
struct labelStat {
    // The set of distinct source and target vertices in the label collection.
    std::unordered_set<uint32_t> distinctSources;
    std::unordered_set<uint32_t> distinctTargets;

    // The total number of edges bearing this label.
    uint32_t noEdges;
};


class SimpleEstimator : public Estimator {

    std::shared_ptr<SimpleGraph> graph;

    // Variables containing graph data.
    std::vector<vertexStat> vertexData;
    std::vector<labelStat> labelData;

public:
    explicit SimpleEstimator(std::shared_ptr<SimpleGraph> &g);
    ~SimpleEstimator() = default;

    void prepare() override ;
    cardStat estimate(RPQTree *q) override ;

    std::vector<std::pair<int, bool>> parseTreeToList(RPQTree *q);
};


#endif //QS_SIMPLEESTIMATOR_H
