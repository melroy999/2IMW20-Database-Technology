//
// Created by Nikolay Yakovets on 2018-02-01.
//

#ifndef QS_SIMPLEESTIMATOR_H
#define QS_SIMPLEESTIMATOR_H

#include "Estimator.h"
#include "SimpleGraph.h"

// A data structure holding information about vertices in the graph.
struct vertexStat {
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

    /*
     * The in and out degrees of the source/target nodes of the label, represented as a frequency map.
     * In these mappings, the key is the in/out degree and the value is the amount of vertices having the in/out degree.
     */
    std::map<uint32_t, uint32_t> sourceOutFrequencies;
    std::map<uint32_t, uint32_t> targetInFrequencies;

    /*
     * Count how many of the target vertices are followed by an edge having the specified label.
     * The key -1 denotes the number of vertices not followed by an edge.
     */
    std::map<int, std::unordered_set<uint32_t>> distinctTargetNodesFollowedByLabel;
    std::map<uint32_t, uint32_t> noEdgesFollowingTargetNodesByLabel;

    /*
     * Count how many of the source vertices are followed by an edge having the specified label,
     * after the application of the current label.
     */
    std::map<int, std::unordered_set<uint32_t>> distinctSourceNodesFollowedByLabel;
    std::map<uint32_t, uint32_t> noEdgesFollowingSourceNodesByLabel;

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

    void printDebugData();
};


#endif //QS_SIMPLEESTIMATOR_H
