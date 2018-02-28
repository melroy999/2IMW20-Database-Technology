//
// Created by Nikolay Yakovets on 2018-02-01.
//

#ifndef QS_SIMPLEESTIMATOR_H
#define QS_SIMPLEESTIMATOR_H

#include "Estimator.h"
#include "SimpleGraph.h"

// A data structure holding information about vertices in the graph.
struct vertexStat {

    // Keep the in and out degree of each vertex.
    std::map<std::pair<uint32_t, bool>, uint32_t> labelDegrees;
};


// A data structure holding information about labels in the graph.
class labelStat {

    // Whether this stat is a twin.
    bool isTwin;

    // The twin of the current label.
    labelStat* twin;

    // The distinct source and target nodes of the label.
    std::unordered_set<uint32_t> distinctSources;
    std::unordered_set<uint32_t> distinctTargets;

    // The number of edges using the label.
    uint32_t edges;

public:

    /**
     * Add the given vertex to the list of distinct source vertices
     *
     * @param v The id of the source vertex
     */
    void addSourceVertex(uint32_t v) {
        distinctSources.insert(v);
    }

    /**
     * Add the given vertex to the list of distinct target vertices
     *
     * @param v The id of the target vertex
     */
    void addTargetVertex(uint32_t v) {
        distinctTargets.insert(v);
    }

    /**
     * Increment the number of edges of the label by the given amount
     *
     * @param n The number of newly discovered edges
     */
    void addEdges(uint32_t n) {
        edges += n;
    }

    /**
     * Get the set of distinct source vertices of the label
     *
     * @return The source vertices if this label is a + label, otherwise, the target vertices of the + label
     */
    const std::unordered_set<uint32_t> &getDistinctSources() const {
        return isTwin ? twin -> getDistinctTargets() : distinctSources;
    }

    /**
     * Get the set of distinct source vertices of the label
     *
     * @return The target vertices if this label is a + label, otherwise, the source vertices of the + label
     */
    const std::unordered_set<uint32_t> &getDistinctTargets() const {
        return isTwin ? twin -> getDistinctSources() : distinctTargets;
    }

    /**
     * Get the number of edges that use the label
     *
     * @return The number of edges of the + label variant of this label
     */
    uint32_t getEdges() const {
        return isTwin ? twin -> edges : edges;
    }

    /**
     * Set the given + label statistic collection to be the twin of this label, converting this label to a twin label
     *
     * @param _twin The + label to make a twin pair with
     */
    void setTwin(labelStat* _twin) {
        twin = _twin;
        isTwin = true;

        // Make sure that the twin is aware that this object is its twin now.
        _twin->twin = this;
    }
};



class SimpleEstimator : public Estimator {

    std::shared_ptr<SimpleGraph> graph;

    // Variables containing graph data.
    std::vector<vertexStat> vertexData;
    std::map<std::pair<uint32_t, bool>, labelStat> labelData;

public:
    explicit SimpleEstimator(std::shared_ptr<SimpleGraph> &g);
    ~SimpleEstimator() = default;

    void prepare() override ;
    cardStat estimate(RPQTree *q) override ;

    std::vector<std::pair<int, bool>> parseTreeToList(RPQTree *q);

    void printDebugData();

};


#endif //QS_SIMPLEESTIMATOR_H
