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
    uint32_t noEdges;

    /*
     * The in and out degrees of the source/target nodes of the label, represented as a frequency map.
     * In these mappings, the key is the in/out degree and the value is the set of vertices having the in/out degree.
     */
    std::map<uint32_t, std::unordered_set<uint32_t>> sourceOutFrequencies;
    std::map<uint32_t, std::unordered_set<uint32_t>> targetInFrequencies;

    /*
     * Count how many of the target vertices are followed by an edge having the specified label.
     */
    std::map<std::pair<uint32_t, bool>, std::unordered_set<uint32_t>> distinctTargetNodesFollowedByLabel;
    std::map<std::pair<uint32_t, bool>, uint32_t> noEdgesFollowingTargetNodesByLabel;

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
        noEdges += n;
    }

    /**
     * Increment the number of encountered source nodes with out degree 'd'
     *
     * @param d The degree that occurred
     * @param v The node the degree occurred for
     */
    void updateSourceOutFrequency(uint32_t d, uint32_t v) {
        sourceOutFrequencies[d].insert(v);
    }

    /**
     * Increment the number of encountered target nodes with in degree 'd'
     *
     * @param d The degree that occurred
     * @param v The node the degree occurred for
     */
    void updateTargetInFrequency(uint32_t d, uint32_t v) {
        targetInFrequencies[d].insert(v);
    }

    /**
     * Add a target node to the set of distinct vertices that are followed by the given label
     *
     * @param label The label for which v has an outgoing edge.
     * @param v The target vertex that is the source of the given label.
     */
    void updateDistinctTargetNodesFollowedByLabel(std::pair<uint32_t , bool> label, uint32_t v) {
        distinctTargetNodesFollowedByLabel[label].insert(v);
    }

    /**
     * Update the number of edges with the given label that originate from the target vertices.
     *
     * @param label The label which the edge belongs to.
     * @param n The number of extra edges that have been found with the given label.
     */
    void updateNoEdgesFollowingTargetNodesByLabel(std::pair<uint32_t , bool> label, uint32_t n) {
        noEdgesFollowingTargetNodesByLabel[label] += n;
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
    uint32_t getNoEdges() const {
        return isTwin ? twin -> noEdges : noEdges;
    }

    /**
     * Get the out degrees of the source nodes expressed as a frequency mapping,
     * where the degree is mapped to the multiplicity.
     *
     * @return The source out frequency if this label is a + label, otherwise, the target in frequencies of the + label
     */
    const std::map<uint32_t, std::unordered_set<uint32_t>> &getSourceOutFrequencies() const {
        return isTwin ? twin -> getTargetInFrequencies() : sourceOutFrequencies;
    }

    /**
     * Get the in degrees of the target nodes expressed as a frequency mapping,
     * where the degree is mapped to the multiplicity.
     *
     * @return The target in frequency if this label is a + label, otherwise, the source out frequencies of the + label
     */
    const std::map<uint32_t, std::unordered_set<uint32_t>> &getTargetInFrequencies() const {
        return isTwin ? twin -> getSourceOutFrequencies() : targetInFrequencies;
    }

    const std::map<std::pair<uint32_t, bool>, std::unordered_set<uint32_t>> &
    getDistinctTargetNodesFollowedByLabel() const {
        return distinctTargetNodesFollowedByLabel;
    }

    const std::map<std::pair<uint32_t, bool>, uint32_t> &getNoEdgesFollowingTargetNodesByLabel() const {
        return noEdgesFollowingTargetNodesByLabel;
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

    std::vector<std::pair<uint32_t, bool>> parseTreeToList(RPQTree *q);

    void printDebugData();

};


#endif //QS_SIMPLEESTIMATOR_H
