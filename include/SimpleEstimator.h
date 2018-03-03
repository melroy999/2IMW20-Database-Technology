//
// Created by Nikolay Yakovets on 2018-02-01.
//

#ifndef QS_SIMPLEESTIMATOR_H
#define QS_SIMPLEESTIMATOR_H

#include "Estimator.h"
#include "SimpleGraph.h"
#include <cmath>

struct pairHasher {
    inline std::size_t operator()(const std::pair<int,int> & v) const {
        return static_cast<size_t>(v.first * 31 + v.second);
    }
};

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
    uint32_t numberOfEdges;

    /*
     * The in and out degrees of the source/target nodes of the label, represented as a frequency map.
     * In these mappings, the key is the in/out degree and the value is the set of vertices having the in/out degree.
     */
    std::map<uint32_t, std::unordered_set<uint32_t>> sourceOutFrequencies;
    std::map<uint32_t, std::unordered_set<uint32_t>> targetInFrequencies;

    /*
     * Next, we evaluate the outcomes of joins with all the other labels, to determine:
     *      - the number of target nodes that can be matched to source nodes in the join
     *      - the number of edges in the next label that are candidates to be in the path
     *      - the distinct number of (s, t) pairs that are in the join
     */
    std::map<std::pair<uint32_t, bool>, uint32_t> distinctTargetsInJoin;
    std::map<std::pair<uint32_t, bool>, uint32_t> numberOfFollowUpEdgesInJoin;
    double squareDegreeEstimation;

public:
    uint32_t sourceOutFrequencyMedian;
    uint32_t targetInFrequencyMedian;

    uint32_t sourceOutFrequencyMin;
    uint32_t sourceOutFrequencyMax;

    uint32_t targetInFrequencyMin;
    uint32_t targetInFrequencyMax;

    void analyzeFrequencyData() {
        if(!getSourceOutFrequencies().empty() && !getTargetInFrequencies().empty()) {
            sourceOutFrequencyMin = getSourceOutFrequencies().begin()->first;
            targetInFrequencyMax = getTargetInFrequencies().begin()->first;

            sourceOutFrequencyMin = getSourceOutFrequencies().begin()->first;
            targetInFrequencyMin = getTargetInFrequencies().begin()->first;

            sourceOutFrequencyMax = getSourceOutFrequencies().rbegin()->first;
            targetInFrequencyMax = getTargetInFrequencies().rbegin()->first;

            unsigned long nodes = getNumberOfDistinctSources();
            for(const auto &f : getSourceOutFrequencies()) {
                nodes -= f.second.size();
                if(nodes <= 0.5 * getNumberOfDistinctSources()) {
                    sourceOutFrequencyMedian = f.first;
                    break;
                }
            }

            nodes = getNumberOfDistinctTargets();
            for(const auto &f : getTargetInFrequencies()) {
                nodes -= f.second.size();
                if(nodes <= 0.5 * getNumberOfDistinctTargets()) {
                    targetInFrequencyMedian = f.first;
                    break;
                }
            }
        }
    }

    void calculateSquareDegreeEstimation() {
        squareDegreeEstimation = 0;
        for(const auto &f : getSourceOutFrequencies()) {
            squareDegreeEstimation += f.second.size() * pow(f.first, 2);
        }

//        for(const auto &f : getTargetInFrequencies()) {
//            squareDegreeEstimation -= f.second.size() * (f.first - 1);
//        }
    }

    //region Getters and setters

    double getSquareDegreeEstimation() const {
        return squareDegreeEstimation;
    }

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
        numberOfEdges += n;
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

    void incrementDistinctTargetsInJoin(const std::pair<uint32_t, bool> label) {
        distinctTargetsInJoin[label]++;
    }


    void updateNumberOfFollowUpEdgesInJoin(const std::pair<uint32_t, bool> label, uint32_t add) {
        numberOfFollowUpEdgesInJoin[label] += add;
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
     * Get the set of distinct source vertices of the label
     *
     * @return The source vertices if this label is a + label, otherwise, the target vertices of the + label
     */
    unsigned long getNumberOfDistinctSources() {
        return isTwin ? twin -> getNumberOfDistinctTargets() : distinctSources.size();
    }

    /**
     * Get the set of distinct source vertices of the label
     *
     * @return The target vertices if this label is a + label, otherwise, the source vertices of the + label
     */
    unsigned long getNumberOfDistinctTargets() {
        return isTwin ? twin -> getNumberOfDistinctSources() : distinctTargets.size();
    }

    /**
     * Get the number of edges that use the label
     *
     * @return The number of edges of the + label variant of this label
     */
    uint32_t getNumberOfEdges() const {
        return isTwin ? twin -> numberOfEdges : numberOfEdges;
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

    const std::map<std::pair<uint32_t, bool>, uint32_t> &getDistinctTargetsInJoin() const {
        return distinctTargetsInJoin;
    }

    uint32_t getDistinctTargetsInJoin(const std::pair<uint32_t, bool> label) {
        return distinctTargetsInJoin[label];
    }

    uint32_t getDistinctSourcesInJoin(const std::pair<uint32_t, bool> label) {
        return twin -> distinctTargetsInJoin[{label.first, !label.second}];
    }

    const std::map<std::pair<uint32_t, bool>, uint32_t> &getNumberOfFollowUpEdgesInJoin() const {
        return numberOfFollowUpEdgesInJoin;
    }

    uint32_t getNumberOfFollowUpEdgesInJoin(const std::pair<uint32_t, bool> label) {
        return numberOfFollowUpEdgesInJoin[label];
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
    //endregion
    
    void printData() {
        std::cout << "\t- is used by " << getNumberOfEdges() << " edges." << std::endl;
        std::cout << "\t- has " << getDistinctSources().size() << " distinct sources." << std::endl;
        std::cout << "\t- has " << getDistinctTargets().size() << " distinct targets." << std::endl;
        std::cout << "\t- source out degree: [min]=" << sourceOutFrequencyMin
                  << ", [average]=" << (double) getNumberOfEdges() / getNumberOfDistinctSources()
                  << ", [median]=" << sourceOutFrequencyMedian
                  << ", [max]=" << sourceOutFrequencyMax << "." << std::endl;
        std::cout << "\t- tatget in degree: [min]=" << targetInFrequencyMin
                  << ", [average]=" << (double) getNumberOfEdges() / getNumberOfDistinctTargets()
                  << ", [median]=" << targetInFrequencyMedian
                  << ", [max]=" << targetInFrequencyMax << "." << std::endl;


        std::cout << "\t- Source-out frequencies (frequency*{id}):" << std::endl << "\t\t[";
        int sum = 0;
        for (auto iter = getSourceOutFrequencies().begin(); iter != getSourceOutFrequencies().end(); iter++) {
            if (iter != getSourceOutFrequencies().begin()) std::cout << ", ";
            std::cout << iter -> second.size() << "*{" << iter -> first << "}";
            sum += iter -> first * iter -> second.size();
        }
        if(sum != getNumberOfEdges()) {
            std::cout << "]"  << std::endl << "\t\tWARNING: sum != " << sum << std::endl;
        } else {
            std::cout << "]"  << std::endl;
        }

        sum = 0;
        std::cout << "\t- Target-in frequencies (frequency*{id}):" << std::endl << "\t\t[";
        for (auto iter = getTargetInFrequencies().begin(); iter != getTargetInFrequencies().end(); iter++) {
            if (iter != getTargetInFrequencies().begin()) std::cout << ", ";
            std::cout << iter -> second.size() << "*{" << iter -> first << "}";
            sum += iter -> first * iter -> second.size();
        }
        if(sum != getNumberOfEdges()) {
            std::cout << "]"  << std::endl << "\t\tWARNING: sum != " << sum << std::endl;
        } else {
            std::cout << "]"  << std::endl;
        }

        std::cout << "\t- Distinct target nodes in join with another label:" << std::endl << "\t\t[";
        for (auto iter = getDistinctTargetsInJoin().begin(); iter != getDistinctTargetsInJoin().end(); iter++) {
            if (iter != getDistinctTargetsInJoin().begin()) std::cout << ", ";
            std::cout << iter -> second << "*{" << iter -> first.first << ", " << (iter -> first.second ? "true" : "false") << "}";
        }
        std::cout << "]"  << std::endl;

        std::cout << "\t- Number of edges in join with another label:" << std::endl << "\t\t[";
        for (auto iter = getNumberOfFollowUpEdgesInJoin().begin(); iter != getNumberOfFollowUpEdgesInJoin().end(); iter++) {
            if (iter != getNumberOfFollowUpEdgesInJoin().begin()) std::cout << ", ";
            std::cout << iter -> second << "*{" << iter -> first.first << ", " << (iter -> first.second ? "true" : "false") << "}";
        }
        std::cout << "]"  << std::endl;
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
