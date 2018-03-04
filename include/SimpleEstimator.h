//
// Created by Nikolay Yakovets on 2018-02-01.
//

#ifndef QS_SIMPLEESTIMATOR_H
#define QS_SIMPLEESTIMATOR_H

#include "Estimator.h"
#include "SimpleGraph.h"
#include <cmath>

// A data structure holding information about vertices in the graph.
struct vertexStat {

    // Keep the in and out degree of each vertex.
    std::map<std::pair<uint32_t, bool>, uint32_t> labelDegrees;
};

struct exCardStat {
    double noOut;
    double noPaths;
    double noIn;

    std::pair<uint32_t, bool> leftLabel;
    std::pair<uint32_t, bool> rightLabel;

    void print() {
        std::cout << "(" << noOut << ", " << noPaths << ", " << noIn << ")" << std::endl;
    }

    explicit operator cardStat() {
        return cardStat {static_cast<uint32_t>(noOut), static_cast<uint32_t>(noPaths), static_cast<uint32_t>(noIn)};
    }
};


// A data structure holding information about labels in the graph.
class labelStat {

    // Whether this stat is a twin.
    bool isTwin;

    // The twin of the current label.
    labelStat* twin;

    // The distinct source and target nodes of the label.
    uint32_t numberOfDistinctSources;
    uint32_t numberOfDistinctTargets;

    /*
     * The in and out degrees of the source/target nodes of the label, represented as a frequency map.
     * In these mappings, the key is the in/out degree and the value is the set of vertices having the in/out degree.
     */
    std::map<uint32_t, uint32_t> sourceOutFrequencies;
    std::map<uint32_t, uint32_t> targetInFrequencies;
    std::map<uint32_t, uint32_t> sourceOutFrequenciesSummation;
    std::map<uint32_t, uint32_t> targetInFrequenciesSummation;

    // The number of edges using the label.
    uint32_t numberOfEdges;

    /*
     * Next, we evaluate the outcomes of joins with all the other labels, to determine:
     *      - the number of target nodes of this label that can be matched to source nodes of another label
     *      - the maximal number of paths that can be created using the two labels.
     *      - The sum of the degrees separately for both sides of the join.
     */
    std::map<std::pair<uint32_t, bool>, uint32_t> numberOfCommonNodesInJoin;
    std::map<std::pair<uint32_t, bool>, uint32_t> numberOfPathsInJoin;
    std::map<std::pair<uint32_t, bool>, uint32_t> degreeSumInLeftJoin;
    std::map<std::pair<uint32_t, bool>, uint32_t> degreeSumInRightJoin;


public:
    void analyze() {
        uint32_t sum = 0;
        for(const auto &labelEntry : sourceOutFrequencies) {
            sourceOutFrequenciesSummation[labelEntry.first] = sum;
            sum += labelEntry.second;
        }

        sum = 0;
        for(const auto &labelEntry : targetInFrequencies) {
            targetInFrequenciesSummation[labelEntry.first] = sum;
            sum += labelEntry.second;
        }
    }


    //region Getters and Setters

    /**
     * Increment the number of edges of the label by the given amount
     *
     * @param n The number of newly discovered edges
     */
    void addEdges(int n) {
        numberOfEdges += n;
    }

    /**
     * Increment the counter of the source vertices.
     */
    void incrementNumberOfDistinctSources() {
        numberOfDistinctSources++;
    }

    /**
     * Increment the counter of the target vertices.
     */
    void incrementNumberOfDistinctTargets() {
        numberOfDistinctTargets++;
    }

    /**
     * Increment the number of encountered source nodes with out degree 'd'
     *
     * @param d The degree that occurred
     */
    void updateSourceOutFrequency(uint32_t d) {
        ++sourceOutFrequencies[d];
    }

    /**
     * Increment the number of encountered target nodes with in degree 'd'
     *
     * @param d The degree that occurred
     */
    void updateTargetInFrequency(uint32_t d) {
        ++targetInFrequencies[d];
    }

    /**
     * Increment the number of nodes this label and the given label have in common.
     *
     * @param label The label with which this label has a node in common.
     */
    void incrementNumberOfCommonNodesInJoin(const std::pair<uint32_t, bool>* label) {
        numberOfCommonNodesInJoin[*label]++;
    }

    /**
     * Add a number of edges that exist between the source nodes of this label and the target nodes of the given label.
     *
     * @param label The label that the join is made with.
     * @param n The number of edges to add.
     */
    void addNumberOfPathsInJoin(const std::pair<uint32_t, bool>* label, int n) {
        numberOfPathsInJoin[*label] += n;
    }

    void addDegreeToLeftSum(const std::pair<uint32_t, bool>* label, int n) {
        degreeSumInLeftJoin[*label] += n;
    }

    void addDegreeToRightSum(const std::pair<uint32_t, bool>* label, int n) {
        degreeSumInRightJoin[*label] += n;
    }

    /**
     * Get the set of distinct source vertices of the label
     *
     * @return The source vertices if this label is a + label, otherwise, the target vertices of the + label
     */
    unsigned long getNumberOfDistinctSources() {
        return isTwin ? twin -> getNumberOfDistinctTargets() : numberOfDistinctSources;
    }

    /**
     * Get the set of distinct source vertices of the label
     *
     * @return The target vertices if this label is a + label, otherwise, the source vertices of the + label
     */
    unsigned long getNumberOfDistinctTargets() {
        return isTwin ? twin -> getNumberOfDistinctSources() : numberOfDistinctTargets;
    }

    /**
     * Get the out degrees of the source nodes expressed as a frequency mapping,
     * where the degree is mapped to the multiplicity.
     *
     * @return The source out frequency if this label is a + label, otherwise, the target in frequencies of the + label
     */
    const std::map<uint32_t, uint32_t> &getSourceOutFrequencies() const {
        return isTwin ? twin -> getTargetInFrequencies() : sourceOutFrequencies;
    }

    /**
     * Get the in degrees of the target nodes expressed as a frequency mapping,
     * where the degree is mapped to the multiplicity.
     *
     * @return The target in frequency if this label is a + label, otherwise, the source out frequencies of the + label
     */
    const std::map<uint32_t, uint32_t> &getTargetInFrequencies() const {
        return isTwin ? twin -> getSourceOutFrequencies() : targetInFrequencies;
    }

    uint32_t getSourceOutFrequencies(uint32_t i) {
        return isTwin ? twin -> getTargetInFrequencies(i) : sourceOutFrequencies[i];
    }

    uint32_t getTargetInFrequencies(uint32_t i) {
        return isTwin ? twin -> getSourceOutFrequencies(i) : targetInFrequencies[i];
    }

    const std::map<uint32_t, uint32_t> &getSourceOutFrequenciesSummation() const {
        return isTwin ? twin -> getTargetInFrequenciesSummation() : sourceOutFrequenciesSummation;
    }

    const std::map<uint32_t, uint32_t> &getTargetInFrequenciesSummation() const {
        return isTwin ? twin -> getSourceOutFrequenciesSummation() : targetInFrequenciesSummation;
    }

    uint32_t getSourceOutFrequenciesSummation(uint32_t i) {
        return isTwin ? twin -> getTargetInFrequenciesSummation(i) : sourceOutFrequenciesSummation[i];
    }

    uint32_t getTargetInFrequenciesSummation(uint32_t i) {
        return isTwin ? twin -> getSourceOutFrequenciesSummation(i) : targetInFrequenciesSummation[i];
    }

    double estimateNumberOfRemovedSourceNodes(double degreeReduction) {
        auto floor = static_cast<uint32_t>(std::ceil(degreeReduction));

        return getSourceOutFrequenciesSummation(floor) + (degreeReduction - floor + 1) * getSourceOutFrequencies(floor);
    }

    double estimateNumberOfRemovedTargetNodes(double degreeReduction) {
        auto floor = static_cast<uint32_t>(std::ceil(degreeReduction));


        auto w = (degreeReduction - floor + 1) * getTargetInFrequencies(floor);
        auto v = getTargetInFrequenciesSummation(floor) + (degreeReduction - floor + 1) * getTargetInFrequencies(floor);

        return getTargetInFrequenciesSummation(floor) + (degreeReduction - floor + 1) * getTargetInFrequencies(floor);
    }

    /**
     * Get the number of edges that use the label
     *
     * @return The number of edges of the + label variant of this label
     */
    int getNumberOfEdges() const {
        return isTwin ? twin -> numberOfEdges : numberOfEdges;
    }

    uint32_t getNumberOfCommonNodesInJoin(const std::pair<uint32_t, bool>* label) {
        return numberOfCommonNodesInJoin[*label];
    }

    uint32_t getNumberOfPathsInJoin(const std::pair<uint32_t, bool>* label) {
        return numberOfPathsInJoin[*label];
    }

    uint32_t getDegreeSumInLeftJoin(const std::pair<uint32_t, bool>* label) {
        return degreeSumInLeftJoin[*label];
    }

    uint32_t getDegreeSumInRightJoin(const std::pair<uint32_t, bool>* label) {
        return degreeSumInRightJoin[*label];
    }

    double getAverageInDegreeLeftJoin(const std::pair<uint32_t, bool>* label) {
        return (double) degreeSumInLeftJoin[*label] / numberOfCommonNodesInJoin[*label];
    }

    double getAverageOutDegreeRightJoin(const std::pair<uint32_t, bool>* label) {
        return (double) degreeSumInRightJoin[*label] / numberOfCommonNodesInJoin[*label];
    }

    const std::map<std::pair<uint32_t, bool>, uint32_t> &getDegreeSumInLeftJoin() const {
        return degreeSumInLeftJoin;
    }

    const std::map<std::pair<uint32_t, bool>, uint32_t> &getDegreeSumInRightJoin() const {
        return degreeSumInRightJoin;
    }

    const std::map<std::pair<uint32_t, bool>, uint32_t> &getNumberOfCommonNodesInJoin() const {
        return numberOfCommonNodesInJoin;
    }

    const std::map<std::pair<uint32_t, bool>, uint32_t> &getNumberOfPathsInJoin() const {
        return numberOfPathsInJoin;
    }

    double getAverageOutDegree() {
        return (double) getNumberOfEdges() / getNumberOfDistinctSources();
    }

    double getAverageInDegree() {
        return (double) getNumberOfEdges() / getNumberOfDistinctTargets();
    }
    //endregion

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
    
    void printData() {
        std::cout << "\t- is used by " << getNumberOfEdges() << " edges" << std::endl;
        std::cout << "\t- has " << getNumberOfDistinctSources() << " distinct sources" << std::endl;
        std::cout << "\t- has " << getNumberOfDistinctTargets() << " distinct targets" << std::endl;
        std::cout << "\t- source out degree: [average]=" << getAverageOutDegree() << std::endl;
        std::cout << "\t- target in degree: [average]=" << getAverageInDegree() << std::endl;

        std::cout << "\t- The number of nodes this label and the given label have in common:" << std::endl << "\t\t[";
        for (auto iter = getNumberOfCommonNodesInJoin().begin(); iter != getNumberOfCommonNodesInJoin().end(); iter++) {
            if (iter != getNumberOfCommonNodesInJoin().begin()) std::cout << ", ";
            std::cout << iter -> second << "*{" << iter -> first.first << ", " << (iter -> first.second ? "true" : "false") << "}";
        }
        std::cout << "]"  << std::endl;

        std::cout << "\t- The total number of paths present between this label and the given label:" << std::endl << "\t\t[";
        for (auto iter = getNumberOfPathsInJoin().begin(); iter != getNumberOfPathsInJoin().end(); iter++) {
            if (iter != getNumberOfPathsInJoin().begin()) std::cout << ", ";
            std::cout << iter -> second << "*{" << iter -> first.first << ", " << (iter -> first.second ? "true" : "false") << "}";
        }
        std::cout << "]"  << std::endl;

        std::cout << "\t- The sum of all degrees of the target nodes in the join:" << std::endl << "\t\t[";
        for (auto iter = getDegreeSumInLeftJoin().begin(); iter != getDegreeSumInLeftJoin().end(); iter++) {
            if (iter != getDegreeSumInLeftJoin().begin()) std::cout << ", ";
            std::cout << iter -> second << "*{" << iter -> first.first << ", " << (iter -> first.second ? "true" : "false") << "}";
        }
        std::cout << "]"  << std::endl;

        std::cout << "\t- The sum of all degrees of the source nodes in the join:" << std::endl << "\t\t[";
        for (auto iter = getDegreeSumInRightJoin().begin(); iter != getDegreeSumInRightJoin().end(); iter++) {
            if (iter != getDegreeSumInRightJoin().begin()) std::cout << ", ";
            std::cout << iter -> second << "*{" << iter -> first.first << ", " << (iter -> first.second ? "true" : "false") << "}";
        }
        std::cout << "]"  << std::endl;

        std::cout << "\t- Source-out frequencies (frequency*{id}):" << std::endl << "\t\t[";
        int sum = 0;
        for (auto iter = getSourceOutFrequencies().begin(); iter != getSourceOutFrequencies().end(); iter++) {
            if (iter != getSourceOutFrequencies().begin()) std::cout << ", ";
            std::cout << iter -> second << "*{" << iter -> first << "}";
            sum += iter -> first * iter -> second;
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
            std::cout << iter -> second << "*{" << iter -> first << "}";
            sum += iter -> first * iter -> second;
        }
        if(sum != getNumberOfEdges()) {
            std::cout << "]"  << std::endl << "\t\tWARNING: sum != " << sum << std::endl;
        } else {
            std::cout << "]"  << std::endl;
        }

        std::cout << "\t- Source-out frequencies sum (sum*{id}):" << std::endl << "\t\t[";
        for (auto iter = getSourceOutFrequenciesSummation().begin(); iter != getSourceOutFrequenciesSummation().end(); iter++) {
            if (iter != getSourceOutFrequenciesSummation().begin()) std::cout << ", ";
            std::cout << iter -> second << "*{" << iter -> first << "}";
        }
        std::cout << "]"  << std::endl;

        std::cout << "\t- Target-in frequencies sum (sum*{id}):" << std::endl << "\t\t[";
        for (auto iter = getTargetInFrequenciesSummation().begin(); iter != getTargetInFrequenciesSummation().end(); iter++) {
            if (iter != getTargetInFrequenciesSummation().begin()) std::cout << ", ";
            std::cout << iter -> second << "*{" << iter -> first << "}";
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

    exCardStat doEstimation(RPQTree *q);
};


#endif //QS_SIMPLEESTIMATOR_H
