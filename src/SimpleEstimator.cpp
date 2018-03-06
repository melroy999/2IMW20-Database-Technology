//
// Created by Nikolay Yakovets on 2018-02-01.
//

#include "SimpleGraph.h"
#include "SimpleEstimator.h"


SimpleEstimator::SimpleEstimator(std::shared_ptr<SimpleGraph> &g){

    // works only with SimpleGraph
    graph = g;
}

// sort on the second item in the pair, then on the first (ascending order)
bool sortPairs2(const std::pair<uint32_t,uint32_t> &a, const std::pair<uint32_t,uint32_t> &b) {
    if (a.second < b.second) return true;
    if (a.second == b.second) return a.first < b.first;
    return false;
}

void SimpleEstimator::prepare() {

    // Initialise vectors with lengths corresponding to num vertices or labels.
    vertexData = std::vector<vertexStat>(graph -> getNoVertices());
    labelData = std::map<std::pair<uint32_t, bool>, labelStat>();

    // Initialize vectors separately for both the positive and negative labels.
    for(uint32_t i = 0; i < graph -> getNoLabels(); i++) {
        labelData[{i, true}] = labelStat();
        labelData[{i, false}] = labelStat();
        labelData[{i, false}].setTwin(&labelData[{i, true}]);
    }

    /*
     * Per label, gather the set of distinct source and target vertices,
     * and the non-distinct count of source and target vertices.
     */
    for(uint32_t i = 0; i < graph -> getNoVertices(); i++) {

        std::sort(graph -> adj[i].begin(), graph -> adj[i].end(), sortPairs2);

        uint32_t prevTarget = 0;
        uint32_t prevLabel = 0;
        bool first = true;

        // Now use the entries in the adjacency matrix to calculate the label data.
        for(const auto &labelTgtPair : graph -> adj[i]) {

            if (first || !(prevTarget == labelTgtPair.second && prevLabel == labelTgtPair.first)) {

                // Update the appropriate out degree counter in the label degree mapping.
                ++vertexData[i].labelDegrees[{labelTgtPair.first, true}];

                first = false;
                prevTarget = labelTgtPair.second;
                prevLabel = labelTgtPair.first;
            }
        }

        std::sort(graph -> reverse_adj[i].begin(), graph -> reverse_adj[i].end(), sortPairs2);

        prevTarget = 0;
        prevLabel = 0;
        first = true;

        // Do the same for the reverse adjacency matrix.
        for(const auto &labelTgtPair : graph -> reverse_adj[i]) {

            if (first || !(prevTarget == labelTgtPair.second && prevLabel == labelTgtPair.first)) {

                // Update the appropriate in degree counter in the label degree mapping.
                ++vertexData[i].labelDegrees[{labelTgtPair.first, false}];

                first = false;
                prevTarget = labelTgtPair.second;
                prevLabel = labelTgtPair.first;
            }
        }
    }

    /*
     * Start gathering data for all the positive labels.
     *
     * First, we want to gather the distinct source and target vertices per label.
     */
    for(uint32_t i = 0; i < graph -> getNoVertices(); i++) {

        for(const auto &labelEntry : vertexData[i].labelDegrees) {

            // Extract the label for readability.
            const std::pair<uint32_t, bool> *label = &labelEntry.first;
            const std::pair<uint32_t, bool> iLabel = {label->first, !label->second};

            /*
             * In the pair:
             *  - first.first corresponds to the label id.
             *  - first.second corresponds to the label direction, true is forwards, false is backwards.
             *  - second corresponds to the degree.
             */

            if(label->second) {
                labelData[{label->first, true}].addEdges(labelEntry.second);
                labelData[{label->first, true}].incrementNumberOfDistinctSources();
                labelData[{label->first, true}].updateSourceOutFrequency(labelEntry.second);
            } else {
                labelData[{label->first, true}].incrementNumberOfDistinctTargets();
                labelData[{label->first, true}].updateTargetInFrequency(labelEntry.second);
            }

            /*
             * Use the vertex data collection to determine join statistics.
             * Here, we should consider the label above as the in label, with others being out labels.
             * In other words, the label's direction should be inverted.
             */
            for(const auto &labelEntry2 : vertexData[i].labelDegrees) {
                labelData[iLabel].incrementNumberOfCommonNodesInJoin(&labelEntry2.first);
                labelData[iLabel].addNumberOfPathsInJoin(&labelEntry2.first, labelEntry.second * labelEntry2.second);
                labelData[iLabel].addDegreeToLeftSum(&labelEntry2.first, vertexData[i].labelDegrees[*label]);
                labelData[iLabel].addDegreeToRightSum(&labelEntry2.first, vertexData[i].labelDegrees[labelEntry2.first]);
            }
        }
    }

    for(uint32_t i = 0; i < graph -> getNoLabels(); i++) {
        labelData[{i, true}].analyze();
    }

    // Print the debug data.
    printDebugData();
}


void estimateInAndOutCardinality(exCardStat* leftStat, exCardStat* rightStat, labelStat* leftData, labelStat* rightData, double* noOut, double* noIn) {
    /*
         * During the join, use the fact that the difference between the two sets of
         * target/source vertices we merge on is not necessarily the empty set.
         *
         * We know for example that certain target vertices in the left subtree
         * are not succeeded by the starting label of the right tree.
         *
         * On the other hand, we have that certain vertices in the right subtree
         * are not proceeded by the end label of the left tree.
         */
    // What is the number of edges that have been terminated and initialized during the merge?
    uint32_t terminatedEdges = leftData->getNumberOfEdges() - leftData->getDegreeSumInLeftJoin(&rightStat->leftLabel);
    uint32_t initializedEdges = rightData->getNumberOfEdges() - leftData->getDegreeSumInRightJoin(&rightStat->leftLabel);

    // What is the probability that a vertex ends up with no connected edges, when an edge gets removed?
    // In other words, what is the chance that a vertex of degree one has its edge removed?
    double pSourceRemove = (double) leftData->getSourceOutFrequencies(1) / leftData->getNumberOfEdges();
    double pTargetRemove = (double) rightData->getTargetInFrequencies(1) / rightData->getNumberOfEdges();

    // What about higher degree edges? In certain cases we have a very low amount of degree 1 edges,
    // and should thus compensate.

    // So, what is the expected number of source and target vertices given the above probabilities?
    // Here we do not simply subtract, as we would risk getting negative numbers.
    *noOut = leftStat->noOut * (1 - (pSourceRemove * terminatedEdges) / leftData->getNumberOfDistinctSources());
    *noIn = rightStat->noIn * (1 - (pTargetRemove * initializedEdges) / rightData->getNumberOfDistinctTargets());
}

exCardStat SimpleEstimator::doEstimation(RPQTree *q) {
    // Estimate only if we are a leaf or want to join.
    if(q -> isLeaf()) {

        // Estimating a leaf is simple, just do a lookup in the data.
        std::string data = q->data;
        std::pair<uint32_t, bool> label =
                {static_cast<const uint32_t &>(std::stoi(data.substr(0, data.size() - 1))), data.back() == '+'};

        return exCardStat {static_cast<double>(labelData[label].getNumberOfDistinctSources()),
                           static_cast<double>(labelData[label].getNumberOfEdges()),
                           static_cast<double>(labelData[label].getNumberOfDistinctTargets()),
                           label,
                           label
        };

    } else if(q -> isConcat()) {

        exCardStat leftStat = doEstimation(q -> left);
        exCardStat rightStat = doEstimation(q -> right);
//        leftStat.print();
//        rightStat.print();
        
        auto leftData = &labelData[leftStat.rightLabel];
        auto rightData = &labelData[rightStat.leftLabel];

        double noOut, noIn = 0;
        estimateInAndOutCardinality(&leftStat, &rightStat, leftData, rightData, &noOut, &noIn);

        // First of all, we have exact data on the number of (non-distinct) paths in the join,
        // and the average expected number of paths in the join using the average in and out degree.
        uint32_t maxPaths = leftData->getNumberOfPathsInJoin(&rightStat.leftLabel);

        // Using the above, we can determine the expected number of follow up edges.
        // Note here that we use the raw path stats, since the max paths measure combined with
        // the degree sum for a specific label already covers for the terminated/initialized edges.
        double leftPathEstimation = leftStat.noPaths * (double) maxPaths / leftData->getDegreeSumInLeftJoin(&rightStat.leftLabel);
        double rightPathEstimation = rightStat.noPaths * (double) maxPaths / leftData->getDegreeSumInRightJoin(&rightStat.leftLabel);
//        double leftPathEstimation = leftStat.noPaths * (double) maxPaths / leftData->getNumberOfEdges();
//        double rightPathEstimation = rightStat.noPaths * (double) maxPaths / rightData->getNumberOfEdges();


        // We can put a better bound on noOut and noIn by observing the join data.
        if(leftStat.leftLabel == leftStat.rightLabel) {
            noOut = std::min((double) leftData->getDegreeSumInLeftJoin(&rightStat.leftLabel), noOut);
        }

        if(rightStat.leftLabel == rightStat.rightLabel) {
            noIn = std::min((double) leftData->getDegreeSumInRightJoin(&rightStat.leftLabel), noIn);
        }

//        std::cout << leftPathEstimation << ", " << rightPathEstimation << std::endl;


        // How many of the given paths is a duplicate?
        // For now, we just take the probability that the source and the target are equal:
//        double pDuplicate = 1.0 / (leftData->getNumberOfDistinctSources() - pSourceRemove * terminatedEdges + rightData->getNumberOfDistinctTargets() - pTargetRemove * initializedEdges);

        double noPaths = ((leftPathEstimation + rightPathEstimation) / 2);

        if(leftStat.rightLabel.first == rightStat.leftLabel.first && leftStat.rightLabel.second != rightStat.leftLabel.second) {
            // We know that at least the number of edges minus the number of source nodes are duplicates.
            unsigned long minimumDuplicateCount = leftData->getNumberOfEdges() - leftData->getNumberOfDistinctSources();
            noPaths -= minimumDuplicateCount;
        }


        exCardStat result = exCardStat {noOut,
                                        std::min(noPaths, noOut * noIn),
                                        noIn,
                                        leftStat.leftLabel,
                                        rightStat.rightLabel
        };

//        q->print();
//        std::cout << std::endl;
//        result.print();


        // Try to merge the two estimates with the join cardinality formula, with uniformity assumptions.
        return result;
    }

    // Return the estimate in the form {#outNodes, #paths, #inNodes}
    return exCardStat {0, 0, 0};
}

cardStat SimpleEstimator::estimate(RPQTree *q) {
//    std::cout << std::endl;

//    auto estimation = doEstimation(q);
//    estimation.print();

    return static_cast<cardStat>(doEstimation(q));
}

/**
 * Display the estimator's saved data in readable form.
 */
void SimpleEstimator::printDebugData() {
    std::cout << std::endl;

    for(uint32_t i = 0; i < graph -> getNoLabels(); i++) {
        for(bool b : {true, false}) {
            auto v = labelData[{i, b}];

            std::cout << "Label {" << i << ", " << (b ? "true" : "false") << "}:" << std::endl;
            v.printData();
            std::cout << std::endl;
        }
    }
}
