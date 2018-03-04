//
// Created by Nikolay Yakovets on 2018-02-01.
//

#include "SimpleGraph.h"
#include "SimpleEstimator.h"


SimpleEstimator::SimpleEstimator(std::shared_ptr<SimpleGraph> &g){

    // works only with SimpleGraph
    graph = g;
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

exCardStat SimpleEstimator::doEstimation(RPQTree *q) {
    // Estimate only if we are a leaf or want to join.
    if(q -> isLeaf()) {

        // Estimating a leaf is simple, just do a lookup in the data.
        std::string data = q->data;
        std::pair<uint32_t, bool> label =
                {static_cast<const uint32_t &>(std::stoi(data.substr(0, data.size() - 1))), data.back() == '+'};

        exCardStat result = exCardStat {static_cast<double>(labelData[label].getNumberOfDistinctSources()),
                                        static_cast<double>(labelData[label].getNumberOfEdges()),
                                        static_cast<double>(labelData[label].getNumberOfDistinctTargets()),
                                        label,
                                        label
        };

//        q->print();
//        std::cout << std::endl;
//        result.print();

        return result;

    } else if(q -> isConcat()) {

        exCardStat leftStat = doEstimation(q -> left);
        exCardStat rightStat = doEstimation(q -> right);
        
        auto leftData = &labelData[leftStat.rightLabel];
        auto rightData = &labelData[rightStat.leftLabel];

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
        uint32_t commonVertices = leftData->getNumberOfCommonNodesInJoin(&rightStat.leftLabel);

        // What is the number of edges that have been terminated and initialized during the merge?
        double terminatedEdges = leftData->getNumberOfEdges() - leftData->getDegreeSumInLeftJoin(&rightStat.leftLabel);
        double initializedEdges = rightData->getNumberOfEdges() - leftData->getDegreeSumInRightJoin(&rightStat.leftLabel);

        // What is the probability that a vertex ends up with no connected edges, when an edge gets removed?
        // In other words, what is the chance that a vertex of degree one has its edge removed?
        double pSourceRemove = (double) leftData->getSourceOutFrequencies(1) / leftData->getNumberOfEdges();
        double pTargetRemove = (double) rightData->getTargetInFrequencies(1) / rightData->getNumberOfEdges();

        // So, what is the expected number of removed vertices?
        // Here we do not simply subtract, as we would risk getting negative numbers.
        double noOut = leftStat.noOut * (1 - (pSourceRemove * terminatedEdges) / leftData->getNumberOfDistinctSources());
        double noIn = rightStat.noIn * (1 - (pTargetRemove * initializedEdges) / rightData->getNumberOfDistinctTargets());

        // What is the effect of the above on the total number of distinct paths?
        // First of all, we have exact data on the number of (non-distinct) paths in the join.
        uint32_t maxPaths = leftData->getNumberOfPathsInJoin(&rightStat.leftLabel);


        double noPaths = rightStat.noPaths * leftStat.noPaths
                         / std::max(leftStat.noIn, rightStat.noOut);


//        if(true || leftStat.rightLabel.first == rightStat.leftLabel.first && leftStat.rightLabel.second != rightStat.leftLabel.second) {
//            noPaths = rightStat.noPaths * leftStat.noPaths * maxPaths / (leftData->getNumberOfEdges() * rightData->getNumberOfEdges());
//        }


//        double noPaths = (leftStat.noPaths * (1 - terminatedEdges / leftData->getNumberOfEdges())) * (rightStat.noPaths * (1 - initializedEdges / rightData->getNumberOfEdges()) / std::max(noIn, noOut));



        std::cout << maxPaths << ", " << maxPaths / noPaths << ", " << (double) maxPaths / noPaths << std::endl;

        /*
         * We know that when we join leftStat.rightLabel and rightStat.leftLabel, that we have at most
         * leftStat.rightLabel.numberOfPathsInJoin[rightStat.leftLabel] not necessarily distinct paths.
         */

        exCardStat result = exCardStat {noOut,
                                        noPaths,
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
    std::cout << std::endl;

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
