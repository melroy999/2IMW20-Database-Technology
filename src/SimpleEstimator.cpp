//
// Created by Nikolay Yakovets on 2018-02-01.
//

#include "SimpleGraph.h"
#include "SimpleEstimator.h"
#include <cmath>


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

        // Now use the entries in the adjacency matrix to calculate the label data.
        for(auto v : graph -> adj[i]) {

            // Update the appropriate out degree counter in the label degree mapping.
            ++vertexData[i].labelDegrees[{v.first, true}];

        }

        // Do the same for the reverse adjacency matrix.
        for(auto v : graph -> reverse_adj[i]) {

            // Update the appropriate in degree counter in the label degree mapping.
            ++vertexData[i].labelDegrees[{v.first, false}];
        }
    }

    /*
     * Start gathering data for all the positive labels.
     */
    for(uint32_t i = 0; i < graph -> getNoVertices(); i++) {
        for(auto pair : vertexData[i].labelDegrees) {
            if(pair.second > 0) {
                /*
                 * In the pair:
                 *  - first.first corresponds to the label id.
                 *  - first.second corresponds to the label direction, true is forwards, false is backwards.
                 *  - second corresponds to the degree.
                 */
                if(pair.first.second) {
                    labelData[{pair.first.first, true}].addEdges(pair.second);
                    labelData[{pair.first.first, true}].addSourceVertex(i);
                } else {
                    labelData[{pair.first.first, true}].addTargetVertex(i);
                }
            }
        }
    }

//    for(uint32_t i = 0; i < graph -> getNoLabels(); i++) {
//        std::cout << labelData[{i, true}].getDistinctSources().size() << "," << labelData[{i, false}].getDistinctTargets().size() << std::endl;
//        std::cout << labelData[{i, true}].getDistinctTargets().size() << "," << labelData[{i, false}].getDistinctSources().size() << std::endl;
//    }

    // Print the debug data.
    printDebugData();
}



cardStat SimpleEstimator::estimate(RPQTree *q) {

    // Convert the parse tree to a more straightforward form to work with i.e. a list.
    std::vector<std::pair<int, bool>> result = parseTreeToList(q);

    /*
     * We can move in two directions: forward and backwards. We start by estimating s and t by
     * getting the number of distinct sources/targets of the first and last label in the chain respectively.
     */
//    auto s = (result.front().second ? labelData[result.front().first].distinctSources : labelData[result.front().first].distinctTargets).size();
//    auto t = (result.back().second ? labelData[result.back().first].distinctTargets : labelData[result.back().first].distinctSources).size();
//
//    // To estimate the number of paths, use the average degree of source vertices bearing the label.
//    float noPaths = labelData[result[0].first].noEdges;
//
//    for(int i = 1; i < result.size(); i++) {
//
//        auto l = result[i];
//        noPaths *= ((float) labelData[l.first].noEdges) / (l.second ? labelData[l.first].distinctTargets : labelData[l.first].distinctSources).size();
//    }
//
//    // Return the estimate in the form {#outNodes, #paths, #inNodes}
//    return cardStat {static_cast<uint32_t>(s), static_cast<uint32_t>(noPaths), static_cast<uint32_t>(t)};

    return cardStat {0, 0, 0};
}

/**
 * Convert the RPQTree to a vector representation of pairs of labels and directions
 * @param q - the tree
 * @return the vector list of labels
 */
std::vector<std::pair<int, bool>> SimpleEstimator::parseTreeToList(RPQTree *q) {

    // Create a vector that will contain the result.
    std::vector<std::pair<int, bool>> result;

    // We keep a queue of nodes to visit.
    std::stack<RPQTree*> visited;
    visited.push(q);

    // Do breadth first search.
    while(!visited.empty()) {

        // Pop the current top node.
        RPQTree* current = visited.top();
        visited.pop();

        if(current->isLeaf()) {

            std::string data = current->data;
            result.emplace_back(std::stoi(data.substr(0, data.size() - 1)), data.back() == '+');
        } else {

            if(current->right != nullptr) {

                visited.push(current->right);
            }

            if(current->left != nullptr) {

                visited.push(current->left);
            }
        }
    }

    // Return the result.
    return result;
}

/**
 * Display the estimator's saved data in readable form.
 */
void SimpleEstimator::printDebugData() {
    std::cout << std::endl;
//
//    for(uint32_t i = 0; i < 100; i++) {
//        auto v = vertexData[i];
//
//        std::cout << "Vertex " << i << ":" << std::endl;
//        std::cout << "\t- Degrees per label:" << std::endl << "\t\t[";
//        for (auto iter = v.labelDegrees.begin(); iter != v.labelDegrees.end(); iter++) {
//            if (iter != v.labelDegrees.begin()) std::cout << ", ";
//            std::cout << "\"" << iter -> first << "\" * " << iter -> second;
//        }
//        std::cout << "]"  << std::endl << std::endl;
//    }



//    for(uint32_t i = 0; i < graph -> getNoLabels(); i++) {
//        auto v = labelData[i];
//
//        std::cout << "Label " << i << ":" << std::endl;
//        std::cout << "\t- is used by " << v.noEdges << " edges." << std::endl;
//        std::cout << "\t- has " << v.distinctSources.size() << " distinct sources." << std::endl;
//        std::cout << "\t- has " << v.distinctTargets.size() << " distinct targets." << std::endl;
//
//        std::cout << "\t- Source-out frequencies (frequency*{id}):" << std::endl << "\t\t[";
//        int sum = 0;
//        for (auto iter = v.sourceOutFrequencies.begin(); iter != v.sourceOutFrequencies.end(); iter++) {
//            if (iter != v.sourceOutFrequencies.begin()) std::cout << ", ";
//            std::cout << iter -> second << "*{" << iter -> first << "}";
//            sum += iter -> first * iter -> second;
//        }
//        std::cout << "]"  << std::endl << "\t\tsum = " << sum << std::endl;
//
//        sum = 0;
//        std::cout << "\t- Target-in frequencies (frequency*{id}):" << std::endl << "\t\t[";
//        for (auto iter = v.targetInFrequencies.begin(); iter != v.targetInFrequencies.end(); iter++) {
//            if (iter != v.targetInFrequencies.begin()) std::cout << ", ";
//            std::cout << iter -> second << "*{" << iter -> first << "}";
//            sum += iter -> first * iter -> second;
//        }
//        std::cout << "]"  << std::endl << "\t\tsum = " << sum << std::endl;
//
//        std::cout << "\t- Distinct target nodes followed by an edge with the given label (distinct targets*{label}):" << std::endl << "\t\t[";
//        for (auto iter = v.distinctTargetNodesFollowedByLabel.begin(); iter != v.distinctTargetNodesFollowedByLabel.end(); iter++) {
//            if (iter != v.distinctTargetNodesFollowedByLabel.begin()) std::cout << ", ";
//            std::cout << iter -> second.size() << "*{" << iter -> first << "}";
//        }
//        std::cout << "]"  << std::endl;
//
//        std::cout << "\t- Number of edges with the given label following the target nodes (noEdges*{label}):" << std::endl << "\t\t[";
//        for (auto iter = v.noEdgesFollowingTargetNodesByLabel.begin(); iter != v.noEdgesFollowingTargetNodesByLabel.end(); iter++) {
//            if (iter != v.noEdgesFollowingTargetNodesByLabel.begin()) std::cout << ", ";
//            std::cout << iter -> second << "*{" << iter -> first << "}";
//        }
//        std::cout << "]"  << std::endl;
//
//
//        std::cout << "\t- Distinct source nodes followed by the current label, "
//                "followed by an edge with the given label (distinct sources*{label}):" << std::endl << "\t\t[";
//        for (auto iter = v.distinctSourceNodesFollowedByLabel.begin(); iter != v.distinctSourceNodesFollowedByLabel.end(); iter++) {
//            if (iter != v.distinctSourceNodesFollowedByLabel.begin()) std::cout << ", ";
//            std::cout << iter -> second.size() << "*{" << iter -> first << "}";
//        }
//        std::cout << "]"  << std::endl;
//
//        std::cout << "\t- Number of edges with the given label following the source nodes followed by the current label (noEdges*{label}):" << std::endl << "\t\t[";
//        for (auto iter = v.noEdgesFollowingSourceNodesByLabel.begin(); iter != v.noEdgesFollowingSourceNodesByLabel.end(); iter++) {
//            if (iter != v.noEdgesFollowingSourceNodesByLabel.begin()) std::cout << ", ";
//            std::cout << iter -> second << "*{" << iter -> first << "}";
//        }
//        std::cout << "]"  << std::endl << std::endl;
//    }
}



