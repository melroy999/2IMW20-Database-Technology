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
    labelData = std::vector<labelStat>(graph -> getNoLabels());

    /*
     * Gather the in and out degrees of the vertices in the graph.
     * Next to that, per label, gather the set of distinct source and target vertices,
     * and the non-distinct count of source and target vertices.
     */
    for(uint32_t i = 0; i < graph -> getNoVertices(); i++) {

        // The in and out degree corresponds to the lengths of the results given in the adjacency matrices.
        vertexData[i].outDegree = static_cast<uint32_t>(graph -> adj[i].size());
        vertexData[i].inDegree = static_cast<uint32_t>(graph -> reverse_adj[i].size());

        // Now use the entries in the adjacency matrix to calculate the label data.
        for(auto v : graph -> adj[i]) {

            // In the adjacency matrix, 'first' is the edge label, 'second' is the target vertex.
            labelData[v.first].noEdges++;
            labelData[v.first].distinctSources.insert(i);
            labelData[v.first].distinctTargets.insert(v.second);

            // Also update the appropriate in and out degree counters in the vertex data.
            ++vertexData[i].labelOutDegrees[v.first];
            ++vertexData[v.second].labelInDegrees[v.first];
        }
    }

    /*
     * For each of the labels in the dataset, count the number of edges that succeed or proceed edges using the label,
     * subdivided by the label used by the succeeding or proceeding edges.
     */
    for(uint32_t i = 0; i < graph -> getNoLabels(); i++) {

        // For each of the target vertices, count the number of outgoing edges per label type.
        for(auto v : labelData[i].distinctTargets) {
            for(auto w : vertexData[v].labelOutDegrees) {
                labelData[i].noSuccessorEdgesPerLabel[w.first] += w.second;
                labelData[i].distinctSourcesPerSuccessorLabel[w.first].insert(v);
            }
        }

        // For each of the source vertices, count the number of ingoing edges per label type.
        for(auto v : labelData[i].distinctSources) {
            for(auto w : vertexData[v].labelInDegrees) {
                labelData[i].noPredecessorEdgesPerLabel[w.first] += w.second;
                labelData[i].distinctTargetsPerPredecessorLabel[w.first].insert(v);
            }
        }
    }

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
    auto s = (result.front().second ? labelData[result.front().first].distinctSources : labelData[result.front().first].distinctTargets).size();
    auto t = (result.back().second ? labelData[result.back().first].distinctTargets : labelData[result.back().first].distinctSources).size();

    // To estimate the number of paths, use the average degree of source vertices bearing the label.
    float noPaths = labelData[result[0].first].noEdges;

    for(int i = 1; i < result.size(); i++) {

        auto l = result[i];
        noPaths *= ((float) labelData[l.first].noEdges) / (l.second ? labelData[l.first].distinctTargets : labelData[l.first].distinctSources).size();

//        if(l.second) {
//            /*
//             * We traverse the edge in expected direction, so we should use the successor value.
//             *
//             * Instead of using the rough distinct targets of the label,
//             * we check how many successor edges the predecessor label has with the desired label.
//             */
//            noPaths *= ((float) labelData[result[i-1].first].noSuccessorEdgesPerLabel[l.first]) / labelData[result[i-1].first].distinctSourcesPerSuccessorLabel[l.first].size();
//
//        } else {
//            /*
//             * We traverse the edge in opposite direction, so we should use the predecessor value.
//             *
//             * Instead of using the rough distinct sources of the label,
//             * we check how many predecessor edges the successor label has with the desired label.
//             */
//            noPaths *= ((float) labelData[result[i-1].first].noPredecessorEdgesPerLabel[l.first]) / labelData[result[i-1].first].distinctTargetsPerPredecessorLabel[l.first].size();
//        }
    }

    // Return the estimate in the form {#outNodes, #paths, #inNodes}
    return cardStat {static_cast<uint32_t>(s), static_cast<uint32_t>(noPaths), static_cast<uint32_t>(t)};
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
    for(uint32_t i = 0; i < graph -> getNoLabels(); i++) {
        auto v = labelData[i];
        std::cout << "Label " << i << ":" << std::endl;
        std::cout << "\t- is used by " << v.noEdges << " edges." << std::endl;
        std::cout << "\t- has " << v.distinctSources.size() << " distinct sources." << std::endl;
        std::cout << "\t- has " << v.distinctTargets.size() << " distinct targets." << std::endl;

        std::cout << "\t- is succeeded by:" << std::endl;
        for(auto w : v.noSuccessorEdgesPerLabel) {
            std::cout << "\t\t- " << w.second << " edges with label " << w.first << ", originating from "
                      << v.distinctSourcesPerSuccessorLabel[w.first].size() << " distinct vertices, with "
                      << labelData[w.first].distinctSources.size() << " available vertices." << std::endl;
        }

        std::cout << "\t- is proceeded by:" << std::endl;
        for(auto w : v.noPredecessorEdgesPerLabel) {
            std::cout << "\t\t- " << w.second << " edges with label " << w.first << ", terminating in "
                      << v.distinctTargetsPerPredecessorLabel[w.first].size() << " distinct vertices, with "
                      << labelData[w.first].distinctTargets.size() << " available vertices." << std::endl;
        }
        std::cout << std::endl;
    }
}

