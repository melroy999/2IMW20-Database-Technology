//
// Created by Nikolay Yakovets on 2018-02-01.
//

#include "SimpleGraph.h"
#include "SimpleEstimator.h"
#include <math.h>

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

//    for(uint32_t i = 0; i < 13; i++) {
//        std::cout << "id:" << i << " in:" << vertexData[i].inDegree << " out:" << vertexData[i].outDegree << std::endl;
//    }
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

    std::cout << std::endl;
    for(int i = 1; i < result.size(); i++) {

        auto l = result[i];
        noPaths *= ((float) labelData[l.first].noEdges) / (l.second ? labelData[l.first].distinctTargets : labelData[l.first].distinctSources).size();
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

