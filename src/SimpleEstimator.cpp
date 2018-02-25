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
    // Initialize the data structures containing the desired data.
    vertexData = std::vector<vertexStat>(graph -> getNoVertices());
    labelData = std::vector<labelStat>(graph -> getNoLabels());

    // We want to gather the in and out degrees of the vertices in the graph.
    // Next to that, per label, we want to gather the set of distinct source and target vertices,
    // and the non-distinct count of source and target vertices.
    for(uint32_t i = 0; i < graph -> getNoVertices(); i++) {
        // The in and out degree correspond with the lengths of the results given in the adjacency matrices.
        vertexData[i].inDegree = static_cast<uint32_t>(graph -> adj[i].size());
        vertexData[i].outDegree = static_cast<uint32_t>(graph -> reverse_adj[i].size());

        // Now use the entries in the adjacency matrix to calculate the label data.
        for(auto v : graph -> adj[i]) {
            labelData[v.first].noEdges++;
            labelData[v.first].distinctSources.insert(i);
            labelData[v.first].distinctTargets.insert(v.second);

            // Also update the appropriate in and out degree counters in the vertex data.
            ++vertexData[i].labelOutDegrees[v.first];
            ++vertexData[v.second].labelInDegrees[v.first];
        }
    }
}

cardStat SimpleEstimator::estimate(RPQTree *q) {
    // Convert the parse tree to a more straightforward form to work with, a list.
    std::vector<std::pair<int, bool>> result = parseTreeToList(q);

    // We have to keep in mind that we can also move in opposite direction.

    // Start with estimating s and t.
    auto s = (result.front().second ? labelData[result.front().first].distinctSources : labelData[result.front().first].distinctTargets).size();
    auto t = (result.back().second ? labelData[result.back().first].distinctTargets : labelData[result.back().first].distinctSources).size();

    // To estimate the number of paths, use the average degree of source vertices bearing the label.
    float noPaths = labelData[result[0].first].noEdges;

    std::cout << std::endl;
    for(int i = 1; i < result.size(); i++) {
        auto l = result[i];
        noPaths *= ((float) labelData[l.first].noEdges) / (l.second ? labelData[l.first].distinctTargets : labelData[l.first].distinctSources).size();
    }

    // perform your estimation here
    return cardStat {static_cast<uint32_t>(s), static_cast<uint32_t>(noPaths), static_cast<uint32_t>(t)};
}

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

