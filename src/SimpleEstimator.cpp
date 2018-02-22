//
// Created by Nikolay Yakovets on 2018-02-01.
//

#include "SimpleGraph.h"
#include "SimpleEstimator.h"

SimpleEstimator::SimpleEstimator(std::shared_ptr<SimpleGraph> &g){

    // works only with SimpleGraph
    graph = g;

    std::cout << "complexity: vertices: " << g->getNoVertices() << ", edges: " << g->getNoEdges() << ", labels: " << g->getNoLabels() << std::endl;
    std::cout << std::endl;
}

void SimpleEstimator::prepare() {

    // do your prep here

    // Notes: good to know, in the pair, first is the label name and second is the target.
    distinctSourceVerticesPerLabel = std::vector<std::unordered_set<int>>(graph -> getNoLabels());
    distinctTargetVerticesPerLabel = std::vector<std::unordered_set<int>>(graph -> getNoLabels());
    totalEdgesPerLabel = std::vector<int>(graph -> getNoLabels());

    // Also note down the distinct vertices that concatenate two labels.
    distinctVerticesPerLabelPair = std::vector<std::vector<std::unordered_set<int>>>(graph -> getNoLabels(), std::vector<std::unordered_set<int>>(graph -> getNoLabels()));

    // Now, for each vertex, get the vertices it is connected to with the label.
    for(int i = 0; i < graph -> getNoVertices(); i++) {
        // Loop over all edges starting at node i.
        for(auto v : graph -> adj[i]) {
            // It starts in i, and ends in v.second with label v.first.
            distinctSourceVerticesPerLabel[v.first].insert(i);
            distinctTargetVerticesPerLabel[v.first].insert(v.second);
            totalEdgesPerLabel[v.first] += 1;

            // We also want to iterate over all edges ending in node i.
            for(auto w : graph -> reverse_adj[i]) {
                // so we know that i is at the end of w, and at the start of v. So add it to the correct set.
                distinctVerticesPerLabelPair[w.first][v.first].insert(i);
            }
        }
    }

    // Print the lengths of all the sets.
    for(int i = 0; i < graph -> getNoLabels(); i++) {
        std::cout << "Label " << i << " has "
                  << distinctSourceVerticesPerLabel[i].size() << " unique starting points." << std::endl;
        std::cout << "Label " << i << " has "
                  << distinctTargetVerticesPerLabel[i].size() << " unique end points." << std::endl;
        std::cout << "Label " << i << " is used by "
                  << totalEdgesPerLabel[i] << " edges." << std::endl;
        std::cout << std::endl;
    }

    // Now, print a table containing all relay data.
    std::cout << "Table containing the amount of vertices that may connect label 1 to label 2." << std::endl;
    for(int i = 0; i < graph -> getNoLabels(); i++) {
        for(int j = 0; j < graph -> getNoLabels(); j++) {
            std::cout << distinctVerticesPerLabelPair[i][j].size() << " \t\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

cardStat SimpleEstimator::estimate(RPQTree *q) {
    std::vector<std::pair<int, bool>> result = parseTreeToList(q);

    // Set the s and t cardinality, using the flag inwards and outwards data.
    auto s = result.front().second ?
                      distinctSourceVerticesPerLabel[result.front().first].size() :
                      distinctTargetVerticesPerLabel[result.front().first].size();
    auto t = !result.back().second ?
                      distinctSourceVerticesPerLabel[result.back().first].size() :
                      distinctTargetVerticesPerLabel[result.back().first].size();

    std::cout << std::endl << "Linear list representation of query: ";
    for (const auto &i: result)
        std::cout << i.first << "-" << i.second << ' ';
    std::cout << std::endl;

    // perform your estimation here
    return cardStat {static_cast<uint32_t>(s), 0, static_cast<uint32_t>(t)};
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

