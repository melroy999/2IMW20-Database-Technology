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
    distinctTargetVerticesPerLabel = std::vector<std::unordered_set<int>>(graph -> getNoLabels());
    distinctSourceVerticesPerLabel = std::vector<std::unordered_set<int>>(graph -> getNoLabels());

    // Now, for each vertex, get the vertices it is connected to with the label.
    for(int i = 0; i < graph -> getNoVertices(); i++) {
        // Loop over all edges starting at node i.
        for(auto v : graph -> adj[i]) {
            // It starts in i, and ends in v.second with label v.first.
            distinctTargetVerticesPerLabel[v.first].insert(i);
            distinctSourceVerticesPerLabel[v.first].insert(v.second);
        }
    }

    // Print the lengths of all the sets.
    for(int i = 0; i < graph -> getNoLabels(); i++) {
        std::cout << "Label " << i << " has "
                  << distinctTargetVerticesPerLabel[i].size() << " unique starting points." << std::endl;
        std::cout << "Label " << i << " has "
                  << distinctSourceVerticesPerLabel[i].size() << " unique end points." << std::endl;
        std::cout << std::endl;
    }
}

cardStat SimpleEstimator::estimate(RPQTree *q) {
    std::vector<std::pair<int, bool>> result = parseTreeToList(q);

    // Set the s and t cardinality, using the flag inwards and outwards data.
    unsigned long s = result.front().second ?
                      distinctTargetVerticesPerLabel[result.front().first].size() :
                      distinctSourceVerticesPerLabel[result.front().first].size();
    unsigned long t = !result.back().second ?
                      distinctTargetVerticesPerLabel[result.back().first].size() :
                      distinctSourceVerticesPerLabel[result.back().first].size();

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

