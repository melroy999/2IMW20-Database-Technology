//
// Created by Nikolay Yakovets on 2018-02-01.
//

#include "SimpleGraph.h"
#include "SimpleEstimator.h"
#include <cmath>
#include <set>

SimpleEstimator::SimpleEstimator(std::shared_ptr<SimpleGraph> &g){

    // works only with SimpleGraph
    graph = g;
}

void SimpleEstimator::prepare() {

    // Initialise vectors with lengths corresponding to num vertices or labels.
    vertexData = std::vector<vertexStat>(graph -> getNoVertices());
    labelData = std::vector<labelStat>(graph -> getNoLabels());

    auto nodeTypes = std::vector<nodeType>(graph -> getNoVertices());

    // We loop over every resource (node) in the graph to obtain the type of every node

    for(uint32_t i = 0; i < graph -> getNoVertices(); i++) {

        for(auto v : graph -> adj[i]) {

            nodeTypes[i].outNodes.insert(v.second);
            nodeTypes[i].outPredicates.insert(v.first);
            nodeTypes[v.second].inPredicates.insert(v.first);
        }
    }

    // Merge all resources with the same type into buckets

    std::vector<std::set<uint32_t>> buckets;
    std::vector<nodeType> bucketDefinitions;



    for(uint32_t i = 0; i < graph -> getNoVertices(); i++) {

        for(nodeType classData : nodeTypes){

            
        }

        auto pred = [](const nodeType & item) {
            return item.outNodes == ;
        };

        ptrdiff_t pos = std::find(bucketDefinitions.begin(), bucketDefinitions.end(), nodeTypes[i]) - bucketDefinitions.begin();

        if(pos >= bucketDefinitions.size()) {
            //old_name_ not found
            pos = bucketDefinitions.size();
//            bucketDefinitions[bucketDefinitions.size()] = nodeTypes[i];
            buckets[pos] = {};
        }

        buckets[pos].insert({i});
    }


    /*
     * Per label, gather the set of distinct source and target vertices,
     * and the non-distinct count of source and target vertices.
     */
    for(uint32_t i = 0; i < graph -> getNoVertices(); i++) {

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
     * With the vertex data complete, we can process the remaining label data.
     * First, we want to gather the in and out degree of the target and source nodes respectively,
     * in the form of a frequency table.
     */
    for(auto v : vertexData) {

        /*
         * For each key-value pair in the in/out degree data, update the appropriate label frequency measurements.
         * In the degree mapping, 'first' is the label id, and 'second' is the degree.
         */
        for(auto w : v.labelOutDegrees) {
            labelData[w.first].sourceOutFrequencies[w.second] += 1;
        }
        for(auto w : v.labelInDegrees) {
            labelData[w.first].targetInFrequencies[w.second] += 1;
        }
    }

    /*
     * Gather the distinct target labels that are followed by an edge with a given label.
     * Next to that, also note down the number of edges with the given label originate from the target vertices.
     * Finally, we note down the distinct source nodes that,
     * after the application of this label, lead to edges with a given label.
     */
    for(uint32_t i = 0; i < graph -> getNoLabels(); i++) {
        for(auto v : labelData[i].distinctTargets) {
            if(vertexData[v].labelOutDegrees.empty()) {
                labelData[i].distinctTargetNodesFollowedByLabel[-1].insert(v);
            }
            for(auto w: vertexData[v].labelOutDegrees) {
                labelData[i].distinctTargetNodesFollowedByLabel[w.first].insert(v);
                labelData[i].noEdgesFollowingTargetNodesByLabel[w.first] += w.second;
            }
        }

        for(auto v : labelData[i].distinctSources) {
            for(auto w : graph -> adj[v]) {
                for(auto u: vertexData[w.second].labelOutDegrees) {
                    labelData[i].distinctSourceNodesFollowedByLabel[u.first].insert(v);
                    labelData[i].noEdgesFollowingSourceNodesByLabel[u.first] += u.second;
                }
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

        std::cout << "\t- Source-out frequencies (frequency*{id}):" << std::endl << "\t\t[";
        int sum = 0;
        for (auto iter = v.sourceOutFrequencies.begin(); iter != v.sourceOutFrequencies.end(); iter++) {
            if (iter != v.sourceOutFrequencies.begin()) std::cout << ", ";
            std::cout << iter -> second << "*{" << iter -> first << "}";
            sum += iter -> first * iter -> second;
        }
        std::cout << "]"  << std::endl << "\t\tsum = " << sum << std::endl;

        sum = 0;
        std::cout << "\t- Target-in frequencies (frequency*{id}):" << std::endl << "\t\t[";
        for (auto iter = v.targetInFrequencies.begin(); iter != v.targetInFrequencies.end(); iter++) {
            if (iter != v.targetInFrequencies.begin()) std::cout << ", ";
            std::cout << iter -> second << "*{" << iter -> first << "}";
            sum += iter -> first * iter -> second;
        }
        std::cout << "]"  << std::endl << "\t\tsum = " << sum << std::endl;

        std::cout << "\t- Distinct target nodes followed by an edge with the given label (distinct targets*{label}):" << std::endl << "\t\t[";
        for (auto iter = v.distinctTargetNodesFollowedByLabel.begin(); iter != v.distinctTargetNodesFollowedByLabel.end(); iter++) {
            if (iter != v.distinctTargetNodesFollowedByLabel.begin()) std::cout << ", ";
            std::cout << iter -> second.size() << "*{" << iter -> first << "}";
        }
        std::cout << "]"  << std::endl;

        std::cout << "\t- Number of edges with the given label following the target nodes (noEdges*{label}):" << std::endl << "\t\t[";
        for (auto iter = v.noEdgesFollowingTargetNodesByLabel.begin(); iter != v.noEdgesFollowingTargetNodesByLabel.end(); iter++) {
            if (iter != v.noEdgesFollowingTargetNodesByLabel.begin()) std::cout << ", ";
            std::cout << iter -> second << "*{" << iter -> first << "}";
        }
        std::cout << "]"  << std::endl;


        std::cout << "\t- Distinct source nodes followed by the current label, "
                "followed by an edge with the given label (distinct sources*{label}):" << std::endl << "\t\t[";
        for (auto iter = v.distinctSourceNodesFollowedByLabel.begin(); iter != v.distinctSourceNodesFollowedByLabel.end(); iter++) {
            if (iter != v.distinctSourceNodesFollowedByLabel.begin()) std::cout << ", ";
            std::cout << iter -> second.size() << "*{" << iter -> first << "}";
        }
        std::cout << "]"  << std::endl;

        std::cout << "\t- Number of edges with the given label following the source nodes followed by the current label (noEdges*{label}):" << std::endl << "\t\t[";
        for (auto iter = v.noEdgesFollowingSourceNodesByLabel.begin(); iter != v.noEdgesFollowingSourceNodesByLabel.end(); iter++) {
            if (iter != v.noEdgesFollowingSourceNodesByLabel.begin()) std::cout << ", ";
            std::cout << iter -> second << "*{" << iter -> first << "}";
        }
        std::cout << "]"  << std::endl << std::endl;
    }
}

