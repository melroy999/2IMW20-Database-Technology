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
        std::unordered_set<std::pair<uint32_t, uint32_t>, pairHasher> encounteredEdges;

        // Now use the entries in the adjacency matrix to calculate the label data.
        for(std::pair<uint32_t,uint32_t> v : graph -> adj[i]) {

            // If the edge has been added already previously, do not add it again as we want distinct paths.
            if(!encounteredEdges.insert(v).second) continue;

            // Update the appropriate out degree counter in the label degree mapping.
            ++vertexData[i].labelDegrees[{v.first, true}];

        }

        // Do the same for the reverse adjacency matrix.
        for(auto v : graph -> reverse_adj[i]) {

            // If the edge has been added already previously, do not add it again as we want distinct paths.
            if(!encounteredEdges.insert(v).second) continue;

            // Update the appropriate in degree counter in the label degree mapping.
            ++vertexData[i].labelDegrees[{v.first, false}];
        }
    }

    /*
     * Start gathering data for all the positive labels.
     *
     * First, we want to gather the distinct source and target vertices per label.
     *
     * Next to that, we want to gather the in and out degree of the target and source nodes respectively,
     * in the form of a frequency table.
     */

    for(uint32_t i = 0; i < graph -> getNoVertices(); i++) {

        for(auto pair : vertexData[i].labelDegrees) {
            /*
             * In the pair:
             *  - first.first corresponds to the label id.
             *  - first.second corresponds to the label direction, true is forwards, false is backwards.
             *  - second corresponds to the degree.
             */
            if(pair.first.second) {
                labelData[{pair.first.first, true}].addEdges(static_cast<uint32_t>(pair.second));
                labelData[{pair.first.first, true}].addSourceVertex(i);
                labelData[{pair.first.first, true}].updateSourceOutFrequency(static_cast<uint32_t>(pair.second), i);
            } else {
                labelData[{pair.first.first, true}].addTargetVertex(i);
                labelData[{pair.first.first, true}].updateTargetInFrequency(static_cast<uint32_t>(pair.second), i);
            }
        }
    }

    for(uint32_t i = 0; i < graph -> getNoLabels(); i++) {
        for (bool b : {true, false}) {
            for (auto v : labelData[{i, b}].getDistinctTargets()) {
                for (auto w: vertexData[v].labelDegrees) {
                    labelData[{i, b}].incrementDistinctTargetsInJoin(w.first);
                    labelData[{i, b}].updateNumberOfFollowUpEdgesInJoin(w.first, static_cast<uint32_t>(w.second));
                }
            }

            labelData[{i, b}].calculateSquareDegreeEstimation();
            labelData[{i, b}].analyzeFrequencyData();
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
    std::cout << std::endl;

    // Convert the parse tree to a more straightforward form to work with i.e. a list.
    std::vector<std::pair<uint32_t, bool>> result = parseTreeToList(q);

    /*
     * We start by estimating s and t by getting the number of distinct sources/targets
     * of the first and last label in the chain respectively.
     */
    auto s = static_cast<uint32_t>(labelData[result.front()].getNumberOfDistinctSources());
    auto t = static_cast<uint32_t>(labelData[result.back()].getNumberOfDistinctTargets());

    // To estimate the number of paths, use the average degree of source vertices bearing the label.
    uint32_t noPaths = labelData[result.front()].getNumberOfEdges();

    // The maximum termination factor.


    for(int i = 1; i < result.size(); i++) {

        // For our estimations, we always require the data of the current and previous label.
        auto l = result[i];
        auto l_prev = result[i - 1];

        // First, get the number of target vertices of l_prev that are connected to the source vertices of l.
        long terminatedTargetNodes = labelData[l_prev].getNumberOfDistinctTargets() -
                                     labelData[l_prev].getDistinctTargetsInJoin(l);

        // Calculate the number of edges we would have to remove.
        float terminationFactor = (float) terminatedTargetNodes / labelData[l_prev].getNumberOfEdges();

        // Remove that factor of paths.
        noPaths = static_cast<uint32_t>(ceilf(noPaths * (1 - terminationFactor)));


        if(l.first == l_prev.first && l.second != l_prev.second) {
            noPaths *= (float) labelData[l].getSquareDegreeEstimation() / labelData[l].getNumberOfEdges();
        } else {
            noPaths *= (float) labelData[l_prev].getNumberOfFollowUpEdgesInJoin(l) /
                       labelData[l_prev].getDistinctTargetsInJoin(l);
        }


//        noPaths *= (float) labelData[l_prev].getNumberOfPairsInJoin(l) /
//                (labelData[l_prev].getNumberOfEdges());




        // We also have edges that just started, and are not connected to any of the previous.
        long instantiatedSourceNodes = labelData[l].getNumberOfDistinctSources() -
                                       labelData[l].getDistinctSourcesInJoin(l_prev);

        // Calculate the number of edges we would have to remove.
        float instantiationFactor = (float) instantiatedSourceNodes / labelData[l].getNumberOfEdges();

        // Remove that factor of paths.
        noPaths = static_cast<uint32_t>(ceilf(noPaths * (1 - instantiationFactor)));
    }

    // Return the estimate in the form {#outNodes, #paths, #inNodes}
    return cardStat {s, noPaths, t};
}

/**
 * Convert the RPQTree to a vector representation of pairs of labels and directions
 * @param q - the tree
 * @return the vector list of labels
 */
std::vector<std::pair<uint32_t, bool>> SimpleEstimator::parseTreeToList(RPQTree *q) {

    // Create a vector that will contain the result.
    std::vector<std::pair<uint32_t, bool>> result;

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
        for(bool b : {true, false}) {
            auto v = labelData[{i, b}];

            std::cout << "Label {" << i << ", " << (b ? "true" : "false") << "}:" << std::endl;
            v.printData();
            std::cout << std::endl;
        }
    }
}



