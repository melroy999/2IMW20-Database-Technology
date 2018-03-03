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

struct pairHasher {
    inline std::size_t operator()(const std::pair<int,int> & v) const {
        return static_cast<size_t>(v.first * 31 + v.second);
    }
};

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
                labelData[{pair.first.first, true}].addEdges(pair.second);
                labelData[{pair.first.first, true}].addSourceVertex(i);
                labelData[{pair.first.first, true}].updateSourceOutFrequency(pair.second, i);
            } else {
                labelData[{pair.first.first, true}].addTargetVertex(i);
                labelData[{pair.first.first, true}].updateTargetInFrequency(pair.second, i);
            }
        }
    }

    /*
     * Gather the distinct target labels that are followed by an edge with a given label.
     * Next to that, also note down the number of edges with the given label originate from the target vertices.
     */
    for(uint32_t i = 0; i < graph -> getNoLabels(); i++) {
        for (bool b : {true, false}) {
            for (auto v : labelData[{i, b}].getDistinctTargets()) {
                for (auto w: vertexData[v].labelDegrees) {
                    labelData[{i, b}].updateDistinctTargetNodesFollowedByLabel(w.first, v);
                    labelData[{i, b}].updateNoEdgesFollowingTargetNodesByLabel(w.first, w.second);
                }
            }

            // Calculate the target frequency square metric for all the labels.
            labelData[{i, b}].calculateSumOfTargetFrequencySquares();
        }
    }

//    for(uint32_t i = 0; i < graph -> getNoLabels(); i++) {
//        std::cout << labelData[{i, true}].getDistinctSources().size() << "," << labelData[{i, false}].getDistinctTargets().size() << std::endl;
//        std::cout << labelData[{i, true}].getDistinctTargets().size() << "," << labelData[{i, false}].getDistinctSources().size() << std::endl;
//    }

    // Print the debug data.
//    printDebugData();
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
    uint32_t noPaths = labelData[result.front()].getNoEdges();

    for(int i = 1; i < result.size(); i++) {

        auto l = result[i];
        auto l_prev = result[i - 1];

        // If the labels are the same, and the direction is opposite, we should use the degree squaring approach.
        if(l.first == l_prev.first && l.second != l_prev.second) {
            uint32_t squareSum = labelData[l_prev].getSumOfTargetFrequencySquares();

            float degree_out = (float) labelData[l].getNoEdges() / labelData[l].getNumberOfDistinctSources();
            float degree_in = (float) labelData[l].getNoEdges() / labelData[l].getNumberOfDistinctTargets();

//            std::cout << "square sum: " << squareSum << std::endl;
//            std::cout << "degree out: " << degree_out << std::endl;
//            std::cout << "degree in: " << degree_in << std::endl << std::endl;

            noPaths *= ((float) squareSum / labelData[l].getNoEdges());
        } else {
            // First, get the number of target vertices of l_prev that are connected to the source vertices of l.
            long terminatedTargetNodes = labelData[l_prev].getNumberOfDistinctTargets() -
                                         labelData[l_prev].getNumberOfDistinctTargetNodesFollowedByLabel(l);

            // Calculate the number of edges we would have to remove.
            float terminationFactor = (float) terminatedTargetNodes / labelData[l_prev].getNoEdges();

            // Remove that factor of paths.
            noPaths = static_cast<uint32_t>(ceilf(noPaths * (1 - terminationFactor)));

            // Simple version:
            noPaths *= (float) labelData[l_prev].getNoEdgesFollowingTargetNodesByLabel(l) /
                       labelData[l_prev].getNumberOfDistinctTargetNodesFollowedByLabel(l);



            // We also have edges that just started, and are not connected to any of the previous.
            long instantiatedSourceNodes = labelData[l].getNumberOfDistinctSources() -
                                           labelData[l].getNumberOfDistinctSourceNodesProceededByLabel(l_prev);

            // Calculate the number of edges we would have to remove.
            float instantiationFactor = (float) instantiatedSourceNodes / labelData[l].getNoEdges();

            // Remove that factor of paths.
            noPaths = static_cast<uint32_t>(ceilf(noPaths * (1 - instantiationFactor)));

        }



        // Using previous vertex data:
//        noPaths = static_cast<uint32_t>(ceilf(noPaths * ((float) labelData[l_prev].getNoEdgesFollowingTargetNodesByLabel(l)) / labelData[l].getDistinctTargets().size()));
//
//        float v = ((float) labelData[l_prev].getNoEdgesFollowingTargetNodesByLabel(l)) / labelData[l].getDistinctTargets().size();
////
////        // The above, with a modifier that puts emphasis on the node with the largest degree.
//        noPaths = static_cast<uint32_t>(ceil(noPaths * (v + (labelData[l].getSourceOutFrequencies().rbegin() -> first - ((float) labelData[l].getNoEdges() / labelData[l].getDistinctTargets().size())) / 16.0)));


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
            std::cout << "\t- is used by " << v.getNoEdges() << " edges." << std::endl;
            std::cout << "\t- has " << v.getDistinctSources().size() << " distinct sources." << std::endl;
            std::cout << "\t- has " << v.getDistinctTargets().size() << " distinct targets." << std::endl;

            std::cout << "\t- Source-out frequencies (frequency*{id}):" << std::endl << "\t\t[";
            int sum = 0;
            for (auto iter = v.getSourceOutFrequencies().begin(); iter != v.getSourceOutFrequencies().end(); iter++) {
                if (iter != v.getSourceOutFrequencies().begin()) std::cout << ", ";
                std::cout << iter -> second.size() << "*{" << iter -> first << "}";
                sum += iter -> first * iter -> second.size();
            }
            if(sum != v.getNoEdges()) {
                std::cout << "]"  << std::endl << "\t\tWARNING: sum != " << sum << std::endl;
            } else {
                std::cout << "]"  << std::endl;
            }

            sum = 0;
            std::cout << "\t- Target-in frequencies (frequency*{id}):" << std::endl << "\t\t[";
            for (auto iter = v.getTargetInFrequencies().begin(); iter != v.getTargetInFrequencies().end(); iter++) {
                if (iter != v.getTargetInFrequencies().begin()) std::cout << ", ";
                std::cout << iter -> second.size() << "*{" << iter -> first << "}";
                sum += iter -> first * iter -> second.size();
            }
            if(sum != v.getNoEdges()) {
                std::cout << "]"  << std::endl << "\t\tWARNING: sum != " << sum << std::endl;
            } else {
                std::cout << "]"  << std::endl;
            }

            std::cout << "\t- Distinct target nodes followed by an edge with the given label (distinct targets*{label}):" << std::endl << "\t\t[";
            for (auto iter = v.getDistinctTargetNodesFollowedByLabel().begin(); iter != v.getDistinctTargetNodesFollowedByLabel().end(); iter++) {
                if (iter != v.getDistinctTargetNodesFollowedByLabel().begin()) std::cout << ", ";
                std::cout << iter -> second.size() << "*{" << iter -> first.first << ", " << (iter -> first.second ? "true" : "false") << "}";
            }
            std::cout << "]"  << std::endl;

            std::cout << "\t- Number of edges with the given label following the target nodes (noEdges*{label}):" << std::endl << "\t\t[";
            for (auto iter = v.getNoEdgesFollowingTargetNodesByLabel().begin(); iter != v.getNoEdgesFollowingTargetNodesByLabel().end(); iter++) {
                if (iter != v.getNoEdgesFollowingTargetNodesByLabel().begin()) std::cout << ", ";
                std::cout << iter -> second << "*{" << iter -> first.first << ", " << (iter -> first.second ? "true" : "false") << "}";
            }
            std::cout << "]"  << std::endl;
//
//
//            std::cout << "\t- Distinct source nodes followed by the current label, "
//                    "followed by an edge with the given label (distinct sources*{label}):" << std::endl << "\t\t[";
//            for (auto iter = v.distinctSourceNodesFollowedByLabel.begin(); iter != v.distinctSourceNodesFollowedByLabel.end(); iter++) {
//                if (iter != v.distinctSourceNodesFollowedByLabel.begin()) std::cout << ", ";
//                std::cout << iter -> second.size() << "*{" << iter -> first << "}";
//            }
//            std::cout << "]"  << std::endl;
//
//            std::cout << "\t- Number of edges with the given label following the source nodes followed by the current label (noEdges*{label}):" << std::endl << "\t\t[";
//            for (auto iter = v.noEdgesFollowingSourceNodesByLabel.begin(); iter != v.noEdgesFollowingSourceNodesByLabel.end(); iter++) {
//                if (iter != v.noEdgesFollowingSourceNodesByLabel.begin()) std::cout << ", ";
//                std::cout << iter -> second << "*{" << iter -> first << "}";
//            }

            std::cout << std::endl;
        }
    }
}



