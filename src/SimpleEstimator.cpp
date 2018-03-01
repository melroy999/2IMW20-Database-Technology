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
        }
    }

    // TODO there currently is a bug where the number of nodes following a target node is double of what it should be.

//    for(uint32_t i = 0; i < graph -> getNoLabels(); i++) {
//        std::cout << labelData[{i, true}].getDistinctSources().size() << "," << labelData[{i, false}].getDistinctTargets().size() << std::endl;
//        std::cout << labelData[{i, true}].getDistinctTargets().size() << "," << labelData[{i, false}].getDistinctSources().size() << std::endl;
//    }

    // Print the debug data.
    printDebugData();
}



cardStat SimpleEstimator::estimate(RPQTree *q) {

    // Convert the parse tree to a more straightforward form to work with i.e. a list.
    std::vector<std::pair<uint32_t, bool>> result = parseTreeToList(q);

    /*
     * We start by estimating s and t by getting the number of distinct sources/targets
     * of the first and last label in the chain respectively.
     */
    auto s = static_cast<uint32_t>(labelData[result.front()].getDistinctSources().size());
    auto t = static_cast<uint32_t>(labelData[result.back()].getDistinctTargets().size());

    // To estimate the number of paths, use the average degree of source vertices bearing the label.
    float noPaths = labelData[result.front()].getNoEdges();

    for(int i = 1; i < result.size(); i++) {

        std::pair<uint32_t, bool> l = result[i];
        std::pair<uint32_t, bool> l_prev = result[i - 1];

        // Simple version:
//        noPaths *= ((float) labelData[l].getNoEdges()) / labelData[l].getDistinctTargets().size();

        // Using previous vertex data:
        noPaths = ceilf(noPaths * ((float) labelData[l_prev].getNoEdgesFollowingTargetNodesByLabel(l)) / labelData[l].getDistinctTargets().size());


        /*
         * Suppose that S are our source nodes, and T are our target nodes for an arbitrary label (l, b)
         *
         * Given the out degrees of the nodes in S, and the in degrees of the nodes in T:
         *      - Can we efficiently estimate the number of edges between S and T?
         *          Let d_out_avg(S) be the average out degree of the source nodes:
         *              sum_{v in S} out-degree(v) / |S| = #edges / |S|
         *          Let d_in_avg(T) be the average in degree of the target nodes:
         *              sum_{v in T} in-degree(v) / |T| = #edges / |T|
         *
         *        Suppose that we are now only interested in the edges originating from a subset of S: S_sub.
         *        Obviously, we could estimate the number of edges as d_out_avg(S_sub) = |S_sub| * d_avg_out(S).
         *        However, the out-degree between vertices might differ greatly, with very large outliers.
         *        Suppose that all these outliers are in S_sub, then d_out_avg(S_sub) >> d_out_avg(S).
         *        Suppose that none of these outliers are in S_sub, then d_out_avg(S_sub) << d_out_avg(S).
         *
         *        So can we estimate the number of edges for S_sub more accurately, keeping in mind that we have outliers?
         */


        /*
         * The increase of path depends on the following factors.
         * Let Q be the query, and (l, b) in Q is the current label we are observing.
         *      - What are the number of source vertices of (l, b) that are not proceeded by label prev(Q, (l, b))?
         *          These node are dangling, and should not be considered part of the path.
         *
         *          How do we remove these dangling paths from the total path count?
         *          First option is to not include them as paths initially.
         *          I.e., valid_edges = total_edges - sum_{v in dangling_vertices} degree(v, (l, b))
         *
         *          How do we estimate: sum_{v in dangling_vertices} degree(v, (l, b))?
         *          First of all, we could estimate it as: average_degree * |dangling_vertices|
         *          However, the average degree of the dangling_vertices might be different for a non random data set.
         *          We could pre-calculate sum_{v in dangling_vertices} degree(v, (l, b)) in our data structure.         *
         *
         *      - What are the number of target vertices of next(Q, (l, b)) that are not succeeded by label (l, b)?
         *          These nodes have a dead end, and should not be considered part of the path.
         */



        // This one does not work out, as the degree is applied to ALL paths, while some have no follow ups.
//        noPaths *= ((float) labelData[l_prev].getNoEdgesFollowingTargetNodesByLabel(l)) / labelData[l_prev].getDistinctTargetNodesFollowedByLabel(l).size();



        // Of the source nodes of the current label, how many do connect to the previous label?
//        float f1 = (float) labelData[l].getDistinctSourceNodesProceededByLabel(l_prev).size() /
//                labelData[l].getDistinctSources().size();

        // Of the target nodes of the previous label, how many do connect to the current label?
//        float f2 = (float) labelData[l_prev].getDistinctTargetNodesFollowedByLabel(l).size() /
//                labelData[l_prev].getDistinctSources().size();

        // How many of these two measurements overlap?


//        noPaths *= f1 * ((float) labelData[l].getNoEdges()) / labelData[l].getDistinctTargets().size();


//        noPaths *= f1 * ((float) labelData[l_prev].getNoEdgesFollowingTargetNodesByLabel(l)) / labelData[l].getDistinctTargetNodesFollowedByLabel(l).size();

//        noPaths *= ((float) labelData[l_prev].getNoEdgesFollowingTargetNodesByLabel()[l].getNoEdges()) / labelData[l].getDistinctTargets().size();
    }

    // Return the estimate in the form {#outNodes, #paths, #inNodes}
    return cardStat {s, static_cast<uint32_t>(noPaths), t};
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



