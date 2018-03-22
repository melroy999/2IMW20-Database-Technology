//
// Created by Nikolay Yakovets on 2018-02-02.
//

#include "SimpleEstimator.h"
#include "SimpleEvaluator.h"

#include <regex>

SimpleEvaluator::SimpleEvaluator(std::shared_ptr<SimpleGraph> &g) {

    // works only with SimpleGraph
    graph = g;
    est = nullptr; // estimator not attached by default
}

void SimpleEvaluator::attachEstimator(std::shared_ptr<SimpleEstimator> &e) {
    est = e;
}

void SimpleEvaluator::prepare() {

    // if attached, prepare the estimator
    if(est != nullptr) est->prepare();

    // prepare other things here.., if necessary

}

cardStat SimpleEvaluator::computeStats(std::shared_ptr<SimpleGraph> &g) {

    cardStat stats {};

    auto _adj = g->adj_ptr ? g->adj_ptr: &g->adj[0];
    auto _reverse_adj = g->reverse_adj_ptr ? g->reverse_adj_ptr: &g->reverse_adj[0];

    for(int source = 0; source < g->getNoVertices(); source++) {
        if(!(*_adj)[source].empty()) stats.noOut++;
    }

    stats.noPaths = g->getNoDistinctEdges();

    for(int target = 0; target < g->getNoVertices(); target++) {
        if(!(*_reverse_adj)[target].empty()) stats.noIn++;
    }

    return stats;
}

std::shared_ptr<SimpleGraph> SimpleEvaluator::project(uint32_t projectLabel, bool inverse, std::shared_ptr<SimpleGraph> &in) {

    auto out = std::make_shared<SimpleGraph>(in->getNoVertices());

    // Why loop here over all the data, while we can just extract the data in one sweep?
    out->addEdges(in, projectLabel, inverse);

    return out;
}

std::shared_ptr<SimpleGraph> SimpleEvaluator::join(std::shared_ptr<SimpleGraph> &left, std::shared_ptr<SimpleGraph> &right) {

    auto out = std::make_shared<SimpleGraph>(left->getNoVertices());
    out->setNoLabels(1);
    out->setDataStructureSizes();

    auto leftMatrix = left->adj_ptr ? left->adj_ptr: &left->adj[0];
    auto rightMatrix = right->adj_ptr ? right->adj_ptr: &right->adj[0];

    for(uint32_t leftSource = 0; leftSource < left->getNoVertices(); leftSource++) {
        for (auto leftTarget : (*leftMatrix)[leftSource]) {

            // try to join the left target with right source
            for (auto rightTarget : (*rightMatrix)[leftTarget]) {

                out->addEdge(leftSource, rightTarget, 0);

            }
        }
    }

    return out;
}

std::shared_ptr<SimpleGraph> SimpleEvaluator::evaluate_aux(RPQTree *q) {

    // evaluate according to the AST bottom-up

    if(q->isLeaf()) {
        // project out the label in the AST
        std::regex directLabel (R"((\d+)\+)");
        std::regex inverseLabel (R"((\d+)\-)");

        std::smatch matches;

        uint32_t label;
        bool inverse;

        if(std::regex_search(q->data, matches, directLabel)) {
            label = (uint32_t) std::stoul(matches[1]);
            inverse = false;
        } else if(std::regex_search(q->data, matches, inverseLabel)) {
            label = (uint32_t) std::stoul(matches[1]);
            inverse = true;
        } else {
            std::cerr << "Label parsing failed!" << std::endl;
            return nullptr;
        }

        return SimpleEvaluator::project(label, inverse, graph);
    }

    if(q->isConcat()) {

        // evaluate the children
        auto leftGraph = SimpleEvaluator::evaluate_aux(q->left);
        auto rightGraph = SimpleEvaluator::evaluate_aux(q->right);

        // join left with right
        return SimpleEvaluator::join(leftGraph, rightGraph);

    }

    return nullptr;
}

/**
 * Get the normalised form of a query, as a vector of nodes
 * @param q
 * @return
 */
std::vector<std::string> SimpleEvaluator::get_query_graph(RPQTree *q, std::vector<std::string> nodes) {

    // evaluate according to the AST bottom-up

    if(q->isLeaf()) {

        nodes.emplace_back(q->data);
        return nodes;
    }

    if(q->isConcat()) {

        auto leftNodes = SimpleEvaluator::get_query_graph(q->left, nodes);
        auto rightNodes = SimpleEvaluator::get_query_graph(q->right, nodes);
        leftNodes.insert(leftNodes.end(), rightNodes.begin(), rightNodes.end());
        return leftNodes;
    }

    return nodes;
}

/**
 * Rewrite a given query tree to an optimal (according to estimates) variant
 * @param query
 * @return
 */
RPQTree * SimpleEvaluator::rewrite_query_tree(RPQTree *query) {

    // First, normalise the query into a vector of nodes
    std::vector<std::string> nodes;
    std::vector<std::string> queryN = get_query_graph(query, nodes);

    // Estimate the cardinalities for each pair
    while(queryN.size() > 2){

        std::vector<uint32_t> estimates;

        for(int i = 0; i < queryN.size() - 1; i++){

            std::string queryPair = "(" + queryN[i] + "/" + queryN[i+1] + ")";

            // Check if this query pair exists in the estimator cache

            if(estCache.find(queryPair) != estCache.end()){

                estimates.push_back(estCache[queryPair]);
            }
            else{

                auto queryTree = RPQTree::strToTree(queryPair);
                cardStat estimate = est->estimate(queryTree);
                estimates.push_back(estimate.noPaths);
                estCache.insert(std::make_pair(queryPair, estimate.noPaths));
            }
        }

        long min_index = std::min_element(estimates.begin(), estimates.end()) - estimates.begin();
        std::string minPair = "(" + queryN[min_index] + "/" + queryN[min_index + 1] + ")";
        queryN[min_index] = minPair;
        queryN.erase(queryN.begin()+min_index+1);
    }

    if(queryN.size() == 2){

        // Now concat the remaining two labels to produce an optimal tree which can be evaluated

        std::string finalQuery = "(" + queryN[0] + "/" + queryN[1] + ")";
        query = RPQTree::strToTree(finalQuery);
    }

    return query;
}

cardStat SimpleEvaluator::evaluate(RPQTree *query) {

    if(est){

        query = rewrite_query_tree(query);
    }

    auto res = evaluate_aux(query);
    return SimpleEvaluator::computeStats(res);
}