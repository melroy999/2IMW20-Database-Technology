//
// Created by Nikolay Yakovets on 2018-02-02.
//

#include "SimpleEstimator.h"
#include "SimpleEvaluator.h"

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

cardStat SimpleEvaluator::computeStats(std::shared_ptr<JoinGraph> &g) {

    cardStat stats {};

    stats.noOut = g->noSources;
    stats.noPaths = g->noEdges;
    stats.noIn = g->noTargets;

    return stats;
}

std::shared_ptr<JoinGraph> SimpleEvaluator::project(uint32_t projectLabel, bool inverse, std::shared_ptr<SimpleGraph> &in) {

    // Why loop here over all the data, while we can just extract the data in one sweep?
    return std::make_shared<JoinGraph>(in->graphs[projectLabel], inverse);
}

std::shared_ptr<JoinGraph> SimpleEvaluator::join(std::shared_ptr<JoinGraph> &left, std::shared_ptr<JoinGraph> &right) {

    return std::make_shared<JoinGraph>(*left.get(), *right.get());
}

// project out the label in the AST
std::regex directLabel (R"((\d+)\+)");
std::regex inverseLabel (R"((\d+)\-)");

std::shared_ptr<JoinGraph> SimpleEvaluator::evaluate_aux(RPQTree *q) {

    // evaluate according to the AST bottom-up

    if(q->isLeaf()) {
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

        std::string queryString = get_query_as_string(q);

        #ifdef RESULTS_CACHE
        if(resultsCache.find(queryString) != resultsCache.end()){

            return resultsCache[queryString];
        }
        #endif

        // evaluate the children
        auto leftGraph = SimpleEvaluator::evaluate_aux(q->left);
        auto rightGraph = SimpleEvaluator::evaluate_aux(q->right);

        // join left with right

        std::shared_ptr<JoinGraph> result = SimpleEvaluator::join(leftGraph, rightGraph);

        #ifdef RESULTS_CACHE
        resultsCache.insert(std::make_pair(queryString, result));
        #endif

        return result;
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
 * Get a query tree as a string
 * @param q
 * @return
 */
std::string SimpleEvaluator::get_query_as_string(RPQTree *q) {

    if(q->isLeaf()) {

        return q->data;
    }

    if(q->isConcat()) {

        std::string leftNodes = SimpleEvaluator::get_query_as_string(q->left);
        std::string rightNodes = SimpleEvaluator::get_query_as_string(q->right);
        return "(" + leftNodes + "/" + rightNodes +")";
    }

    return "";
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

            #ifdef EST_CACHE
            if(estCache.find(queryPair) != estCache.end()){

                estimates.push_back(estCache[queryPair]);
            }
            else{
            #endif
                auto queryTree = RPQTree::strToTree(queryPair);
                cardStat estimate = est->estimate(queryTree);
                delete(queryTree);
                estimates.push_back(estimate.noPaths);

            #ifdef EST_CACHE
                estCache.insert(std::make_pair(queryPair, estimate.noPaths));
            }
            #endif
        }

        long min_index = std::min_element(estimates.begin(), estimates.end()) - estimates.begin();
        std::string minPair = "(" + queryN[min_index] + "/" + queryN[min_index + 1] + ")";
        queryN[min_index] = minPair;
        queryN.erase(queryN.begin()+min_index+1);
    }

    if(queryN.size() == 2){

        // Now concat the remaining two labels to produce an optimal tree which can be evaluated

        std::string finalQuery = "(" + queryN[0] + "/" + queryN[1] + ")";
        return RPQTree::strToTree(finalQuery);
    }

    return query;
}

cardStat SimpleEvaluator::evaluate(RPQTree *query) {

    std::string queryString = get_query_as_string(query);

    #ifdef FINAL_CACHE
    if(finalCache.find(queryString) != finalCache.end()){

        return finalCache[queryString];
    }
    #endif

    std::shared_ptr<JoinGraph> res;
    if(est){
        auto optimized_query = rewrite_query_tree(query);
//        std::cout << std::endl << "Optimized query tree: ";
//        optimized_query->print();

        res = evaluate_aux(optimized_query);

        if(query->isConcat()) {
            delete(optimized_query);
        }
    } else {
        res = evaluate_aux(query);
    }

    auto stats = SimpleEvaluator::computeStats(res);

    #ifdef FINAL_CACHE
    finalCache.insert(std::make_pair(queryString, stats));
    #endif

    return stats;
}