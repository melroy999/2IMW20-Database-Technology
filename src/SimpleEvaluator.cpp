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

    for(int source = 0; source < g->getNoVertices(); source++) {
        if(!g->adj[source].empty()) stats.noOut++;
    }

    stats.noPaths = g->getNoDistinctEdges();

    for(int target = 0; target < g->getNoVertices(); target++) {
        if(!g->reverse_adj[target].empty()) stats.noIn++;
    }

    return stats;
}

std::shared_ptr<SimpleGraph> SimpleEvaluator::project(uint32_t projectLabel, bool inverse, std::shared_ptr<SimpleGraph> &in) {

    auto out = std::make_shared<SimpleGraph>(in->getNoVertices());
    out->setNoLabels(in->getNoLabels());

    if(!inverse) {
        // going forward
        for(uint32_t source = 0; source < in->getNoVertices(); source++) {
            for (auto labelTarget : in->adj[source]) {

                auto label = labelTarget.first;
                auto target = labelTarget.second;

                if (label == projectLabel)
                    out->addEdge(source, target, label);
            }
        }
    } else {
        // going backward
        for(uint32_t source = 0; source < in->getNoVertices(); source++) {
            for (auto labelTarget : in->reverse_adj[source]) {

                auto label = labelTarget.first;
                auto target = labelTarget.second;

                if (label == projectLabel)
                    out->addEdge(source, target, label);
            }
        }
    }

    return out;
}

std::shared_ptr<SimpleGraph> SimpleEvaluator::join(std::shared_ptr<SimpleGraph> &left, std::shared_ptr<SimpleGraph> &right) {

    auto out = std::make_shared<SimpleGraph>(left->getNoVertices());
    out->setNoLabels(1);

    for(uint32_t leftSource = 0; leftSource < left->getNoVertices(); leftSource++) {
        for (auto labelTarget : left->adj[leftSource]) {

            int leftTarget = labelTarget.second;
            // try to join the left target with right source
            for (auto rightLabelTarget : right->adj[leftTarget]) {

                auto rightTarget = rightLabelTarget.second;
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
 * Perform post-order tree traversal to get all possible pairs (that is left to right)
 * @param q
 * @return
 */
RPQTree * SimpleEvaluator::get_pairs(RPQTree *q, RPQTree *currentTree, RPQTree *bestTree) {

    if(q->isLeaf()) {

        if(currentTree){

            std::string tree = get_string_query(currentTree);

            std::string full = "(" + tree + "/" + q->data + ")";
            RPQTree * newTree = RPQTree::strToTree(full);
            labelPairs.emplace_back(newTree);

            cardStat estimate = est->estimate(newTree);
            estimates.emplace_back(estimate.noPaths);
        }

        return q;
    }

    if(q->isConcat()) {

        // Check if this node is a "leaf" node

        if(bestTree){

            std::string qString = get_string_query(q);
            std::string bestTreeString = get_string_query(bestTree);

            if(qString == bestTreeString){

                if(currentTree){

                    std::string tree = get_string_query(currentTree);

                    std::string full = "(" + tree + "/" + bestTreeString + ")";
                    RPQTree * newTree = RPQTree::strToTree(full);
                    labelPairs.emplace_back(newTree);

                    cardStat estimate = est->estimate(newTree);
                    estimates.emplace_back(estimate.noPaths);
                }

                return q;
            }
        }

        // Traverse the left subtree recursively
        RPQTree *leftTree = SimpleEvaluator::get_pairs(q->left, currentTree, bestTree);

        // Traverse the right tree recursively
        RPQTree *rightTree = SimpleEvaluator::get_pairs(q->right, leftTree, bestTree);

        return rightTree;
    }

    return nullptr;
}

std::string SimpleEvaluator::get_string_query(RPQTree *q) {

    // evaluate according to the AST bottom-up

    if(q->isLeaf()) {

        std::string label = q->data;
        return label;
    }

    if(q->isConcat()) {

        std::string rightLabel = SimpleEvaluator::get_string_query(q->right);
        std::string leftLabel = SimpleEvaluator::get_string_query(q->left);
        std::string full = "(" + leftLabel + "/" + rightLabel + ")";
        return full;
    }

    return nullptr;
}

RPQTree * SimpleEvaluator::rewrite_query_tree(RPQTree *query) {

    RPQTree *bestTree = nullptr;
    std::string oldQuery = get_string_query(query);
    std::string newQuery = "";

    do{
        labelPairs.clear();
        estimates.clear();
        oldQuery = get_string_query(query);

        // Get the estimates for all possible pairs in the tree
        get_pairs(query, nullptr, bestTree);

        // Use the lowest estimate to rewrite the query so that pair is a leaf pair

        if(!estimates.empty()){

            long min_index = std::min_element(estimates.begin(), estimates.end()) - estimates.begin();
            RPQTree * minPair = labelPairs.at(min_index);

            // replace all occurrences of + or - with escaped variants
            std::regex rgx_Meta (R"(([\^\$\\\.\+\?\(\)\[\]\{\}\|]))");
            std::string leftPair = std::regex_replace(get_string_query(minPair->left), rgx_Meta, R"(\$1)");
            std::string rightPair = std::regex_replace(get_string_query(minPair->right), rgx_Meta, R"(\$1)");

            std::string leftBracket = leftPair + "\\)+/" + rightPair;
            std::string rightBracket = leftPair + "/\\(+" + rightPair;
            std::string bothBrackets = leftPair + "\\)+/\\(+" + rightPair;

            std::regex lr(leftBracket);
            std::regex rr(rightBracket);
            std::regex br(bothBrackets);

            std::string minPairString = get_string_query(minPair);

            if (std::regex_search (oldQuery, rr)){

                newQuery = std::regex_replace(oldQuery, rr, "(" + minPairString);
            }
            else if(std::regex_search (oldQuery, lr)){

                newQuery = std::regex_replace(oldQuery, lr, minPairString + ")");
            }
            else if(std::regex_search (oldQuery, br)){

                newQuery = std::regex_replace(oldQuery, br, minPairString);
            }
            else{
                newQuery = oldQuery;
            }

            query = RPQTree::strToTree(newQuery);

            bestTree = minPair;
        }
    }
    while(!estimates.empty() && newQuery != oldQuery);

    return query;
}

cardStat SimpleEvaluator::evaluate(RPQTree *query) {

    if(est){

        query = rewrite_query_tree(query);
    }

    auto res = evaluate_aux(query);
    return SimpleEvaluator::computeStats(res);
}