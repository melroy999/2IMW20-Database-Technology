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
std::string SimpleEvaluator::rewrite_query_tree(RPQTree *q, std::string label) {

    if(q->isLeaf()) {

        if(!label.empty()){

            std::string full = "(" + label + "/" + q->data + ")";
            RPQTree * newTree = RPQTree::strToTree(full);
            labelPairs.emplace_back(newTree);

            cardStat estimate = est->estimate(newTree);
            estimates.emplace_back(estimate.noPaths);
        }

        return q->data;
    }

    if(q->isConcat()) {

        // Traverse the left subtree recursively
        std::string leftLabel = SimpleEvaluator::rewrite_query_tree(q->left, label);

        // Traverse the right tree recursively
        std::string rightLabel = SimpleEvaluator::rewrite_query_tree(q->right, leftLabel);

        return rightLabel;
    }

    return "";
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

cardStat SimpleEvaluator::evaluate(RPQTree *query) {

    // Get the estimates for all possible pairs in the tree
    rewrite_query_tree(query, "");
    std::string queryAsString = get_string_query(query);

    // Use the lowest estimate to rewrite the query so that pair is a leaf pair
    long min_index = std::min_element(estimates.begin(), estimates.end()) - estimates.begin();
    RPQTree * minPair = labelPairs.at(min_index);

    std::string leftBracket = minPair->left->data + ")/" + minPair->right->data;
    std::string rightBracket = minPair->left->data + "/(" + minPair->right->data;

    if (queryAsString.find(leftBracket) != std::string::npos) {

        queryAsString.replace(queryAsString.find(leftBracket),leftBracket.length(), "(" + minPair->left->data + "/" + minPair->right->data + "))");
    }
    else if(queryAsString.find(rightBracket) != std::string::npos){

        queryAsString.replace(queryAsString.find(rightBracket),rightBracket.length(), "((" + minPair->left->data + "/" + minPair->right->data + ")");
    }

    query = RPQTree::strToTree(queryAsString);
    auto res = evaluate_aux(query);
    return SimpleEvaluator::computeStats(res);
}