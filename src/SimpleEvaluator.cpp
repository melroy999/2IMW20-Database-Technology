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

bool found_left = false;

RPQTree * SimpleEvaluator::rewrite_tree(RPQTree *q, RPQTree *newTree) {

    // If we encounter a node which is equal to the tree being looked for, return that tree

    std::string queryString = get_string_query(q);
    std::string newString = get_string_query(newTree);

    if(queryString == newString){

        return q;
    }

    // Traverse the tree all the way to the left hand side

    if(q->isConcat()){

        RPQTree *leftTree = SimpleEvaluator::rewrite_tree(q->left, newTree);
        q->left = leftTree;

        // Once there is nothing on the left of the tree anymore,
        // check if the right of the node, is the same as the left hand part of the new tree

        std::string leftSide = get_string_query(newTree->left);
        std::string rightNodeText = get_string_query(q->right);

        if(leftSide == rightNodeText && !found_left){

            // If so, remove the right hand side and bump up the left hand side of the tree

            q = q->left;
            found_left = true;
            return q;
        }


        // Check if we are at a node where it is the same as the right hand side of the new tree,
        // and if so, replace the node with the contents of the new tree

        std::string rightSide = get_string_query(newTree->right);

        if(queryString == rightSide && found_left){

            return newTree;
        }

        // If none of this holds, keep traversing right until we find something

        RPQTree *rightTree = SimpleEvaluator::rewrite_tree(q->right, newTree);
        q->right = rightTree;

        return q;
    }

    if(q->isLeaf()){

        // Check if we are at a node where it is the same as the right hand side of the new tree,
        // and if so, replace the node with the contents of the new tree

        std::string rightSide = get_string_query(newTree->right);

        if(queryString == rightSide && found_left){

            return newTree;
        }

        return q;
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

cardStat SimpleEvaluator::evaluate(RPQTree *query) {

    RPQTree *bestTree = nullptr;
    std::string oldQuery = get_string_query(query);
    std::string newQuery = "";

    do{
        labelPairs.clear();
        estimates.clear();
        found_left = false;
        oldQuery = get_string_query(query);

        // Get the estimates for all possible pairs in the tree
        get_pairs(query, nullptr, bestTree);

        // Use the lowest estimate to rewrite the query so that pair is a leaf pair

        if(!estimates.empty()){

            long min_index = std::min_element(estimates.begin(), estimates.end()) - estimates.begin();
            RPQTree * minPair = labelPairs.at(min_index);

            query = rewrite_tree(query, minPair);
            newQuery = get_string_query(query);
            bestTree = minPair;
        }
    }
    while(!estimates.empty() && newQuery != oldQuery);

    std::string queryAsString = get_string_query(query);
    auto res = evaluate_aux(query);
    return SimpleEvaluator::computeStats(res);
}