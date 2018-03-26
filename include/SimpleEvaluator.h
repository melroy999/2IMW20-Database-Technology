//
// Created by Nikolay Yakovets on 2018-02-02.
//

#ifndef QS_SIMPLEEVALUATOR_H
#define QS_SIMPLEEVALUATOR_H


#include <memory>
#include <cmath>
#include "SimpleGraph.h"
#include "RPQTree.h"
#include "Evaluator.h"
#include "Graph.h"

class SimpleEvaluator : public Evaluator {

    std::shared_ptr<SimpleGraph> graph;
    std::shared_ptr<SimpleEstimator> est;
    std::map<std::string, uint32_t> estCache;
    std::map<std::string, cardStat> finalCache;
    std::map<std::string, std::shared_ptr<SimpleGraph>> resultsCache;

public:

    explicit SimpleEvaluator(std::shared_ptr<SimpleGraph> &g);
    ~SimpleEvaluator() = default;

    void prepare() override ;
    cardStat evaluate(RPQTree *query) override ;

    void attachEstimator(std::shared_ptr<SimpleEstimator> &e);

    std::shared_ptr<SimpleGraph> evaluate_aux(RPQTree *q);
    RPQTree * rewrite_query_tree(RPQTree *q);
    std::vector<std::string> get_query_graph(RPQTree *q, std::vector<std::string> nodes);
    std::string get_query_as_string(RPQTree *q);
    static std::shared_ptr<SimpleGraph> project(uint32_t label, bool inverse, std::shared_ptr<SimpleGraph> &g);
    std::shared_ptr<SimpleGraph> join(std::shared_ptr<SimpleGraph> &left, std::shared_ptr<SimpleGraph> &right);

    static cardStat computeStats(std::shared_ptr<SimpleGraph> &g);

};


#endif //QS_SIMPLEEVALUATOR_H
