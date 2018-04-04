//
// Created by Nikolay Yakovets on 2018-02-02.
//

#ifndef QS_SIMPLEEVALUATOR_H
#define QS_SIMPLEEVALUATOR_H

// #define EST_CACHE
// #define FINAL_CACHE
// #define RESULTS_CACHE

#include <memory>
#include <cmath>
#include "SimpleGraph.h"
#include "RPQTree.h"
#include "Evaluator.h"
#include "Graph.h"

class SimpleEvaluator : public Evaluator {

    std::shared_ptr<SimpleGraph> graph;
    std::shared_ptr<SimpleEstimator> est;

    #ifdef EST_CACHE
    std::map<std::string, uint32_t> estCache;
    #endif
    #ifdef FINAL_CACHE
    std::map<std::string, cardStat> finalCache;
    #endif
    #ifdef RESULTS_CACHE
    std::map<std::string, std::shared_ptr<JoinGraph>> resultsCache;
    #endif

public:

    explicit SimpleEvaluator(std::shared_ptr<SimpleGraph> &g);
    ~SimpleEvaluator() = default;

    void prepare() override ;
    cardStat evaluate(RPQTree *query) override ;

    void attachEstimator(std::shared_ptr<SimpleEstimator> &e);

    std::shared_ptr<JoinGraph> evaluate_aux(RPQTree *q);
    RPQTree * rewrite_query_tree(RPQTree *q);
    std::vector<std::string> get_query_graph(RPQTree *q, std::vector<std::string> nodes);
    std::string get_query_as_string(RPQTree *q);
    static std::shared_ptr<JoinGraph> project(uint32_t label, bool inverse, std::shared_ptr<SimpleGraph> &g);
    std::shared_ptr<JoinGraph> join(std::shared_ptr<JoinGraph> &left, std::shared_ptr<JoinGraph> &right);
    static cardStat computeStats(std::shared_ptr<JoinGraph> &g);

};


#endif //QS_SIMPLEEVALUATOR_H
