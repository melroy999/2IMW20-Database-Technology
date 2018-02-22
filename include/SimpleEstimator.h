//
// Created by Nikolay Yakovets on 2018-02-01.
//

#ifndef QS_SIMPLEESTIMATOR_H
#define QS_SIMPLEESTIMATOR_H

#include "Estimator.h"
#include "SimpleGraph.h"

class SimpleEstimator : public Estimator {

    std::shared_ptr<SimpleGraph> graph;
    std::vector<std::unordered_set<int>> distinctTargetVerticesPerLabel;
    std::vector<std::unordered_set<int>> distinctSourceVerticesPerLabel;

public:
    explicit SimpleEstimator(std::shared_ptr<SimpleGraph> &g);
    ~SimpleEstimator() = default;

    void prepare() override ;
    cardStat estimate(RPQTree *q) override ;
    std::vector<std::pair<int, bool>> parseTreeToList(RPQTree *q);
};


#endif //QS_SIMPLEESTIMATOR_H
