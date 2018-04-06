//
// Created by Nikolay Yakovets on 2018-02-03.
//

#ifndef QS_EVALUATOR_H
#define QS_EVALUATOR_H

#include <memory>
#include "Graph.h"
#include "Estimator.h"

class Evaluator {

public:
    virtual void prepare() = 0;
    virtual cardStat evaluate(RPQTree *query) = 0;

};


#endif //QS_EVALUATOR_H
