//
// Created by Nikolay Yakovets on 2018-02-01.
//

#include <cmath>
#include <SimpleEstimator.h>

#define CHECK_BIT(var,pos) ((var) & (1ULL<<(pos)))

// sort on the second item in the pair, then on the first (ascending order)
bool sortEdges(const std::pair<uint32_t,uint32_t> &a, const std::pair<uint32_t,uint32_t> &b) {
    if (a.second < b.second) return true;
    if (a.second == b.second) return a.first < b.first;
    return false;
}

SimpleEstimator::SimpleEstimator(std::shared_ptr<SimpleGraph> &g){

    // works only with SimpleGraph
    graph = g;
}

uint32_t SimpleEstimator::countBitsSet(std::vector<uint64_t>* result) {
    uint32_t sum = 0;
    for(auto v : *result) {
        if(v == 0) continue;
        sum += __builtin_popcountll(v);
    }
    return sum;
}

std::vector<uint64_t> SimpleEstimator::doAnd(const std::vector<uint64_t> *t, const std::vector<uint64_t> *s) {
    std::vector<uint64_t> result(t->size());
    for(unsigned long i = t->size() ; i -- > 0 ; ) {
        result[i] = (*t)[i] & (*s)[i];
    }
    return result;
}

void SimpleEstimator::prepare() {

    // The number of bits we need to store all the vertices.
    auto noBins = static_cast<unsigned long>(ceil((double) graph->getNoVertices() / 64));

    // Resize the join data matrix.
    joinData.resize(2 * graph->getNoLabels(), std::vector<joinStat>(2 * graph->getNoLabels()));

    // Create a collection of label data.
    labelData = std::vector<labelStat>();
    for(uint32_t i = 0; i < 2 * graph->getNoLabels(); i++) {
        labelData.emplace_back(i, i >= graph->getNoLabels(), noBins);
//        if(i >= graph->getNoLabels())
//            labelData[i].setTwin(&labelData[i % graph->getNoLabels()]);
//
//        if(i == 2 * graph->getNoLabels() - 7) {
//            std::cout << "A" << std::endl;
//        }
    }

    // Do this outside of the loop, as it somehow dereferences pointers if we don't.
    for(uint32_t i = 0; i < graph->getNoLabels(); i++) {
        labelData[i + graph->getNoLabels()].setTwin(&labelData[i]);
    }


    for(uint32_t i = 0; i < graph->getNoVertices(); i++) {
        // Presort both the adjacency and reverse adjacency graphs.
        std::sort(graph -> adj[i].begin(), graph -> adj[i].end(), sortEdges);
        std::sort(graph -> reverse_adj[i].begin(), graph -> reverse_adj[i].end(), sortEdges);

        // Filter out duplicates.
        uint32_t prevTarget = 0;
        uint32_t prevLabel = 0;
        bool first = true;

        for(const auto &edge : graph->adj[i]) {
            if (first || !(prevTarget == edge.second && prevLabel == edge.first)) {
                labelData[edge.first].sources[i / 64] |= 1ULL << (i % 64);
                labelData[edge.first].targets[edge.second / 64] |= 1ULL << (edge.second % 64);
                labelData[edge.first].numEdges++;

                first = false;
                prevTarget = edge.second;
                prevLabel = edge.first;
            }
        }
    }

    for(uint32_t i = 0; i < graph->getNoLabels(); i++) {
        labelData[i].numSources = countBitsSet(&labelData[i].sources);
        labelData[i].numTargets = countBitsSet(&labelData[i].targets);
//        std::cout << "Label " << i << ": #sources=" << labelData[i].numSources << ", #targets=" << labelData[i].numTargets << std::endl;
    }
//    std::cout << std::endl;

    // Find all the join information.
    // Here, joining the labels i and j should correspond to inv(j), inv(i).
    for(uint32_t i = 0; i < 2 * graph->getNoLabels(); i++) {
        for(uint32_t j = 0; j < 2 * graph->getNoLabels(); j++) {
            joinData[i][j].join(&labelData[i], &labelData[j]);

//            std::cout << "(" << i << ", " << j << "): " << joinData[i][j].numCommonNodes << std::endl;
        }
    }
//    std::cout << std::endl;

    // Calculate the source and target nodes per label.
    for(auto v : joinData) {
        for(auto w : v) {

            uint32_t paths = 0;

            // We dont want to waste time on joins without any edges.
            if(w.numCommonNodes != 0) {

                // Change the bucket size of the source and target collections to fit all the vertices.
                w.changeBucketSize(noBins);

                for(uint32_t i = 0; i < w.commonNodes.size(); i++) {
                    auto bucket = w.commonNodes[i];

                    if(bucket != 0ULL) {
                        for(uint32_t j = 0; j < 64; j++) {
                            if(CHECK_BIT(bucket, j)) {
                                int sources = 0;

                                // Filter out duplicates.
                                uint32_t prevTarget = 0;
                                uint32_t prevLabel = 0;
                                bool first = true;

                                for(const auto &edge : w.source < graph->getNoLabels() ? graph->reverse_adj[i * 64 + j] : graph->adj[i * 64 + j]) {
                                    if(edge.first == w.source % graph->getNoLabels()) {
                                        if (first || !(prevTarget == edge.second && prevLabel == edge.first)) {
                                            w.sourceNodes[edge.second / 64] |= 1ULL << (edge.second % 64);
                                            sources++;

                                            first = false;
                                            prevTarget = edge.second;
                                            prevLabel = edge.first;
                                        }
                                    }
                                }

                                prevTarget = 0;
                                prevLabel = 0;
                                first = true;

                                int targets = 0;
                                for(const auto &edge : w.target < graph->getNoLabels() ? graph->adj[i * 64 + j] : graph->reverse_adj[i * 64 + j]) {
                                    if(edge.first == w.target % graph->getNoLabels()) {

                                        if (first || !(prevTarget == edge.second && prevLabel == edge.first)) {
                                            w.targetNodes[edge.second / 64] |= 1ULL << (edge.second % 64);
                                            targets++;

                                            first = false;
                                            prevTarget = edge.second;
                                            prevLabel = edge.first;
                                        }
                                    }
                                }

                                paths += sources * targets;
                            }
                        }
                    }
                }
            }

            w.numPaths = paths;
            w.numSourceNodes = countBitsSet(&w.sourceNodes);
            w.numTargetNodes = countBitsSet(&w.targetNodes);

//            std::cout << "(" << w.source << ", " << w.target << "): " << w.numSourceNodes << ", " << w.numTargetNodes << std::endl;

        }
    }
//    std::cout << std::endl;


    // do your prep here

}



cardStat SimpleEstimator::estimate(RPQTree *q) {

    // perform your estimation here

    // XOR the targets of label 0+ with the sources of label 1+
    auto result = doAnd(&labelData[0].targets, &labelData[1].sources);
//
//
//
//    std::vector<uint64_t> result(labelData[0].sources.size());
//    for(int i = 0; i < graph->getNoVertices(); i++) {
//        if((labelData[0].sources[i / 64] >> i % 64) & 1U) {
//            for(const auto &loc : labelData[0].edges[i]) {
//                result[loc.first] |= 1ULL << loc.first;
//            }
//        }
//    }



    return cardStat {0, 0, 0};
}