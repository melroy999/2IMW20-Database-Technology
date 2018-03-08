//
// Created by Nikolay Yakovets on 2018-02-01.
//

#include <SimpleEstimator.h>

#define CHECK_BIT(var,pos) ((var) & (1ULL<<(pos)))

// sort on the second item in the pair, then on the first (ascending order)
bool SimpleEstimator::sortEdges(const std::pair<uint32_t,uint32_t> &a, const std::pair<uint32_t,uint32_t> &b) {
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
        labelData.emplace_back(i, i % graph->getNoLabels(), i >= graph->getNoLabels(), noBins);
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
                labelData[edge.first].insertEdge(i, edge.second);

                first = false;
                prevTarget = edge.second;
                prevLabel = edge.first;
            }
        }
    }

    for(uint32_t i = 0; i < graph->getNoLabels(); i++) {
        labelData[i].calculateSize();
    }

    // Find all the join information.
    // Here, joining the labels i and j should correspond to inv(j), inv(i).
    for(uint32_t i = 0; i < 2 * graph->getNoLabels(); i++) {
        for(uint32_t j = 0; j < 2 * graph->getNoLabels(); j++) {
            joinData[i][j].join(&labelData[i], &labelData[j]);
        }
    }

    // Calculate the source and target nodes per label.
    for(auto &v : joinData) {
        for(auto &w : v) {


            // The total number of paths.
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
                                            w.numSourceEdges++;

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
                                            w.numTargetEdges++;

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

    // Calculate the distinct s,t pairs.
//    for(auto &v : joinData) {
//        for (auto &w : v) {
//
//            // The total number of s, t pairs.
//            uint32_t paths = 0;
//
//            // We dont want to waste time on joins without any edges.
//            if(w.numSourceNodes != 0) {
//
//                for(uint32_t i = 0; i < w.sourceNodes.size(); i++) {
//                    auto bucket = w.sourceNodes[i];
//
//                    if(bucket != 0ULL) {
//                        for(uint32_t j = 0; j < 64; j++) {
//                            if(CHECK_BIT(bucket, j)) {
//
//                                std::unordered_set<uint32_t> pairs;
//
//                                // Check which vertices we can reach starting in the source node i * 64 + j.
//                                for(const auto &edge : w.source < graph->getNoLabels() ? graph->adj[i * 64 + j] : graph->reverse_adj[i * 64 + j]) {
//                                    if (edge.first == w.source % graph->getNoLabels()) {
//
//                                        for(const auto &edge2 : w.target < graph->getNoLabels() ? graph->adj[edge.second] : graph->reverse_adj[edge.second]) {
//                                            if (edge.first == w.target % graph->getNoLabels()) {
//
//                                                pairs.insert(i * 64 + j + graph -> getNoVertices() * edge2.second);
//                                            }
//                                        }
//                                    }
//                                }
//
//                                paths += pairs.size();
//                            }
//                        }
//                    }
//                }
//            }
//
//            w.uniquePaths = paths;
//        }
//    }

//    for(const auto &data : labelData) {
//        data.print();
//    }
//
//    for(const auto &data : joinData) {
//        for(const auto &data2 : data) {
//            if(data2.numCommonNodes != 0) data2.print();
//        }
//    }

    // do your prep here

}



cardStat SimpleEstimator::estimate(RPQTree *q) {
    std::cout << std::endl;

    // perform your estimation here

    // XOR the targets of label 0+ with the sources of label 1+
//    auto result = doAnd(&labelData[0].targets, &labelData[1].sources);
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

    exCardStat evaluation = doEstimation(q);

    return static_cast<cardStat>(evaluation);
}

exCardStat SimpleEstimator::doEstimation(RPQTree *q) {
    // Estimate only if we are a leaf or want to join.
    if(q -> isLeaf()) {

        // Estimating a leaf is simple, just do a lookup in the data.
        std::string data = q->data;
        uint32_t label = std::stoi(data.substr(0, data.size() - 1)) + (data.back() == '+' ? 0 : graph->getNoLabels());

        exCardStat result = exCardStat {static_cast<double>(labelData[label].getNumSources()),
                                        static_cast<double>(labelData[label].getNumEdges()),
                                        static_cast<double>(labelData[label].getNumTargets())
        };
        result.vertices.push_back(label);

        return result;

    } else if(q -> isConcat()) {

        exCardStat leftStat = doEstimation(q -> left);
        exCardStat rightStat = doEstimation(q -> right);

        auto leftData = &labelData[leftStat.vertices.back()];
        auto rightData = &labelData[rightStat.vertices.front()];

        // Get the join data of the join between the left and right labels.
        auto join = &joinData[leftStat.vertices.back()][rightStat.vertices.front()];

        // Estimate the number of out vertices and the number of in vertices.
        double noOut, noIn;

        // If the number of vertices in our cardinality stat is larger than one, we have done at least one join.
        if(leftStat.vertices.size() > 1) {
            // Recall the join made between the two left neighboring labels:
            auto previousJoin = &joinData[leftStat.vertices.end()[-2]][leftStat.vertices.end()[-1]];

            // Find out which of the source nodes in the join are part of the previous join.
            auto usableSources = doAnd(&join->sourceNodes, &previousJoin->commonNodes);
            auto numberOfUsableSources = countBitsSet(&usableSources);

            // Compared to the interface, how many vertices are in the terminated nodes?
            unsigned long terminatedNodes = previousJoin->numCommonNodes - numberOfUsableSources;

            // Which number of source nodes can actually be used?
            auto averageDegree = (double) previousJoin->numSourceEdges / previousJoin->numCommonNodes;
            auto n = std::min(((averageDegree * numberOfUsableSources) / join->numSourceNodes), (double) 1);

            noOut = leftStat.noOut - terminatedNodes;
            noOut *= n;
        } else {
            // We can determine the maximum number of source nodes exactly, since we have no predecessors.
            noOut = join->numSourceNodes;
        }

        if(rightStat.vertices.size() > 1) {
            double pTargetTermination = (double) (rightData->getNumTargets() - join->numTargetNodes) / rightData->getNumTargets();

            // Estimate the number of noIn.
            noIn = rightStat.noIn * (1 - pTargetTermination);
        } else {
            // We can determine the number of target vertices exactly.
            noIn = join->numTargetNodes;
        }

        // We can put a better bound on noOut and noIn by observing the join data.
        if(leftStat.vertices.front() == leftStat.vertices.back()) {
            noOut = std::min((double) join->numSourceNodes, noOut);
        }

        if(rightStat.vertices.front() == rightStat.vertices.back()) {
            noIn = std::min((double) join->numTargetNodes, noIn);
        }

        // First of all, we have exact data on the number of (non-distinct) paths in the join,
        // and the average expected number of paths in the join using the average in and out degree.
        uint32_t maxPaths = join->numPaths;

        // Using the above, we can determine the expected number of follow up edges.
        // Note here that we use the raw path stats, since the max paths measure combined with
        // the degree sum for a specific label already covers for the terminated/initialized edges.
        double leftPathEstimation = leftStat.noPaths * (double) maxPaths / join->numSourceEdges;
        double rightPathEstimation = rightStat.noPaths * (double) maxPaths / join->numTargetEdges;

        double noPaths = (leftPathEstimation + rightPathEstimation) / 2;

        if(abs(leftStat.vertices.back() - rightStat.vertices.front()) == graph -> getNoLabels()) {
            // We know that at least the number of edges minus the number of source nodes are duplicates.
            unsigned long minimumDuplicateCount = leftData->getNumEdges() - leftData->getNumSources();
            noPaths *= (1 - (double) minimumDuplicateCount / join->numPaths);
        }

        // How many of the given paths is a duplicate?
        // For now, we just take the probability that the source and the target are equal:
        double pDuplicate = 1.0 / (leftData->getNumSources() + rightData->getNumTargets());

        noPaths *= (1 - pDuplicate);

        // Try to merge the two estimates with the join cardinality formula, with uniformity assumptions.
        exCardStat result = exCardStat {
                noOut,
                noPaths,
                noIn
        };

        // Update the cardinality stat list of vertices.
        result.vertices.insert(result.vertices.end(), leftStat.vertices.begin(), leftStat.vertices.end());
        result.vertices.insert(result.vertices.end(), rightStat.vertices.begin(), rightStat.vertices.end());

        return result;
    }

    // Return the estimate in the form {#outNodes, #paths, #inNodes}
    return exCardStat {0, 0, 0};
}
