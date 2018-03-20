//
// Created by Nikolay Yakovets on 2018-02-01.
//

#include <SimpleEstimator.h>

#define CHECK_BIT(var,pos) ((var) & (1ULL<<(pos)))
#define SET_BIT(pos) (1ULL << (pos % 64))

// sort on the second item in the pair, then on the first (ascending order)
bool SimpleEstimator::sortEdges(const std::pair<uint32_t,uint32_t> &a, const std::pair<uint32_t,uint32_t> &b) {
    // Slight alteration, we want to sort on label instead.
    if (a.first < b.first) return true;
    if (a.first == b.first) return a.second < b.second;
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

    // If the sizes do not correspond, return an empty join.
    std::vector<uint64_t> result(t->size());
    if(t->size() != s->size()) {
        return result;
    }

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
        uint32_t degreeCount = 0;

        for(const auto &edge : graph->adj[i]) {
            if (first || !(prevTarget == edge.second && prevLabel == edge.first)) {

                if(!first && prevLabel != edge.first) {
                    if(degreeCount == 1) labelData[prevLabel].incrementSourceOut1();
                    degreeCount = 0;
                }

                labelData[edge.first].insertEdge(i, edge.second);

                first = false;
                prevTarget = edge.second;
                prevLabel = edge.first;
                degreeCount++;
            }
        }
        if(degreeCount == 1) labelData[prevLabel].incrementSourceOut1();

        prevTarget = 0;
        prevLabel = 0;
        first = true;
        degreeCount = 0;

        // Count which target vertices have an in degree of one.
        for(const auto &edge : graph->reverse_adj[i]) {
            if (first || !(prevTarget == edge.second && prevLabel == edge.first)) {

                if(!first && prevLabel != edge.first) {
                    if(degreeCount == 1) labelData[prevLabel].incrementTargetIn1();
                    degreeCount = 0;
                }

                first = false;
                prevTarget = edge.second;
                prevLabel = edge.first;
                degreeCount++;
            }
        }
        if(degreeCount == 1) labelData[prevLabel].incrementTargetIn1();
    }

    for(uint32_t i = 0; i < graph->getNoLabels(); i++) {
        labelData[i].calculateSize();
//        labelData[i].print();
//        labelData[i + graph->getNoLabels()].print();
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

            // The total number of paths, and the distinct number of pairs.
            uint32_t paths = 0;
//            uint32_t pairs = 0;

            // A scratch list to keep all the found edges.
//            auto pairCollection = std::vector<uint64_t>(2 * noBins);

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
                                            w.sourceNodes[edge.second / 64] |= SET_BIT(edge.second);
                                            sources++;
                                            w.numSourceEdges++;

                                            // Works incorrectly.
//                                            if(!CHECK_BIT(pairCollection[edge.second / 64], edge.second % 64) || !CHECK_BIT(pairCollection[noBins + i / 64], i % 64)) {
//                                                pairs++;
//                                                pairCollection[edge.second / 64] |= SET_BIT(edge.second);
//                                                pairCollection[noBins + i / 64] |= SET_BIT(i);
//                                            }

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
                                            w.targetNodes[edge.second / 64] |= SET_BIT(edge.second);
                                            targets++;
                                            w.numTargetEdges++;

                                            // Works incorrectly.
//                                            if(!CHECK_BIT(pairCollection[noBins + edge.second / 64], edge.second % 64) || !CHECK_BIT(pairCollection[i / 64], i % 64)) {
//                                                pairs++;
//                                                pairCollection[noBins + edge.second / 64] |= SET_BIT(edge.second);
//                                                pairCollection[i / 64] |= SET_BIT(i);
//                                            }

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
//            w.uniquePaths = pairs;
            w.numSourceNodes = countBitsSet(&w.sourceNodes);
            w.numTargetNodes = countBitsSet(&w.targetNodes);
        }
    }
}

cardStat SimpleEstimator::estimate(RPQTree *q) {

    return static_cast<cardStat>(doEstimation(q));
}

/**
 * Estimate the cardinality of a leaf node, which corresponds to looking up values in our data structure.
 *
 * @param q The query parse tree.
 * @return The cardinality data of the leaf node.
 */
exCardStat SimpleEstimator::estimateLeafNode(RPQTree *q) {

    // Estimating a leaf's cardinality is simple, just do a lookup in the data.
    std::string data = q->data;

    // Parse the data.
    uint32_t label = std::stoi(data.substr(0, data.size() - 1)) + (data.back() == '+' ? 0 : graph->getNoLabels());

    // Report the result.
    exCardStat result = exCardStat {
            static_cast<double>(labelData[label].getNumSources()),
            static_cast<double>(labelData[label].getNumEdges()),
            static_cast<double>(labelData[label].getNumTargets())
    };
    result.vertices.push_back(label);

    return result;
}

/**
 * Estimate the cardinality of a join between two leaf nodes, which can be achieved through a lookup.
 *
 * @param leftStat The cardinality estimation of the left subtree.
 * @param rightStat The cardinality estimation of the right subtree.
 * @return The cardinality data of the join between two leaf nodes.
 */
exCardStat SimpleEstimator::estimateSimpleJoin(exCardStat *leftStat, exCardStat *rightStat) {

    // Get the join data of the join between the left and right labels.
    auto join = &joinData[leftStat->vertices.back()][rightStat->vertices.front()];

    // De-duplication factor.
    uint32_t dedup = 0;
    if(abs(leftStat->vertices.back() - rightStat->vertices.front()) == graph -> getNoLabels()) {
        // We know that at least the number of edges minus the number of source nodes are duplicates.
        dedup = join->numSourceEdges - join->numSourceNodes;
    }

    // Copy over the data.
    return exCardStat {
            static_cast<double>(join->numSourceNodes),
            static_cast<double>(join->numPaths - dedup),
            static_cast<double>(join->numTargetNodes)
    };
}

void estimateInAndOutCardinality(exCardStat* leftStat, exCardStat* rightStat, labelStat* leftData, labelStat* rightData,
                                 joinStat* joinStat, double* noOut, double* noIn) {
    /*
         * During the join, use the fact that the difference between the two sets of
         * target/source vertices we merge on is not necessarily the empty set.
         *
         * We know for example that certain target vertices in the left subtree
         * are not succeeded by the starting label of the right tree.
         *
         * On the other hand, we have that certain vertices in the right subtree
         * are not proceeded by the end label of the left tree.
         */
    // What is the number of edges that have been terminated and initialized during the merge?
    uint32_t terminatedEdges = leftData->getNumEdges() - joinStat->numSourceEdges;
    uint32_t initializedEdges = rightData->getNumEdges() - joinStat->numTargetEdges;

    // What is the probability that a vertex ends up with no connected edges, when an edge gets removed?
    // In other words, what is the chance that a vertex of degree one has its edge removed?
    double pSourceRemove = (double) leftData->getSourceOut1() / leftData->getNumEdges();
    double pTargetRemove = (double) rightData->getTargetIn1() / rightData->getNumEdges();

    // What about higher degree edges? In certain cases we have a very low amount of degree 1 edges,
    // and should thus compensate.

    // So, what is the expected number of source and target vertices given the above probabilities?
    // Here we do not simply subtract, as we would risk getting negative numbers.
    *noOut = leftStat->noOut * (1 - (pSourceRemove * terminatedEdges) / leftData->getNumSources());
    *noIn = rightStat->noIn * (1 - (pTargetRemove * initializedEdges) / rightData->getNumTargets());
}

/**
 * Estimate the cardinality of a join between two sub-queries.
 *
 * @param leftStat The cardinality estimation of the left subtree.
 * @param rightStat The cardinality estimation of the right subtree.
 * @return The cardinality data of the join between two sub-queries.
 */
exCardStat SimpleEstimator::estimateJoin(exCardStat *leftStat, exCardStat *rightStat) {

    auto leftData = &labelData[leftStat->vertices.back()];
    auto rightData = &labelData[rightStat->vertices.front()];

    // Get the join data of the join between the left and right labels.
    auto join = &joinData[leftStat->vertices.back()][rightStat->vertices.front()];

    double noOut = 0;
    double noPaths = 0;
    double noIn = 0;

    // Estimate the noOut and noIn metric the old way.
    estimateInAndOutCardinality(leftStat, rightStat, leftData, rightData, join, &noOut, &noIn);

    // First of all, we have exact data on the number of (non-distinct) paths in the join,
    // and the average expected number of paths in the join using the average in and out degree.
    uint32_t maxPaths = join->numPaths;

    // Using the above, we can determine the expected number of follow up edges.
    // Note here that we use the raw path stats, since the max paths measure combined with
    // the degree sum for a specific label already covers for the terminated/initialized edges.
    double leftPathEstimation = leftStat->noPaths * (double) maxPaths / join->numSourceEdges;
    double rightPathEstimation = rightStat->noPaths * (double) maxPaths / join->numTargetEdges;

    // We can put a better bound on noOut and noIn by observing the join data.
    if(leftStat->vertices.front() == leftStat->vertices.back()) {
        noOut = std::min((double) join->numSourceEdges, noOut);
    }

    if(rightStat->vertices.front() == rightStat->vertices.back()) {
        noIn = std::min((double) join->numTargetEdges, noIn);
    }

    noPaths = (leftPathEstimation + rightPathEstimation) / 2;

    if(abs(leftStat->vertices.back() - rightStat->vertices.front()) == graph -> getNoLabels()) {
        // We know that at least the number of edges minus the number of source nodes are duplicates.
        unsigned long minimumDuplicateCount = leftData->getNumEdges() - leftData->getNumSources();
        noPaths -= minimumDuplicateCount;
    }

    // How many of the given paths is a duplicate?
    // For now, we just take the probability that the source and the target are equal:
//        double pDuplicate = 1.0 / (leftData->getNumSources() + rightData->getNumTargets());
//        noPaths *= (1 - pDuplicate);

    // Try to merge the two estimates with the join cardinality formula, with uniformity assumptions.
    return exCardStat {
            noOut,
            std::min(noPaths, noOut * noIn),
            noIn
    };
}



exCardStat SimpleEstimator::doEstimation(RPQTree *q) {
    // Estimate only if we are a leaf or want to join.
    if(q -> isLeaf()) {
        return estimateLeafNode(q);
    } else if(q -> isConcat()) {

        exCardStat leftStat = doEstimation(q -> left);
        exCardStat rightStat = doEstimation(q -> right);
        exCardStat result;

        if(leftStat.vertices.size() == 1 && rightStat.vertices.size() == 1) {
            result = estimateSimpleJoin(&leftStat, &rightStat);
        } else {
            result = estimateJoin(&leftStat, &rightStat);
        }

        // Update the cardinality stat list of vertices.
        result.vertices.insert(result.vertices.end(), leftStat.vertices.begin(), leftStat.vertices.end());
        result.vertices.insert(result.vertices.end(), rightStat.vertices.begin(), rightStat.vertices.end());

        return result;
    }

    // Return the estimate in the form {#outNodes, #paths, #inNodes}
    return exCardStat {0, 0, 0};
}
