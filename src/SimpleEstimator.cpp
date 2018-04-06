//
// Created by Nikolay Yakovets on 2018-02-01.
//

#include <SimpleEstimator.h>

SimpleEstimator::SimpleEstimator(std::shared_ptr<SimpleGraph> &g){

    // works only with SimpleGraph
    graph = g;
}

void SimpleEstimator::prepare() {

    // The number of bits we need to store all the vertices.
    auto noBins = static_cast<unsigned long>(ceil((double) graph->getNoVertices() / 64));

    // Resize the join data matrix.
    joinData.resize(2 * graph->getNoLabels(), std::vector<joinStat>(2 * graph->getNoLabels()));

    // Create a collection of label data.
    labelData = std::vector<labelStat>();
    for(uint32_t i = 0; i < 2 * graph->getNoLabels(); i++) {
        labelData.emplace_back(i, i % graph->getNoLabels(), i >= graph->getNoLabels(), graph->graphs[i % graph->getNoLabels()]);
    }

    // Do this outside of the loop, as it somehow dereferences pointers if we don't.
    for(uint32_t i = 0; i < graph->getNoLabels(); i++) {
        labelData[i + graph->getNoLabels()].setTwin(&labelData[i]);
    }

    // Find all the join information.
    // Here, joining the labels i and j should correspond to inv(j), inv(i).
    for(uint32_t i = 0; i < 2 * graph->getNoLabels(); i++) {
        for(uint32_t j = 0; j < 2 * graph->getNoLabels(); j++) {
            joinData[i][j].join(&labelData[i], &labelData[j]);
        }
    }


    #define CHECK_BIT(var,pos) ((var) & (1ULL<<(pos)))

    // Calculate the required 2-gram join statistics.
    for(auto &v : joinData) {
        for(auto &w : v) {

            // We dont want to waste time on joins without any edges.
            if(w.numCommonNodes != 0) {
                // The total number of paths, and the distinct number of pairs.
                uint32_t paths = 0;

                // Change the bucket size of the source and target collections to fit all the vertices.
                w.changeBucketSize(noBins);

                // Create pointers to the appropriate adjacency lists.
                std::vector<std::vector<Entry>*>* sourceAdj = w.source < graph->getNoLabels() ?
                                                &graph->graphs[w.source % graph->getNoLabels()].reverse_adj:
                                                &graph->graphs[w.source % graph->getNoLabels()].adj;

                std::vector<std::vector<Entry>*>* targetAdj = w.target < graph->getNoLabels() ?
                                                &graph->graphs[w.target % graph->getNoLabels()].adj:
                                                &graph->graphs[w.target % graph->getNoLabels()].reverse_adj;

                // For each of the common blocks
                for(uint32_t i = 0; i < w.commonNodes.size(); i++) {
                    auto bucket = w.commonNodes[i];

                    if(bucket != 0ULL) {
                        // Now, decompress each block in the common range j to j + 8.
                        for(uint32_t j = 0; j < 8; j++) {
                            // Check if the block has any values in the bucket.
                            if((bucket & 255) == 0) {
                                // If not, skip.
                                bucket >>= 8;
                                continue;
                            }

                            bucket >>= 8;

                            std::vector<Entry>* sourceBlocks = (*sourceAdj)[(i << 3) + j];
                            std::vector<Entry>* targetBlocks = (*targetAdj)[(i << 3) + j];

                            // We want that both collections are non-null, otherwise there are no join results.
                            if(sourceBlocks && targetBlocks) {

                                // Since both are non null, we can decompress the data.
                                std::vector<std::vector<uint32_t>> sources(8);
                                std::vector<std::vector<uint32_t>> targets(8);

                                for(const Entry &e : *sourceBlocks) {
                                    // Check each bit position.
                                    for(uint32_t u = 0; u < 64; u++) {
                                        if(CHECK_BIT(e.v, u)) {
                                            sources[u >> 3].push_back((e.i << 3) + (u & 7));
                                        }
                                    }
                                }

                                for(const Entry &e : *targetBlocks) {
                                    // Check each bit position.
                                    for(uint32_t u = 0; u < 64; u++) {
                                        if(CHECK_BIT(e.v, u)) {
                                            targets[u >> 3].push_back((e.i << 3) + (u & 7));
                                        }
                                    }
                                }

                                // Set the source and target flags, and calculate the number of paths.
                                for(uint32_t u = 0; u < 8; u++) {
                                    if(sources[u].empty() || targets[u].empty()) {
                                        continue;
                                    }

                                    w.numSourceEdges += sources[u].size();
                                    w.numTargetEdges += targets[u].size();
                                    paths += sources[u].size() * targets[u].size();

                                    for(const auto &vertex : sources[u]) {
                                        w.sourceNodes[vertex >> 6] |= SET_BIT(vertex & 63);
                                    }

                                    for(const auto &vertex : targets[u]) {
                                        w.targetNodes[vertex >> 6] |= SET_BIT(vertex & 63);
                                    }
                                }
                            }
                        }
                    }
                }

                w.numPaths = paths;
                w.numSourceNodes = countBitsSet(w.sourceNodes);
                w.numTargetNodes = countBitsSet(w.targetNodes);
            }
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
