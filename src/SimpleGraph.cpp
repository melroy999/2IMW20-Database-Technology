//
// Created by Nikolay Yakovets on 2018-01-31.
//

#include <cmath>
#include "SimpleGraph.h"

SimpleGraph::SimpleGraph(uint32_t n)   {
    setNoVertices(n);
}

uint32_t SimpleGraph::getNoVertices() const {
    return V;
}

void SimpleGraph::setNoVertices(uint32_t n) {
    V = n;
}

uint32_t SimpleGraph::getNoEdges() const {
    return E;
}

uint32_t SimpleGraph::getNoLabels() const {
    return L;
}

void SimpleGraph::setNoLabels(uint32_t noLabels) {
    L = noLabels;
    graphs.resize(noLabels);
}

void SimpleGraph::readFromContiguousFile(const std::string &fileName) {

    std::fstream file(fileName, std::ios_base::in);

    std::string graphStructure;
    file >> graphStructure;

    std::istringstream iss(graphStructure);
    std::string line;
    uint32_t V = 0;
    uint32_t E = 0;
    uint32_t L = 0;
    int j = 0;
    while(std::getline(iss, line, ',')) {
        if(j == 0) {
            V = (uint32_t) std::stoul(line);
        }  else if(j == 1) {
            E = (uint32_t) std::stoul(line);
        } else if(j == 2) {
            L = (uint32_t) std::stoul(line);
        }
        j++;
    }

    setNoLabels(L);
    setNoVertices(V);

    std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> edges;
    edges.reserve(E);

    for(int i = 0; i < E; i++) {
        uint32_t subject, predicate, object;

        file >> subject >> predicate >> object;
        file.ignore(INT32_MAX, '\n');

        // We place them in the order of importance. We want the order predicate > subject > object.
        edges.emplace_back(predicate, subject, object);
    }

    // Sort the vector.
    sort(edges.begin(), edges.end(),
         [](
                 const std::tuple<uint32_t, uint32_t, uint32_t> &a,
                 const std::tuple<uint32_t, uint32_t, uint32_t> &b
         ) -> bool {
             // The order in the tuple is label, source, target.
             if(std::get<0>(a) < std::get<0>(b)) {
                 return true;
             } else if(std::get<0>(a) == std::get<0>(b)) {
                 // Find which block each vertex belongs to.
                 if((std::get<1>(a) >> 3) < (std::get<1>(b) >> 3)) {
                     // A is in a previous row.
                     return true;
                 } else if((std::get<1>(a) >> 3) == (std::get<1>(b) >> 3)) {
                     // Check if a is in a earlier column than b.
                     return (std::get<2>(a) >> 3) < (std::get<2>(b) >> 3);
                 }
             }
             return false;
         });

    // Remove the duplicates.
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());

    // The starting point of our iterations, and the end point.
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t>>::iterator start;
    auto end = edges.begin();

    // For each of the labels, do the insertions.
    for(uint32_t i = 0; i < L; i++) {
        // Get the set of iterators.
        start = end;
        end = std::find_if(start, edges.end(), [&i] (const std::tuple<uint32_t, uint32_t, uint32_t> &a) { return std::get<0>(a) == i + 1; });

        // The current block graph.
        BlockGraph* graph = &graphs[i];
        graph->setDataStructureSizes(V);

        // The current block x and y starting positions.
        uint32_t px = 0;
        uint32_t py = 0;
        uint64_t v = 0ULL;
        uint64_t w = 0ULL;
        uint32_t entryHelper = 0;
        bool first = true;

        // Continue until we end up at another label.
        while(start != end) {

            // Which block is the pair in?
            uint32_t x = (std::get<1>(*start) >> 3);
            uint32_t y = (std::get<2>(*start) >> 3);

            if(first || !(x == px && y == py)) {

//                std::cout << std::bitset<8>(entryHelper) << std::endl;
//                for(uint32_t k = 0; k < 64; k += 8) {
//                    std::cout << std::bitset<8>(v >> k) << std::endl;
//                }
//                std::cout << std::endl;
//
//                std::cout << std::bitset<8>(entryHelper >> 8) << std::endl;
//                for(uint32_t k = 0; k < 64; k += 8) {
//                    std::cout << std::bitset<8>(w >> k) << std::endl;
//                }
//                std::cout << std::endl;

                // Add the block to the graph.
                graph->addEdge(px, py, {py, entryHelper, v}, {px, (entryHelper >> 8) | ((entryHelper & 255) << 8), w});

                // Set new previous values.
                px = x;
                py = y;
                first = false;

                // Reset the data.
                v = 0ULL;
                w = 0ULL;
                entryHelper = 0;
            }

            // A set of bits denoting whether a row or column has values.
            entryHelper |= SET_BIT(std::get<1>(*start) & 7);
            entryHelper |= SET_BIT(8 + (std::get<2>(*start) & 7));

            // Set the bit in v.
            v |= SET_BIT(8 * (std::get<1>(*start) & 7) + (std::get<2>(*start) & 7));
            w |= SET_BIT(8 * (std::get<2>(*start) & 7) + (std::get<1>(*start) & 7));

            start++;
        }

        // Flush the remaining data.
        graph->addEdge(px, py, {py, entryHelper, v}, {px, (entryHelper << 8) | ((entryHelper & 255) >> 8), w});

        // Finalize the graph.
        graph->finalize();
    }

    // Sort each of the blocks in the reverse adjacency list on i coordinate.
    for(auto &g : graphs) {
        for(auto &vertices : g.reverse_adj) {
            if(vertices) {
                std::sort(vertices->begin(), vertices->end(), [](const Entry &a, const Entry b) -> bool {
                    return a.i < b.i;
                });
            }
        }
    }

    for(uint32_t i = 0; i < L; i++) {
        graphs[i].reportData(i);
    }

    // Close the file, as we are done.
    file.close();
}
