#include <iostream>
#include <chrono>
#include <SimpleGraph.h>
#include <Estimator.h>
#include <SimpleEstimator.h>
#include <SimpleEvaluator.h>


struct query {
    std::string s;
    std::string path;
    std::string t;

    void print() {
        std::cout << s << ", " << path << ", " << t << std::endl;
    }
};

std::vector<query> parseQueries(std::string &fileName) {

    std::vector<query> queries {};

    std::string line;
    std::ifstream graphFile { fileName };

    std::regex edgePat (R"((.+),(.+),(.+))");

    while(std::getline(graphFile, line)) {
        std::smatch matches;

        // match edge data
        if(std::regex_search(line, matches, edgePat)) {
            auto s = matches[1];
            auto path = matches[2];
            auto t = matches[3];

            queries.emplace_back(query{s, path, t});
        }
    }

    graphFile.close();

    if(queries.size() == 0) std::cout << "Did not parse any queries... Check query file." << std::endl;

    return queries;
}

int estimatorBench(std::string &graphFile, std::string &queriesFile) {

    std::cout << "\n(1) Reading the graph into memory and preparing the estimator...\n" << std::endl;

    // read the graph
    auto g = std::make_shared<SimpleGraph>();

    auto start = std::chrono::steady_clock::now();
    try {
        g->readFromContiguousFile(graphFile);
    } catch (std::runtime_error &e) {
        std::cerr << e.what() << std::endl;
        return 0;
    }

    auto end = std::chrono::steady_clock::now();
    std::cout << "Time to read the graph into memory: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms" << std::endl;

    // prepare the estimator
    auto est = std::make_unique<SimpleEstimator>(g);
    start = std::chrono::steady_clock::now();
    est->prepare();
    end = std::chrono::steady_clock::now();
    std::cout << "Time to prepare the estimator: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms" << std::endl;

    std::cout << "\n(2) Running the query workload..." << std::endl;

    for(auto query : parseQueries(queriesFile)) {

        // perform estimation
        // parse the query into an AST
        std::cout << "\nProcessing query: ";
        query.print();
        RPQTree *queryTree = RPQTree::strToTree(query.path);
        std::cout << "Parsed query tree: ";
        queryTree->print();

        start = std::chrono::steady_clock::now();
        auto estimate = est->estimate(queryTree);
        end = std::chrono::steady_clock::now();

        std::cout << "\nEstimation (noOut, noPaths, noIn) : ";
        estimate.print();
        std::cout << "Time to estimate: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms" << std::endl;

        // perform evaluation
        auto ev = std::make_unique<SimpleEvaluator>(g);
        ev->prepare();
        start = std::chrono::steady_clock::now();
        auto actual = ev->evaluate(queryTree);
        end = std::chrono::steady_clock::now();

        std::cout << "Actual (noOut, noPaths, noIn) : ";
        actual.print();
        std::cout << "Time to evaluate: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms" << std::endl;

        // clean-up
        delete(queryTree);

    }

    return 0;
}

int evaluatorBench(std::string &graphFile, std::string &queriesFile) {

    std::cout << "\n(1) Reading the graph into memory and preparing the evaluator...\n" << std::endl;

    // read the graph
    auto g = std::make_shared<SimpleGraph>();

    auto start = std::chrono::steady_clock::now();
    try {
        g->readFromContiguousFile(graphFile);
    } catch (std::runtime_error &e) {
        std::cerr << e.what() << std::endl;
        return 0;
    }

    auto end = std::chrono::steady_clock::now();
    std::cout << "Time to read the graph into memory: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms" << std::endl;

    // prepare the evaluator
    auto est = std::make_shared<SimpleEstimator>(g);
    auto ev = std::make_unique<SimpleEvaluator>(g);
    ev->attachEstimator(est);

    start = std::chrono::steady_clock::now();
    ev->prepare();
    end = std::chrono::steady_clock::now();
    std::cout << "Time to prepare the evaluator: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms" << std::endl;

    std::cout << "\n(2) Running the query workload..." << std::endl;

    for(auto query : parseQueries(queriesFile)) {

        // perform estimation
        // parse the query into an AST
        std::cout << "\nProcessing query: ";
        query.print();
        RPQTree *queryTree = RPQTree::strToTree(query.path);
        std::cout << "Parsed query tree: ";
        queryTree->print();

        // perform the evaluation
        start = std::chrono::steady_clock::now();
        auto actual = ev->evaluate(queryTree);
        end = std::chrono::steady_clock::now();

        std::cout << "\nActual (noOut, noPaths, noIn) : ";
        actual.print();
        std::cout << "Time to evaluate: " << std::chrono::duration<double, std::milli>(end - start).count() << " ms" << std::endl;

        // clean-up
        delete(queryTree);

    }

    return 0;
}


int main(int argc, char *argv[]) {

    if(argc < 3) {
        std::cout << "Usage: quicksilver <graphFile> <queriesFile>" << std::endl;
        return 0;
    }

    // args
    std::string graphFile {argv[1]};
    std::string queriesFile {argv[2]};

//    estimatorBench(graphFile, queriesFile);
    evaluatorBench(graphFile, queriesFile);

    return 0;
}



