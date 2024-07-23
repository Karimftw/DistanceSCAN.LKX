#include "../include/graph.h"
Graph graph;

void Graph::init_d_neighbors() {
    Timer timer(D_NEIGHBOR_TIME);
    path_weight = vector<unordered_map<int, double >>(n, unordered_map<int, double>{});
    d_neighbors = vector<vector<int >>(n, vector<int>());
    for (int n_id = 0; n_id < n; ++n_id) {
        d_neighbors[n_id].emplace_back(n_id);
        priority_queue<idpair, vector<idpair>, cmp_idpair> pq;
        pq.push(make_pair(n_id, 0));
        set_path_weight(n_id, n_id, 0);
        for (int nei_id: d_neighbors[n_id]) {
            pq.push(make_pair(nei_id, get_path_weight(n_id, nei_id)));
        }
        while (!pq.empty()) {
            auto cur_node = pq.top();
            pq.pop();
            if (cmp_double(cur_node.second, get_path_weight(n_id, cur_node.first)) == 1)continue;
            if (cur_node.first > n_id) {
                d_neighbors[n_id].emplace_back(cur_node.first);
                d_neighbors[cur_node.first].emplace_back(n_id);
            }

            for (int nei_id: adj_list[cur_node.first]) {
                if (nei_id < n_id) continue;
                double nei_dis = cur_node.second + edge_weight[cur_node.first][nei_id];
                if (cmp_double(nei_dis, config.distance) < 1) {
                    double old_dis = get_path_weight(n_id, nei_id);
                    if ((cmp_double(old_dis, -1) == 0 || cmp_double(old_dis, nei_dis) > 0)) {
                        set_path_weight(n_id, nei_id, nei_dis);
                        pq.push(make_pair(nei_id, nei_dis));
                    }
                }
            }
        }
    }
}

void Graph::init(const string &graph_path) {
    Timer timer(READ_GRAPH_TIME);
    INFO("Reading graph ...");
    this->data_folder = graph_path;
    init_nm();
    adj_list = vector<vector<int >>(n, vector<int>());
    edge_weight = vector<unordered_map<int, double >>(n, unordered_map<int, double>());
    string graph_file = data_folder;// + FILESEP;
    if (directed) {
        graph_file += "graph.txt";
    } else if (config.operation == CONVERT_GRAPH) {
        graph_file += "undirect_graph.txt";
        weighted = false;
    } else {
        if (data_folder.find("graphs_coauthors") != data_folder.npos) {
            graph_file += "undirect_graph2.txt";
        } else if (data_folder.find("LFR") != data_folder.npos || config.operation == EXPONLFR) {
            graph_file += "undirect_graph.txt";
        } else {
            //graph_file += "uniform_weighted_graph.txt";
            graph_file += "jac_graph.txt";
        }
    }

    FILE *fin = fopen(graph_file.c_str(), "r");
    if (weighted) {
        int t1, t2;
        double w;
        while (fscanf(fin, "%d%d%lf", &t1, &t2, &w) != EOF) {
            if (t1 == t2)continue;
            adj_list[t1].push_back(t2);
            edge_weight[t1][t2] = w;
        }
    } else {
        int t1, t2;
        while (fscanf(fin, "%d%d", &t1, &t2) != EOF) {
            if (t1 == t2)continue;
            adj_list[t1].push_back(t2);
        }
    }
    fclose(fin);
    result.n = this->n;
    result.m = this->m;
    cout << "init graph graph n: " << this->n << " m: " << this->m << endl;

    // Reading words from a text file
    string words_file = data_folder + "words.txt";  // Update the path as needed
    ifstream infile(words_file);
    if (!infile) {
        cerr << "Could not open the words file: " << words_file << endl;
        exit(1);
    }

    vector<string> words;
    string word;
    while (infile >> word) {
        words.push_back(word);
    }
    infile.close();

    cout << "Read " << words.size() << " words from " << words_file << endl;
    for (const auto &w : words) {
        cout << w << endl;
    }

    clusterID = vector<int>(n, -1);
    is_core = vector<int>(n, -1);
    similarity = vector<unordered_map<int, bool >>(n, unordered_map<int, bool>{});
    if (config.algo == W_SCAN) {
        weight_degree = vector<double>(n, 0);
        whole_weight = 0;
        if (config.similarityType == Config::cos) {
            for (int i = 0; i < n; ++i) {
                for (int nei: adj_list[i]) {
                    weight_degree[i] += edge_weight[i][nei] * edge_weight[i][nei];
                }
                weight_degree[i] = sqrt(weight_degree[i] + 1);
                whole_weight += weight_degree[i] + 1;
            }
        }
    }
    if (data_folder.find("graphs_coauthors") != data_folder.npos || data_folder.find("LFR") != data_folder.npos) {
        reweighted(config.type);
    }
    if (config.operation == CLUSTER_VALIDATION) {
        jac_res = vector<unordered_map<int, double >>(n, unordered_map<int, double>{});
    }
}

// Rest of your existing code...
