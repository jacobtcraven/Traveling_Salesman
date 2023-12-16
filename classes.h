#ifndef CLASSES_H
#define CLASSES_H

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <utility>
#include <limits>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <chrono>

using namespace std;

int finc = 1;
int finc2 = 1;
int finc3 = 1;

class Node{
    public:
        int ID;
        double key;
        int parent;
        double x;
        double y;

        Node(int ID, double key)
            : ID{ID}, key{key}
        {}
        Node()
            : ID{-1}, key{-1} , parent{-1} {}
        Node(int ID, double key, double x, double y)
            : ID{ID}, key{key} , parent{-1}, x{x}, y{y} {}
};



bool operator > (const Node& A, const Node& B) {
    return A.key > B.key;
}

bool operator < (const Node& A, const Node& B) {
    return A.key < B.key;
}

bool operator == (const Node& A, const Node& B) {
    return A.key == B.key;
}

class Heap{
    friend class PriorityQueue;
    private:
        vector<Node> heap;
        int heapsize;
        vector<int> location;
    public:
        Heap(){
            heap.resize(20);
            heapsize = 20;
            location.resize(20);
        }
        Heap(int n){
            heap.resize(n);
            for (int i = 0; i < n; i++){
                heap[i].ID = i;
                heap[i].key = -1;
            }
            heapsize = 0;
            location.resize(n, -1);
        }
        Heap(Node A[]){
            int size = 8;
            for(int i=0; i < size; i++){
                heap[i] = A[i];
                location[i] = A[i].ID;
            }
            heapsize = size;
            buildHeap();
        }
        int parent(int i){
            return (i-1) / 2;
        }
        int left(int i){
            return heap[2 * i + 1].key;
        }
        int right(int i){
            return heap[2 * i + 2].key;
        }
        void heapify(int i){
            int smallest = -1;
            int left = 2 * i + 1;
            int right = 2 * i + 2;

            if (left < heapsize && heap[left].key < heap[i].key){
                smallest = left;
            }
            else{
                smallest = i;
            }
            
            if (right < heapsize && heap[right].key < heap[smallest].key)
                smallest = right;

            if (smallest != i) {
                Node hold(-1,-1);
                hold = heap[i];
                heap[i] = heap[smallest];
                heap[smallest] = hold;
                heapify(smallest);
            }
        }
        void buildHeap(){
            for (int i = heapsize / 2 - 1; i >= 0; i--){
                heapify(i);
            }
        }
        void increaseHeapSize(){
            heapsize += 1;
        }
};

class PriorityQueue : public Heap{
    public:
        vector<bool> inQueue;  




        PriorityQueue(int maxSize){
            inQueue.resize(maxSize, false);
            heap.resize(maxSize);
            heapsize = 0;
            location.resize(maxSize, -1);

        }
        void insert(Node newNode){
            heapsize += 1;
            inQueue[newNode.ID] = true;  
            Node A(heapsize - 1, 10000000);
            heap[heapsize - 1] = A;
            location[A.ID] = heapsize - 1;
            changeKey(heapsize - 1, newNode.key);
        }
        Node get(int i){
            return heap[i];
        }
        void changeKey(int i, int key){
            int rID = location[i];
            if (key > heap[rID].key){
                return;
            }
            Node hold(-1, -1);
            heap[rID].key = key;
            while (rID > 0 && heap[parent(rID)].key > heap[rID].key){
                hold = heap[rID];
                
                heap[rID] = heap[parent(rID)];
                heap[parent(rID)] = hold;
                location[parent(rID)] = rID;
                location[rID] = parent(rID);
                rID = parent(rID);
                
            }
        }
        Node peek() {
            return heap[0];
        }
        Node extract() {
            if (heapsize < 1){
                cerr << "heap underflow";
            }
            Node min = heap[0];
            int indx = 0;
            for (int i = 0; i < heapsize; i++){
                if(heap[i] < min){
                    min = heap[i];
                    indx = i;
                }
            }
            inQueue[min.ID] = false; 
            heap[indx] = heap[heapsize - 1];
            heapsize -= 1;
            heapify(0);
            return min;
        }
        bool isEmpty(){
            if (heapsize == 0){
                return true;
            }
            else{
                return false;
            }
        }
        bool contains(int ID) {
            if (ID >= 0 && ID < inQueue.size()) {
                return inQueue[ID];
            }
            else{
                return false;
            }
        }


};

class Graph {
public:
    Graph(int V);
    void addE(int u, int v, double weight);
    void addVert(int u, int z, double w);
    int vCount() const;
    const vector<pair<int, double>> adjVert(int u) const;
    vector<Node> vertices;
    vector<vector<double>> wTrack;
private:
    int verts;
    vector<list<pair<int, double>>> adjList;
};

Graph::Graph(int V) : verts(V), adjList(V) {
    wTrack = vector<vector<double>>(V + 1, vector<double>(V + 1, 0));
}

void Graph::addE(int u, int v, double weight) {
    adjList[u].emplace_back(v, weight);
    vertices.push_back(Node(-1,-1,1,1));
    wTrack[u][v] = weight;
}

void Graph::addVert(int u, int z, double w){
    vertices.push_back(Node(-1,w, u, z));
    wTrack[u][z] = w;
}

int Graph::vCount() const {
    return verts;
}

const vector<pair<int, double>> Graph::adjVert(int u) const {
    vector<pair<int, double>> result;
    

    for (const auto& pair : adjList[u]) {
        result.push_back(pair);
    }

    return result;
}

struct inn{
    int V1 = -1;
    int V2 = -1;
    double weight = -100.0;
    inn(int V1, int V2, double w) : V1{V1}, V2{V2}, weight{w} {}
};

void Prim(Graph& G, int r) {
    std::chrono::steady_clock::time_point start = std::chrono::steady_clock::now();
    int V = G.vCount();
    vector<int> parent(V, -1);
    vector<double> key(V, numeric_limits<double>::max());
    vector<bool> inMST(V, false);
    vector<inn> innmst;
    double totDist = 0.0;
    
    for (int i = 0; i < V; i++){
        innmst.push_back(inn(-1,-1,-1));
    }

    PriorityQueue Q(V);

    for (int u = 0; u < V; ++u) {
        Q.insert(Node(u, 100000));
    }

    key[r] = 0;



    while (!Q.isEmpty()) {
        Node u = Q.extract();
        inMST[u.ID] = true;

        const vector<pair<int, double>> adjVertices = G.adjVert(u.ID);
        for (const pair<int, double>& edge : adjVertices) {
            int v = edge.first;
            double weight = edge.second;

            if (!inMST[v] && weight < key[v]) {
                parent[v] = u.ID;
                key[v] = weight;
                innmst[v - 1] = inn(u.ID,v,weight);
                Q.changeKey(v, weight);
            }
        }
    }

    Graph G2(V);
    for (int i = 0; i < innmst.size() - 1; i++){
        G2.addE(innmst[i].V1,innmst[i].V2,innmst[i].weight);
        G2.addE(innmst[i].V2,innmst[i].V1,innmst[i].weight);
        totDist += innmst[i].weight;
    }
    std::chrono::steady_clock::time_point finish = std::chrono::steady_clock::now();
    ofstream OFile("primOut" + to_string(finc) +  ".txt");
    if (OFile.is_open()) {
        OFile << fixed << setprecision(1);
        vector<pair<int, double>> edges;
        for (int i = 0; i < V; i++){
            OFile << i << " ";
            edges = G2.adjVert(i);
            for (const pair<int, double>& edge : edges) {
                OFile << edge.first << " " << edge.second << " ";
            }
            OFile << endl;
        }
    OFile << totDist << endl;
    OFile << "Prim Time " + to_string(finc) + " = " << std::chrono::duration_cast<std::chrono::microseconds>(finish - start).count() << "[µs]" << std::endl;
    OFile << "Prim Time " + to_string(finc) + " = " << std::chrono::duration_cast<std::chrono::seconds> (finish - start).count() << "[s]" << std::endl;
    OFile.close();
    finc += 1;
    } 
    else {
        cerr << "Error creating file" << endl;
        return;
    }

}

void ApproxTour(Graph& G, int r) {
    std::chrono::steady_clock::time_point startt = std::chrono::steady_clock::now();
    int V = G.vCount();
    vector<int> parent(V, -1);
    vector<double> key(V, numeric_limits<double>::max());
    vector<bool> inMST(V, false);
    vector<inn> innmst;
    vector<int> tour;
    
    for (int i = 0; i < V; i++){
        innmst.push_back(inn(-1,-1,-1));
    }

    key[r] = 0;
    int cur = 0;
    int next = 0;
    tour.push_back(0);
    double totpd = 0.0;
    vector<vector<double>> wmst;
    wmst = vector<vector<double>>(V + 1, vector<double>(V + 1, 0));
    inMST[0] = true;
    vector<pair<int, double>> adjVertices;

    for(int i = 0; i < V - 1; i++){
        double min = numeric_limits<double>::max();


        adjVertices = G.adjVert(cur);
        for (const pair<int, double>& edge : adjVertices) {
            if(edge.second < min & !inMST[edge.first]){
                min = edge.second;
                next = edge.first;
            }
        }
        wmst[cur][next] = min;
        wmst[next][cur] = min;
        inMST[next] = true;
        cur = next;
        tour.push_back(cur);
        totpd += min;
    }
    adjVertices = G.adjVert(cur);
    bool checkp = false;
    int checkpc = 0;
    while(checkp == false){
        if(adjVertices[checkpc].first == 0){
            totpd += adjVertices[checkpc].second;
            checkp = true;
        }
    }

    for (int i = 0; i < V ; i++){
        cout << tour[i] << " ";
    }
    cout << 0 << endl << totpd << endl;

    std::chrono::steady_clock::time_point finisht = std::chrono::steady_clock::now();
    ofstream OFile("approxTour" + to_string(finc3) +  ".txt");
    if (OFile.is_open()) {
        OFile << fixed << setprecision(1);
        OFile << "Approximate Tour: ";
        for (int i = 0; i < V ; i++){
            OFile << tour[i] << " ";
        }
        OFile << 0 << endl << "Distance: " << totpd << endl;
        OFile << "Approximate Tour Time " + to_string(finc) + " = " << std::chrono::duration_cast<std::chrono::microseconds>(finisht - startt).count() << "[µs]" << std::endl;
        OFile << "Approximate Tour Time " + to_string(finc) + " = " << std::chrono::duration_cast<std::chrono::seconds> (finisht - startt).count() << "[s]" << std::endl;
        OFile.close();
    } 
    else {
        cerr << "Error creating file" << endl;
        return;
    }
    finc3 += 1;
}


double calculateDistance(int a, const int b, const vector<vector<double>>& distances) {


    double distance = distances[a][b];
    return distance;
}

double calculateTourDistance(const vector<int>& tour, const vector<Node>& cities, const vector<vector<double>>& distances, double min) {
    double totalDistance = 0.0;

    for (size_t i = 0; i < tour.size() - 1; ++i) {
        totalDistance += distances[tour[i]][tour[i+1]];
        if (totalDistance > min){
            return totalDistance;
        }
    }

    totalDistance += calculateDistance(tour.back(), tour.front(), distances);

    return totalDistance;
}

void bruteForceTSP(string graphName) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
    struct eHold{
        int V1 = -1;
        int V2 = -1;
        double weight = - 100.0;
        int inner = 0;

        eHold(int V1,int V2,double weight) : V1(V1), V2(V2), weight(weight) {}
    };

    ifstream input_file(graphName);
    if (!input_file.is_open()) {
        cerr << "Error: Unable to open file" << endl;
    }

    string line;
    int counter = 0;
    int tot = 0;
    vector<eHold> E1;
    while (getline(input_file, line)) {
        istringstream iss(line);
        int V1 = -1;
        int V2 = -1;
        double weight = -100.0;

        iss >> V1;

        int eCount =  0;

        while (iss >> V2 >> weight){
            E1.push_back(eHold(V1,V2,weight));
            tot++;
            eCount++;
        }
        E1[counter].inner = eCount;
        counter++;
    }
    input_file.close();

    Graph graph(counter); 

    cout << fixed << setprecision(1);
    for (int i = 0; i < tot; i++) {
        graph.addE(E1[i].V1,E1[i].V2, E1[i].weight);
    }

    int numCities = graph.vCount();
    vector<int> currentTour(numCities);
    vector<int> bestTour(numCities);
    double minDistance = numeric_limits<double>::infinity();
    vector<vector<double>> distances(numCities, vector<double>(numCities, 0.0));

    for (int i = 0; i < numCities; ++i) {
        currentTour[i] = i;
    }
    for (int i = 0; i < numCities; ++i) {
        for (int j = i + 1; j < numCities; ++j) {
            distances[i][j] = graph.wTrack[i][j];
            distances[j][i] = graph.wTrack[j][i];
        }
    }

    do {
        double currentDistance = calculateTourDistance(currentTour, graph.vertices, distances, minDistance);

        if (currentDistance < minDistance) {
            minDistance = currentDistance;
            bestTour = currentTour;
        }
    } while (next_permutation(currentTour.begin() + 1, currentTour.end())); 

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();




    ofstream OFile("bruteOut" + to_string(finc2) + ".txt");
    if (OFile.is_open()) {
        OFile << "Best Tour: ";
        for (int city : bestTour) {
            OFile << city << " ";
        }
        OFile << "\nDistance: " << minDistance << endl;
        OFile << "Brute Time " + to_string(finc2) + " = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[µs]" << std::endl;
        OFile << "Brute Time " + to_string(finc2) + " = " << std::chrono::duration_cast<std::chrono::seconds> (end - begin).count() << "[s]" << std::endl;
        finc2 += 1;
    } 
    else {
        cerr << "Error creating file" << endl;
        return;
    }

    Prim(graph, 0);
    ApproxTour(graph, 0);
    cout << "Finished Graph " + to_string(finc2 - 1) << endl;
}

#endif

