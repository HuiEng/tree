// There are two kinds of ktree nodes- branch nodes and leaf nodes
// Both contain a signature matrix, plus their own signature
// (the root node signature does not matter and can be blank)
// Branch nodes then contain 'order' links to other nodes
// Leaf nodes do not.
// However, as leaf nodes may become branch nodes, we allocate
// the space anyway.
// As the space to be used is determined at runtime, we use
// parallel arrays, not structs

#ifndef INCLUDE_TREAP2_HPP
#define INCLUDE_TREAP2_HPP

#include <omp.h>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include "bloom_filter.hpp"
#include "read.hpp"

using namespace std;

// typedef unsigned char cell_type;
bloom_parameters parameters;
size_t signatureSize; // Signature size (depends on element in BF, obtained while read binary)
size_t maxWindow = 0;
size_t partree_order = 10;
size_t HD_threshold = 0;
size_t HD_t = 90;
vector<vector<cell_type>> rep;
double match_t = 0.85;

void toBinaryIdx(FILE *stream, cell_type letter, int offset)
{
    int binary[bits_per_char];
    for (int n = 0; n < bits_per_char; n++)
        binary[bits_per_char - 1 - n] = (letter >> n) & 1;

    for (int n = 0; n < bits_per_char; n++)
    {
        if (binary[n] > 0)
        {
            fprintf(stream, "%d,", n + offset);
        }
    }
}

void printBF(FILE *stream, vector<cell_type> bf, bool idxOnly = true)
{
    for (int i = 0; i < bf.size(); i++)
    {
        toBinaryIdx(stream, bf[i], i * bits_per_char);
    }
    fprintf(stream, "\n");
}

// aggregate list of windows (BFs) to a set of minimisers within the seq (regardless of windows)
vector<cell_type> do_union(vector<cell_type> a, vector<cell_type> b)
{
    vector<cell_type> minimisers;
    for (int i = i; i < min(a.size(), b.size()); i++)
    {
        minimisers.push_back(a[i] | b[i]);
    }
    return minimisers;
}

// aggregate list of windows (BFs) to a set of minimisers within the seq (regardless of windows)
vector<cell_type> getMinimiserSet(vector<vector<cell_type>> a)
{
    vector<cell_type> minimisers = a[0];
    for (int i = i; i < a.size(); i++)
    {
        for (int j = 0; j < a[i].size(); j++)
        {
            minimisers[j] |= a[i][j];
        }
    }
    return minimisers;
}

size_t calcHD(vector<cell_type> a, vector<cell_type> b)
{
    size_t c = 0;
    for (size_t i = 0; i < a.size(); i++)
    {
        c += __builtin_popcountll(a[i] ^ b[i]);
    }
    return c;
}

size_t calcInter(vector<cell_type> a, vector<cell_type> b)
{
    size_t c = 0;
    for (size_t i = 0; i < a.size(); i++)
    {
        c += __builtin_popcountll(a[i] & b[i]);
    }
    return c;
}

size_t countSetBits(vector<cell_type> a)
{
    size_t c = 0;
    for (size_t i = 0; i < a.size(); i++)
    {
        c += __builtin_popcountll(a[i]);
    }
    return c;
}

size_t compareSeqsHD(vector<vector<cell_type>> a, vector<vector<cell_type>> b)
{
    size_t distance = 0;
    for (int i = 0; i < min(a.size(), b.size()); i++)
    {
        distance += calcHD(a[i], b[i]);
    }
    return distance;
}

size_t compareSeqsWindows(vector<vector<cell_type>> mean, vector<vector<cell_type>> x, size_t hd_threshold = HD_threshold)
{
    size_t matchingWindows = 0;
    size_t minSize = min(mean.size(), x.size());
    if (minSize == 0)
        return 0;

    for (int w = 0; w < minSize; w++) // compare window by window
    {
        size_t setBits = countSetBits(x[w]);
        size_t threshold = min(hd_threshold, setBits);
        size_t distance = calcInter(mean[w], x[w]);
        if (distance > threshold)
        {
            matchingWindows++;
        }
    }

    // return matchingWindows;
    return matchingWindows;
}

int calcPriority(vector<vector<cell_type>> a, vector<vector<cell_type>> b)
{
    // return rand() % 100; //-compareSeqsHD(a, b);
    return 0;
}

// return true if at least half the windows match
bool matchThreshold(vector<vector<cell_type>> mean, vector<vector<cell_type>> x, size_t hd_thres = HD_threshold, double threshold = 0.5)
{
    size_t minSize = min(mean.size(), x.size());
    // fprintf(stderr, "(%zu,%zu,", mean.size(), mean.size());
    if (minSize == 0)
        return false;

    size_t matchingWindows = compareSeqsWindows(x, mean, hd_thres);

    // fprintf(stderr, "%zu, %zu >= %f\n", hd_thres, matchingWindows, minSize * threshold);

    // return matchingWindows;
    return matchingWindows >= minSize * threshold;
}

bool firstMatch(vector<vector<cell_type>> mean, vector<vector<cell_type>> x, size_t w = 0)
{
    size_t setBits = countSetBits(x[w]);
    // size_t threshold = min(HD_threshold, setBits);
    size_t distance = calcInter(mean[w], x[w]);
    return distance >= setBits;
}

// return true if at least half the windows match
size_t compareSeqsInter(vector<vector<cell_type>> mean, vector<vector<cell_type>> x, double threshold = 0.5)
{
    size_t minSize = min(mean.size(), x.size()); // min(countSetBits(mean), countSetBits(x));
    if (minSize == 0)
        return false;

    size_t distance = 0;
    size_t mean_bits = 0;
    size_t x_bits = 0;

    for (int w = 0; w < minSize; w++) // compare window by window
    {
        distance += calcInter(mean[w], x[w]);
        mean_bits += countSetBits(mean[w]);
        x_bits += countSetBits(x[w]);
    }

    // return matchingWindows;
    return distance;
}

// return true if at least half the windows match
bool goLeft(vector<vector<cell_type>> mean, vector<vector<cell_type>> x)
{
    // for (int i = 0; i < x.size(); i++) // compare window by window
    // {
    //     for (int j = 0; j < mean.size(); j++)
    //     {
    //         if (x[i] == mean[j])
    //         // if (calcInter(x[i],mean[j])>3)
    //         {
    //             return true;
    //         }
    //     }
    // }
    // // return matchingWindows;
    // return false;

    vector<cell_type> a = getMinimiserSet(mean);
    vector<cell_type> b = getMinimiserSet(x);
    size_t inter = calcInter(a, b);
    size_t max = countSetBits(b);

    // fprintf(stderr, "~%zu,%zu\n", inter, max);

    return inter > max * match_t;
}

// return true if at least half the windows match
bool interThreshold(vector<vector<cell_type>> mean, vector<vector<cell_type>> x, double threshold = 0.5)
{
    size_t minSize = min(mean.size(), x.size()); // min(countSetBits(mean), countSetBits(x));
    if (minSize == 0)
        return false;

    size_t distance = 0;
    size_t mean_bits = 0;
    size_t x_bits = 0;

    for (int w = 0; w < minSize; w++) // compare window by window
    {
        distance += calcInter(mean[w], x[w]);
        mean_bits += countSetBits(mean[w]);
        x_bits += countSetBits(x[w]);
    }

    // return matchingWindows;
    return distance >= min(mean_bits, x_bits) * threshold;
}

// https://www.tutorialspoint.com/cplusplus-program-to-implement-treap
struct TreapNod2
{ // node declaration
    vector<vector<cell_type>> data;
    vector<cell_type> flatten; // union
    int priority;
    size_t id;
    TreapNod2 *l, *r, *parent;
    TreapNod2(vector<vector<cell_type>> d, size_t i)
    { // constructor
        this->data = d;
        this->flatten = getMinimiserSet(d);
        this->priority = rand() % 100; // calcPriority(parent->data, d); //compareSeqsWindows(rep, d); // compareSeqsInter(rep, d); //-i; // rand() % 100;
        // this->priority = 0;
        this->l = this->r = this->parent = nullptr;
        this->id = i;
    }
};
// void rotLeft(TreapNod2 *&root)
// { // left rotation
//     TreapNod2 *R = root->r;
//     TreapNod2 *X = root->r->l;

//     //?
//     TreapNod2 *P = root->parent;
//     root->parent = R;
//     R->l = root;
//     root->r = X;
//     root = R;
//     root->parent = P; //?
// }
// void rotRight(TreapNod2 *&root)
// { // right rotation
//     TreapNod2 *L = root->l;
//     TreapNod2 *Y = root->l->r;
//     //?
//     TreapNod2 *P = root->parent;
//     root->parent = L;
//     // Y->parent->data = L->data;
//     // root->parent->data = Y->data;
//     // fprintf(stderr,"???%zu,%zu,%zu---", root->parent->id, root->id, root->l->id);
//     //
//     L->r = root;
//     // L->r->parent = L;
//     // Y->parent = root->parent;
//     root->l = Y;
//     root = L;
//     root->parent = P; //?
//     // fprintf(stderr,"???%zu,%zu,%zu\n", root->parent->id, root->id, root->r->id);
// }
void rotLeft(TreapNod2 *&root)
{ // left rotation
    TreapNod2 *R = root->r;
    TreapNod2 *X = root->r->l;
    TreapNod2 *parent = root->parent;
    root->parent = R;
    if (X != nullptr)
    {
        X->parent = root;
        // X->priority = calcPriority(X->parent->data, X->data);
    }
    // root->priority = calcPriority(parent->data, root->data);
    R->l = root;
    root->r = X;
    root = R;
    root->parent = parent;
    // root->priority = calcPriority(parent->data, root->data);
}
void rotRight(TreapNod2 *&root)
{ // right rotation
    TreapNod2 *L = root->l;
    TreapNod2 *Y = root->l->r;
    TreapNod2 *parent = root->parent;
    root->parent = L;
    if (Y != nullptr)
    {
        Y->parent = root;
        // Y->priority = calcPriority(Y->parent->data, Y->data);
    }
    // root->priority = calcPriority(parent->data, root->data);
    L->r = root;
    root->l = Y;
    root = L;
    root->parent = parent;
    // root->priority = calcPriority(parent->data, root->data);
}
void insertNod(TreapNod2 *&root, vector<vector<cell_type>> d, size_t idx)
{ // insertion
    if (root == nullptr)
    {
        root = new TreapNod2(d, idx);
        return;
    }
    if (matchThreshold(d, root->data))
    {
        // if similar to exiting, do nothing
        return;
    }
    else if (root->l != nullptr && root->r != nullptr)
    {
        // if got children, go to the closer one
        size_t l = interThreshold(d, root->l->data, match_t);
        size_t r = interThreshold(d, root->r->data, match_t);

        if (l >= r)
        {
            insertNod(root->l, d, idx);
            if (root->l != nullptr && root->l->priority > root->priority)
                rotRight(root);
        }
        else
        {
            insertNod(root->r, d, idx);
            if (root->r != nullptr && root->r->priority > root->priority)
                rotLeft(root);
        }
    }
    // else if (root->l != nullptr || root->r != nullptr){
    //     // if only 1 child, decide whether go down or spawn another child

    // }

    // size_t hd = compareSeqsHD(d, root->data);
    // fprintf(stderr, "???%zu,%zu\n", root->id, hd);
    // if (hd < HD_t)
    else if (interThreshold(d, root->data, match_t))
    {
        insertNod(root->l, d, idx);
        if (root->l != nullptr && root->l->priority > root->priority)
            rotRight(root);
    }
    else
    {
        insertNod(root->r, d, idx);
        if (root->r != nullptr && root->r->priority > root->priority)
            rotLeft(root);
    }
}

void insertNod2(TreapNod2 *&root, TreapNod2 *&parent, vector<vector<cell_type>> d, size_t idx)
{ // insertion
    if (root == nullptr)
    {
        root = new TreapNod2(d, idx);
        root->parent = parent;
        // fprintf(stderr, "***%zu", root->id);
        return;
    }
    if (interThreshold(d, root->data))
    {
        insertNod2(root->l, root, d, idx);
        if (root->l != nullptr && root->l->priority > root->priority)
            rotRight(root);
    }
    else
    {
        insertNod2(root->r, root, d, idx);
        if (root->r != nullptr && root->r->priority > root->priority)
            rotLeft(root);
    }
}

void insertNod3(TreapNod2 *&root, TreapNod2 *&parent, vector<vector<cell_type>> d, size_t idx)
{ // insertion
    if (root == nullptr)
    {
        root = new TreapNod2(d, idx);
        root->parent = parent;
        // fprintf(stderr, "***%zu", root->id);
        // fprintf(stderr, "\n");
        if (parent != nullptr)
        {
            root->priority = calcPriority(parent->data, root->data);
        }
        return;
    }
    // fprintf(stderr, "%zu,", root->id);
    if (matchThreshold(d, root->data, HD_threshold, match_t))
    {
        return;
    }
    else if (goLeft(d, root->data))
    {
        insertNod3(root->l, root, d, idx);
        // if (root->l != nullptr && root->l->priority < calcPriority(root->parent->data, root->l->data))
        if (root->l != nullptr && root->l->priority > root->priority)
        {
            fprintf(stderr, "<%zu,", root->id);
            rotRight(root);
        }
    }
    else
    {
        insertNod3(root->r, root, d, idx);
        // // fprintf(stderr,"<%d,%d>",root->r->priority,calcPriority(root->parent->data,root->r->data));
        //?
        // if (root->r != nullptr && root->r->priority < calcPriority(root->parent->data, root->r->data))
        if (root->r != nullptr && root->r->priority > root->priority)
        {
            fprintf(stderr, ">%zu,", root->id);
            rotLeft(root);
        }
    }
}


// if there are 2 children, return true if the left child is closer to query
bool goLeft2(TreapNod2 *node, vector<vector<cell_type>> x)
{
    if (node->l != nullptr && node->r != nullptr)
    { // if 2 children, return true if left is closer
        vector<cell_type> n = getMinimiserSet(x);
        vector<cell_type> l = getMinimiserSet(node->l->data);
        vector<cell_type> r = getMinimiserSet(node->r->data);
        // vector<cell_type> l = node->l->flatten;
        // vector<cell_type> r = node->r->flatten;
        size_t inter_l = calcInter(n, l);
        size_t inter_r = calcInter(n, r);
        return inter_l >= inter_r;
    }
    else if (node->l == nullptr)
    { 
        return false;
    }
    else{
        return true;
    }




    // if (node->l == nullptr)
    // { // if no child, fill the left first
    //     return true;
    // }
    // else if (node->r != nullptr)
    // { // if 2 children, return true if left is closer
    //     vector<cell_type> n = getMinimiserSet(x);
    //     // vector<cell_type> l = getMinimiserSet(node->l->data);
    //     // vector<cell_type> r = getMinimiserSet(node->r->data);
    //     vector<cell_type> l = node->l->flatten;
    //     vector<cell_type> r = node->r->flatten;
    //     size_t inter_l = calcInter(n, l);
    //     size_t inter_r = calcInter(n, r);
    //     return inter_l >= inter_r;
    // }
    // else
    // { 
    //     // // if 1 child, proceed with left if it is "close enough"
    //     // return goLeft(node->l->data, x);

    //     // if 1 child, fill the right
    //     return false;
    // }
}


void insertNod4(TreapNod2 *&root, TreapNod2 *&parent, vector<vector<cell_type>> d, size_t idx, int temp_p)
{ // insertion
    if (root == nullptr)
    {
        root = new TreapNod2(d, idx);
        root->parent = parent;
        if (parent != nullptr)
        {
            // root->priority = calcPriority(parent->data, root->data);
            root->priority = temp_p;
        }
        return;
    }

    root->flatten = do_union(root->flatten, getMinimiserSet(d));
    // root->priority -= compareSeqsHD(root->data,d);

    // fprintf(stderr, "%zu,", root->id);
    if (matchThreshold(d, root->data, HD_threshold, match_t))
    {
        return;
    }
    else if (goLeft(d, root->data))
    // else if (goLeft2(root, d))
    {
        // if (root->priority <= temp_p)
        // // if (root->l!=nullptr && !goLeft(d, root->l->data))
        // {
        //     fprintf(stderr, "<<<%zu,", root->id);

        //     TreapNod2 *node = new TreapNod2(d, idx);
        //     node->priority = temp_p;
        //     TreapNod2 *T = root;
        //     node->parent = T->parent;
        //     node->l = T;
        //     T->parent = node;
        //     root = node;
        //     return;
        // }
        insertNod4(root->l, root, d, idx, temp_p);
    }
    else
    {
        // if (root->priority <= temp_p)
        // // if (root->r!=nullptr && !goLeft(d, root->r->data))
        // {
        //     fprintf(stderr, "<<<%zu,", root->id);

        //     TreapNod2 *node = new TreapNod2(d, idx);
        //     node->priority = temp_p;
        //     TreapNod2 *T = root;
        //     node->parent = T->parent;
        //     node->r = T;
        //     T->parent = node;
        //     root = node;
        //     return;
        // }
        insertNod4(root->r, root, d, idx, temp_p);
    }
}

void deleteNod(TreapNod2 *&root, vector<vector<cell_type>> key)
{
    if (root == nullptr)
        return;
    if (root->data == key)
    {
        if (root->l == nullptr && root->r == nullptr)
        {
            delete root;
            root = nullptr;
        }
        else if (root->l && root->r)
        {
            if (root->l->priority < root->r->priority)
            {
                rotLeft(root);
                deleteNod(root->l, key);
            }
            else
            {
                rotRight(root);
                deleteNod(root->r, key);
            }
        }
        // node to be deleted has only one child
        else
        {
            TreapNod2 *child = (root->l) ? root->l : root->r;
            TreapNod2 *curr = root;
            root = child;
            delete curr;
        }
        return;
    }
    if (matchThreshold(key, root->data))
        deleteNod(root->l, key);
    deleteNod(root->r, key);
}

void displayTreap(TreapNod2 *root, int space = 0, int height = 10)
{ // display treap
    if (root == nullptr)
        return;
    space += height;
    displayTreap(root->l, space);
    cout << endl;
    for (int i = height; i < space; i++)
        cout << ' ';
    // cout << root->data << "(" << root->priority << ")\n";
    cout << endl;
    displayTreap(root->r, space);
}

void printTreapNode(TreapNod2 *root)
{
    fprintf(stdout, "{\"node\":\"%zu(%d)", root->id, root->priority);
    // fprintf(stdout, "{\"node\":\"%zu(%zu;%d)", root->id, root->parent->id, -root->priority);
    // fprintf(stdout, "{\"node\":\"%zu", root->id);
    fprintf(stdout, "\",\"children\":[");
}

void printTreap(TreapNod2 *root)
{
    if (root == nullptr)
    {
        fprintf(stdout, "{}");
        return;
    }

    printTreapNode(root);
    printTreap(root->l);
    fprintf(stdout, ",");
    printTreap(root->r);
    fprintf(stdout, "]}");
}

void findNod(TreapNod2 *root, vector<vector<cell_type>> d, size_t &best, size_t score_m, size_t score_hd)
{
    if (root == nullptr)
    {
        // fprintf(stderr, "---%zu\n", best->id);
        return;
    }
    // fprintf(stderr, "---%zu\n", root->id);
    size_t matching = compareSeqsWindows(d, root->data);
    size_t hd = compareSeqsHD(d, root->data);
    if (matching >= score_m)
    {
        // best.push_back(root->id);

        score_m = matching;
        if (hd <= score_hd)
        {
            best = root->id;
            score_hd = hd;
        }
    }

    findNod(root->l, d, best, score_m, score_hd);
    findNod(root->r, d, best, score_m, score_hd);
}

size_t findNod(TreapNod2 *root, vector<vector<cell_type>> d)
{
    if (root == nullptr)
    {
        fprintf(stderr, "***NA\n");
        return 55555;
    }
    if (matchThreshold(d, root->data, HD_threshold, match_t))
    {
        // fprintf(stderr, "f\n");
        return root->id;
    }
    else if (goLeft(d, root->data))
    //else if (goLeft2(root,d))
    {
        fprintf(stderr, "<%zu,",root->id);
        return findNod(root->l, d);
    }
    else
    {
        fprintf(stderr, ">%zu,",root->id);
        return findNod(root->r, d);
    }
}

size_t findNod2(TreapNod2 *root, vector<vector<cell_type>> d)
{
    if (root == nullptr)
    {
        // fprintf(stderr, "***NA\n");
        return 55555;
    }
    if (matchThreshold(d, root->data, HD_threshold, match_t))
    {
        // fprintf(stderr, "f\n");
        return root->id;
    }
    else if (goLeft2(root, d))
    {
        // fprintf(stderr, "<,");
        return findNod2(root->l, d);
    }
    else
    {
        // fprintf(stderr, ">,");
        return findNod2(root->r, d);
    }
}

void findNodCond(TreapNod2 *root, vector<vector<cell_type>> d, int &best, size_t score_m)
{
    if (root == nullptr)
    {
        // fprintf(stderr, "---%zu\n", best->id);
        return;
    }
    // fprintf(stderr, "---%zu\n", root->id);
    size_t matching = compareSeqsWindows(d, root->data);
    if (matching >= score_m)
    {
        // best.push_back(root->id);
        best = root->id;
        score_m = matching;
    }

    findNodCond(root->l, d, best, score_m);
    findNodCond(root->r, d, best, score_m);
}

// int treaptest()
// {
//     int nums[] = {1, 7, 6, 4, 3, 2, 8, 9, 10};
//     int a = sizeof(nums) / sizeof(int);
//     TreapNod2 *root = nullptr;
//     srand(time(nullptr));
//     for (int n : nums)
//         insertNod(root, n);
//     cout << "Constructed Treap:\n\n";
//     displayTreap(root);
//     cout << "\nDeleting node 8:\n\n";
//     deleteNod(root, 8);
//     displayTreap(root);
//     cout << "\nDeleting node 3:\n\n";
//     deleteNod(root, 3);
//     displayTreap(root);

//     fprintf(stdout, "\nvar treeData=");
//     printTreap(root);
//     return 0;
// }

class treap
{
public:
    TreapNod2 *root;
    size_t capacity = 0; // Set during construction, currently can't change

    void reserve(size_t capacity)
    {
        // For safety, only call this at startup currently
        if (this->capacity != 0)
        {
            fprintf(stderr, "Reserve can only be called from 0 capacity\n");
            exit(1);
        }
        this->capacity = capacity;

        // //#pragma omp parallel
        // {
        //     //#pragma omp single
        //     {
        //         childCounts.resize(capacity);
        //     }

        // }
    }

    treap(size_t capacity, size_t threshold, vector<vector<cell_type>> key)
    {
        reserve(capacity);
        HD_threshold = threshold;
        root = nullptr;
        // root = new TreapNod2(key, 0);
        // root->priority = 1000;
        // root->parent = new TreapNod2(key, 0);
        rep = key;
    }

    void insert(vector<vector<cell_type>> d, size_t idx)
    {
        // int best = -1;
        // findNodCond(root, n, best, n.size() * match_t);
        // if (best == -1)
        // {
        //     insertNod2(root, root, n, idx);
        // }
        // insertNod3(root, root, n, idx);
        size_t dest = findNod(root, d);
        if (dest == 55555)
        {
            TreapNod2 *parent = nullptr;
            insertNod4(root, parent, d, idx, rand() % 100);
        }
    }

    void printTreeJson()
    {
        fprintf(stdout, "\nvar treeData=");
        printTreap(root);
    }
};

#endif