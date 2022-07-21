// There are two kinds of ktree nodes- branch nodes and leaf nodes
// Both contain a signature matrix, plus their own signature
// (the root node signature does not matter and can be blank)
// Branch nodes then contain 'order' links to other nodes
// Leaf nodes do not.
// However, as leaf nodes may become branch nodes, we allocate
// the space anyway.
// As the space to be used is determined at runtime, we use
// parallel arrays, not structs

#ifndef INCLUDE_TREAP_HPP
#define INCLUDE_TREAP_HPP

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

size_t compareSeqsWindows(vector<vector<cell_type>> mean, vector<vector<cell_type>> x)
{
    size_t matchingWindows = 0;
    size_t minSize = min(mean.size(), x.size());
    if (minSize == 0)
        return 0;

    for (int w = 0; w < minSize; w++) // compare window by window
    {
        size_t setBits = countSetBits(x[w]);
        size_t threshold = min(HD_threshold, setBits);
        size_t distance = calcInter(mean[w], x[w]);
        if (distance > threshold)
        {
            matchingWindows++;
        }
    }

    // return matchingWindows;
    return matchingWindows;
}

// return true if at least half the windows match
bool matchThreshold(vector<vector<cell_type>> mean, vector<vector<cell_type>> x, double threshold = 0.5)
{
    size_t minSize = min(mean.size(), x.size());
    if (minSize == 0)
        return false;

    size_t matchingWindows = compareSeqsWindows(mean, x);

    // return matchingWindows;
    return matchingWindows >= minSize * threshold;
}

// https://www.tutorialspoint.com/cplusplus-program-to-implement-treap
struct TreapNod
{ // node declaration
    vector<vector<cell_type>> data;
    int priority;
    size_t id;
    TreapNod *l, *r;
    TreapNod(vector<vector<cell_type>> d, size_t i)
    { // constructor
        this->data = d;
        this->priority = -i; // rand() % 100;
        this->l = this->r = nullptr;
        this->id = i;
    }
};
void rotLeft(TreapNod *&root)
{ // left rotation
    TreapNod *R = root->r;
    TreapNod *X = root->r->l;
    R->l = root;
    root->r = X;
    root = R;
}
void rotRight(TreapNod *&root)
{ // right rotation
    TreapNod *L = root->l;
    TreapNod *Y = root->l->r;
    L->r = root;
    root->l = Y;
    root = L;
}
void insertNod(TreapNod *&root, vector<vector<cell_type>> d, size_t idx)
{ // insertion
    if (root == nullptr)
    {
        root = new TreapNod(d, idx);
        return;
    }
    // if (matchThreshold(d, root->data))
    size_t hd = compareSeqsHD(d, root->data);
    fprintf(stderr,"???%zu,%zu\n", root->id, hd);
    if (hd < HD_t)
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

bool searchExactNod(TreapNod *root, vector<vector<cell_type>> key)
{
    if (root == nullptr)
        return false;
    if (root->data == key)
        return true;
    if (matchThreshold(key, root->data))
        return searchExactNod(root->l, key);
    return searchExactNod(root->r, key);
}

TreapNod *searchNod(TreapNod *root, TreapNod *node, vector<vector<cell_type>> key)
{
    if (node == nullptr)
        return nullptr;
    fprintf(stderr, "%zu,", node->id);
    if (matchThreshold(key, node->data, 0.5))
        return node;
    if (compareSeqsHD(key, node->data) < HD_t)
        // if (matchThreshold(key, node->data, 0.3))
        return searchNod(root, node->l, key);
    return searchNod(root, node->r, key);
}

TreapNod *findBest(size_t best, TreapNod *final, TreapNod *node, vector<vector<cell_type>> key)
{
    if (node == nullptr)
    {
        // fprintf(stderr, "\n*** %zu\n", final->id);
        // return nullptr;
        return final;
    }

    size_t match = compareSeqsWindows(key, node->data);
    size_t max = min(key.size(), node->data.size());

    if (match > best)
    {
        // fprintf(stderr, "---%zu,%zu,%zu\n", match, best, node->id);
        best = match;
        final = node;
    }

    if (match == max) // found match
    {
        // fprintf(stderr, "here\n");
        return node;
    }
    else if (compareSeqsHD(key, node->data) < HD_t)
    // else if (match<max*0.3)
    {
        // fprintf(stderr, ">%zu,%zu,%zu,%zu\n", match, best, node->id, final->id);
        return findBest(best, final, node->r, key);
    }

    else
    {
        // fprintf(stderr, "*%zu,%zu,%zu\n", match, best, node->id);
        return findBest(best, final, node->l, key);
    }
}

TreapNod *findBestHD(size_t best, TreapNod *final, TreapNod *node, vector<vector<cell_type>> key)
{
    if (node == nullptr)
    {
        // fprintf(stderr, "\n*** %zu\n", final->id);
        // return nullptr;
        return final;
    }

    size_t HD = compareSeqsHD(key, node->data);
    size_t max = min(key.size(), node->data.size());

    if (HD < best)
    {
        // fprintf(stderr, "---%zu,%zu,%zu\n", match, best, node->id);
        best = HD;
        final = node;
    }

    if (HD == 0) // found match
    {
        // fprintf(stderr, "here\n");
        return node;
    }
    else
    {
        // fprintf(stderr, ">%zu,%zu,%zu,%zu\n", match, best, node->id, final->id);
        return findBestHD(best, final, node->r, key);
    }

    // else
    // {
    //     fprintf(stderr, "*%zu,%zu,%zu\n", match, best, node->id);
    //     return findBest(best, final, node->l, key);
    // }
}

TreapNod *searchNodTest(TreapNod *root, vector<vector<cell_type>> key)
{
    if (root == nullptr)
        return nullptr;
    if (matchThreshold(key, root->data, 0.5))
        return root;
    if (matchThreshold(key, root->data, 0.25))
        return searchNodTest(root->l, key);
    return searchNodTest(root->r, key);
}

TreapNod *searchNodPrint(TreapNod *root, vector<vector<cell_type>> key)
{
    if (root == nullptr)
        return nullptr;
    fprintf(stderr, "%zu,", root->id);
    if (matchThreshold(key, root->data, 0.5))
    {
        // fprintf(stderr, "%zu,", root->id);
        return root;
    }

    if (matchThreshold(key, root->data, 0.4))
    {
        fprintf(stderr, "R");
        return searchNodPrint(root->l, key);
    }
    fprintf(stderr, "L");
    return searchNodPrint(root->r, key);
}

// // depth first search
// void searchAll(TreapNod *root)
// {
//     if (root == nullptr)
//     {
//         return;
//     }

//     fprintf(stderr, "%zu,", root->id);

//     searchAll(root->l);
//     searchAll(root->r);
// }

// depth first search
TreapNod *searchAll(size_t best, TreapNod *final, TreapNod *node, vector<vector<cell_type>> key)
{
    if (node == nullptr)
    {
        return final;
    }

    size_t match = compareSeqsWindows(key, node->data);
    size_t max = min(key.size(), node->data.size());
    fprintf(stderr, "---%zu,%zu,%zu\n", match, best, node->id);

    if (match >= best)
    {
        // fprintf(stderr, "---%zu,%zu,%zu\n", match, best, node->id);
        best = match;
        final = node;
    }
    return searchAll(best, final, node->l, key);
    return searchAll(best, final, node->r, key);

    // if (match == max) // found match
    // {
    //     // fprintf(stderr, "here\n");
    //     return node;
    // }
    // else if (compareSeqsHD(key, node->data) < HD_t)
    // // else if (match<max*0.3)
    // {
    //     // fprintf(stderr, ">%zu,%zu,%zu,%zu\n", match, best, node->id, final->id);
    //     return searchAll(best, final, node->r, key);
    // }

    // else
    // {
    //     // fprintf(stderr, "*%zu,%zu,%zu\n", match, best, node->id);
    //     return searchAll(best, final, node->l, key);
    // }
}

void deleteNod(TreapNod *&root, vector<vector<cell_type>> key)
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
            TreapNod *child = (root->l) ? root->l : root->r;
            TreapNod *curr = root;
            root = child;
            delete curr;
        }
        return;
    }
    if (matchThreshold(key, root->data))
        deleteNod(root->l, key);
    deleteNod(root->r, key);
}

void displayTreap(TreapNod *root, int space = 0, int height = 10)
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

void printTreapNode(TreapNod *root)
{
    fprintf(stdout, "{\"node\":\"%zu(%d)", root->id, root->priority);
    fprintf(stdout, "\",\"children\":[");
}

void printTreap(TreapNod *root)
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

// int treaptest()
// {
//     int nums[] = {1, 7, 6, 4, 3, 2, 8, 9, 10};
//     int a = sizeof(nums) / sizeof(int);
//     TreapNod *root = nullptr;
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
    TreapNod *root;
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

    treap(size_t capacity, size_t threshold)
    {
        reserve(capacity);
        HD_threshold = threshold;
        root = nullptr;
    }

    void insert(vector<vector<cell_type>> n, size_t idx)
    {
        // TreapNod *dest = searchNod(root, root, n);
        TreapNod *dest = searchAll(n.size() / 2, nullptr, root, n);

        if (dest == nullptr)
        {
            insertNod(root, n, idx);
            fprintf(stderr, "\n>>> %zu\n", idx);
        }
    }

    void printTreeJson()
    {
        fprintf(stdout, "\nvar treeData=");
        printTreap(root);
    }
};

#endif