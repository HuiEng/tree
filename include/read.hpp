// There are two kinds of ktree nodes- branch nodes and leaf nodes
// Both contain a signature matrix, plus their own signature
// (the root node signature does not matter and can be blank)
// Branch nodes then contain 'order' links to other nodes
// Leaf nodes do not.
// However, as leaf nodes may become branch nodes, we allocate
// the space anyway.
// As the space to be used is determined at runtime, we use
// parallel arrays, not structs

#ifndef INCLUDE_READ_HPP
#define INCLUDE_READ_HPP

#include <iostream>
#include <vector>
#include <string>
#include <experimental/filesystem>
#include <fstream>

using namespace std;
namespace fs = std::experimental::filesystem;
typedef size_t sig_type;

unsigned long long int readSignatures(const string file, vector<cell_type> &sigs)
{
    ifstream rf(file, ios::out | ios::binary);

    // get length of file:
    rf.seekg(0, rf.end);
    int len = rf.tellg();
    len = len - (int)sizeof(unsigned long long int); // / sizeof(uint64_t);
    rf.seekg(0, rf.beg);

    unsigned long long int length;
    if (rf)
        rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));

    sigs.resize(len);
    size_t i = 0;
    while (rf)
    {
        rf.read((char *)&sigs[i], sizeof(cell_type));
        i++;
    }
    rf.close();

    return length;
}

// read all files in given folder
unsigned long long int readSignaturesMultiple(const string folder, vector<cell_type> &sigs)
{
    unsigned long long int length;
    // keep track of sequence index
    FILE *pFile = fopen((folder + "-seqs.txt").c_str(), "w");

    for (const auto &entry : fs::directory_iterator(folder))
    {
        if (entry.path().extension() == ".bin")
        {
            fprintf(pFile, "%s\n", entry.path().c_str());
            ifstream rf(entry.path(), ios::out | ios::binary);

            // get length of file:
            rf.seekg(0, rf.end);
            int len = rf.tellg();
            len = len - (int)sizeof(unsigned long long int); // / sizeof(uint64_t);
            rf.seekg(0, rf.beg);

            if (rf)
                rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));

            auto old_size = sigs.size();
            sigs.resize(old_size + len);
            size_t i = 0;
            while (rf)
            {
                rf.read((char *)&sigs[old_size + i], sizeof(cell_type));
                i++;
            }
            rf.close();
        }
    }
    return length;
}

// read all files from the input txt, each line is the full path of the file to read
unsigned long long int readSignaturesSpecific(const string ref_file_path, vector<cell_type> &sigs)
{
    unsigned long long int length;
    // keep track of sequence index
    std::ifstream ref_file(ref_file_path);
    string file;
    while (getline(ref_file, file))
    { // read data from file object and put it into string.
        ifstream rf(file.c_str(), ios::out | ios::binary);
        fprintf(stderr, "reading %s\n", file.c_str()); // << file << "\n"; //print the data of the string

        // get length of file:
        rf.seekg(0, rf.end);
        int len = rf.tellg();
        len = len - (int)sizeof(unsigned long long int); // / sizeof(uint64_t);
        rf.seekg(0, rf.beg);
        fprintf(stderr, "%d\n", len);

        if (rf)
            rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));
        fprintf(stderr, "here22\n");
        auto old_size = sigs.size();

        sigs.resize(old_size + len);
        size_t i = 0;
        fprintf(stderr, "here\n");
        while (rf)
        {
            rf.read((char *)&sigs[old_size + i], sizeof(cell_type));
            i++;
        }
        rf.close();
    }

    return 1;
}

vector<vector<sig_type>> readPartition(const string file_path)
{
    ifstream rf(file_path, ios::out | ios::binary);

    // get length of file:
    rf.seekg(0, rf.end);
    int len = rf.tellg();
    len = len - (int)sizeof(unsigned long long int); // / sizeof(uint64_t);
    rf.seekg(0, rf.beg);

    vector<vector<sig_type>> seqs;
    vector<sig_type> tseq;
    sig_type temp = 0;
    while (rf)
    {
        rf.read((char *)&temp, sizeof(sig_type));
        if (temp == -1)
        {
            seqs.push_back(tseq);
            tseq.clear();
        }
        else
        {
            tseq.push_back(temp);
        }
    }
    rf.close();

    return seqs;
}

bool isEmpty(vector<cell_type> bf)
{
    for (cell_type val : bf)
    {
        if (val != 0)
            return false;
    }
    return true;
}

vector<vector<vector<cell_type>>> readPartitionBF(const string file_path)
{
    ifstream rf(file_path, ios::out | ios::binary);

    // get length of file:
    rf.seekg(0, rf.end);
    int len = rf.tellg();
    len = len - (int)sizeof(unsigned long long int); // / sizeof(uint64_t);
    rf.seekg(0, rf.beg);

    unsigned long long int length;
    if (rf)
        rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));

    vector<vector<vector<cell_type>>> seqs;
    vector<vector<cell_type>> tseq; // list of BFs for a seq

    //? 1 window
    // cout << "length: " << length << "\n";
    vector<cell_type> bf(length);
    cell_type temp = 0;
    size_t i = 0;

    while (rf)
    {
        rf.read((char *)&bf[i], sizeof(cell_type));
        i++;
        if (i == length)
        {
            if (isEmpty(bf))
            {
                seqs.push_back(tseq);
                tseq.clear();
            }
            else
            {
                tseq.push_back(bf);
            }
            fill(bf.begin(), bf.end(), 0);
            i = 0;
        }
    }
    rf.close();

    return seqs;
}

#endif