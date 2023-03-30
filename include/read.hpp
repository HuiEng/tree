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
#include "bloom_filter.hpp"

using namespace std;
namespace fs = std::experimental::filesystem;
typedef size_t sig_type;
size_t max_fileSize = 1000000000;

// output signature size, if file too big => return empty "sigs", use readSignaturesBatch
unsigned long long int readSignatures(const string file, vector<cell_type> &sigs)
{
    ifstream rfSize(file, ios::binary | ios::ate);
    if (!rfSize.is_open())
    {
        fprintf(stderr, "Invalid File. Please try again\n");
        exit(0);
    }
    if (rfSize.tellg() > max_fileSize)
    {
        fprintf(stderr, "File too big, please read in batches\n");
        return 0;
    }
    ifstream rf(file, ios::out | ios::binary);

    // get length of file:
    rf.seekg(0, rf.end);
    unsigned long long int len = rf.tellg();
    len = len - sizeof(unsigned long long int); // / sizeof(uint64_t);
    rf.seekg(0, rf.beg);

    unsigned long long int length;
    if (rf)
        rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));

    try
    {
        sigs.resize(len);
    }
    catch (const std::exception &e) // caught by reference to base
    {
        // std::cout << " a standard exception was caught, with message '"
        //           << e.what() << "'\n";
        return length;
    }
    size_t i = 0;
    while (rf)
    {
        rf.read((char *)&sigs[i], sizeof(cell_type));
        i++;
    }
    rf.close();

    return length;
}

// if too large to use readSignatures, use batch, signatureSize should be outputed from the prev function
vector<vector<cell_type>> readSignaturesBatch(const string file, size_t size, size_t &signatureSize)
{
    vector<vector<cell_type>> sigs;
    ifstream rf(file, ios::out | ios::binary);
    if (!rf.is_open())
    {
        fprintf(stderr, "Invalid File. Please try again\n");
        exit(0);
    }

    unsigned long long int length;
    if (rf)
        rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));
    signatureSize = length;
    size = size * length;
    vector<cell_type> temp(size);
    size_t i = 0;

    while (rf)
    {
        rf.read((char *)&temp[i], sizeof(cell_type));
        i++;

        if (i == size)
        {
            sigs.push_back(temp);
            temp.clear();
            temp.resize(size);
            i = 0;
        }
    }
    temp.resize(i);
    sigs.push_back(temp);
    rf.close();

    return sigs;
}

vector<size_t> getIndices(size_t seqCount, size_t sampleSize)
{
    vector<size_t> indices(seqCount);
    iota(indices.begin(), indices.end(), 0);
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    shuffle(indices.begin(), indices.end(), default_random_engine(seed));
    indices.resize(sampleSize);
    sort(indices.begin(), indices.end(), std::greater<>());
    fprintf(stderr, "Random sampling %zu seqs...\n", sampleSize);

    return indices;
}

unsigned long long int readSignaturesSample(const string file, vector<cell_type> &sigs, size_t sampleSize)
{
    ifstream rf(file, ios::out | ios::binary);
    if (!rf.is_open())
    {
        fprintf(stderr, "Invalid File. Please try again\n");
        exit(0);
    }

    // get length of file:
    rf.seekg(0, rf.end);
    unsigned long long int len = rf.tellg();
    len = len - sizeof(unsigned long long int); // / sizeof(uint64_t);
    rf.seekg(0, rf.beg);

    unsigned long long int length;
    if (rf)
        rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));
    size_t seqCount = len / length;
    fprintf(stderr, "File contains %zu seqs...\n", seqCount);

    if (sampleSize > seqCount)
    {
        sigs.resize(len);
        size_t i = 0;
        while (rf)
        {
            rf.read((char *)&sigs[i], sizeof(cell_type));
            i++;
        }
    }
    else
    {
        vector<size_t> indices = getIndices(seqCount, sampleSize);
        // for (size_t i : indices)
        // {
        //     fprintf(stderr, "%zu,", i);
        // }
        // fprintf(stderr, "\n");

        size_t size = sampleSize * length;
        sigs.resize(size);
        size_t idx = indices.back() * length;
        indices.pop_back();
        size_t i = 0;
        size_t c = 0;
        cell_type t;
        while (rf)
        {
            if (i == idx)
            {
                // fprintf(stderr, "%zu,", idx / length);
                for (size_t n = 0; n < length; n++)
                {
                    rf.read((char *)&sigs[c], sizeof(cell_type));
                    c++;
                }
                i += length;
                idx = indices.back() * length;
                indices.pop_back();
                if (c == size)
                {
                    // fprintf(stderr, "done\n");
                    break;
                }
            }
            else
            {
                rf.read((char *)&t, sizeof(cell_type));
                i++;
            }
        }
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

vector<vector<vector<cell_type>>> readPartitionBF(const string file_path, size_t &signatureSize)
{
    vector<vector<vector<cell_type>>> seqs;
    vector<vector<cell_type>> tseq; // list of BFs for a seq

    ifstream rfSize(file_path, ios::binary | ios::ate);
    if (!rfSize.is_open())
    {
        fprintf(stderr, "Invalid File. Please try again\n");
        exit(0);
    }
    if (rfSize.tellg() > max_fileSize)
    {
        fprintf(stderr, "File too big, please read in batches\n");
        return seqs;
    }
    ifstream rf(file_path, ios::out | ios::binary);

    unsigned long long int length;
    if (rf)
        rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));
    signatureSize = length;

    //? 1 window
    // cout << "length: " << length << "\n";
    vector<cell_type> bf(length);
    cell_type temp = 0;
    size_t i = 0;

    try
    {

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
    }
    catch (const std::exception &e) // caught by reference to base
    {
        std::cout << " a standard exception was caught, with message '"
                  << e.what() << "'\n";
        seqs.clear();
    }
    rf.close();

    return seqs;
}

vector<vector<vector<vector<cell_type>>>> readPartitionBFBatch(const string file_path, size_t size, size_t &signatureSize)
{
    ifstream rf(file_path, ios::out | ios::binary);
    if (!rf.is_open())
    {
        fprintf(stderr, "Invalid File. Please try again\n");
        exit(0);
    }

    unsigned long long int length;
    if (rf)
        rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));
    signatureSize = length;

    vector<vector<vector<vector<cell_type>>>> seqs_batch;
    vector<vector<vector<cell_type>>> seqs;
    vector<vector<cell_type>> tseq; // list of BFs for a seq

    //? 1 window
    // cout << "length: " << length << "\n";
    vector<cell_type> bf(length);
    size_t i = 0;
    size_t count = 0;

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
                count++;
                if (count == size)
                {
                    seqs_batch.push_back(seqs);
                    seqs.clear();
                    count = 0;
                }
            }
            else
            {
                tseq.push_back(bf);
            }
            fill(bf.begin(), bf.end(), 0);
            i = 0;
        }
    }
    seqs_batch.push_back(seqs);
    rf.close();
    return seqs_batch;
}


vector<vector<vector<vector<cell_type>>>> readPartitionBFBatch(const string file_path, size_t size, size_t &signatureSize, size_t cap)
{
    ifstream rf(file_path, ios::out | ios::binary);
    if (!rf.is_open())
    {
        fprintf(stderr, "Invalid File. Please try again\n");
        exit(0);
    }

    unsigned long long int length;
    if (rf)
        rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));
    signatureSize = length;

    vector<vector<vector<vector<cell_type>>>> seqs_batch;
    vector<vector<vector<cell_type>>> seqs;
    vector<vector<cell_type>> tseq; // list of BFs for a seq

    //? 1 window
    // cout << "length: " << length << "\n";
    vector<cell_type> bf(length);
    size_t i = 0;
    size_t count = 0;

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
                count++;

                if (count == cap)
                {
                    fprintf(stderr, "Read seqs cap at %zu seqs...\n", cap);
                    break;
                }
                else if (count % size == 0)
                {
                    seqs_batch.push_back(seqs);
                    seqs.clear();
                }
            }
            else
            {
                tseq.push_back(bf);
            }
            fill(bf.begin(), bf.end(), 0);
            i = 0;
        }
    }
    seqs_batch.push_back(seqs);
    rf.close();
    return seqs_batch;
}


size_t estimateSeqCount(ifstream &rf)
{
    // get length of file:
    rf.seekg(0, rf.end);
    unsigned long long int len = rf.tellg();
    len = len - sizeof(unsigned long long int); // / sizeof(uint64_t);
    rf.seekg(0, rf.beg);

    unsigned long long int length;
    if (rf)
        rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));

    size_t winNum = 0;
    vector<cell_type> bf(length);
    size_t i = 0;
    while (rf)
    {
        rf.read((char *)&bf[i], sizeof(cell_type));
        i++;
        if (i == length)
        {
            if (isEmpty(bf))
            {
                rf.seekg(0, rf.beg);
                return (len * 1.0 / length) / winNum;
            }
            else
            {
                winNum++;
            }
            fill(bf.begin(), bf.end(), 0);
            i = 0;
        }
    }
}

vector<vector<vector<cell_type>>> readPartitionBFSample(const string file_path, size_t &signatureSize, size_t sampleSize)
{
    ifstream rf(file_path, ios::out | ios::binary);
    if (!rf.is_open())
    {
        fprintf(stderr, "Invalid File. Please try again\n");
        exit(0);
    }

    size_t seqCount = estimateSeqCount(rf);

    unsigned long long int length;
    if (rf)
        rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));
    fprintf(stderr, "File contains approximately %zu seqs...\n", seqCount);
    signatureSize = length;

    //? 1 window
    // cout << "length: " << length << "\n";
    vector<cell_type> bf(length);
    cell_type temp = 0;
    size_t i = 0;

    vector<vector<vector<cell_type>>> seqs;
    vector<vector<cell_type>> tseq; // list of BFs for a seq

    if (sampleSize > seqCount)
    {
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
    }
    else
    {
        // actual sample size might be smaller due to the estimated seqCount
        // just ignore for now
        vector<size_t> indices = getIndices(seqCount, sampleSize);
        // for (size_t i : indices)
        // {
        //     fprintf(stderr, "%zu,", i);
        // }
        // fprintf(stderr, "\n");

        size_t idx = indices.back();
        indices.pop_back();
        size_t c = 0;
        size_t last_seen = 0;
        size_t n = sizeof(unsigned long long int);

        while (rf)
        {
            rf.read((char *)&bf[i], sizeof(cell_type));
            i++;
            n++;
            if (i == length)
            {
                if (isEmpty(bf))
                {
                    if (c == idx)
                    {
                        // fprintf(stderr, "%zu,", idx);
                        seqs.push_back(tseq);
                        idx = indices.back();
                        indices.pop_back();
                        last_seen = n;

                        if (seqs.size() == sampleSize)
                        {
                            break;
                        }
                    }
                    c++;
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

        if (seqs.size() < sampleSize)
        {
            fprintf(stderr, "Read all %zu seqs...\n", c);
            fill(bf.begin(), bf.end(), 0);
            tseq.clear();
            i = 0;

            rf.clear();
            rf.seekg(last_seen, rf.beg);
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

                        if (seqs.size() == sampleSize)
                        {
                            break;
                        }
                    }
                    else
                    {
                        tseq.push_back(bf);
                    }
                    fill(bf.begin(), bf.end(), 0);
                    i = 0;
                }
            }
        }
        else
        {
            fprintf(stderr, "Read up to %zu seqs...\n", c);
        }
    }
    rf.close();

    return seqs;
}

#endif