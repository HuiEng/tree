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
#include "bf/all.hpp"

using namespace bf;

using namespace std;
namespace fs = std::experimental::filesystem;

template <typename bf_type>
size_t readSignatures(const string file, vector<bf_type> &sigs)
{
    ifstream rf(file, ios::binary);

    double fp;
    if (rf)
        rf.read(reinterpret_cast<char *>(&fp), sizeof(fp));

    size_t capacity;
    if (rf)
        rf.read(reinterpret_cast<char *>(&capacity), sizeof(capacity));

    // get number of bf:
    rf.seekg(0, rf.end);
    int len = rf.tellg();
    len = len - (int)sizeof(double) - (int)sizeof(size_t); // / sizeof(uint64_t);
    len = len / sizeof(std::vector<bitvector::block_type>);
    rf.seekg(0, rf.beg);
    // sigs.resize(len);

    int i = 0;
    while (i < len)
    {
        bf_type bf(fp, capacity);
        std::vector<bitvector::block_type> temp(bf.storage().size());
        rf.read((char *)&temp[0], sizeof(temp));

        bf.update(temp);
        i++;
        sigs.push_back(bf);
    }
    std::cout << "done\n";
    rf.close();
}

// // read all files in given folder
// unsigned long long int readSignaturesMultiple(const string folder, vector<bf_type> &sigs)
// {
//     unsigned long long int length;
//     // keep track of sequence index
//     FILE *pFile = fopen((folder + "-seqs.txt").c_str(), "w");

//     for (const auto &entry : fs::directory_iterator(folder))
//     {
//         if (entry.path().extension() == ".bin")
//         {
//             fprintf(pFile, "%s\n", entry.path().c_str());
//             ifstream rf(entry.path(), ios::out | ios::binary);

//             // get length of file:
//             rf.seekg(0, rf.end);
//             int len = rf.tellg();
//             len = len - (int)sizeof(unsigned long long int); // / sizeof(uint64_t);
//             rf.seekg(0, rf.beg);

//             if (rf)
//                 rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));

//             auto old_size = sigs.size();
//             sigs.resize(old_size + len);
//             size_t i = 0;
//             while (rf)
//             {
//                 rf.read((char *)&sigs[old_size + i], sizeof(cell_type));
//                 i++;
//             }
//             rf.close();
//         }
//     }
//     return length;
// }

// // read all files from the input txt, each line is the full path of the file to read
// unsigned long long int readSignaturesSpecific(const string ref_file_path, vector<bf_type> &sigs)
// {
//     unsigned long long int length;
//     // keep track of sequence index
//     std::ifstream ref_file(ref_file_path);
//     string file;
//     while (getline(ref_file, file))
//     { //read data from file object and put it into string.
//         ifstream rf(file.c_str(), ios::out | ios::binary);
//         fprintf(stderr, "reading %s\n", file.c_str()); // << file << "\n"; //print the data of the string

//         // get length of file:
//         rf.seekg(0, rf.end);
//         int len = rf.tellg();
//         len = len - (int)sizeof(unsigned long long int); // / sizeof(uint64_t);
//         rf.seekg(0, rf.beg);
//         fprintf(stderr, "%d\n", len);

//         if (rf)
//             rf.read(reinterpret_cast<char *>(&length), sizeof(unsigned long long int));
//         fprintf(stderr, "here22\n");
//         auto old_size = sigs.size();

//         sigs.resize(old_size + len);
//         size_t i = 0;
//         fprintf(stderr, "here\n");
//         while (rf)
//         {
//             rf.read((char *)&sigs[old_size + i], sizeof(cell_type));
//             i++;
//         }
//         rf.close();
//     }

//     return 1;
// }

#endif