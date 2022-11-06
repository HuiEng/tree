// #include "bf_ktree.hpp"
#include <fstream>
#include <random>
#include "merge_main_cmdline.hpp"
#include "bloom_filter.hpp"
#include "read.hpp"

using namespace std;
typedef pair<string, string> config_type;

static merge_main_cmdline args; // Command line switches and arguments

vector<config_type> readConfig(const char *path)
{
    FILE *fp = fopen(path, "r");
    if (!fp)
    {
        fprintf(stderr, "Failed to load %s\n", path);
        exit(1);
    }
    vector<config_type> data;
    for (;;)
    {
        char colname[100];
        char fileBuf[100];
        if (fscanf(fp, "%100[^,],%100[^\n]\n", colname, fileBuf) < 1)
            break;
        // break;
        // fprintf(stderr,"%d,%s\n",isChunk, (char*)&fileBuf[0]);
        data.push_back(make_pair(string(colname), string(fileBuf)));
    }

    // for (config_type val : data)
    // {
    //     fprintf(stderr, "%s---%s\n", val.first.c_str(), val.second.c_str());
    // }
    fclose(fp);

    return data;
}

// read index column for the first file
vector<string> readFirstFile(const char *path)
{
    FILE *fp = fopen(path, "r");
    if (!fp)
    {
        fprintf(stderr, "Failed to load %s\n", path);
        exit(1);
    }

    fprintf(stderr, "reading %s\n", path);
    char t[100];
    fscanf(fp, "%s", t); // skip first line

    vector<string> data;
    for (;;)
    {
        // size_t idx;
        // size_t cnt;
        // if (fscanf(fp, "%zu,%zu\n", &idx, &cnt) < 1)
        char fileBuf[100];
        if (fscanf(fp, "%s\n", fileBuf) < 1)
            break;

        data.push_back(string(fileBuf));
    }
    return data;
}

// use this function to append the other files to the first one
vector<string> readFile(const char *path, vector<string> &data)
{
    FILE *fp = fopen(path, "r");
    if (!fp)
    {
        fprintf(stderr, "Failed to load %s\n", path);
        exit(1);
    }

    fprintf(stderr, "reading %s\n", path);
    char fileBuf[100];
    fscanf(fp, "%s", fileBuf); // skip first line

    size_t i = 0;
    for (;;)
    {
        size_t idx;
        size_t cnt;
        if (fscanf(fp, "%zu,%zu\n", &idx, &cnt) < 1)
            break;

        data[i] = data[i] + "," + to_string(cnt);
        i++;
    }
    return data;
}


// // use this function to append the other files to the first one
// vector<string> readFile(const char *path, vector<string> &data)
// {
//     FILE *fp = fopen(path, "r");
//     if (!fp)
//     {
//         fprintf(stderr, "Failed to load %s\n", path);
//         exit(1);
//     }

//     fprintf(stderr, "reading %s\n", path);
//     char fileBuf[100];
//     fscanf(fp, "%s", fileBuf); // skip first line

//     size_t i = 0;
//     for (;;)
//     {
//         size_t i;
//         size_t j;
//         float cnt;
//         if (fscanf(fp, "%zu,%zu,%f\n", &i,&j, &cnt) < 1)
//             break;

//         data[i] = data[i] + "," + to_string(cnt);
//         i++;
//     }
//     return data;
// }

void doMerge(vector<config_type> config, FILE *pFile)
{
    vector<string> results = readFirstFile(config[0].second.c_str());
    string colname = "i," + config[0].first;

    for (int i = 1; i < config.size(); i++)
    {
        readFile(config[i].second.c_str(), results);
        colname = colname + "," + config[i].first;
    }

    fprintf(pFile, "%s\n", colname.c_str());

    for (string val : results)
    {
        fprintf(pFile, "%s\n", val.c_str());
    }
}

int merge_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    std::ios::sync_with_stdio(false); // No sync with stdio -> faster]
    vector<config_type> config = readConfig(args.bf_input_arg);
    FILE *pFile = fopen(args.output_arg, "w");
    doMerge(config, pFile);
    return 0;
}
