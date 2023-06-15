// #include "bf_ktree.hpp"
#include <fstream>
#include <random>
#include <filesystem>
#include "histo_main_cmdline.hpp"
static histo_main_cmdline args; // Command line switches and arguments
using namespace std;
namespace fs = std::filesystem;

// output signature size, if file too big => return empty "sigs", use readSignaturesBatch
void readSimFile(FILE *pFile, const char *file)
{
    vector<size_t> histo(11);

    FILE *fp = fopen(file, "r");
    if (!fp)
    {
        fprintf(stderr, "Failed to load %s\n", file);
        exit(1);
    }

    fprintf(stderr, "reading %s\n", file);
    char t[100];
    fscanf(fp, "%s", t); // skip first line

    // vector<string> data;
    for (;;)
    {
        size_t i, j;
        double similarity;
        if (fscanf(fp, "%zu,%zu,%lf\n", &i, &j, &similarity) < 1)
            break;

        histo[similarity / 10]++;
        // char fileBuf[100];
        // if (fscanf(fp, "%s\n", fileBuf) < 1)
        //     break;

        // // data.push_back(string(fileBuf));
        // cout << fileBuf << endl;
    }

    for (size_t i = 0; i < histo.size(); i++)
    {
        fprintf(pFile, "%zu,%zu\n", i * 10, histo[i]);
    }
}

int histo_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    std::ios::sync_with_stdio(false); // No sync with stdio -> faster

    string bfIn = args.bf_input_arg;

    size_t firstindex = bfIn.find_last_of("/") + 1;
    size_t lastindex = bfIn.find_last_of("-");
    string rawname = bfIn.substr(firstindex, lastindex - firstindex);

    if (args.output_given)
    {
        rawname = rawname + "-" + args.output_arg;
    }
    FILE *pFile = fopen((rawname + "-histo.txt").c_str(), "w");
    readSimFile(pFile, bfIn.c_str());
    return 0;
}
