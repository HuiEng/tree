#include <filesystem>
#include "build_partition_main_cmdline.hpp"

using namespace std;
namespace fs = std::filesystem;

static build_partition_main_cmdline args; // Command line switches and arguments

int test_main(int argc, char *argv[])
{
    args.parse(argc, argv);
    std::ios::sync_with_stdio(false); // No sync with stdio -> faster

    string inputFile = args.input_arg;
    if (args.folder_arg)
    {

        string delimiter = "/*";
        string folder = inputFile.substr(0, inputFile.find(delimiter));
        string ext = inputFile.substr(inputFile.find(delimiter) + delimiter.size(), inputFile.size() - 1);
        fprintf(stderr, "Reading folder %s\n", folder.c_str());
        for (const auto &entry : fs::directory_iterator(folder))
        {
            if (entry.path().extension() == ext)
            {
                cout << entry.path() << endl;
            }
        }
    }

    return 0;
}
