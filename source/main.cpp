#include <fstream>
#include <random>
#include "minimiser.hpp"
// #include "bf_ktree.hpp"
#include "bloom_filter.hpp"
// #include "bloom_filter3.hpp"
// #include "bf/object.hpp"

// #include "bf/all.hpp"
// using namespace bf;

static size_t kmerLength = 4;   // Kmer length
static size_t windowLength = 8; // window length

typedef int(main_func_t)(int argc, char *argv[]);
main_func_t build_main;
main_func_t build_partition_main;
main_func_t similarity_main;
main_func_t count_main;
main_func_t ktree_main;
main_func_t minimiser_main;
main_func_t partree_main;
main_func_t partree_multi_main;
main_func_t tree_main;
main_func_t sos;

struct cmd_func
{
  const char *cmd;
  //  std::string  cmd;
  main_func_t *func;
};
cmd_func cmd_list[] = {
    {"build", &build_main},
    {"sim", &similarity_main},
    {"count", &count_main},
    {"ktree", &ktree_main},
    {"minimiser", &minimiser_main},
    {"pbuild", &build_partition_main},
    {"partree", &partree_main},
    {"partreeM", &partree_multi_main},
    {"tree", &tree_main},

    /* help in all its form. Must be first non-command */
    {"help", &sos},
    {"-h", &sos},
    {"-help", &sos},
    {"--help", &sos},
    {"-?", &sos},
    {"", 0}};

void __sos(std::ostream *os)
{
  *os << "Usage: sm <cmd> [options] arg..." << std::endl
      << "Where <cmd> is one of: ";
  bool comma = false;
  for (cmd_func *ccmd = cmd_list; ccmd->func != sos; ccmd++)
  {
    *os << (comma ? ", " : "") << ccmd->cmd;
    comma = true;
  }
  *os << "." << std::endl;
  *os << "Options:" << std::endl
      << "  --help           Display this message" << std::endl;
}

int sos(int argc, char *argv[])
{
  __sos(&std::cout);
  return 0;
}


int main(int argc, char *argv[])
{

  std::string error;

  if (argc < 2)
  {
    error = "Too few arguments";
  }
  else
  {
    for (cmd_func *ccmd = cmd_list; ccmd->func != 0; ccmd++)
    {
      if (!strcmp(ccmd->cmd, argv[1]))
        //      if(!ccmd->cmd.compare(argv[1]))
        return ccmd->func(argc - 1, argv + 1);
    }
    error = "Unknown command '";
    error += argv[1];
    error += "'\n";
  }

  std::cerr << error << std::endl;
  __sos(&std::cerr);
  return 1;
}
