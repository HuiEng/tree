#ifndef __primary_tree_MAIN_CMDLINE_HPP__
#define __primary_tree_MAIN_CMDLINE_HPP__

#include <stdint.h>
#include <unistd.h>
#include <stdlib.h>
#include <getopt.h>
#include <errno.h>
#include <string.h>
#include <stdexcept>
#include <string>
#include <limits>
#include <vector>
#include <iostream>
#include <sstream>
#include <memory>

class primary_tree_main_cmdline
{
  // Boiler plate stuff. Conversion from string to other formats
  static bool adjust_double_si_suffix(double &res, const char *suffix)
  {
    if (*suffix == '\0')
      return true;
    if (*(suffix + 1) != '\0')
      return false;

    switch (*suffix)
    {
    case 'a':
      res *= 1e-18;
      break;
    case 'f':
      res *= 1e-15;
      break;
    case 'p':
      res *= 1e-12;
      break;
    case 'n':
      res *= 1e-9;
      break;
    case 'u':
      res *= 1e-6;
      break;
    case 'm':
      res *= 1e-3;
      break;
    case 'k':
      res *= 1e3;
      break;
    case 'M':
      res *= 1e6;
      break;
    case 'G':
      res *= 1e9;
      break;
    case 'T':
      res *= 1e12;
      break;
    case 'P':
      res *= 1e15;
      break;
    case 'E':
      res *= 1e18;
      break;
    default:
      return false;
    }
    return true;
  }

  static double conv_double(const char *str, ::std::string &err, bool si_suffix)
  {
    char *endptr = 0;
    errno = 0;
    double res = strtod(str, &endptr);
    if (endptr == str)
    {
      err.assign("Invalid floating point string");
      return (double)0.0;
    }
    if (errno)
    {
      err.assign(strerror(errno));
      return (double)0.0;
    }
    bool invalid =
        si_suffix ? !adjust_double_si_suffix(res, endptr) : *endptr != '\0';
    if (invalid)
    {
      err.assign("Invalid character");
      return (double)0.0;
    }
    return res;
  }

  static int conv_enum(const char *str, ::std::string &err, const char *const strs[])
  {
    int res = 0;
    for (const char *const *cstr = strs; *cstr; ++cstr, ++res)
      if (!strcmp(*cstr, str))
        return res;
    err += "Invalid constant '";
    err += str;
    err += "'. Expected one of { ";
    for (const char *const *cstr = strs; *cstr; ++cstr)
    {
      if (cstr != strs)
        err += ", ";
      err += *cstr;
    }
    err += " }";
    return -1;
  }

  template <typename T>
  static bool adjust_int_si_suffix(T &res, const char *suffix)
  {
    if (*suffix == '\0')
      return true;
    if (*(suffix + 1) != '\0')
      return false;

    switch (*suffix)
    {
    case 'k':
      res *= (T)1000;
      break;
    case 'M':
      res *= (T)1000000;
      break;
    case 'G':
      res *= (T)1000000000;
      break;
    case 'T':
      res *= (T)1000000000000;
      break;
    case 'P':
      res *= (T)1000000000000000;
      break;
    case 'E':
      res *= (T)1000000000000000000;
      break;
    default:
      return false;
    }
    return true;
  }

  template <typename T>
  static T conv_int(const char *str, ::std::string &err, bool si_suffix)
  {
    char *endptr = 0;
    errno = 0;
    long long int res = strtoll(str, &endptr, 0);
    if (endptr == str)
    {
      err.assign("Invalid signed int string");
      return (T)0;
    }
    if (errno)
    {
      err.assign(strerror(errno));
      return (T)0;
    }
    bool invalid =
        si_suffix ? !adjust_int_si_suffix(res, endptr) : *endptr != '\0';
    if (invalid)
    {
      err.assign("Invalid character");
      return (T)0;
    }
    if (res > ::std::numeric_limits<T>::max() ||
        res < ::std::numeric_limits<T>::min())
    {
      err.assign("Value out of range");
      return (T)0;
    }
    return (T)res;
  }

  template <typename T>
  static T conv_uint(const char *str, ::std::string &err, bool si_suffix)
  {
    char *endptr = 0;
    errno = 0;
    while (isspace(*str))
    {
      ++str;
    }
    if (*str == '-')
    {
      err.assign("Negative value");
      return (T)0;
    }
    unsigned long long int res = strtoull(str, &endptr, 0);
    if (endptr == str)
    {
      err.assign("Invalid unsigned int string");
      return (T)0;
    }
    if (errno)
    {
      err.assign(strerror(errno));
      return (T)0;
    }
    bool invalid =
        si_suffix ? !adjust_int_si_suffix(res, endptr) : *endptr != '\0';
    if (invalid)
    {
      err.assign("Invalid character");
      return (T)0;
    }
    if (res > ::std::numeric_limits<T>::max())
    {
      err.assign("Value out of range");
      return (T)0;
    }
    return (T)res;
  }

  template <typename T>
  static ::std::string vec_str(const std::vector<T> &vec)
  {
    ::std::ostringstream os;
    for (typename ::std::vector<T>::const_iterator it = vec.begin();
         it != vec.end(); ++it)
    {
      if (it != vec.begin())
        os << ",";
      os << *it;
    }
    return os.str();
  }

  class string : public ::std::string
  {
  public:
    string() : ::std::string() {}
    explicit string(const ::std::string &s) : std::string(s) {}
    explicit string(const char *s) : ::std::string(s) {}
    int as_enum(const char *const strs[])
    {
      ::std::string err;
      int res = conv_enum((const char *)this->c_str(), err, strs);
      if (!err.empty())
        throw ::std::runtime_error(err);
      return res;
    }

    uint32_t as_uint32_suffix() const { return as_uint32(true); }
    uint32_t as_uint32(bool si_suffix = false) const
    {
      ::std::string err;
      uint32_t res = conv_uint<uint32_t>((const char *)this->c_str(), err, si_suffix);
      if (!err.empty())
      {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to uint32_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    uint64_t as_uint64_suffix() const { return as_uint64(true); }
    uint64_t as_uint64(bool si_suffix = false) const
    {
      ::std::string err;
      uint64_t res = conv_uint<uint64_t>((const char *)this->c_str(), err, si_suffix);
      if (!err.empty())
      {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to uint64_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    int32_t as_int32_suffix() const { return as_int32(true); }
    int32_t as_int32(bool si_suffix = false) const
    {
      ::std::string err;
      int32_t res = conv_int<int32_t>((const char *)this->c_str(), err, si_suffix);
      if (!err.empty())
      {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to int32_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    int64_t as_int64_suffix() const { return as_int64(true); }
    int64_t as_int64(bool si_suffix = false) const
    {
      ::std::string err;
      int64_t res = conv_int<int64_t>((const char *)this->c_str(), err, si_suffix);
      if (!err.empty())
      {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to int64_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    int as_int_suffix() const { return as_int(true); }
    int as_int(bool si_suffix = false) const
    {
      ::std::string err;
      int res = conv_int<int>((const char *)this->c_str(), err, si_suffix);
      if (!err.empty())
      {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to int_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    long as_long_suffix() const { return as_long(true); }
    long as_long(bool si_suffix = false) const
    {
      ::std::string err;
      long res = conv_int<long>((const char *)this->c_str(), err, si_suffix);
      if (!err.empty())
      {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to long_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
    double as_double_suffix() const { return as_double(true); }
    double as_double(bool si_suffix = false) const
    {
      ::std::string err;
      double res = conv_double((const char *)this->c_str(), err, si_suffix);
      if (!err.empty())
      {
        ::std::string msg("Invalid conversion of '");
        msg += *this;
        msg += "' to double_t: ";
        msg += err;
        throw ::std::runtime_error(msg);
      }
      return res;
    }
  };

public:
  const char *input_arg;
  const char *tag_arg;
  size_t minimiser_match_arg;
  size_t capacity_arg;
  size_t sizeCap_arg;
  const char *topology_in_arg;
  const char *topology_out_arg;
  double split_threshold_arg;
  double stay_threshold_arg;
  size_t iteration_arg;
  size_t method_arg;
  unsigned seed_arg;
  size_t tree_order_arg;

  bool input_given;
  bool tag_given;
  bool minimiser_match_given;
  bool capacity_given;
  bool sizeCap_given;
  bool random_arg;
  bool topology_in_given;
  bool topology_out_given;
  bool split_threshold_given;
  bool stay_threshold_given;
  bool debug_arg;
  bool print_arg;
  bool force_split_arg;
  bool tree_order_given;
  bool iteration_given;
  bool single_arg;
  bool multiple_arg;

  enum
  {
    START_OPT = 1000,
    USAGE_OPT,
    DEBUG_OPT = 'd'
  };

  primary_tree_main_cmdline() : input_arg(""), tag_arg(""), tag_given(false), tree_order_arg(0),
                                minimiser_match_arg(0), capacity_arg(0), sizeCap_arg(0), method_arg(0),
                                input_given(false),
                                minimiser_match_given(false), capacity_given(false), sizeCap_given(false),
                                random_arg(false), iteration_arg(0), seed_arg(0),
                                topology_in_given(false), topology_in_arg(""),
                                topology_out_given(false), topology_out_arg(""),
                                split_threshold_given(false), split_threshold_arg(0),
                                stay_threshold_given(false), stay_threshold_arg(0),
                                debug_arg(false), print_arg(false), force_split_arg(false), tree_order_given(false),
                                iteration_given(false), single_arg(false), multiple_arg(false)
  {
  }

  primary_tree_main_cmdline(int argc, char *argv[]) : input_arg(""), tag_arg(""), tag_given(false), tree_order_arg(0),
                                                      minimiser_match_arg(0), capacity_arg(0), sizeCap_arg(0), method_arg(0),
                                                      input_given(false),
                                                      minimiser_match_given(false), capacity_given(false), sizeCap_given(false),
                                                      random_arg(false), iteration_arg(0), seed_arg(0),
                                                      topology_in_given(false), topology_in_arg(""),
                                                      topology_out_given(false), topology_out_arg(""),
                                                      split_threshold_given(false), split_threshold_arg(0),
                                                      stay_threshold_given(false), stay_threshold_arg(0),
                                                      debug_arg(false), print_arg(false), force_split_arg(false), tree_order_given(false),
                                                      iteration_given(false), single_arg(false), multiple_arg(false)
  {
    parse(argc, argv);
  }

  void parse(int argc, char *argv[])
  {
    static struct option long_options[] = {
        {"input", 1, 0, 'i'},
        {"debug", 0, 0, DEBUG_OPT},
        {"print", 0, 0, 'p'},
        {"tag", 1, 0, 'T'},
        {"single", 0, 0, 'x'},
        {"capacity", 1, 0, 'C'},
        {"size_cap", 1, 0, 'c'},
        {"minimiser_match", 1, 0, 'm'},
        {"split", 1, 0, 'L'},
        {"stay", 1, 0, 'S'},
        {"random", 1, 0, 'R'},
        {"iteration", 0, 0, 'I'},
        {"topology_in", 1, 0, 'q'},
        {"topology_out", 1, 0, 'w'},
        {"help", 0, 0, 'h'},
        {"usage", 0, 0, USAGE_OPT},
        {"version", 0, 0, 'V'},
        {0, 0, 0, 0}};
    static const char *short_options = "hVi:o:c:C:R:dfI:q:w:S:L:m:T:pxfF:M";

    ::std::string err;
#define CHECK_ERR(type, val, which)                                                      \
  if (!err.empty())                                                                      \
  {                                                                                      \
    ::std::cerr << "Invalid " #type " '" << val << "' for [" which "]: " << err << "\n"; \
    exit(1);                                                                             \
  }
    while (true)
    {
      int index = -1;
      int c = getopt_long(argc, argv, short_options, long_options, &index);
      if (c == -1)
        break;
      switch (c)
      {
      case ':':
        ::std::cerr << "Missing required argument for "
                    << (index == -1 ? ::std::string(1, (char)optopt) : std::string(long_options[index].name))
                    << ::std::endl;
        exit(1);
      case 'h':
        ::std::cout << usage() << "\n\n"
                    << help() << std::endl;
        exit(0);
      case USAGE_OPT:
        ::std::cout << usage() << "\nUse --help for more information." << std::endl;
        exit(0);
      case DEBUG_OPT:
        debug_arg = true;
        break;
      case 'p':
        print_arg = true;
        break;
      case 'f':
        force_split_arg = true;
        break;
      case 'F':
        force_split_arg = true;
        tree_order_given = true;
        tree_order_arg = conv_uint<size_t>((const char *)optarg, err, false);
        CHECK_ERR(size_t, optarg, "-f=size_t")
        break;
      case 'x':
        single_arg = true;
        break;
      case 'M':
        multiple_arg = true;
        break;
      case 'V':
        print_version();
        exit(0);
      case '?':
        ::std::cerr << "Use --usage or --help for some help\n";
        exit(1);
      case 'i':
        input_given = true;
        input_arg = optarg;
        break;
      case 'T':
        tag_given = true;
        tag_arg = optarg;
        break;
      case 'q':
        topology_in_given = true;
        topology_in_arg = optarg;
        break;
      case 'w':
        topology_out_given = true;
        topology_out_arg = optarg;
        break;
      case 'L':
        split_threshold_given = true;
        split_threshold_arg = conv_double((const char *)optarg, err, false);
        CHECK_ERR(double, optarg, "-t, --split=double")
        // split_threshold_arg = conv_uint<size_t>((const char *)optarg, err, false);
        // CHECK_ERR(size_t, optarg, "-t, --split=size_t")
        break;
      case 'S':
        stay_threshold_given = true;
        stay_threshold_arg = conv_double((const char *)optarg, err, false);
        CHECK_ERR(double, optarg, "-t, --stay=double")
        break;
      case 'm':
        minimiser_match_given = true;
        minimiser_match_arg = conv_uint<size_t>((const char *)optarg, err, false);
        CHECK_ERR(size_t, optarg, "-m, --minimiser_match=size_t")
        break;
      case 'C':
        capacity_given = true;
        capacity_arg = conv_uint<size_t>((const char *)optarg, err, false);
        CHECK_ERR(size_t, optarg, "-C, --capacity=size_t")
        break;
      case 'c':
        sizeCap_given = true;
        sizeCap_arg = conv_uint<size_t>((const char *)optarg, err, false);
        CHECK_ERR(size_t, optarg, "-c, --size_cap=size_t")
        break;
      case 'o':
        method_arg = conv_uint<size_t>((const char *)optarg, err, false);
        CHECK_ERR(size_t, optarg, "-m=size_t")
        break;
      case 'R':
        random_arg = true;
        seed_arg = conv_uint<unsigned>((const char *)optarg, err, false);
        CHECK_ERR(size_t, optarg, "-R=unsigned")
        break;
      case 'I':
        iteration_given = true;
        iteration_arg = conv_uint<size_t>((const char *)optarg, err, false);
        CHECK_ERR(size_t, optarg, "-I, --iteration=size_t")
        break;
      }
    }

    // // Parse arguments
    // if(argc - optind < 1)
    //   error("Requires at least 1 argument.");
  }
  static const char *usage() { return "Usage: primary_tree primary_tree [options] "; }
  class error
  {
    int code_;
    std::ostringstream msg_;

    // Select the correct version (GNU or XSI) version of
    // strerror_r. strerror_ behaves like the GNU version of strerror_r,
    // regardless of which version is provided by the system.
    static const char *strerror__(char *buf, int res)
    {
      return res != -1 ? buf : "Invalid error";
    }
    static const char *strerror__(char *buf, char *res)
    {
      return res;
    }
    static const char *strerror_(int err, char *buf, size_t buflen)
    {
      return strerror__(buf, strerror_r(err, buf, buflen));
    }
    struct no_t
    {
    };

  public:
    static no_t no;
    error(int code = EXIT_FAILURE) : code_(code) {}
    explicit error(const char *msg, int code = EXIT_FAILURE) : code_(code)
    {
      msg_ << msg;
    }
    error(const std::string &msg, int code = EXIT_FAILURE) : code_(code)
    {
      msg_ << msg;
    }
    error &operator<<(no_t)
    {
      char buf[1024];
      msg_ << ": " << strerror_(errno, buf, sizeof(buf));
      return *this;
    }
    template <typename T>
    error &operator<<(const T &x)
    {
      msg_ << x;
      return (*this);
    }
    ~error()
    {
      ::std::cerr << "Error: " << msg_.str() << "\n"
                  << usage() << "\n"
                  << "Use --help for more information"
                  << ::std::endl;
      exit(code_);
    }
  };
  static const char *help()
  {
    return "Read minimisers BF and get all-against-all Jaccard primary_tree\n\n"
           "Options (default value in (), *required):\n"
           " -d, --debug                              print debug files\n"
           " -p, --print                              print debug msg to stderr [default = false]\n"
           " -x, --single                             input is singleBF [default = false]\n"
           " -o,                                      output path\n"
           " -F,                                      force split root [default=false]\n"
           " -f,                                      tree order force split root [default=5]\n"
           " -i, --input                              input file or folder\n"
           " -M                                       input is a folder"
           " -T, --tag                                [optional] tag for output cluster file =>${tag}${computeFilename}\n"
           " -m, --minimiser_match                    minimiser_match_threshold [default=4], you need to know the number of minimiser per window\n"
           " -C, --capacity                           primary_tree capacity [default=0]\n"
           " -c, --size_cap                           Only cluster the first c seqs [default=10000]\n"
           " -q, --topology_in                        Read tree topology from this folder [must end with '/']\n"
           " -w, --topology_out                       Print tree topology into this folder [must end with '/']\n"
           " -L, --split                              split node if HD >= signature length * split_threshold [default=1]\n"
           " -S, --stay                               stay in node if HD <= signature length * stay_threshold [default=0]\n"
           " -R,                                      inserting seqs in random order [defalut seed=0]\n"
           " -I, --iteration                          number or reinsertion [default=0]\n"
           "     --usage                              Usage\n"
           " -h, --help                               This message\n"
           " -V, --version                            Version";
  }
  static const char *hidden() { return ""; }
  void print_version(::std::ostream &os = std::cout) const
  {
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.0"
#endif
    os << PACKAGE_VERSION << "\n";
  }
  void primary_tree(::std::ostream &os = std::cout)
  {
    os << " input_given:" << input_given << "\t"
       << " input_arg:" << input_arg << "\n";
    os << " minimiser_match_given:" << minimiser_match_given << "\t"
       << " minimiser_match_arg:" << minimiser_match_arg << "\n";
    os << " capacity_given:" << capacity_given << "\t"
       << " capacity_arg:" << capacity_arg << "\n";
    os << " random_arg:" << random_arg << "\n"
       << " seed_arg:" << seed_arg << "\n";
    os << " iteration_arg:" << iteration_arg << "\n";
    os << " topology_in_given:" << topology_in_given << "\t"
       << " topology_in_arg:" << topology_in_arg << "\n";
    os << " split_threshold_given:" << split_threshold_given << "\t"
       << " split_threshold_arg:" << split_threshold_arg << "\n";
    os << " stay_threshold_given:" << stay_threshold_given << "\t"
       << " stay_threshold_arg:" << stay_threshold_arg << "\n";
    os << " tag_given:" << tag_given << "\t"
       << " tag_arg:" << tag_arg << "\n";
  }
};
#endif // __primary_tree_MAIN_CMDLINE_HPP__"
