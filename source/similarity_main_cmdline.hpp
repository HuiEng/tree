#ifndef __SIMILARITY_MAIN_CMDLINE_HPP__
#define __SIMILARITY_MAIN_CMDLINE_HPP__

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

class similarity_main_cmdline
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
  const char *bf_input_arg;
  const char *output_arg;
  size_t threshold_arg;
  size_t window_arg;
  size_t element_arg;
  size_t batch_arg;
  double min_sim_arg;

  bool bf_input_given;
  bool output_given;
  bool threshold_given;
  bool window_given;
  bool element_given;
  bool multiple_arg;
  bool all_kmer_arg;
  bool global_arg;
  bool local_arg;
  bool skip_arg;
  bool min_sim_given;
  bool batch_given;

  enum
  {
    START_OPT = 1000,
    USAGE_OPT
  };

  similarity_main_cmdline() : bf_input_arg(""), output_arg(""),
                              threshold_arg(0), window_arg(0),
                              element_arg(0), batch_arg(0), min_sim_arg(0),
                              bf_input_given(false), output_given(false), min_sim_given(false),
                              threshold_given(false), window_given(false),
                              element_given(false), multiple_arg(false),
                              all_kmer_arg(false), global_arg(false),
                              skip_arg(false), local_arg(false), batch_given(false)
  {
  }

  similarity_main_cmdline(int argc, char *argv[]) : bf_input_arg(""), output_arg(""),
                                                    threshold_arg(0), window_arg(0),
                                                    element_arg(0), batch_arg(0), min_sim_arg(0),
                                                    bf_input_given(false), output_given(false), min_sim_given(false),
                                                    threshold_given(false), window_given(false),
                                                    element_given(false), multiple_arg(false),
                                                    all_kmer_arg(false), global_arg(false),
                                                    skip_arg(false), local_arg(false), batch_given(false)
  {
    parse(argc, argv);
  }

  void parse(int argc, char *argv[])
  {
    static struct option long_options[] = {
        {"input", 1, 0, 'i'},
        {"batch", 1, 0, 'b'},
        {"output", 1, 0, 'o'},
        {"capacity", 1, 0, 'c'},
        {"threshold", 1, 0, 't'},
        {"min", 1, 0, 'M'},
        {"multiple", 0, 0, 'm'},
        {"all", 0, 0, 'A'},
        {"skip", 0, 0, 'S'},
        {"global", 0, 0, 'G'},
        {"local", 0, 0, 'L'},
        {"help", 0, 0, 'h'},
        {"usage", 0, 0, USAGE_OPT},
        {"version", 0, 0, 'V'},
        {0, 0, 0, 0}};
    static const char *short_options = "hVi:o:t:w:M:mAGLSb:";

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
      case 'V':
        print_version();
        exit(0);
      case '?':
        ::std::cerr << "Use --usage or --help for some help\n";
        exit(1);
      case 'i':
        bf_input_given = true;
        bf_input_arg = optarg;
        break;
      case 'o':
        output_given = true;
        output_arg = optarg;
        break;
      case 'b':
        batch_given = true;
        batch_arg = conv_uint<size_t>((const char *)optarg, err, false);
        CHECK_ERR(size_t, optarg, "-b, --batch=size_t")
        break;
      case 't':
        threshold_given = true;
        threshold_arg = conv_uint<size_t>((const char *)optarg, err, false);
        CHECK_ERR(size_t, optarg, "-t, --threshold=size_t")
        break;
      case 'w':
        window_given = true;
        window_arg = conv_uint<size_t>((const char *)optarg, err, false);
        CHECK_ERR(size_t, optarg, "-w, --windowLength=size_t")
        break;
      case 'M':
        min_sim_given = true;
        min_sim_arg = conv_double((const char *)optarg, err, false);
        CHECK_ERR(double, optarg, "-s, --BF-elementCount=double")
        break;
      case 'm':
        multiple_arg = true;
        break;
      case 'A':
        all_kmer_arg = true;
        break;
      case 'S':
        skip_arg = true;
        break;
      case 'G':
        global_arg = true;
        break;
      case 'L':
        local_arg = true;
        break;
      }
    }

    // // Parse arguments
    // if(argc - optind < 1)
    //   error("Requires at least 1 argument.");
  }
  static const char *usage() { return "Usage: tree sim [options] "; }
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
    return "Read minimisers BF and get all-against-all similarity, default is with matching window\n\n"
           "Options (default value in (), *required):\n"
           " -i, --input                              BF input path\n"
           " -b, --batch                              recommended for processing large input, do in batch size [default=false]\n"
           " -o, --output                             output path is {input}{-output}_sim.txt\n"
           " -t, --threshold                          matching threshold per window  [default=4]\n"
           " -M, --min                                only print entry greater than this value, and do skip [default=0.5]\n"
           " -m,                                      output one binary file per seq [default=FALSE], give folder name with -b\n"
           " -A, --all                                input is all kmers\n"
           " -G, --global                             count with global minimiser set\n"
           " -L, --local                              count with local jaccard\n"
           " {no flag}                                count % matching minimsers by window\n"
           " -S, --skip                               skip printing half the output\n"
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
  void similarity(::std::ostream &os = std::cout)
  {
    os << " bf_input_given:" << bf_input_given << "\t"
       << " bf_input_arg:" << bf_input_arg << "\n";
    os << " output_given:" << output_given << "\t"
       << " output_arg:" << output_arg << "\n";
    os << " threshold_given:" << threshold_given << "\t"
       << " threshold_arg:" << threshold_arg << "\n";
    os << " window_given:" << window_given << "\t"
       << " window_arg:" << window_arg << "\n";
    os << " element_given:" << element_given << "\t"
       << " element_arg:" << element_arg << "\n";
    os << " multiple_arg:" << multiple_arg << "\n";
  }
};
#endif // __SIMILARITY_MAIN_CMDLINE_HPP__"
