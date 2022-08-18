#ifndef __BUILD_PARTITION_MAIN_CMDLINE_HPP__
#define __BUILD_PARTITION_MAIN_CMDLINE_HPP__

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

class build_partition_main_cmdline
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
  const char *output_arg;
  uint8_t kmer_arg;
  uint32_t window_arg;
  size_t size_arg;

  bool input_given;
  bool output_given;
  bool kmer_given;
  bool window_given;
  bool element_given;
  bool multiple_arg;
  bool canonical_arg;
  bool size_given;
  bool debug;

  enum
  {
    START_OPT = 1000,
    USAGE_OPT
  };

  build_partition_main_cmdline() : input_arg(""), output_arg(""),
                                   kmer_arg(0), window_arg(0),
                                   input_given(false), output_given(false),
                                   kmer_given(false), window_given(false),
                                   element_given(false), multiple_arg(false),
                                   canonical_arg(false), size_given(false),
                                   size_arg(0), debug(false)
  {
  }

  build_partition_main_cmdline(int argc, char *argv[]) : input_arg(""), output_arg(""),
                                                         kmer_arg(0), window_arg(0),
                                                         input_given(false), output_given(false),
                                                         kmer_given(false), window_given(false),
                                                         element_given(false), multiple_arg(false),
                                                         canonical_arg(false), size_given(false),
                                                         size_arg(0), debug(false)
  {
    parse(argc, argv);
  }

  void parse(int argc, char *argv[])
  {
    static struct option long_options[] = {
        {"input", 1, 0, 'i'},
        {"output", 1, 0, 'o'},
        {"multiple", 0, 0, 'm'},
        {"canonical", 0, 0, 'C'},
        {"partition", 0, 0, 'p'},
        {"size", 0, 0, 's'},
        {"help", 0, 0, 'h'},
        {"usage", 0, 0, USAGE_OPT},
        {"version", 0, 0, 'V'},
        {"debug", 0, 0, 'd'},
        {0, 0, 0, 0}};
    static const char *short_options = "hVi:o:k:w:mCs:d";

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
        input_given = true;
        input_arg = optarg;
        break;
      case 'o':
        output_given = true;
        output_arg = optarg;
        break;
      case 'k':
        kmer_given = true;
        kmer_arg = conv_uint<uint8_t>((const char *)optarg, err, false);
        CHECK_ERR(uint8_t, optarg, "-k, --kmerLength=uint8_t")
        break;
      case 'w':
        window_given = true;
        window_arg = conv_uint<uint32_t>((const char *)optarg, err, false);
        CHECK_ERR(uint8_t, optarg, "-w, --windowLength=uint32_t")
        break;
      case 'm':
        multiple_arg = true;
        break;
      case 'C':
        canonical_arg = true;
        break;
      case 's':
        size_given = true;
        size_arg = conv_uint<size_t>((const char *)optarg, err, false);
        CHECK_ERR(size_t, optarg, "-s, --size=size_t")
        break;
      case 'd':
        debug = true;
        break;
      }
    }

    // // Parse arguments
    // if(argc - optind < 1)
    //   error("Requires at least 1 argument.");
  }
  static const char *usage() { return "Usage: sm pbuild [options] "; }
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
    return "partition input seq(s) into window of length w and select s minimisers per window and store in a bf\n\n"
           "Options (default value in (), *required):\n"
           " -i, --input                              input fasta file\n"
           " -o, --output                             output path\n"
           " -C, --canonical                          canonical [default=FALSE]\n"
           " -s, --size                               number of minimiser per window [default=3]\n"
           " -k,                                      kmer length [default=4]\n"
           " -w,                                      window length [default=8]\n"
           " -m,                                      output one binary file per seq [default=FALSE], give folder name with -b\n"
           "     --usage                              Usage\n"
           " -h, --help                               This message\n"
           " -V, --version                            Version\n"
           " -d, --debug                              print mer in string in stdout";
  }
  static const char *hidden() { return ""; }
  void print_version(::std::ostream &os = std::cout) const
  {
#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.0"
#endif
    os << PACKAGE_VERSION << "\n";
  }
  void build_partition(::std::ostream &os = std::cout)
  {
    os << " input_given:" << input_given << "\t"
       << " input_arg:" << input_arg << "\n";
    os << " output_given:" << output_given << "\t"
       << " output_arg:" << output_arg << "\n";
    os << " kmer_given:" << kmer_given << "\t"
       << " kmer_arg:" << kmer_arg << "\n";
    os << " window_given:" << window_given << "\t"
       << " window_arg:" << window_arg << "\n";
    os << " multiple_arg:" << multiple_arg << "\n";
    os << " canonical_arg:" << canonical_arg << "\n";
    os << " size_given:" << size_given << "\n";
    os << " size_arg:" << size_arg << "\n";
  }
};
#endif // __BUILD_PARTITION_MAIN_CMDLINE_HPP__"
