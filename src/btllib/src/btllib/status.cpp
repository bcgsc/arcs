#include "btllib/status.hpp"

#include <string>

namespace btllib {

std::string
get_time()
{
  time_t now;
  time(&now);
  char buf[sizeof("2011-10-08T07:07:09Z")];
  std::tm tm_result = {};
  localtime_r(&now, &tm_result);
  std::strftime(buf, sizeof buf, "%F %T", &tm_result);
  return std::string(buf);
}

void
log_info(const std::string& msg)
{
  std::cerr << ('[' + get_time() + "]" + PRINT_COLOR_INFO + "[INFO] " +
                PRINT_COLOR_END + msg + '\n')
            << std::flush;
}

void
log_warning(const std::string& msg)
{
  std::cerr << ('[' + get_time() + "]" + PRINT_COLOR_WARNING + "[WARNING] " +
                PRINT_COLOR_END + msg + '\n')
            << std::flush;
}

void
log_error(const std::string& msg)
{
  std::cerr << ('[' + get_time() + "]" + PRINT_COLOR_ERROR + "[ERROR] " +
                PRINT_COLOR_END + msg + '\n')
            << std::flush;
}

void
check_info(bool condition, const std::string& msg)
{
  if (condition) {
    log_info(msg);
  }
}

void
check_warning(bool condition, const std::string& msg)
{
  if (condition) {
    log_warning(msg);
  }
}

void
check_error(bool condition, const std::string& msg)
{
  if (condition) {
    log_error(msg);
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }
}

std::string
get_strerror()
{
  static const size_t buflen = 1024;
  char buf[buflen];
// POSIX and GNU implementation of strerror_r differ, even in function signature
// and so we need to check which one is used
#if __APPLE__ ||                                                               \
  ((_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600) && !_GNU_SOURCE)
  strerror_r(errno, buf, buflen);
  return buf;
#else
  return strerror_r(errno, buf, buflen);
#endif
}

void
check_stream(const std::ios& stream, const std::string& name)
{
  if (!stream.good()) {
    log_error("'" + name + "' stream error: " + get_strerror());
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  }
}

} // namespace btllib