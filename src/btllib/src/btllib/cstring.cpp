#include "btllib/cstring.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <string>

namespace btllib {

CString::CString(const CString& cstr)
{
  if (cstr.s_size + 1 > s_cap) {
    change_cap(cstr.s_size + 1);
  }
  s_size = cstr.s_size;
  std::memcpy(s, cstr.s, s_size + 1);
}

CString::CString(CString&& cstr) noexcept
  : s_size(cstr.s_size)
{
  std::swap(s, cstr.s);
  cstr.clear();
  std::swap(s_cap, cstr.s_cap);
}

CString::CString(const std::string& str)
{
  if (str.size() + 1 > s_cap) {
    change_cap(str.size() + 1);
  }
  s_size = str.size();
  std::memcpy(s, str.c_str(), s_size + 1);
}

CString&
CString::operator=(const CString& cstr)
{
  if (this == &cstr) {
    return *this;
  }
  if (cstr.s_size + 1 > s_cap) {
    change_cap(cstr.s_size + 1);
  }
  s_size = cstr.s_size;
  std::memcpy(s, cstr.s, s_size + 1);
  return *this;
}

CString&
CString::operator=(CString&& cstr) noexcept
{
  std::swap(s, cstr.s);
  s_size = cstr.s_size;
  cstr.clear();
  std::swap(s_cap, cstr.s_cap);
  return *this;
}

CString&
CString::operator=(const std::string& str)
{
  if (str.size() + 1 > s_cap) {
    change_cap(str.size() + 1);
  }
  s_size = str.size();
  std::memcpy(s, str.c_str(), s_size + 1);
  return *this;
}

CString&
CString::operator+=(const CString& cstr)
{
  const auto new_size = s_size + cstr.s_size;
  if (new_size + 1 > s_cap) {
    const auto factor = size_t(std::pow(
      2,
      std::ceil(std::log2(double(new_size + 1)) - std::log2(double(s_size)))));
    change_cap(s_size * factor);
  }
  std::memcpy(s + s_size, cstr.s, cstr.s_size);
  s_size = new_size;
  return *this;
}

CString&
CString::operator+=(const std::string& str)
{
  const auto new_size = s_size + str.size();
  if (new_size + 1 > s_cap) {
    const auto factor = size_t(std::pow(
      2,
      std::ceil(std::log2(double(new_size + 1)) - std::log2(double(s_size)))));
    change_cap(s_size * factor);
  }
  std::memcpy(s + s_size, str.c_str(), str.size());
  s[new_size] = '\0';
  s_size = new_size;
  return *this;
}

CString&
CString::operator+=(const char c)
{
  const auto new_size = s_size + 1;
  if (new_size + 1 > s_cap) {
    const auto factor = size_t(std::pow(
      2,
      std::ceil(std::log2(double(new_size + 1)) - std::log2(double(s_size)))));
    change_cap(s_size * factor);
  }
  s[new_size - 1] = c;
  s[new_size] = '\0';
  s_size = new_size;
  return *this;
}

void
CString::clear()
{
  s_size = 0;
  s[0] = '\0';
}

void
CString::pop_back()
{
  s_size -= 1;
  s[s_size] = '\0';
}

void
CString::change_cap(const size_t new_cap)
{
  s_cap = new_cap;
  s = (char*)std::realloc(s, new_cap); // NOLINT
}

void
CString::resize(const size_t n, const char c)
{
  if (n > s_size) {
    change_cap(n + 1);
    for (size_t i = s_size; i < n; i++) {
      s[i] = c;
    }
  }
  s_size = n;
  s[s_size] = '\0';
}

CString&
CString::erase(const size_t pos, size_t len)
{
  if (pos + len > size()) {
    len = size() - pos;
  }
  const ssize_t to_move = ssize_t(size()) - ssize_t(pos) - ssize_t(len);
  if (to_move > 0 && to_move < ssize_t(size())) {
    std::memmove(s + pos, s + pos + len, to_move);
  }
  resize(size() - len);
  return *this;
}

} // namespace btllib