#include <algorithm>
#include <cstring>
#include <iostream>

class String {
private:
  size_t size_;
  size_t capacity_;
  char* string_;
public:
  String() : size_(0), capacity_(2), string_(new char[capacity_]) {
    string_[0] = '\0';
  }

  String(int size, int capacity, char* string) : size_(size), capacity_(capacity), string_(string) {
    string_[size_ + 1] = '\0';
  }

  String(const String& str) : size_(str.size_), capacity_(str.capacity_), string_(new char[capacity_]) {
    std::copy(str.string_, str.string_ + size_ + 1, string_);
  }

  String& operator=(const String& str) {
    if (&str == this) {
      return *this;
    }
    delete[] string_;
    size_ = str.size_;
    capacity_ = str.capacity_;
    string_ = new char[capacity_];
    memcpy(string_, str.string_, size_ + 1);
    return *this;
  }

  String(const char* const_string) : size_(strlen(const_string)), capacity_(size_ + 1), string_(new char[capacity_]) {
    memcpy(string_, const_string, size_ + 1);
  }

  String(size_t num, const char ch) : size_(num), capacity_(num + 1), string_(new char[capacity_]) {
    std::fill(string_, string_ + num, ch);
    string_[size_] = '\0';
  }

  String(char ch) : size_(1), capacity_(2), string_(new char[capacity_]) {
    string_[0] = ch;
    string_[1] = '\0';
  }

  const char& operator[](int ind) const {
    return string_[ind];
  }

  char& operator[](int ind) {
    return string_[ind];
  }

  size_t length() const {
    return size_;
  }

  size_t size() const {
    return size_;
  }

  size_t capacity() const {
    return capacity_ - 1;
  }

  void push_back(char ch) {
    if (size_ + 1 >= capacity_) {
      capacity_ *= 2;
      char* string1 = new char[capacity_];
      memcpy(string1, string_, size_ + 1);
      string1[size_] = ch;
      string1[size_ + 1] = '\0';
      ++size_;
      delete[] string_;
      string_ = string1;
    } else {
      string_[size_] = ch;
      string_[size_ + 1] = '\0';
      ++size_;
    }
  }

  void pop_back() {
    string_[size_ - 1] = '\0';
    --size_;
  }

  const char& front() const {
    return *(string_);
  }

  char& front() {
    return *(string_);
  }

  const char& back() const {
    return *(string_ + size_ - 1);
  }

  char& back() {
    return *(string_ + size_ - 1);
  }

  String& operator+=(char ch) {
    push_back(ch);
    return *this;
  }

  String& operator+=(const String& str) {
    while (size_ + str.size_ >= capacity_) {
      capacity_ *= 2;
    }
    char* string1 = new char[capacity_];
    memcpy(string1, string_, size_);
    memcpy(string1 + size_, str.string_, str.size_ + 1);
    size_ += str.size_;
    delete[] string_;
    string_ = string1;
    return *this;
  }
  size_t Find(const String& str, bool reverse) const {
    if (size_ < str.size_ || str.size_ == 0) {
      return size_;
    }
    ssize_t start = reverse ? size_ - str.size_ : 0;
    ssize_t end = reverse ? 0 : size_ - str.size_;
    ssize_t step = reverse ? -1 : 1;

    for (ssize_t i = start; i != end; i += step) {
      size_t j = 0;
      for (j = 0; string_[i + j] == str.string_[j] && j < str.size_; ++j);
      if (j == str.size_) {
        return i;
      }
    }
    return size_;
  }
  size_t find(const String& str) const {
    return Find(str, false);
  }

  size_t rfind(const String& str) const {
    return Find(str, true);
  }

  String substr(size_t start, size_t count) const {
    int size = count;
    int capacity = count + 1;
    char* string = new char[count + 1];
    memcpy(string, string_ + start, count);
    string_[count] = '\0';
    return String(size, capacity, string);
  }

  bool empty() const {
    return size_ == 0;
  }

  void clear() {
    capacity_ = 2;
    size_ = 0;
    string_[0] = '\0';
  }

  void shrink_to_fit() {
    if (size_ < 2) {
      capacity_ = 2;
    }
    capacity_ = size_ + 1;
    char* str = new char[capacity_];
    memcpy(str, string_, size_ + 1);
    delete[] string_;
    string_ = str;
  }

  const char* data() const {
    return string_;
  }

  char* data() {
    return string_;
  }

  ~String() {
    delete[] string_;
  }
};

std::ostream& operator<<(std::ostream& out, const String& str) {
  return out << str.data();
}

std::istream& operator>>(std::istream& in, String& str) {
  str.clear();
  char ch;
  in.get(ch);
  while (!isspace(ch) && !in.eof()) {
    str.push_back(ch);
    in.get(ch);
  }
  return in;
}

String operator+(const String& ch, const String& str) {
  String new_str;
  new_str += ch;
  new_str += str;
  return new_str;
}

bool operator==(const String& str1, const String& str2) {
  if (str1.size() != str2.size()) return false;
  return !strcmp(str1.data(), str2.data());
}

bool operator!=(const String& str1, const String& str2) {
  return !(str1 == str2);
}

bool operator<(const String& str1, const String& str2) {
  if (str1.size() < str2.size()) return true;
  if (str1.size() > str2.size()) return true;
  return strcmp(str1.data(), str2.data()) < 0;
}

bool operator<=(const String& str1, const String& str2) {
  return (str1 < str2 || str1 == str2);
}

bool operator>(const String& str1, const String& str2) {
  return !(str1 <= str2);
}

bool operator>=(const String& str1, const String& str2) {
  return !(str1 < str2);
}
