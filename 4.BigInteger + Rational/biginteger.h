#include <iostream>
#include <vector>

class BigInteger {
private:
  static const int base_ = 10e4;
  int sign_ = 0;
  std::vector<int> num_;
public:
  BigInteger() : sign_(0) {
    num_.push_back(0);
  }

  BigInteger(int num) : sign_(num > 0 ? 1 : -1) {
    if (num == 0) {
      num_.push_back(0);
      sign_ = 0;
      return;
    }
    num *= sign_;
    while (num > 0) {
      num_.push_back(num % base_);
      num /= base_;
    }
  }

  friend bool operator ==(const BigInteger&, const BigInteger&);

  friend bool operator <(const BigInteger&, const BigInteger&);

  BigInteger& operator +=(const BigInteger& bi) {
    if (*this == 0) {
      *this = bi;
      return *this;
    }
    if (bi == 0) {
      return *this;
    }
    if (sign_ == -1) {
      if (sign_ == -1) {
        sign_ = 1;
        *this += -(bi);
      }
      else return const_cast<BigInteger&>(bi) - (-*this);
    }
    else if (bi.sign_ == -1) return *this - (-bi);

    int in_mind = 0;
    for (size_t i = 0; i < std::max(num_.size(), bi.num_.size()) || in_mind; ++i) {
      if (i == num_.size())
        num_.push_back(0);
      num_[i] += in_mind + (i < bi.num_.size() ? bi.num_[i] : 0);
      in_mind = num_[i] >= base_;
      if (in_mind)
        num_[i] -= base_;
    }
    return *this;
  }

  BigInteger& operator -=(const BigInteger& bi) {
    //???
  }

  BigInteger& operator *=(const BigInteger& bi) {
    //???
  }

  BigInteger& operator /=(const BigInteger& bi) {
    //???
  }

  BigInteger& operator %=(const BigInteger& bi) {
    //???
  }

  BigInteger operator -() const {
    BigInteger bi(*this);
    bi.sign_ *= -1;
    return bi;
  }

  BigInteger& operator ++() {
    return *this += 1;
  }

  BigInteger& operator --() {
    return *this -= 1;
  }

  BigInteger operator ++(int) {
    BigInteger tmp = *this;
    ++(*this);
    return tmp;
  }

  BigInteger operator --(int) {
    BigInteger tmp = *this;
    --(*this);
    return tmp;
  }

  BigInteger string_to_bi(const std::string& str) {
    num_.clear();
    sign_ = 1;
    int k = 0;
    if (str[0] == '-') {
      sign_ = -1;
      k = 1;
    }
    for (int i = str.size() - 1; i >= k; i -= 4) {
      if (i < 4) {
        num_.push_back(std::stoi(str.substr(k, i + 1 - k)));
      } else {
        num_.push_back(std::stoi(str.substr(i - 3, 4)));
      }
    }
    return *this;
  }

  std::string toString() const {
    std::string str;
    if (sign_ == -1) str += "-";
    for (int i = num_.size() - 1; i >= 0; --i) {
      if (i == num_.size() - 1) {
        str += std::to_string(num_[i]);
      } else {
        if (num_[i] == 0)
          str += "0000";
        else if (num_[i] > 999)
          str += std::to_string(num_[i]);
        else if (num_[i] > 99)
          str += "0" + std::to_string(num_[i]);
        else if (num_[i] > 9)
          str += "00" + std::to_string(num_[i]);
        else if (num_[i] > 0)
          str += "000" + std::to_string(num_[i]);
      }
    }
    return str;
  }

  explicit operator bool() const {
    return num_.back() == 0;
  }


};

bool operator ==(const BigInteger& bi1, const BigInteger& bi2) {
  if (bi1.num_.size() != bi2.num_.size() || bi1.sign_ != bi2.sign_) return false;
  for (size_t i = 0; i < bi1.num_.size(); ++i) {
    if (bi1.num_[i] != bi2.num_[i]) {
      return false;
    }
  }
  return true;
}

bool operator !=(const BigInteger& bi1, const BigInteger& bi2) {
  return !(bi1 == bi2);
}

bool operator <(const BigInteger& bi1, const BigInteger& bi2) {
  //???
}

bool operator <=(const BigInteger& bi1, const BigInteger& bi2) {
  return (bi1 < bi2 || bi1 == bi2);
}

bool operator >(const BigInteger& bi1, const BigInteger& bi2) {
  return !(bi1 <= bi2);
}

bool operator >=(const BigInteger& bi1, const BigInteger& bi2) {
  return !(bi1 < bi2);
}

std::istream& operator >>(std::istream& in, BigInteger& bi) {
  std::string str;
  in >> str;
  bi.string_to_bi(str);
  return in;
}

std::ostream& operator <<(std::ostream& out, const BigInteger& bi) {
  out << bi.toString();
  return out;
}

BigInteger operator +(const BigInteger& bi1, const BigInteger& bi2) {
  BigInteger bi = bi1;
  bi += bi2;
  return bi;
}

BigInteger operator -(const BigInteger& bi1, const BigInteger& bi2) {
  BigInteger bi = bi1;
  bi -= bi2;
  return bi;
}
