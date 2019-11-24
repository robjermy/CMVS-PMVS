#include <stdexcept>
#include <sstream>

class ExceptionFormatter {
public:
  ExceptionFormatter() {}
  ~ExceptionFormatter() {}

  template<typename T>
  ExceptionFormatter& operator<<(const T& value) {
    _stream << value;
    return *this;
  }

  std::string str() const {
    return _stream.str();
  }

  operator std::string() const {
    return _stream.str();
  }

private:
  std::stringstream _stream;
};
