
#pragma once

#include <iostream>
#include <string>

#include <geogram/basic/memory.h>

namespace GEO {

struct Logger {

  static void initialize() {};
  static void terminate() {};

  static std::ostream& out(const std::string& name)
  {
    return std::cout << " [" << name << "]";
  }

  static std::ostream& err(const std::string& name)
  {
    return std::cerr << "E[" << name << "]";
  }

  static std::ostream& warn(const std::string& name)
  {
    return std::cerr << "W[" << name << "]";
  }

  static void div(const std::string& name){}

  static Logger* instance() { return &sLogger; }

  bool is_quiet() const { return true; }

private:
  static Logger sLogger;
};

}
