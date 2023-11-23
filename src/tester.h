#pragma once

#include <fcntl.h>
#include <math.h>
#include <unistd.h>

#include <chrono>
#include <cstdio>
#include <exception>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#ifndef STDOUT_REDIRECT
#define STDOUT_REDIRECT "unit_tests_output.txt"
#endif

class TestSuite;

class TestResult {
 public:
  TestResult() : failed_(0), passed_(0), asserted_(0), exceptions_(0) {}

  void summary() {
    std::cout << "Summary" << std::endl
              << passed_ << " assertions passed out of " << asserted_
              << " with " << failed_ << " and " << exceptions_ << " exceptions"
              << std::endl;
  }

  bool successful() const { return (failed_ == 0 && exceptions_ == 0); }
  unsigned long n_failed() const { return failed_; }
  unsigned long n_assert() const { return asserted_; }
  unsigned long n_exceptions() const { return exceptions_; }

  void failed() { failed_++; }
  void passed() { passed_++; }
  void asserted() { asserted_++; }
  void exception() { exceptions_++; }

 private:
  unsigned long failed_;
  unsigned long passed_;
  unsigned long asserted_;
  unsigned long exceptions_;
};

class TestCase {
 public:
  TestCase(const char* name, TestSuite* suite) : name_(name), suite_(suite) {}

  virtual ~TestCase() {}

  virtual void run(TestResult& __result__) = 0;

  const std::string& name() { return name_; }

 protected:
  std::string name_;
  TestSuite* suite_;
};

class TestDriver;

class TestSuite {
 public:
  TestSuite(const char* name) : name_(name) {}

  TestSuite(const char* name, TestDriver* driver);

  const std::string& name() const { return name_; }

  void add(TestCase* t) { cases_.push_back(t); }

  std::size_t ntests() const { return cases_.size(); }

  int run(TestResult& __result__) {
    unsigned long n_failed0 = __result__.n_failed();
    (void)(n_failed0);  // suppresses the warning of unused n_failed0
#ifndef STANDALONE

#ifndef STDOUT_REDIRECT
#error "must define where stdout is redirected"
#endif

    // redirection for stdout
    fflush(stdout);
    int fd = dup(fileno(stdout));
    FILE* x = freopen(STDOUT_REDIRECT, "a", stdout);
    (void)(x);  // suppresses the warning of unused x

    // redirection for cout
    std::stringstream out;
    std::streambuf* coutbuf = std::cout.rdbuf();  // save old buf
    std::cout.rdbuf(out.rdbuf());                 // redirect std::cout to out

    std::cout << "stdout output:" << std::endl;

#endif

    std::cout << "\n\n==== Test suite: " << name_ << " ====\n\n";
    std::cout << "running " << cases_.size() << " tests..." << std::endl;

    auto t0 = std::chrono::high_resolution_clock::now();
    for (std::size_t i = 0; i < cases_.size(); i++) {
      try {
        std::cout << "\n--> running case " << cases_[i]->name() << "\n\n";
        cases_[i]->run(__result__);
      } catch (std::exception& e) {
        std::cout << "unexpected exception!!" << std::endl;
        std::cout << e.what() << std::endl;
        __result__.exception();
      }
    }
    auto t1 = std::chrono::high_resolution_clock::now();
    auto t =
        std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cout << "\n ==== Total time = " << t / 1000.0 << " seconds. ====\n";

#ifndef STANDALONE

    // write the cout re-direction to the file too
    std::flush(std::cout);
    fprintf(stdout, "\nstd::cout output:\n\n%s", out.str().c_str());

    fflush(stdout);
    dup2(fd, fileno(stdout));
    close(fd);

    // reset to standard output
    std::cout.rdbuf(coutbuf);

    if (__result__.n_failed() > n_failed0)
      std::cout << "suite " << name_ << " failed with "
                << __result__.n_failed() - n_failed0 << " new errors!"
                << std::endl;

#endif

    if (__result__.successful()) return 0;
    return -1;
  }

 private:
  std::string name_;
  std::vector<TestCase*> cases_;
};

class TestDriver {
 public:
  void add(TestSuite* suite) { suites_.push_back(suite); }

  int run(TestResult& __result__) {
    int result = 0;
    for (std::size_t i = 0; i < suites_.size(); i++) {
      int n_assert0 = __result__.n_assert();
      std::cout << "running " << std::setw(32) << std::left
                << suites_[i]->name() << " with " << std::setw(3)
                << suites_[i]->ntests() << " tests ... " << std::endl;
      ;
      auto t0 = std::chrono::high_resolution_clock::now();
      try {
        if (suites_[i]->run(__result__) < 0) result = -1;
      } catch (...) {
        result = -1;
      }
      int n_assert = __result__.n_assert() - n_assert0;
      auto t1 = std::chrono::high_resolution_clock::now();
      auto t = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0)
                   .count();

      std::cout << std::right << std::setw(16) << "--> done: " << std::setw(6)
                << n_assert << " assertions in " << std::setw(6) << t / 1000.0
                << " seconds." << std::endl;
    }
    return result;
  }

 private:
  std::vector<TestSuite*> suites_;
};
extern TestDriver* __driver__;

inline TestSuite::TestSuite(const char* name, TestDriver* driver)
    : TestSuite(name) {
  driver->add(this);
}

#define UT_ASSERT_EQUALS(X, Y)                                                 \
  do {                                                                         \
    __result__.asserted();                                                     \
    if (X != Y) {                                                              \
      __result__.failed();                                                     \
      std::cout << suite_->name() << "::" << name_ << " assertion " << #X      \
                << " = " << #Y << " failed on line " << __LINE__               \
                << " of file " << __FILE__ << std::endl;                       \
      std::cout << "--> expected " << Y << " but received " << X << std::endl; \
    } else                                                                     \
      __result__.passed();                                                     \
  } while (0)

#define UT_ASSERT(X)                                                         \
  do {                                                                       \
    __result__.asserted();                                                   \
    if (!(X)) {                                                              \
      __result__.failed();                                                   \
      std::cout << suite_->name() << "::" << name_ << " assertion " << #X    \
                << " failed on line " << __LINE__ << " of file " << __FILE__ \
                << std::endl;                                                \
    } else                                                                   \
      __result__.passed();                                                   \
  } while (0)

#define UT_ASSERT_NEAR(X, Y, Z)                                             \
  do {                                                                      \
    __result__.asserted();                                                  \
    if (fabs(Y - X) > Z) {                                                  \
      __result__.failed();                                                  \
      std::cout << suite_->name() << "::" << __FUNCTION__ << ": assertion " \
                << #X << " (" << X << ") == " << #Y << "(" << Y << ") +/- " \
                << Z << "failed on line " << __LINE__ << " of file "        \
                << __FILE__ << std::endl;                                   \
    } else                                                                  \
      __result__.passed();                                                  \
  } while (0)

#define UT_ASSERT_CLOSE(X, Y, Z, Z0) UT_ASSERT_NEAR(X, Y, Z)
#define UT_ASSERT_SMALL(X, Z) UT_ASSERT_NEAR(X, 0, Z)

#define UT_CATCH_EXCEPTION(X)                                                \
  do {                                                                       \
    __result__.asserted();                                                   \
    try {                                                                    \
      (X);                                                                   \
      __result__.failed();                                                   \
      std::cout << suite_->name() << "::" << name_                           \
                << " expected exception on line " << __LINE__ << " of file " \
                << __FILE__ << std::endl;                                    \
    } catch (...) {                                                          \
      __result__.passed();                                                   \
      break;                                                                 \
    }                                                                        \
  } while (0)

#ifdef STANDALONE
#define UT_TEST_SUITE(X) \
  namespace suite_##X {  \
    TestSuite __suite__(#X);
#else
#define UT_TEST_SUITE(X) \
  namespace suite_##X {  \
    static TestSuite __suite__(#X, __driver__);
#endif

#define UT_TEST_CASE(X)                                   \
  class X : public TestCase {                             \
   public:                                                \
    X() : TestCase(#X, &__suite__) { suite_->add(this); } \
    void run(TestResult& __result__)

#define UT_TEST_CASE_SKIP(X)                              \
  class X : public TestCase {                             \
   public:                                                \
    X() : TestCase(#X, &__suite__) { suite_->add(this); } \
    void run(TestResult& __result__) {                    \
      std::cout << "skipping test!" << std::endl;         \
    }                                                     \
    void run_skip(TestResult& __result__)

#define UT_TEST_CASE_END(X) \
  }                         \
  ;                         \
  static X instance_##X;

#ifdef STANDALONE
#define UT_TEST_SUITE_END(X)                                          \
  }                                                                   \
  int main(int argc, char* argv[]) {                                  \
    TestResult __result__;                                            \
    suite_##X::__suite__.run(__result__);                             \
    __result__.summary();                                             \
    if (__result__.n_failed() != 0 || __result__.n_exceptions() != 0) \
      return 1;                                                       \
    return 0;                                                         \
  }
#else
#define UT_TEST_SUITE_END(X) }
#endif
