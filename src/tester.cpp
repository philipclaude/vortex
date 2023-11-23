#ifdef STANDALONE
#error "this should be used for entire unit tests"
#endif

#include "tester.h"

TestDriver* __driver__ = new TestDriver;

int main() {
#ifdef STDOUT_REDIRECT
  FILE* fid = fopen(STDOUT_REDIRECT, "w");
  fclose(fid);
#endif
  TestResult result;
  __driver__->run(result);
  result.summary();
  delete __driver__;
  if (!result.successful()) return 1;
  return 0;
}
