// fake quickjs.h for libbf
#pragma once

#include <assert.h>
#define QJS_ASSERT(something) { if (!(something)) printf("%s:%d - %s: assert\n", __FILE__, __LINE__, __FUNCTION__); assert(something); }
#define QJS_ABORT() { printf("%s:%d - %s: abort\n", __FILE__, __LINE__, __FUNCTION__); assert(0); }
static int IsDebuggerPresent() { return 0; }
static void DebugBreak() { }
