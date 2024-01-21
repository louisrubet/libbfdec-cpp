// fake quickjs.h for libbf
#pragma once

#include <assert.h>
static void QJS_ASSERT(int something) { assert(something); }
static void QJS_ABORT() { assert(0); }
static int IsDebuggerPresent() { return 0; }
static void DebugBreak() { }
