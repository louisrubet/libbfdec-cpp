#include "libbfdec++.h"

bf_context_t Bfdec::bf_ctx{.realloc_func = Bfdec::bf_realloc};
