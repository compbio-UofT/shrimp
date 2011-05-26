#include "my-alloc.h"

size_t max_mem = 0;
size_t crt_mem = 0;
size_t alert_mem = 0;

bool my_alloc_initialized = false;
bool warned_max = false;
bool warned_fail = false;
