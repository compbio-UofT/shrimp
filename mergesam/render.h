#ifndef __RENDER_H__
#define __REDNER_H__

#include "sam2pretty_lib.h"

size_t render_sam_unaligned_bound(pretty * pa );
size_t render_sam_unaligned_string(pretty * pa, char * buffer, size_t buffer_size);
void render_sam_unaligned(pretty * pa,bool inplace);

size_t render_fastx_bound(pretty * pa );
size_t render_fastx_string(pretty * pa, char * buffer, size_t buffer_size);
void render_fastx(pretty* pa, bool inplace );

size_t render_sam_bound(pretty * pa );
size_t render_sam_string(pretty * pa, char * buffer, size_t buffer_size);
void render_sam(pretty * pa, bool inplace);

#endif
