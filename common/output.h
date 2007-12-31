/*	$Id$	*/

char *readtostr(uint32_t *read, u_int len, bool use_colours, int initbp);
void output_pretty(FILE *, const char *, const char *, const struct
    sw_full_results *, const char *, const char *, uint32_t *, uint32_t,
    uint32_t, bool, uint32_t *, u_int, int, bool);
void output_normal(FILE *, const char *, const char *, const struct
    sw_full_results *, uint32_t, uint32_t, bool, uint32_t *, u_int, int, bool,
    bool);
