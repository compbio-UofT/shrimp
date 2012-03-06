/*	$Id: output.h,v 1.5 2008/06/06 18:24:09 rumble Exp $	*/

char *readtostr(const uint32_t *, u_int, bool, int);
char *output_format_line(bool);
char *output_pretty(const char *, const char *, const struct
    sw_full_results *, uint32_t *, uint32_t, bool, uint32_t *, u_int,int, bool);
char *output_normal(const char *, const char *, const struct
    sw_full_results *, uint32_t, bool, uint32_t *, u_int, int, bool, bool);
char *alignment_edit_string(char const *, char const *);
