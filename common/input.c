/*	$Id$	*/

#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

#include "../common/input.h"
#include "../common/util.h"

/* 
 * Trim leading/trailing whitespace.
 */
static char *
strtrim(char *str)
{
	char *ret;

	assert(str != NULL);

	while (isspace((int)*str) && *str != '\0')
		str++;

	ret = str;

	while (!isspace((int)*str) && *str != '\0')
		str++;
	
	*str = '\0';

	return (ret);
}

/*
 * Extract the 'key' and 'value' from 'key=value' in 'str'. Accomodates for
 * leading/trailing whitespace and removes any single outer double-quote pairs.
 *
 * 'str' is modified in place.
 *
 * Returns true if input was well-formatted and 'key' and 'val' are the
 * extracted values.
 */
static bool
parse_keyvalue(char *str, char **key, char **val)
{

	assert(str != NULL && key != NULL && val != NULL);

	while (isspace((int)*str) && *str != '\0')
		str++;

	if (*str == '\0')
		return (false);

	*key = str;

	while (*str != '=' && *str != '\0')
		str++;

	if (*str == '\0')
		return (false);

	*str++ = '\0';

	while (isspace((int)*str) && *str != '\0')
		str++;

	if (*str == '"') {
		char last;

		*val = str + 1;

		last = *str++;
		while (true) {
			while (*str != '"' && *str != '\0')
				last = *str++;

			/* Don't bail early on escaped quotes... */
			if (*str == '"' && last != '\\')
				break;
			else if (*str == '\0')
				break;

			last = *str++;
		}

		if (*str == '\0')
			return (false);

		*str = '\0';
	} else {
		*val = str;

		while (!isspace((int)*str) && *str != '\0')
			str++;

		*str = '\0';
		*val = strtrim(*val);
	}

	*key = strtrim(*key);

	return (true);
}

/* Remove "\=" and "\"" sequences added to escape whacky input files. */
static const char *
unescapestr(const char *str)
{
	static char *buf;
	static u_int buflen;

	int i, j;
	u_int l;

	l = strlen(str) + 1;
	if (buf == NULL || buflen < l) {
		if (buf != NULL)
			free(buf);
		buf = xmalloc(l);
		buflen = l;
	}

	for (i = j = 0; str[i] != '\0'; i++, j++) {
		if (str[i] == '\\' && (str[i+1] == '=' || str[i+1] == '"'))
			j--;
		else
			buf[j] = str[i];
	}
	buf[j] = '\0';

	assert(j < buflen);

	return (buf);
}

static void
handle_keyvalue(struct input *inp, char *key, char *val)
{
	bool error = false;

	assert(inp != NULL && key != NULL && val != NULL);

	switch (key[0]) {
	case 'd':
		if (strcmp(key, "dels") == 0)
			inp->deletions = strtoul(val, NULL, 0);
		else
			error = true;
		break;

	case 'g':
		if (strcmp(key, "g") == 0)
			inp->genome = xstrdup(unescapestr(val));
		else if (strcmp(key, "g_str") == 0) {
			if (*val == '-')
				inp->flags |= INPUT_FLAG_IS_REVCMPL;
		/* NB: internally 0 is first position, in output 1 is. adjust.*/
		} else if (strcmp(key, "g_start") == 0)
			inp->genome_start = strtoul(val, NULL, 0) - 1;
		else if (strcmp(key, "g_end") == 0)
			inp->genome_end = strtoul(val, NULL, 0) - 1;
		else
			error = true;
		break;

	case 'i':
		if (strcmp(key, "ins") == 0)
			inp->insertions = strtoul(val, NULL, 0);
		else
			error = true;
		break;

	case 'm':
		if (strcmp(key, "match") == 0)
			inp->matches = strtoul(val, NULL, 0);
		else
			error = true;
		break;

	case 'n':
		if (strcmp(key, "normodds") == 0) {
			inp->normodds = atof(val);
			inp->flags |= INPUT_FLAG_HAS_NORMODDS;
		} else
			error = true;
		break;

	case 'p':
		if (strcmp(key, "pgenome") == 0) {
			inp->pgenome = atof(val);
			inp->flags |= INPUT_FLAG_HAS_PGENOME;
		} else if (strcmp(key, "pchance") == 0) {
			inp->pchance = atof(val);
			inp->flags |= INPUT_FLAG_HAS_PCHANCE;
		} else
			error = true;
		break;

	case 'r':	
		if (strcmp(key, "r") == 0)
			inp->read = xstrdup(unescapestr(val));
		else if (strcmp(key, "r_seq") == 0)
			inp->read_seq = xstrdup(val);
		else if (strcmp(key, "r_start") == 0)
			inp->read_start = strtoul(val, NULL, 0) - 1;
		else if (strcmp(key, "r_end") == 0)
			inp->read_end = strtoul(val, NULL, 0) - 1;
		else if (strcmp(key, "r_len") == 0)
			inp->read_length = strtoul(val, NULL, 0);
		else
			error = true;
		break;

	case 's':
		if (strcmp(key, "score") == 0)
			inp->score = strtoul(val, NULL, 0);
		else if (strcmp(key, "subs") == 0)
			inp->mismatches = strtoul(val, NULL, 0);
		else
			error = true;
		break;

	case 'x':
		if (strcmp(key, "xovers") == 0)
			inp->crossovers = strtoul(val, NULL, 0);
		else
			error = true;
		break;

	default:
		error = true;
	}

	if (error)
		fprintf(stderr, "warning: unknown key in input (%s)\n", key);
}

/*
 * Parse the key-value paired output created by output_normal in output.c.
 * This will only parse lines beginning with '>', so it'll work just fine
 * with both rmapper's standard and pretty printed output formats.
 *
 * Returns false on EOF.
 */
bool
input_parseline(FILE *fp, struct input *inp)
{
	char buf[8192];

	assert(fp != NULL && inp != NULL);

	memset(inp, 0, sizeof(*inp));

	while (true) {
		/* XXX: Assume sizeof(buf) is always big enough. */
		if (fgets(buf, sizeof(buf), fp) == NULL)
			return (false);

		if (buf[0] == '>') {
			int i;
			char *key, *val;
			
			i = strlen(buf);

			while (--i >= 0) {
				/* Be careful with escapes... */
				if (buf[i] == '=' && i > 0 &&
				    buf[i-1] != '\\') {
					while (--i >= 0) {
						if (!isspace((int)buf[i]))
							break;
					}
					while (--i >= 0) {
						if (isspace((int)buf[i]))
							break;
					}
					if (i < 0)
						break;
					if (!parse_keyvalue(&buf[i],&key,&val))
						return (false);

					if (i == 0 && *key == '>')
						key++;
					handle_keyvalue(inp, key, val);
				}
			}

			if (i < 0)
				break;
		}
	}

	return (true);
}
