{
  idx = match($1, SUFFIX)
  if (idx == 0)
    next
  base = substr($1, 1, idx-1)
  if ($1 == prevprev && $1 != prev && base == prev_base ) {
    print prevprev, prev, $1
    exit 1;
  }
  prevprev = prev;
  prev = $1;
  prev_base = base;
}
