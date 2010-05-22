{
  if ($1 == prevprev && $1 == prev && and(prevprev_flag,0x40) != 0 && and(prev_flag,0x40) == 0 && and($2,0x40) != 0) {
    print $1, prevprev_flag, prev_flag, $2
    exit 1;
  }
  prevprev = prev;
  prevprev_flag = prev_flag;
  prev = $1;
  prev_flag = $2;
}
