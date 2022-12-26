#!/usr/bin/awk -f

## usage: set_cols(array)
##   sets the following values in "array" with tput. printing them will format
##   any text afterwards. colors and formats are:
##     bold - bold text (can be combined with a color)
##     black - black text
##     red - red text
##     green - green text
##     yellow - yellow text
##     blue - blue text
##     magenta - magenta text
##     cyan - cyan text
##     white - white text
##     reset - resets to default settings
function set_cols(array) {
  # bold
  cmd = "tput bold";
  cmd | getline array["bold"];
  close(cmd);
  # black
  cmd = "tput setaf 0";
  cmd | getline array["black"];
  close(cmd);
  # red
  cmd = "tput setaf 1";
  cmd | getline array["red"];
  close(cmd);
  # green
  cmd = "tput setaf 2";
  cmd | getline array["green"];
  close(cmd);
  # yellow
  cmd = "tput setaf 3";
  cmd | getline array["yellow"];
  close(cmd);
  # blue
  cmd = "tput setaf 4";
  cmd | getline array["blue"];
  close(cmd);
  # magenta
  cmd = "tput setaf 5";
  cmd | getline array["magenta"];
  close(cmd);
  # cyan
  cmd = "tput setaf 6";
  cmd | getline array["cyan"];
  close(cmd);
  # white
  cmd = "tput setaf 7";
  cmd | getline array["white"];
  close(cmd);
  # reset
  cmd = "tput sgr0";
  cmd | getline array["reset"];
  close(cmd);
}



# You can do whatever you want with this stuff, but a thanks is always
# appreciated
#!/usr/bin/awk -f

## usage: create_line(array, max [, sep [, qualifier [, quote_type] ] ])
## Generates an output line in quoted CSV format, from the contents of "array"
## "array" is expected to be an indexed array (1-indexed). "max" is the highest
## index to be used. "sep", if provided, is the field separator. If it is more
## than one character, the first character in the string is used. By default,
## it is a comma. "qualifier", if provided, is the quote character. Like "sep",
## it is one character. The default value is `"'. "quote_type", if provided, is
## used to determine how the output fields are quoted. Valid values are given
## below. For example, the array: a[1]="foo"; a[2]="bar,quux"; a[3]="blah\"baz"
## when called with create_line(a, 3), will return: "foo","bar,quux","blah""baz"
## note: expects a non-sparse array. empty or unset values will become
## empty fields
## Valid values for "quote_type":
##   "t": Quote all strings, do not quote numbers. This is the default
##   "a": Quote all fields
##   "m": Only quote fields with commas or quote characters in them
function create_line(arr, len, sep, q, type,    i, out, c, new) {
  # set "sep" if the arg was provided, using the first char
  if (length(sep)) {
    sep = substr(sep, 1, 1);
  # default
  } else {
    sep = ",";
  }

  # validate "type"
  if (!length(type) || type !~ /^[tam]$/) {
    type = "t";
  }

  # set "q" if the arg was provided, using the first char
  if (length(q)) {
    q = substr(q, 1, 1);
  # default
  } else {
    q = "\"";
  }

  # empty the output string
  out = "";

  # iterate over the array elements
  for (i=1; i<=len; i++) {
    # determine if the output string needs to be quoted
    toquote = 0;
    if (type == "t") {
      if (arr[i] ~ /[^0-9.]/ || index(arr[i], sep) || index(arr[i], q)) {
        toquote = 1;
      }
    } else if (type == "a") {
      toquote = 1;
    } else {
      if (index(arr[i], sep) || index(arr[i], q)) {
        toquote = 1;
      }
    }

    # create output string
    if (toquote) {
      new = "";
      while (c = index(arr[i], q)) {
        new = new substr(arr[i], 1, c - 1) q q;
        arr[i] = substr(arr[i], c + 1);
      }
      new = new arr[i];

      # quote escaped string, add to output with sep
      out = (i > 1) ? out sep q new q : q new q;

      # no quotes needed, just add to output with sep
    } else {
      out = (i > 1) ? out sep arr[i] : arr[i];
    }
  }

  # return output string
  return out;
}

## usage: qsplit(string, array [, sep [, qualifier] ])
## a version of split() designed for CSV-like data. splits "string" on "sep"
## (,) if not provided, into array[1], array[2], ... array[n]. returns "n", or
## "-1 * n" if the line is incomplete (it has an uneven number of quotes). both
## "sep" and "qualifier" will use the first character in the provided string.
## uses "qualifier" (" if not provided) and ignores "sep" within quoted fields.
## doubled qualifiers are considered escaped, and a single qualifier character
## is used in its place. for example, foo,"bar,baz""blah",quux will be split as
## such: array[1] = "foo"; array[2] = "bar,baz\"blah"; array[3] = "quux";
function qsplit(str, arr, sep, q,    a, len, cur, isin, c) {
  delete arr;

  # set "sep" if the argument was provided, using the first char
  if (length(sep)) {
    sep = substr(sep, 1, 1);
  # otherwise, use ","
  } else {
    sep = ",";
  }

  # set "q" if the argument was provided, using the first char
  if (length(q)) {
    q = substr(q, 1, 1);
  # otherwise, use '"'
  } else {
    q = "\"";
  }

  # split the string into the temporary array "a", one element per char
  len = split(str, a, "");

  # "cur" contains the current element of 'arr' the function is assigning to
  cur = 1;
  # boolean, whether or not the iterator is in a quoted string
  isin = 0;
  # iterate over each character
  for (c=1; c<=len; c++) {
    # if the current char is a quote...
    if (a[c] == q) {
      # if the next char is a quote, and the previous character is not a
      # delimiter, it's an escaped literal quote (allows empty fields 
      # that are quoted, such as "foo","","bar")
      if (a[c+1] == q && a[c-1] != sep) {
        arr[cur] = arr[cur] a[c];
        c++;

      # otherwise, it's a qualifier. switch boolean
      } else {
        isin = ! isin;
      }

    # if the current char is the separator, and we're not within quotes
    } else if (a[c] == sep && !isin) {
      # increment array element
      cur++;

    # otherwise, just append to the current element
    } else {
      arr[cur] = arr[cur] a[c];
    }
  }

  # return length
  return cur * (isin ? -1 : 1);
}



# You can do whatever you want with this stuff, but a thanks is always
# appreciated
#!/usr/bin/awk -f

## usage: abs(number)
## returns the absolute value of "number"
function abs(num) {
  return num < 0 ? -num : num;
}

## usage: ceil(number)
## returns "number" rounded UP to the nearest int
function ceil(num) {
  if (num < 0) {
    return int(num);
  } else {
    return int(num) + (num == int(num) ? 0 : 1);
  }
}

## usage: ceiling(multiple, number)
## returns "number" rounded UP to the nearest multiple of "multiple"
function ceiling(mult, num,    r) {
  return (r = num % mult) ? num + (mult - r) : num;
}

## usage: change_base(number, start_base, end_base)
## converts "number" from "start_base" to "end_base"
## bases must be between 2 and 64. the digits greater than 9 are represented
## by the lowercase letters, the uppercase letters, @, and _, in that order.
## if ibase is less than or equal to 36, lowercase and uppercase letters may
## be used interchangeably to represent numbers between 10 and 35.
## returns 0 if any argument is invalid
function change_base(num, ibase, obase,
                     chars, c, l, i, j, cur, b10, f, fin, isneg) {
  # convert number to lowercase if ibase <= 36
  if (ibase <= 36) {
    num = tolower(num);
  }

  # determine if number is negative. if so, set isneg=1 and remove the '-'
  if (sub(/^-/, "", num)) {
    isneg = 1;
  }

  # determine if inputs are valid
  if (num ~ /[^[:xdigit:]]/ || ibase != int(ibase) || obase != int(obase) ||
      ibase < 2 || ibase > 64 || obase < 2 || obase > 64) {
    return 0;
  }

  # set letters to numbers conversion array
  if (ibase > 10 || obase > 10) {
    # set chars[] array to convert letters to numbers
    c = "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ@_";
    l = length(c);

    j = 10;
    for (i=1; i<=l; i++) {
      cur = substr(c, i, 1);
      chars[cur] = j;
      chars[j] = cur;

      j++;
    }
  }
  
  # convert to base 10
  if (ibase != 10) { 
    l = length(num);

    j = b10 = 0;
    for (i=l; i>0; i--) {
      c = substr(num, i, 1);

      # if char is a non-digit convert to dec
      if (c !~ /^[0-9]$/) {
        c = chars[c];
      }

      # check to make sure value isn't too great for base
      if (+c >= ibase) {
        return 0;
      }

      b10 += c * (ibase ^ j++);
    }
  } else {
    # num is already base 10
    b10 = num;
  }
  
  # convert from base 10 to obase
  if (obase != 10) {
    # build number backwards
    j = 0;
    do {
      f[++j] = (c = b10 % obase) > 9 ? chars[c] : c;
      b10 = int(b10 / obase);
    } while (b10);

    # reverse number
    fin = f[j];
    for (i=j-1; i>0; i--) {
      fin = fin f[i];
    }
  } else {
    # num has already been converted to base 10
    fin = b10;
  }

  # add '-' if number was negative
  if (isneg) {
    fin = "-" fin;
  }

  return fin;
}

## usage: format_num(number)
## adds commas to "number" to make it more readable. for example,
## format_num(1000) will return "1,000", and format_num(123456.7890) will
## return "123,456.7890". also trims leading zeroes
## returns 0 if "number" is not a valid number
function format_num(num,    is_float, b, e, i, len, r, out) {
  # trim leading zeroes
  sub(/^0+/, "", num);

  # make sure "num" is a valid number
  if (num ~ /[^0-9.]/ || num ~ /\..*\./) {
    return 0;
  }
  
  # if "num" is not an int, split it into pre and post decimal parts.
  # use sub() instead of int() because int() can be funny for float arithmetic
  # results
  if (num ~ /\./) {
    is_float = 1; # flag "num" as a float
    b = e = num;
    sub(/\..*/, "", b);
    sub(/.*\./, "", e);

  # otherwise, just assign the number to "b"
  } else {
    is_float = 0;
    b = num;
  }

  len = length(b)

  # only do anything if the pre-decimal section is greater than 3 digits
  if (len < 3) {
    return num;
  }

  # start by assigning the last 3 pre-decimal digits to out
  out = substr(b, len - 2);

  # loop backwards over each grouping of 3 numbers after that, prepending
  # each to out (with a comma)
  for (i=len-5; i>0; i-=3) {
    out = substr(b, i, 3) "," out;
  }

  # if the length is not a multiple of 3, prepend the remaining digits
  if (r = len % 3) {
    out = substr(b, 1, r) "," out;
  }

  # if number was a float, add the post-decimal digits back on
  if (is_float) {
    out = out "." e;
  }

  # return the formatted number
  return out;
}

## usage: str_to_num(string)
## examines "string", and returns its numeric value. if "string" begins with a
## leading 0, assumes that "string" is an octal number. if "string" begins with
## a leading "0x" or "0X", assumes that "string" is a hexadecimal number.
## otherwise, decimal is assumed.
function str_to_num(str,    base, isneg, l, i, j, chars, c, num) {
  # convert to all lowercase
  str = tolower(str);

  # determine if number is negative. if so, set isneg=1 and remove the '-'
  if (sub(/^-/, "", num)) {
    isneg = 1;
  }

  # examine the string, to determine the base and trim said base information
  if (sub(/^0x/, "", str)) {
    base = 16;
  } else if (sub(/^0/, "", str)) {
    base = 8;
  } else {
    base = 10;
  }

  # trim everything from the first non-number character to the end
  if (base == 16) {
    sub(/[^[:xdigit:]].*/, "", str);
  } else {
    sub(/[^[:digit:]].*/, "", str);
  }

  # if the base is octal, but there's a number >= 8, set it to decimal instead
  if (base == 8 && str ~ /[89]/) {
    base = 10;
  }

  # don't need to convert if the base is 10
  if (base == 10) {
    return isneg ? -str : +str;
  }

  # set letters for hex
  if (base == 16) {
    chars["a"] = 10; chars["b"] = 11; chars["c"] = 12;
    chars["d"] = 13; chars["e"] = 14; chars["f"] = 15;
  }

  # convert to base 10
  l = length(str);

  j = num = 0;
  for (i=l; i>0; i--) {
    c = substr(str, i, 1);

    # if char is a non-digit convert to dec
    if (c !~ /^[0-9]$/) {
      c = chars[c];
    }

    num += c * (base ^ j++);
  }
  
  # return the number
  return isneg ? -num : +num;
}

## usage: floor(multiple, number)
## returns "number" rounded DOWN to the nearest multiple of "multiple"
function floor(mult, num) {
  return num - (num % mult);
}

## usage: round(multiple, number)
## returns "number" rounded to the nearest multiple of "multiple"
function round(mult, num,    r) {
  if (num % mult < mult / 2) {
    return num - (num % mult);
  } else {
    return (r = num % mult) ? num + (mult - r) : num;
  }
}

## usage: rint(number)
## returns "number" rounded to the nearest integer
function rint(num,    n) {
  if (num < 0) {
    return (num - (n = int(num)) < -.5) ? n - 1 : n;
  } else {
    return (num - (n = int(num)) >= .5) ? n + 1 : n;
  }
}

## usage: isint(string)
## returns 1 if "string" is a valid integer, otherwise 0
function isint(str) {
  if (str !~ /^-?[0-9]+$/) {
    return 0;
  }

  return 1;
}

## usage: isnum(string)
## returns 1 if "string" is a valid number, otherwise 0
function isnum(str) {
  # use a regex comparison because 'str == str + 0' has issues with some floats
  if (str !~ /^-?[0-9.]+$/ || str ~ /\..*\./) {
    return 0;
  }

  return 1;
}

## usage: isprime(number)
## returns 1 if "number" is a prime number, otherwise 0. "number" must be a
## positive integer greater than one
function isprime(num,    i, s) {
  # check to make sure "num" is a valid positive int (and not 1)
  if (num !~ /^[0-9]+$/ || num <= 1) {
    return 0;
  }

  # 1, 2, and 3 are prime
  if (num <= 3) {
    return 1;
  }
  
  # check if even or divisible by 3
  if (!(num % 2) || !(num % 3)) {
    return 0;
  }
  
  # use naive method, fermats little theorem had overflow and did not work
  # for primes larger than 1021
  s = sqrt(num);
  for (i=5; i<=s; i+=2) {
    if (!(num % i)) {
      return 0;
    }
  }

  return 1;
}

## usage: gcd(a, b)
## returns the greatest common denominator (greatest common factor) of a and b.
## both a and b must be positive integers. uses the recursive euclid algorithm.
function gcd(a, b,    f) {
  # check to make sure both numbers are positive ints
  if (!f) {
    if (a !~ /^[0-9]+$/ || !a || b !~ /^[0-9]+$/ || !b) {
      return 0;
    }
  }

  if (b) {
    return gcd(b, a % b, 1);

  } else {
    # return the absolute value
    return a < 0 ? -a : a;
  }
}

## usage: lcm(a, b)
## returns the least common multiple of a and b. both a and b must be positive
## integers.
function lcm(a, b,    m, l) {
  # check to make sure both numbers are positive ints
  if (a !~ /^[0-9]+$/ || !a || b !~ /^[0-9]+$/ || !b) {
    return 0;
  }

  m = 0;
  while ((l = ++m * a) % b);

  return l;
}

## usage: calc_e()
## approximates e by calculating the sumation from k=0 to k=50 of 1/k!
## returns 10 decimal places
function calc_e(lim,    e, k, i, f) {
  for (k=0; k<=50; k++) {
    # calculate factorial
    f = 1;
    for (i=1; i<=k; i++) {
      f = f * i;
    }

    # add to e
    e += 1 / f;
  }

  return sprintf("%0.10f", e);
}


## usage: calc_pi()
## returns pi, with an accuracy of 10 decimal places
function calc_pi() {
  return sprintf("%0.10f", atan2(0, -1));
}

## usage: calc_tau()
## returns tau, with an accuracy of 10 decimal places
function calc_tau() {
  return sprintf("%0.10f", 2 * atan2(0, -1));
}

## usage: deg_to_rad(degrees)
## converts degrees to radians
function deg_to_rad(deg,    tau) {
  tau = 8 * atan2(1,1);

  return (deg/360) * tau;
}

## usage: rad_to_deg(radians)
## converts radians to degrees
function rad_to_deg(rad,    tau) {
  tau = 8 * atan2(1,1);

  return (rad/tau) * 360;
}

## usage: tan(expr)
## returns the tangent of expr, which is in radians
function tan(ang) {
  return sin(ang)/cos(ang);
}

## usage: csc(expr)
## returns the cosecant of expr, which is in radians
function csc(ang) {
  return 1/sin(ang);
}

## usage: sec(expr)
## returns the secant of expr, which is in radians
function sec(ang) {
  return 1/cos(ang);
}

## usage: cot(expr)
## returns the cotangent of expr, which is in radians
function cot(ang) {
  return cos(ang)/sin(ang);
}



# You can do whatever you want with this stuff, but a thanks is always
# appreciated
#!/usr/bin/awk -f

# comparison function
# usage: __compare(a, b, how)
# compares "a" and "b" based on "how", returning 0 for false and 1 for true.
# required for all of the msort() functions below
function __compare(a, b, how) {
  # standard comparisons
  if (how == "std asc") {
    return a < b;
  } else if (how == "std desc") {
    return a > b;

  # force string comps
  } else if (how == "str asc") {
    return "a" a < "a" b;
  } else if (how == "str desc") {
    return "a" a > "a" b;

  # force numeric
  } else if (how == "num asc") {
    return +a < +b;
  } else if (how == "num desc") {
    return +a > +b;
  }
}

# actual sorting function
# usage: __mergesort(array, len, how)
# sorts the values in "array" in-place, from indices 1 to "len", based
# on the comparison mode "how" (see the msort() description).
# required for all of the msort() functions below
function __mergesort(array, len, how,
                     tmpa, alen, a, tmpb, blen, b, half, cur, pos, tmp) {
  # if there are 10 elements or less, use an insertion sort and return
  if (len <= 10) {
    # loop over each item, starting with the second
    for (cur=2; cur<=len; cur++) {
      pos = cur;
      # shift the item down the list into position
      while (pos > 1 && __compare(array[pos], array[pos-1], how)) {
        tmp = array[pos];
        array[pos] = array[pos-1];
        array[pos-1] = tmp;

        pos--;
      }
    }

    # return
    return len;
  }

  # determine the halfway point of the indices
  half = int(len / 2);

  # create temp arrays of the two halves
  a = 0;
  for (i=1; i<=half; i++) {
    tmpa[++a] = array[i];

    # remove the index from the original array
    delete array[i];
  }
  b = 0;
  for (i=half+1; i<=len; i++) {
    tmpb[++b] = array[i];

    # remove the index from the original array
    delete array[i];
  }

  # sort the two halves with recursive calls
  alen = __mergesort(tmpa, a, how);
  blen = __mergesort(tmpb, b, how);

  # merge the two halves
  len = 0;
  a = b = 1;
  # loop while there is still an element in either array
  while (a <= alen || b <= blen) {
    # a sorts first
    if (a <= alen && (b > blen || __compare(tmpa[a], tmpb[b], how))) {
      array[++len] = tmpa[a];
      delete tmpa[a++]; # remove the index from the temp array

    # b sorts first
    } else {
      array[++len] = tmpb[b];
      delete tmpb[b++]; # remove the index from the temp array
    }
  }

  # return the length
  return len;
}

# actual sorting function for the msortv() function
# usage: __mergesortv(array, values, len, how)
# sorts the values in "array" on the original values in "values", from indices
# 1 through "len", based on the comparison mode "how" (see the msortv()
# description). required for all of the msortv() functions below
function __mergesortv(array, values, len, how,
                      tmpa, tmpva, alen, a, tmpb, tmpvb, blen, b,
                      half, cur, pos, tmp) {
  # if there are 10 elements or less, use an insertion sort and return
  if (len <= 10) {
    # loop over each item, starting with the second
    for (cur=2; cur<=len; cur++) {
      pos = cur;
      # shift the item down the list into position
      while (pos > 1 && __compare(values[pos], values[pos-1], how)) {
        tmp = array[pos];
        array[pos] = array[pos-1];
        array[pos-1] = tmp;
        tmp = values[pos];
        values[pos] = values[pos-1];
        values[pos-1] = tmp;

        pos--;
      }
    }

    # return
    return len;
  }

  # determine the halfway point of the indices
  half = int(len / 2);

  # create temp arrays of the two halves
  a = 0;
  for (i=1; i<=half; i++) {
    tmpa[++a] = array[i];
    tmpva[a] = values[i];

    # remove the index from the original array
    delete array[i];
  }
  b = 0;
  for (i=half+1; i<=len; i++) {
    tmpb[++b] = array[i];
    tmpvb[b] = values[i];

    # remove the index from the original array
    delete array[i];
  }

  # sort the two halves with recursive calls
  alen = __mergesortv(tmpa, tmpva, a, how);
  blen = __mergesortv(tmpb, tmpvb, b, how);

  # merge the two halves
  len = 0;
  a = b = 1;
  # loop while there is still an element in either array
  while (a <= alen || b <= blen) {
    # a sorts first
    if (a <= alen && (b > blen || __compare(tmpva[a], tmpvb[b], how))) {
      array[++len] = tmpa[a];
      values[len] = tmpva[a];
      delete tmpva[a];
      delete tmpa[a++]; # remove the index from the temp array

    # b sorts first
    } else {
      array[++len] = tmpb[b];
      values[len] = tmpvb[b];
      delete tmpvb[b];
      delete tmpb[b++]; # remove the index from the temp array
    }
  }

  # return the length
  return len;
}



## usage: msort(s, d [, how])
## sorts the elements in the array "s" using awk's normal rules for comparing
## values, creating a new sorted array "d" indexed with sequential integers
## starting with 1. returns the length, or -1 if an error occurs.. leaves the
## indices of the source array "s" unchanged. the optional string "how" controls
## the direction and the comparison mode. uses the merge sort algorithm, with an
## insertion sort when the list size gets small enough. this is not a stable
## sort. requires the __compare() and __mergesort() functions.
## valid values for "how" are:
##   "std asc"
##     use awk's standard rules for comparison, ascending. this is the default
##   "std desc"
##     use awk's standard rules for comparison, descending.
##   "str asc"
##     force comparison as strings, ascending.
##   "str desc"
##     force comparison as strings, descending.
##   "num asc"
##     force a numeric comparison, ascending.
##   "num desc"
##     force a numeric comparison, descending.
function msort(array, out, how,    count, i) {
  # make sure how is correct
  if (length(how)) {
    if (how !~ /^(st[rd]|num) (a|de)sc$/) {
      return -1;
    }

  # how was not passed, use the default
  } else {
    how = "std asc";
  }
  
  # loop over each index, and generate a new array with the same values and
  # sequential indices
  count = 0;
  for (i in array) {
    out[++count] = array[i];
  }

  # actually sort
  return __mergesort(out, count, how);
}

## usage: imsort(s [, how])
## the bevavior is the same as that of msort(), except that the array "s" is
## sorted in-place. the original indices are destroyed and replaced with
## sequential integers. everything else is described in msort() above.
function imsort(array, how,    tmp, count, i) {
  # make sure how is correct
  if (length(how)) {
    if (how !~ /^(st[rd]|num) (a|de)sc$/) {
      return -1;
    }

  # how was not passed, use the default
  } else {
    how = "std asc";
  }
  
  # loop over each index, and generate a new array with the same values and
  # sequential indices
  count = 0;
  for (i in array) {
    tmp[++count] = array[i];
    delete array[i];
  }

  # copy tmp back over array
  for (i=1; i<=count; i++) {
    array[i] = tmp[i];
    delete tmp[i];
  }

  # actually sort
  return __mergesort(array, count, how);
}

## usage: msorti(s, d [, how])
## the behavior is the same as that of msort(), except that the array indices
## are used for sorting, not the array values. when done, the new array is
## indexed numerically, and the values are those of the original indices.
## everything else is described in msort() above.
function msorti(array, out, how,    count, i) {
  # make sure how is correct
  if (length(how)) {
    if (how !~ /^(st[rd]|num) (a|de)sc$/) {
      return -1;
    }

  # how was not passed, use the default
  } else {
    how = "std asc";
  }

  # loop over each index, and generate a new array with the original indices
  # mapped to new numeric ones
  count = 0;
  for (i in array) {
    out[++count] = i;
  }

  # actually sort
  return __mergesort(out, count, how);
}

## usage: imsorti(s [, how])
## the bevavior is the same as that of msorti(), except that the array "s" is
## sorted in-place. the original indices are destroyed and replaced with
## sequential integers. everything else is described in msort() and msorti()
## above.
function imsorti(array, how,    tmp, count, i) {
  # make sure how is correct
  if (length(how)) {
    if (how !~ /^(st[rd]|num) (a|de)sc$/) {
      return -1;
    }

  # how was not passed, use the default
  } else {
    how = "std asc";
  }

  # loop over each index, and generate a new array with the original indices
  # mapped to new numeric ones
  count = 0;
  for (i in array) {
    tmp[++count] = i;
    delete array[i];
  }

  # copy tmp back over the original array
  for (i=1; i<=count; i++) {
    array[i] = tmp[i];
    delete tmp[i];
  }

  # actually sort
  return __mergesort(array, count, how);
}

## usage: msortv(s, d [, how])
## sorts the indices in the array "s" based on the values, creating a new
## sorted array "d" indexed with sequential integers starting with 1, and the
## values the indices of "s". returns the length, or -1 if an error occurs.
## leaves the source array "s" unchanged. the optional string "how" controls
## the direction and the comparison mode. uses the merge sort algorithm, with
## an insertion sort when the list size gets small enough. this is not a stable
## sort. requires the __compare() and __mergesortv() functions. valid values for
## "how" are explained in the msort() function above.
function msortv(array, out, how,    values, count, i) {
  # make sure how is correct
  if (length(how)) {
    if (how !~ /^(st[rd]|num) (a|de)sc$/) {
      return -1;
    }

  # how was not passed, use the default
  } else {
    how = "std asc";
  }

  # loop over each index, and generate two new arrays: the original indices
  # mapped to numeric ones, and the values mapped to the same indices
  count = 0;
  for (i in array) {
    count++;
    out[count] = i;
    values[count] = array[i];
  }

  # actually sort
  return __mergesortv(out, values, count, how);
}



# You can do whatever you want with this stuff, but a thanks is always
# appreciated
#!/usr/bin/awk -f

## usage: getopts(optstring [, longopt_array ])
## Parses options, and deletes them from ARGV. "optstring" is of the form
## "ab:c". Each letter is a possible option. If the letter is followed by a
## colon (:), then the option requires an argument. If an argument is not
## provided, or an invalid option is given, getopts will print the appropriate
## error message and return "?". Returns each option as it's read, and -1 when
## no options are left. "optind" will be set to the index of the next
## non-option argument when finished. "optarg" will be set to the option's
## argument, when provided. If not provided, "optarg" will be empty. "optname"
## will be set to the current option, as provided. Getopts will delete each
## option and argument that it successfully reads, so awk will be able to treat
## whatever's left as filenames/assignments, as usual. If provided,
## "longopt_array" is the name of an associative array that maps long options to
## the appropriate short option (do not include the hyphens on either).
## Sample usage can be found in the examples dir, with gawk extensions, or in
## the ogrep script for a POSIX example: https://github.com/e36freak/ogrep
function getopts(optstring, longarr,    opt, trimmed, hasarg, repeat) {
  hasarg = repeat = 0;
  optarg = "";
  # increment optind
  optind++;

  # return -1 if the current arg is not an option or there are no args left
  if (ARGV[optind] !~ /^-/ || optind >= ARGC) {
    return -1;
  }

  # if option is "--" (end of options), delete arg and return -1
  if (ARGV[optind] == "--") {
    for (i=1; i<=optind; i++) {
      delete ARGV[i];
    }
    return -1;
  }

  # if the option is a long argument...
  if (ARGV[optind] ~ /^--/) {
    # trim hyphens
    trimmed = substr(ARGV[optind], 3);
    # if of the format --foo=bar, split the two. assign "bar" to optarg and
    # set hasarg to 1
    if (trimmed ~ /=/) {
      optarg = trimmed;
      sub(/=.*/, "", trimmed); sub(/^[^=]*=/, "", optarg);
      hasarg = 1;
    }
    
    # invalid long opt
    if (!(trimmed in longarr)) {
      printf("unrecognized option -- '%s'\n", ARGV[optind]) > "/dev/stderr";
      return "?";
    }

    opt = longarr[trimmed];
    # set optname by prepending dashes to the trimmed argument
    optname = "--" trimmed;

  # otherwise, it is a short option
  } else {
    # remove the hyphen, and get just the option letter
    opt = substr(ARGV[optind], 2, 1);
    # set trimmed to whatevers left
    trimmed = substr(ARGV[optind], 3);

    # invalid option
    if (!index(optstring, opt)) {
      printf("invalid option -- '%s'\n", opt) > "/dev/stderr";
      return "?";
    }

    # if there is more to the argument than just -o
    if (length(trimmed)) {
      # if option requires an argument, set the rest to optarg and hasarg to 1
      if (index(optstring, opt ":")) {
        optarg = trimmed;
        hasarg = 1;

      # otherwise, prepend a hyphen to the rest and set repeat to 1, so the
      # same arg is processed again without the first option
      } else {
        ARGV[optind] = "-" trimmed;
        repeat = 1;
      }
    }

    # set optname by prepending a hypen to opt
    optname = "-" opt;
  }

  # if the option requires an arg and hasarg is 0
  if (index(optstring, opt ":") && !hasarg) {
    # increment optind, check if no arguments are left
    if (++optind >= ARGC) {
      printf("option requires an argument -- '%s'\n", optname) > "/dev/stderr";
      return "?";
    }

    # set optarg
    optarg = ARGV[optind];

  # if repeat is set, decrement optind so we process the same arg again
  # mutually exclusive to needing an argument, otherwise hasarg would be set
  } else if (repeat) {
    optind--;
  }

  # delete all arguments up to this point, just to make sure
  for (i=1; i<=optind; i++) {
    delete ARGV[i];
  }

  # return the option letter
  return opt;
}



# You can do whatever you want with this stuff, but a thanks is always
# appreciated
#!/usr/bin/awk -f

# comparison function for the *psort* functions
# usage: __pcompare(a, b, patterns, max, how)
# compares "a" and "b" based on "patterns" and "how", returning 0 for false and
# 1 for true. "patterns" is an indexed array of regexes, from 1 through "max".
# each regex takes priority over subsequent regexes, followed by non-matching
# values. required for all of the psort() functions below
function __pcompare(a, b, pattens, plen, how,    p) {
  # loop over each regex in order, and check if either value matches
  for (p=1; p<=plen; p++) {
    # if the first matches...
    if (a ~ p) {
      # check if the second also matches. if so, do a normal comparison
      if (b ~ p) {
        # standard comparisons
        if (how == "std asc") {
          return a < b;
        } else if (how == "std desc") {
          return a > b;

        # force string comps
        } else if (how == "str asc") {
          return "a" a < "a" b;
        } else if (how == "str desc") {
          return "a" a > "a" b;

        # force numeric
        } else if (how == "num asc") {
          return +a < +b;
        } else if (how == "num desc") {
          return +a > +b;
        }

      # if the second doesn't match, the first sorts higher
      } else {
        return 1;
      }

    # if the second matches but the first didn't, the second sorts higher
    } else if (b ~ p) {
      return 0;
    }
  }

  # no patterns matched, do a normal comparison
  return __compare(a, b, how);
}


# actual sorting function for the *psort* functions
# sorts the values in "array" in-place, from indices "left" to "right", based
# on "how" and the array "patterns" (see the psort() description)
# required for all of the psort() functions below
function __pquicksort(array, left, right, patterns, plen, how,
                      piv, mid, tmp) {
  # return if array contains one element or less
  if ((right - left) <= 0) {
    return;
  }

  # choose random pivot
  piv = int(rand() * (right - left + 1)) + left;

  # swap left and pivot
  tmp = array[piv];
  array[piv] = array[left];
  array[left] = tmp;
  
  mid = left;
  # iterate over each element from the second to the last, and compare
  for (piv=left+1; piv<=right; piv++) {
    # if the comparison based on "how" is true...
    if (__pcompare(array[piv], array[left], patterns, plen, how)) {
      # increment mid
      mid++;

      # swap mid and pivot
      tmp = array[piv];
      array[piv] = array[mid];
      array[mid] = tmp;
    }
  }

  # swap left and mid
  tmp = array[mid];
  array[mid] = array[left];
  array[left] = tmp;
  
  # recursively sort the two halves
  __pquicksort(array, left, mid - 1, patterns, plen, how);
  __pquicksort(array, mid + 1, right, patterns, plen, how);
}


## usage: psort(s, d, patts, max [, how])
## sorts the values of the array "s", based on the rules below. creates a new
## sorted array "d" indexed with sequential integers starting with 1. "patts"
## is a compact (*non-sparse) 1-indexed array containing regular expressions.
## "max" is the length of the "patts" array. returns the length of the "d"
## array. valid values for "how" are explained below. uses the quicksort
## algorithm, with a random pivot to avoid worst-case behavior on already sorted
## arrays. requires the __pcompare() and __pquicksort() functions.
##
##  Sorting rules:
##  - When sorting, values matching an expression in the "patts" array will
##    take priority over any other values
##  - Each expression in the "patts" array will have priority in ascending
##    order by index. "patts[1]" will have priority over "patts[2]" and
##    "patts[3]", etc
##  - Values both matching the same regex will be compared as usual
##  - All non-matching values will be compared as usual
##
## valid values for "how" are:
##   "std asc"
##     use awk's standard rules for comparison, ascending. this is the default
##   "std desc"
##     use awk's standard rules for comparison, descending.
##   "str asc"
##     force comparison as strings, ascending.
##   "str desc"
##     force comparison as strings, descending.
##   "num asc"
##     force a numeric comparison, ascending.
##   "num desc"
##     force a numeric comparison, descending.
function psort(array, out, patterns, plen, how,    count, i) {
  # make sure how is correct
  if (length(how)) {
    if (how !~ /^(st[rd]|num) (a|de)sc$/) {
      return -1;
    }

  # how was not passed, use the default
  } else {
    how = "std asc";
  }
  
  # loop over each index, and generate a new array with the same values and
  # sequential indices
  count = 0;
  for (i in array) {
    out[++count] = array[i];
  }

  # seed the random number generator
  srand();

  # actually sort
  __pquicksort(out, 1, count, patterns, plen, how);

  # return the length
  return count;
}

## usage: ipsort(s, patts, max [, how])
## the bevavior is the same as that of psort(), except that the array "s" is
## sorted in-place. the original indices are destroyed and replaced with
## sequential integers. everything else is described in psort() above.
function ipsort(array, patterns, plen, how,    tmp, count, i) {
  # make sure how is correct
  if (length(how)) {
    if (how !~ /^(st[rd]|num) (a|de)sc$/) {
      return -1;
    }

  # how was not passed, use the default
  } else {
    how = "std asc";
  }
  
  # loop over each index, and generate a new array with the same values and
  # sequential indices
  count = 0;
  for (i in array) {
    tmp[++count] = array[i];
    delete array[i];
  }

  # copy tmp back over array
  for (i=1; i<=count; i++) {
    array[i] = tmp[i];
    delete tmp[i];
  }

  # seed the random number generator
  srand();

  # actually sort
  __pquicksort(array, 1, count, patterns, plen, how);

  # return the length
  return count;
}

## usage: psorti(s, d, patts, max [, how])
## the behavior is the same as that of psort(), except that the array indices
## are used for sorting, not the array values. when done, the new array is
## indexed numerically, and the values are those of the original indices.
## everything else is described in psort() above.
function psorti(array, out, patterns, plen, how,    count, i) {
  # make sure how is correct
  if (length(how)) {
    if (how !~ /^(st[rd]|num) (a|de)sc$/) {
      return -1;
    }

  # how was not passed, use the default
  } else {
    how = "std asc";
  }

  # loop over each index, and generate a new array with the original indices
  # mapped to new numeric ones
  count = 0;
  for (i in array) {
    out[++count] = i;
  }

  # seed the random number generator
  srand();

  # actually sort
  __pquicksort(out, 1, count, patterns, plen, how);

  # return the length
  return count;
}

## usage: ipsorti(s, patts, max [, how])
## the bevavior is the same as that of psorti(), except that the array "s" is
## sorted in-place. the original indices are destroyed and replaced with
## sequential integers. everything else is described in psort() and psorti()
## above.
function ipsorti(array, patterns, plen, how,    tmp, count, i) {
  # make sure how is correct
  if (length(how)) {
    if (how !~ /^(st[rd]|num) (a|de)sc$/) {
      return -1;
    }

  # how was not passed, use the default
  } else {
    how = "std asc";
  }

  # loop over each index, and generate a new array with the original indices
  # mapped to new numeric ones
  count = 0;
  for (i in array) {
    tmp[++count] = i;
    delete array[i];
  }

  # copy tmp back over the original array
  for (i=1; i<=count; i++) {
    array[i] = tmp[i];
    delete tmp[i];
  }

  # seed the random number generator
  srand();

  # actually sort
  __pquicksort(array, 1, count, patterns, plen, how);

  # return the length
  return count;
}


# You can do whatever you want with this stuff, but a thanks is always
# appreciated
#!/usr/bin/awk -f

# comparison function
# usage: __compare(a, b, how)
# compares "a" and "b" based on "how", returning 0 for false and 1 for true.
# required for all of the qsort() functions below
function __compare(a, b, how) {
  # standard comparisons
  if (how == "std asc") {
    return a < b;
  } else if (how == "std desc") {
    return a > b;

  # force string comps
  } else if (how == "str asc") {
    return "a" a < "a" b;
  } else if (how == "str desc") {
    return "a" a > "a" b;

  # force numeric
  } else if (how == "num asc") {
    return +a < +b;
  } else if (how == "num desc") {
    return +a > +b;
  }
}

# actual sorting function
# sorts the values in "array" in-place, from indices "left" to "right", based
# on the comparison mode "how" (see the qsort() description).
# required for all of the qsort() functions below
function __quicksort(array, left, right, how,    piv, mid, tmp) {
  # return if array contains one element or less
  if ((right - left) <= 0) {
    return;
  }

  # choose random pivot
  piv = int(rand() * (right - left + 1)) + left;

  # swap left and pivot
  tmp = array[piv];
  array[piv] = array[left];
  array[left] = tmp;
  
  mid = left;
  # iterate over each element from the second to the last, and compare
  for (piv=left+1; piv<=right; piv++) {
    # if the comparison based on "how" is true...
    if (__compare(array[piv], array[left], how)) {
      # increment mid
      mid++;

      # swap mid and pivot
      tmp = array[piv];
      array[piv] = array[mid];
      array[mid] = tmp;
    }
  }

  # swap left and mid
  tmp = array[mid];
  array[mid] = array[left];
  array[left] = tmp;
  
  # recursively sort the two halves
  __quicksort(array, left, mid - 1, how);
  __quicksort(array, mid + 1, right, how);
}

# actual sorting function for the qsortv() function
# sorts the indices in "array" on the original values in "values", from indices
# "left" to "right", based on the comparison mode "how" (see the qsortv()
# description)
# required for the qsortv() function below
function __vquicksort(array, values, left, right, how,    piv, mid, tmp) {
  # return if array contains one element or less
  if ((right - left) <= 0) {
    return;
  }

  # choose random pivot
  piv = int(rand() * (right - left + 1)) + left;

  # swap left and pivot
  tmp = array[piv];
  array[piv] = array[left];
  array[left] = tmp;
  tmp = values[piv];
  values[piv] = values[left];
  values[left] = tmp;
  
  mid = left;
  # iterate over each element from the second to the last, and compare
  for (piv=left+1; piv<=right; piv++) {
    # if the comparison based on "how" is true...
    if (__compare(values[piv], values[left], how)) {
      # increment mid
      mid++;

      # swap mid and pivot
      tmp = array[piv];
      array[piv] = array[mid];
      array[mid] = tmp;
      tmp = values[piv];
      values[piv] = values[mid];
      values[mid] = tmp;
    }
  }

  # swap left and mid
  tmp = array[mid];
  array[mid] = array[left];
  array[left] = tmp;
  tmp = values[mid];
  values[mid] = values[left];
  values[left] = tmp;
  
  # recursively sort the two halves
  __vquicksort(array, values, left, mid - 1, how);
  __vquicksort(array, values, mid + 1, right, how);
}



## usage: qsort(s, d [, how])
## sorts the elements in the array "s" using awk's normal rules for comparing
## values, creating a new sorted array "d" indexed with sequential integers
## starting with 1. returns the length, or -1 if an error occurs.. leaves the
## indices of the source array "s" unchanged. the optional string "how" controls
## the direction and the comparison mode. uses the quick sort algorithm, with a
## random pivot to avoid worst-case behavior on already sorted arrays. this is
## not a stable sort. requires the __compare() and __quicksort() functions.
## valid values for "how" are:
##   "std asc"
##     use awk's standard rules for comparison, ascending. this is the default
##   "std desc"
##     use awk's standard rules for comparison, descending.
##   "str asc"
##     force comparison as strings, ascending.
##   "str desc"
##     force comparison as strings, descending.
##   "num asc"
##     force a numeric comparison, ascending.
##   "num desc"
##     force a numeric comparison, descending.
function qsort(array, out, how,    count, i) {
  # make sure how is correct
  if (length(how)) {
    if (how !~ /^(st[rd]|num) (a|de)sc$/) {
      return -1;
    }

  # how was not passed, use the default
  } else {
    how = "std asc";
  }
  
  # loop over each index, and generate a new array with the same values and
  # sequential indices
  count = 0;
  for (i in array) {
    out[++count] = array[i];
  }

  # seed the random number generator
  srand();

  # actually sort
  __quicksort(out, 1, count, how);

  # return the length
  return count;
}

## usage: iqsort(s [, how])
## the bevavior is the same as that of qsort(), except that the array "s" is
## sorted in-place. the original indices are destroyed and replaced with
## sequential integers. everything else is described in qsort() above.
function iqsort(array, how,    tmp, count, i) {
  # make sure how is correct
  if (length(how)) {
    if (how !~ /^(st[rd]|num) (a|de)sc$/) {
      return -1;
    }

  # how was not passed, use the default
  } else {
    how = "std asc";
  }
  
  # loop over each index, and generate a new array with the same values and
  # sequential indices
  count = 0;
  for (i in array) {
    tmp[++count] = array[i];
    delete array[i];
  }

  # copy tmp back over array
  for (i=1; i<=count; i++) {
    array[i] = tmp[i];
    delete tmp[i];
  }

  # seed the random number generator
  srand();

  # actually sort
  __quicksort(array, 1, count, how);

  # return the length
  return count;
}

## usage: qsorti(s, d [, how])
## the behavior is the same as that of qsort(), except that the array indices
## are used for sorting, not the array values. when done, the new array is
## indexed numerically, and the values are those of the original indices.
## everything else is described in qsort() above.
function qsorti(array, out, how,    count, i) {
  # make sure how is correct
  if (length(how)) {
    if (how !~ /^(st[rd]|num) (a|de)sc$/) {
      return -1;
    }

  # how was not passed, use the default
  } else {
    how = "std asc";
  }

  # loop over each index, and generate a new array with the original indices
  # mapped to new numeric ones
  count = 0;
  for (i in array) {
    out[++count] = i;
  }

  # seed the random number generator
  srand();

  # actually sort
  __quicksort(out, 1, count, how);

  # return the length
  return count;
}

## usage: iqsorti(s [, how])
## the bevavior is the same as that of qsorti(), except that the array "s" is
## sorted in-place. the original indices are destroyed and replaced with
## sequential integers. everything else is described in qsort() and qsorti()
## above.
function iqsorti(array, how,    tmp, count, i) {
  # make sure how is correct
  if (length(how)) {
    if (how !~ /^(st[rd]|num) (a|de)sc$/) {
      return -1;
    }

  # how was not passed, use the default
  } else {
    how = "std asc";
  }

  # loop over each index, and generate a new array with the original indices
  # mapped to new numeric ones
  count = 0;
  for (i in array) {
    tmp[++count] = i;
    delete array[i];
  }

  # copy tmp back over the original array
  for (i=1; i<=count; i++) {
    array[i] = tmp[i];
    delete tmp[i];
  }

  # seed the random number generator
  srand();

  # actually sort
  __quicksort(array, 1, count, how);

  # return the length
  return count;
}

## usage: qsortv(s, d [, how])
## sorts the indices in the array "s" based on the values, creating a new
## sorted array "d" indexed with sequential integers starting with 1, and the
## values the indices of "s". returns the length, or -1 if an error occurs.
## leaves the source array "s" unchanged. the optional string "how" controls
## the direction and the comparison mode. uses the quicksort algorithm, with a
## random pivot to avoid worst-case behavior on already sorted arrays. this is
## not a stable sort. requires the __compare() and __vquicksort() functions.
## valid values for "how" are explained in the qsort() function above.
function qsortv(array, out, how,    values, count, i) {
  # make sure how is correct
  if (length(how)) {
    if (how !~ /^(st[rd]|num) (a|de)sc$/) {
      return -1;
    }

  # how was not passed, use the default
  } else {
    how = "std asc";
  }

  # loop over each index, and generate two new arrays: the original indices
  # mapped to numeric ones, and the values mapped to the same indices
  count = 0;
  for (i in array) {
    count++;
    out[count] = i;
    values[count] = array[i];
  }

  # seed the random number generator
  srand();

  # actually sort
  __vquicksort(out, values, 1, count, how);

  # return the length
  return count;
}



# You can do whatever you want with this stuff, but a thanks is always
# appreciated
#!/usr/bin/awk -f

# actual shuffle function
# shuffles the values in "array" in-place, from indices "left" to "right".
# required for all of the shuf() functions below
function __shuffle(array, left, right,    r, i, tmp) {
  # loop backwards over the elements
  for (i=right; i>left; i--) {
    # generate a random number between the start and current element
    r = int(rand() * (i - left + 1)) + left;

    # swap current element and randomly generated one
    tmp = array[i];
    array[i] = array[r];
    array[r] = tmp;
  }
}



## usage: shuf(s, d)
## shuffles the array "s", creating a new shuffled array "d" indexed with
## sequential integers starting with one. returns the length, or -1 if an error
## occurs. leaves the indices of the source array "s" unchanged. uses the knuth-
## fisher-yates algorithm. requires the __shuffle() function.
function shuf(array, out,    count, i) {
  # loop over each index, and generate a new array with the same values and
  # sequential indices
  count = 0;
  for (i in array) {
    out[++count] = array[i];
  }

  # seed the random number generator
  srand();

  # actually shuffle
  __shuffle(out, 1, count);

  # return the length
  return count;
}

## usage: ishuf(s)
## the behavior is the same as that of shuf(), except the array "s" is sorted
## in-place. the original indices are destroyed and replaced with sequential
## integers. everything else is described in shuf() above.
function ishuf(array,    tmp, count, i) {
  # loop over each index, and generate a new array with the same values and
  # sequential indices
  count = 0;
  for (i in array) {
    tmp[++count] = array[i];
    delete array[i];
  }

  # copy tmp back over array
  for (i=1; i<=count; i++) {
    array[i] = tmp[i];
    delete tmp[i];
  }

  # seed the random number generator
  srand();

  # actually shuffle
  __shuffle(array, 1, count);

  # return the length
  return count;
}

## usage: shufi(s, d)
## the bevavior is the same as that of shuf(), except that the array indices
## are shuffled, not the array values. when done, the new array is indexed
## numerically, and the values are those of the original indices. everything
## else is described in shuf() above.
function shufi(array, out,    count, i) {
  # loop over each index, and generate a new array with the original indices
  # mapped to new numeric ones
  count = 0;
  for (i in array) {
    out[++count] = i;
  }

  # seed the random number generator
  srand();

  # actually shuffle
  __shuffle(out, 1, count);

  # return the length
  return count;
}

## usage: ishufi(s)
## the behavior is tha same as that of shufi(), except that the array "s" is
## sorted in-place. the original indices are destroyed and replaced with
## sequential integers. everything else is describmed in shuf() and shufi()
## above.
function ishufi(array,    tmp, count, i) {
  # loop over each index, and generate a new array with the original indices
  # mapped to new numeric ones
  count = 0;
  for (i in array) {
    tmp[++count] = i;
    delete array[i];
  }

  # copy tmp back over the original array
  for (i=1; i<=count; i++) {
    array[i] = tmp[i];
    delete tmp[i];
  }

  # seed the random number generator
  srand();

  # actually shuffle
  __shuffle(array, 1, count);

  # return the length
  return count;
}



# You can do whatever you want with this stuff, but a thanks is always
# appreciated
#!/usr/bin/awk -f

# comparison function
# compares "A" and "b" based on "how", returning 0 for false and 1 for true
# required for all max() and min() functions below
function __mcompare(a, b, how) {
 # standard comparison
  if (how == "std") {
    return a > b;

  # force string comp
  } else if (how == "str") {
    return "a" a > "a" b;

  # force numeric
  } else if (how == "num") {
    return +a > +b;
  }
}



## usage: center(string [, width])
## returns "string" centered based on "width". if "width" is not provided (or 
## is 0), uses the width of the terminal, or 80 if standard output is not open
## on a terminal.
## note: does not check the length of the string. if it's wider than the
## terminal, it will not center lines other than the first. for best results,
## combine with fold() (see the cfold script in the examples directory for a
## script that does exactly this)
function center(str, cols,    off, cmd) {
  if (!cols) {
    # checks if stdout is a tty
    if (system("test -t 1")) {
      cols = 80;
    } else {
      cmd = "tput cols";
      cmd | getline cols;
      close(cmd);
    }
  }

  off = int((cols/2) + (length(str)/2));

  return sprintf("%*s", off, str);
}

## usage: delete_arr(array)
## deletes every element in "array"
function delete_arr(arr) {
  split("", arr);
}

## usage: fold(string, sep [, width])
## returns "string", wrapped, with lines broken on "sep" to "width" columns.
## "sep" is a list of characters to break at, similar to IFS in a POSIX shell.
## if "sep" is empty, wraps at exactly "width" characters. if "width" is not
## provided (or is 0), uses the width of the terminal, or 80 if standard output
## is not open on a terminal.
## note: currently, tabs are squeezed to a single space. this will be fixed
function fold(str, sep, cols,    out, cmd, i, len, chars, c, last, f, first) {
  if (!cols) {
    # checks if stdout is a tty
    if (system("test -t 1")) {
      cols = 80;
    } else {
      cmd = "tput cols";
      cmd | getline cols;
      close(cmd);
    }
  }

  # squeeze tabs and newlines to spaces
  gsub(/[\t\n]/, " ", str);

  # if "sep" is empty, just fold on cols with substr
  if (!length(sep)) {
    len = length(str);

    out = substr(str, 1, cols);
    for (i=cols+1; i<=len; i+=cols) {
      out = out "\n" substr(str, i, cols);
    }

    return out;

  # otherwise, we have to loop over every character (can't split() on sep, it
  # would destroy the existing separators)
  } else {
    # split string into char array
    len = split(str, chars, "");
    # set boolean, used to assign the first line differently
    first = 1;

    for (i=1; i<=len; i+=last) {
      f = 0;
      for (c=i+cols-1; c>=i; c--) {
        if (index(sep, chars[c])) {
          last = c - i + 1;
          f = 1;
          break;
        }
      }

      if (!f) {
        last = cols;
      }

      if (first) {
        out = substr(str, i, last);
        first = 0;
      } else {
        out = out "\n" substr(str, i, last);
      }
    }
  }

  # return the output
  return out;
}

## usage: ssub(ere, repl [, in])
## behaves like sub, except returns the result and doesn't modify the original
function ssub(ere, repl, str) {
  # if "in" is not provided, use $0
  if (!length(str)) {
    str = $0;
  }

  # substitute
  sub(ere, repl, str);
  return str;
}

## usage: sgsub(ere, repl [, in])
## behaves like gsub, except returns the result and doesn't modify the original
function sgsub(ere, repl, str) {
  # if "in" is not provided, use $0
  if (!length(str)) {
    str = $0;
  }

  # substitute
  gsub(ere, repl, str);
  return str;
}

## usage: lsub(str, repl [, in])
## substites the string "repl" in place of the first instance of "str" in the
## string "in" and returns the result. does not modify the original string.
## if "in" is not provided, uses $0.
function lsub(str, rep, val,    len, i) {
  # if "in" is not provided, use $0
  if (!length(val)) {
    val = $0;
  }

  # get the length of val, in order to know how much of the string to remove
  if (!(len = length(str))) {
    # if "str" is empty, just prepend "rep" and return
    val = rep val;
    return val;
  }

  # substitute val for rep
  if (i = index(val, str)) {
    val = substr(val, 1, i - 1) rep substr(val, i + len);
  }

  # return the result
  return val;
}

## usage: glsub(str, repl [, in])
## behaves like lsub, except it replaces all occurances of "str"
function glsub(str, rep, val,    out, len, i, a, l) {
  # if "in" is not provided, use $0
  if (!length(val)) {
    val = $0;
  }
  # empty the output string
  out = "";

  # get the length of val, in order to know how much of the string to remove
  if (!(len = length(str))) {
    # if "str" is empty, adds "rep" between every character and returns
    l = split(val, a, "");
    for (i=1; i<=l; i++) {
      out = out rep a[i];
    }

    return out rep;
  }

  # loop while 'val' is in 'str'
  while (i = index(val, str)) {
    # append everything up to the search string, and the replacement, to out
    out = out substr(val, 1, i - 1) rep;
    # remove everything up to and including the first instance of str from val
    val = substr(val, i + len);
  }

  # append whatever is left in val to out and return
  return out val;
}

## usage: shell_esc(string)
## returns the string escaped so that it can be used in a shell command
function shell_esc(str) {
  gsub(/'/, "'\\''", str);

  return "'" str "'";
}

## usage: str_to_arr(string, array)
## converts string to an array, one char per element, 1-indexed
## returns the array length
function str_to_arr(str, arr) {
  return split(str, arr, "");
}

## usage: extract_range(string, start, stop)
## extracts fields "start" through "stop" from "string", based on FS, with the
## original field separators intact. returns the extracted fields.
function extract_range(str, start, stop,    i, re, out) {
  # if FS is the default, trim leading and trailing spaces from "string" and
  # set "re" to the appropriate regex
  if (FS == " ") {
    gsub(/^[[:space:]]+|[[:space:]]+$/, "", str);
    re = "[[:space:]]+";
  } else {
    re = FS;
  }

  # remove fields 1 through start - 1 from the beginning
  for (i=1; i<start; i++) {
    if (match(str, re)) {
      str = substr(str, RSTART + RLENGTH);

    # there's no FS left, therefore the range is empty
    } else {
      return "";
    }
  }

  # add fields start through stop - 1 to the output var
  for (i=start; i<stop; i++) {
    if (match(str, re)) {
      # append the field to the output
      out = out substr(str, 1, RSTART + RLENGTH - 1);

      # remove the field from the line
      str = substr(str, RSTART + RLENGTH);

    # no FS left, just append the rest of the line and return
    } else {
      return out str;
    }
  }

  # append the last field and return
  if (match(str, re)) {
    return out substr(str, 1, RSTART - 1);
  } else {
    return out str;
  }
}

## usage: fwidths(width_spec [, string])
## extracts substrings from "string" according to "width_spec" from left to
## right and assigns them to $1, $2, etc. also assigns the NF variable. if
## "string" is not supplied, uses $0. "width_spec" is a space separated list of
## numbers that specify field widths, just like GNU awk's FIELDWIDTHS variable.
## if there is data left over after the last width_spec, adds it to a final
## field. returns the value for NF.
function fwidths(wspec, str,    fw, i, len) {
  if (!length(str)) {
    str = $0;
  }

  # turn wspec into the array fw
  len = split(wspec, fw, / /);
  
  # loop over each wspec value, while the string is not exhausted
  for (i=1; i <= len && length(str); i++) {
    # assign the field
    $i = substr(str, 1, fw[i]);

    # chop the value off of the original string
    str = substr(str, fw[i] + 1);
  }

  # if there's anything left, add another field
  if (length(str)) {
    $i = str;
  } else {
    i--;
  }

  # set and return NF
  return NF = i;
}

## usage: fwidths_arr(width_spec, array [, string])
## the behavior is the same as that of fwidths(), except that the values are
## assigned to "array", indexed with sequential integers starting with 1.
## returns the length. everything else is described in fwidths() above.
function fwidths_arr(wspec, arr, str,    fw, i, len) {
  if (!length(str)) {
    str = $0;
  }

  # turn wspec into the array fw
  len = split(wspec, fw, / /);

  # loop over each wspec value, while the string is not exhausted
  for (i=1; i <= len && length(str); i++) {
    # assign the array element
    arr[i] = substr(str, 1, fw[i]);

    # chop the value off of the original string
    str = substr(str, fw[i] + 1);
  }

  # if there's anything left, add another field
  if (length(str)) {
    arr[i] = str;
  } else {
    i--;
  }

  # return the array length
  return i;
}

## usage: lsplit(str, arr, sep)
## splits the string "str" into array elements "arr[1]", "arr[2]", .., "arr[n]",
## and returns "n". all elements of "arr" are deleted before the split is
## performed. the separation is done on the literal string "sep".
function lsplit(str, arr, sep,    len, slen, i) {
  # empty "arr"
  split("", arr);

  # if "sep" is empty, just do a normal split
  if (!(slen = length(sep))) {
    return split(str, arr, "");
  }

  # loop while "sep" is matched
  while (i = index(str, sep)) {
    # append field to array
    arr[++len] = substr(str, 1, i - 1);

    # remove that portion (with the sep) from the string
    str = substr(str, i + slen);
  }

  # append last field to "arr"
  arr[++len] = str;

  # return the length
  return len;
}

## usage: ssplit(str, arr, seps [, ere])
## similar to GNU awk 4's "seps" functionality for split(). splits the string
## "str" into the array "arr" and the separators array "seps" on the regular
## expression "ere", and returns the number of fields. the value of "seps[i]"
## is the separator that appeared in front of "arr[i+1]". if "ere" is omitted or
## empty, FS is used instead. if "ere" is a single space, leading whitespace in
## "str" will go into the extra array element "seps[0]" and trailing whitespace
## will go into the extra array element "seps[len]", where "len" is the return
## value.
## note: /regex/ style quoting cannot be used for "ere".
function ssplit(str, arr, seps, ere,    len, totrim) {
  # if "ere" is unset or empty, use FS
  if (!length(ere)) {
    ere = FS;
  }

  # if "ere" is a single space...
  if (ere == " ") {
    # set it to match all spaces
    ere = "[[:space:]]+";

    # trim leading whitespace and assign it to seps[0]
    if (match(str, /[^[:space:]]/)) {
      seps[0] = substr(str, 1, RSTART - 1);
      str = substr(str, RSTART);

    # no non-space characters in the line, just return
    } else {
      return 0;
    }

    # don't put an empty element after the last separator
    totrim = 1;
  }


  # loop while "ere" is matched 
  while (match(str, ere)) {
    # append field and sep to arrays
    len++;
    arr[len] = substr(str, 1, RSTART - 1);
    seps[len] = substr(str, RSTART, RLENGTH);

    # remove matched portion from the string
    str = substr(str, RSTART + RLENGTH);
  }

  # append last field to "arr" if needed
  if (length(str) || !totrim) {
    arr[++len] = str;
  }

  # return the length
  return len;
}

## usage: ends_with(string, substring)
## returns 1 if "strings" ends with "substring", otherwise 0
function ends_with(string, s) {
  return substr(string, length(string) - length(s) + 1) == s;
}

## usage: trim(string)
## returns "string" with leading and trailing whitespace trimmed
function trim(str) {
  gsub(/^[[:blank:]]+|[[:blank:]]+$/, "", str);

  return str;
}

## usage: rev(string)
## returns "string" backwards
function rev(str,    a, len, i, o) {
  # split string into character array
  len = split(str, a, "");

  # iterate backwards and append to the output string
  for (i=len; i>0; i--) {
    o = o a[i];
  }

  return o;
}

## usage: max(array [, how ])
## returns the maximum value in "array", 0 if the array is empty, or -1 if an
## error occurs. the optional string "how" controls the comparison mode.
## requires the __mcompare() function.
## valid values for "how" are:
##   "std"
##     use awk's standard rules for comparison. this is the default
##   "str"
##     force comparison as strings
##   "num"
##     force a numeric comparison
function max(array, how,    m, i, f) {
  # make sure how is correct
  if (length(how)) {
    if (how !~ /^(st[rd]|num)$/) {
      return -1;
    }

  # how was not passed, use the default
  } else {
    how = "std";
  }

  m = 0;
  f = 1;

  # loop over each array value
  for (i in array) {
    # if this is the first iteration, use the value as m
    if (f) {
      m = array[i];
      f = 0;

      continue;
    }

    # otherwise, if it's greater than "m", reassign it
    if (__mcompare(array[i], m, how)) {
      m = array[i];
    }
  }

  return m;
}

## usage: maxi(array [, how ])
## the behavior is the same as that of max(), except that the array indices are
## used, not the array values. everything else is explained in max() above.
function maxi(array, how,    m, i, f) {
  # make sure how is correct
  if (length(how)) {
    if (how !~ /^(st[rd]|num)$/) {
      return -1;
    }

  # how was not passed, use the default
  } else {
    how = "std";
  }

  m = 0;
  f = 1;

  # loop over each index
  for (i in array) {
    # if this is the first iteration, use the value as m
    if (f) {
      m = i;
      f = 0;

      continue;
    }

    # otherwise, if it's greater than "m", reassign it
    if (__mcompare(i, m, how)) {
      m = i;
    }
  }

  return m;
}

## usage: min(array [, how ])
## the behavior is the same as that of max(), except that the minimum value is
## returned instead of the maximum. everything else is explained in max() above.
function min(array, how,    m, i, f) {
  # make sure how is correct
  if (length(how)) {
    if (how !~ /^(st[rd]|num)$/) {
      return -1;
    }

  # how was not passed, use the default
  } else {
    how = "std";
  }

  m = 0;
  f = 1;

  # loop over each index
  for (i in array) {
    # if this is the first iteration, use the value as m
    if (f) {
      m = array[i];
      f = 0;

      continue;
    }

    # otherwise, if it's less than "m", reassign it
    if (__mcompare(m, array[i], how)) {
      m = array[i];
    }
  }

  return m;
}

## usage: mini(array [, how ])
## the behavior is the same as that of min(), except that the array indices are
## used instead of the array values. everything else is explained in min() and
## max() above.
function mini(array, how,    m, i, f) {
  # make sure how is correct
  if (length(how)) {
    if (how !~ /^(st[rd]|num)$/) {
      return -1;
    }

  # how was not passed, use the default
  } else {
    how = "std";
  }

  m = 0;
  f = 1;

  # loop over each index
  for (i in array) {
    # if this is the first iteration, use the value as m
    if (f) {
      m = i;
      f = 0;

      continue;
    }

    # otherwise, if it's less than "m", reassign it
    if (__mcompare(m, i, how)) {
      m = i;
    }
  }

  return m;
}



# You can do whatever you want with this stuff, but a thanks is always
# appreciated
#!/usr/bin/awk -f

## usage: isatty(fd)
## Checks if "fd" is open on a tty. Returns 1 if so, 0 if not, and -1 if an
## error occurs
function isatty(fd) {
  # make sure fd is an int
  if (fd !~ /^[0-9]+$/) {
    return -1;
  }

  # actually test
  return !system("test -t " fd);
}

## usage: mktemp(template [, type])
## creates a temporary file or directory, safely, and returns its name.
## if template is not a pathname, the file will be created in ENVIRON["TMPDIR"]
## if set, otherwise /tmp. the last six characters of template must be "XXXXXX",
## and these are replaced with a string that makes the filename unique. type, if
## supplied, is either "f", "d", or "u": for file, directory, or dry run (just
## returns the name, doesn't create a file), respectively. If template is not
## provided, uses "tmp.XXXXXX". Files are created u+rw, and directories u+rwx,
## minus umask restrictions. returns -1 if an error occurs.
function mktemp(template, type,
                c, chars, len, dir, dir_esc, rstring, i, out, out_esc, umask,
                cmd) {
  # portable filename characters
  c = "012345689ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz";
  len = split(c, chars, "");

  # make sure template is valid
  if (length(template)) {
    if (template !~ /XXXXXX$/) {
      return -1;
    }

  # template was not supplied, use the default
  } else {
    template = "tmp.XXXXXX";
  }

  # make sure type is valid
  if (length(type)) {
    if (type !~ /^[fdu]$/) {
      return -1;
    }

  # type was not supplied, use the default
  } else {
    type = "f";
  }

  # if template is a path...
  if (template ~ /\//) {
    dir = template;
    sub(/\/[^/]*$/, "", dir);
    sub(/.*\//, "", template);

  # template is not a path, determine base dir
  } else {
    if (length(ENVIRON["TMPDIR"])) {
      dir = ENVIRON["TMPDIR"];
    } else {
      dir = "/tmp";
    }
  }

  # escape dir for shell commands
  esc_dir = dir;
  sub(/'/, "'\\''", esc_dir);
  esc_dir = "'" esc_dir "'";

  # if this is not a dry run, make sure the dir exists
  if (type != "u" && system("test -d " esc_dir)) {
    return -1;
  }

  # get the base of the template, sans Xs
  template = substr(template, 0, length(template) - 6);
  
  # generate the filename
  do {
    rstring = "";
    for (i=0; i<6; i++) {
      c = chars[int(rand() * len) + 1];
      rstring = rstring c;
    }
    
    out_esc = out = dir "/" template rstring;
    sub(/'/, "'\\''", out_esc);
    out_esc = "'" out_esc "'";
  } while (!system("test -e " out_esc));

  # if needed, create the filename
  if (type == "f") {
    system("touch " out_esc);
    cmd = "umask";
    cmd | getline umask;
    close(cmd);
    umask = substr(umask, 2, 1);
    system("chmod 0" 6 - umask "00 " out_esc);
  } else if (type == "d") {
    system("mkdir " out_esc);
    cmd = "umask";
    cmd | getline umask;
    close(cmd);
    umask = substr(umask, 2, 1);
    system("chmod 0" 7 - umask "00 " out_esc);
  }

  # return the filename
  return out;
}



# You can do whatever you want with this stuff, but a thanks is always
# appreciated
#!/usr/bin/awk -f

## usage: month_to_num(month)
## converts human readable month to the decimal representation
## returns the number, -1 if the month doesn't exist
function month_to_num(mon,    months, m) {
  # populate months[] array
  months["january"] =  1; months["february"] =  2; months["march"]     =  3;
  months["april"]   =  4; months["may"]      =  5; months["june"]      =  6;
  months["july"]    =  7; months["august"]   =  8; months["september"] =  9;
  months["october"] = 10; months["november"] = 11; months["december"]  = 12;

  # also populate abbreviations
  for (m in months) {
    months[substr(m, 1, 3)] = months[m];
  }

  # convert month to lowercase
  mon = tolower(mon);

  # check if month exists
  if (mon in months) {
    return months[mon];
  } else {
    return -1;
  }
}

## usage: day_to_num(day)
## converts human readable day to the decimal representation
## returns the number, -1 if the day doesn't exist
## like date +%w, sunday is 0
function day_to_num(day,    days, d) {
  # populate days[] array
    days["sunday"]    = 0; days["monday"]   = 1; days["tuesday"] = 2;
    days["wednesday"] = 3; days["thursday"] = 4; days["friday"]  = 5;
    days["saturday"]  = 6;

  # also populate abbreviations
    days["sun"]   = 0; days["mon"] = 1; days["tues"] = 2; days["wed"] = 3;
    days["thurs"] = 4; days["fri"] = 5; days["sat"]  = 6;

  # convert day to lowercase
    day = tolower(day);

  # check if day exists
  if (day in days) {
    return days[day];
  } else {
    return -1;
  }
}

## usage: hr_to_sec(timestamp)
## converts HH:MM:SS or MM:SS to seconds
## returns -1 if invalid format
function hr_to_sec(time,    t, l, i, j) {
  # check for valid format
  if (time !~ /^[0-9]+(:[0-9][0-9])?:[0-9][0-9]$/) {
    return -1;
  }

  # convert
  l = split(time, t, /:/);
  
  j = time = 0;
  for (i=l; i>0; i--) {
    time += t[i] * (60 ^ j++);
  }

  return time;
}

## usage: sec_to_hr(seconds)
## converts seconds to HH:MM:SS
function sec_to_hr(sec,    m, s) {
  s = sec % 60;
  sec = int(sec / 60);
  m = sec % 60;
  sec = int(sec / 60);

  return sprintf("%02d:%02d:%02d", sec, m, s);
}

## usage: ms_to_hr(milliseconds)
## converts milliseconds to a "time(1)"-similar human readable format, such
## as 1m4.356s
function ms_to_hr(ms,    m, s, ns) {
  ms = ms / 1000;
  s = int(ms);
  m = int(s / 60);
  ns = s % 60;

  return sprintf("%dm%0.3fs", m, ns + (ms - s));
}

## usage: add_day_suff(day_of_month)
## prepends the appropriate suffix to "day_of_month". for example,
## add_day_suff(1) will return "1st", and add_day_suff(22) will return "22nd"
## returns -1 if "day_of_month" is not a positive integer
function add_day_suff(day) {
  # make sure day is a positive int
  if (day !~ /^[0-9]+$/ || day <= 0) {
    return -1;
  }

  # append prefix
  if ((day > 3 && day < 21) || day ~ /[04-9]$/) {
    return day "th";
  } else if (day ~ /1$/) {
    return day "st";
  } else if (day ~ /2$/) {
    return day "nd";
  } else {
    return day "rd";
  }
}




# You can do whatever you want with this stuff, but a thanks is always
# appreciated
