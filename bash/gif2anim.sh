#!/bin/sh -
#
# gif2anim [options] image_anim.gif...
#
# This script uses the ImageMagick 'convert' program to take a gif animation
# (normally ending in '_anim.gif' and NOT just '.gif') and converts it into a
# collection of X Pixmap files (or GIFs) as well as a "animation sequence
# file" which lists all settings needed to re-build the original GIF
# animation.
#
# The animation file uses the string "BASENAME" to represent the basenames for
# the individual frame images that it loads when re-creating the original GIF
# animation.  This allows you to specify a different set of images without
# changling any other settings.
#
# The resulting animation file can be used to re-build the animation either
# using the IM options contained (ignoring comments), or by giving the file to
# the "anim2gif" script which will also do the basename filename substitutions
# specified.
#
# OPTIONS
#      -c             coalesce animation before parsing
#      -t             Add time synchronization comment before each frame
#      -l             just list the anim file to stdout, no images
#      -v             Be verbose in animation conversions
#      -n             no images, just create the '.anim' file
#      -g             use an GIF suffix for frame images (default)
#      -x             use an XPM suffix for frame images
#      -s suffix      use this suffix for the frame images
#      -i initframe   number of the first frame image (def=1)
#      -b framename   basename for the individual frames
#      -o file.anim   output to this .anim file (- for stdout)
#
# NOTE: This script reqires both the ImageMagick 'convert' program and
# 'identify' programs to extract the images and information needed.
#
####
#
# WARNING: Input arguments are NOT tested for correctness.
# This script represents a security risk if used ONLINE.
# I accept no responsiblity for misuse. Use at own risk.
#
# Anthony Thyssen    25 April 1996    <anthony@cit.griffith.edu.au>
#
# 1 November 2005  Expanded with options for basename handling.
#
PROGNAME=$(type "$0" | awk '{print $3}') # search for executable on path
PROGDIR=$(dirname "$PROGNAME")           # extract directory of program
PROGNAME=$(basename "$PROGNAME")         # base name of program
Usage() {                                # output the script comments as docs
  echo >&2 "$PROGNAME:" "$@"
  sed >&2 -n '/^###/q; /^#/!q; s/^#//; s/^ //; 3s/^/Usage: /; 2,$ p' \
    "$PROGDIR/$PROGNAME"
  exit 10
}

sfx="gif"   # Default output suffix of images to produce
init=1      # the count for the first frame in sequence
fmt="%03d"  # The format for frame count
frames=true # output frame images by defult
basename=   # use this name for the individual frame images
animfile=   # no default animfile to output to.
coalesce=   # coalesce flag (off)
VERBOSE=    # be verbose (off)
time_sync=0 # add time sync comments (off)
awk=awk

while [ $# -gt 0 ]; do
  case "$1" in
  --help | --doc*) Usage ;;
  -c) coalesce='-coalesce' ;; # coalesce animation (on)
  -t) time_sync=1 ;;          # add time sync comments (on)
  -l)
    frames=
    animfile=-
    ;;                # List anim file to standard output, no images
  -v) VERBOSE=true ;; # Be verbose
  -n) frames= ;;      # don't output any individual frames
  -s)
    sfx="$2"
    shift
    ;;             # output frame images with this suffix
  -g) sfx="gif" ;; # output frame images as GIF's (default)
  -x) sfx="xpm" ;; # output frame images as XPM's
  -i)
    init="$2"
    shift
    ;; # the count for the first frame in sequence
  -b)
    basename="$2"
    shift
    ;; # use this basename for frames
  -o)
    animfile="$2"
    shift
    ;; # output .anim file here
  --)
    shift
    break
    ;; # end of user options
  -*) Usage "Unknown option \"$1\"" ;;
  *) break ;; # end of user options
  esac
  shift # next option
done

if [ $# -eq 0 ]; then
  Usage "No GIF animations given"
fi

# echo without return (all systems, including old BSD systems)
if [ "X$(echo -n)" = "X-n" ]; then
  echo_n() { echo ${1+"$@"}'\c'; }
else
  echo_n() { echo -n ${1+"$@"}; }
fi

# Is the command available on this machine?
cmd_found() {
  case "$(type $1 2>&1)" in *'not found'*) return 1 ;; esac
  return 0
}
cmd_found convert || Usage "No IM convert command found -- ABORTING"
cmd_found identify || Usage "No IM identify command found -- ABORTING"

# output animation details to stdout
if [ "X$animfile" = "X-" ]; then
  animfile=/dev/stdout
fi

# use gawk if it is present. (mawk has a parse bug under Ubuntu)
if cmd_found gawk; then
  awk=gawk
fi

# -----------------------------------------------------

for i in "$@"; do
  [ "$VERBOSE" ] && echo_n "converting \"$i\" "
  if [ ! -f "$i" ]; then
    echo >&2 "Unable find file \"$i\""
    continue
  fi

  # --- find out the type ---
  name="$(expr "//$i" : '.*/\([^/]*\)')"       # remove path to file
  suffix="$(expr "$name" : '.*\.\([^./]*\)$')" # extract last suffix
  name="$(expr "$name" : '\(.*\)\.[^.]*$')"    # remove last suffix

  case "$name.$suffix" in
  *.gif) ;;
  *)
    [ "$VERBOSE" ] && echo_n "${b}"
    echo >&2 "Skipping non-GIF file: \"$i\""
    continue
    ;;
  esac

  name=$(echo "$name" | sed 's/_anim//') # remove any "_anim" part
  [ "$basename" ] && name=$basename
  anim_output=${animfile:-"$name.anim"}

  # Generate Animation file...
  date=$(date +'%Y-%m-%d %R:%S')
  convert "$i" $coalesce -verbose info:- |
    $awk ' #  Parse IM "identify" output, for almost direct use by "convert"
      BEGIN { print "#";
              print "# Animation Sequence File \"'"$i"'\"";
              print "# using a BASENAME of     \"'"$name"'\"";
              print "# Extracted on    '"$date"'";
              print "#";
              frame='"$init"'
              last_offset="+0+0";
              time_sync='"$time_sync"';
              time=0;
            }
      /^  Iterations:/    { if ( ! loop ) {
                              print "-loop " $NF;
                              loop = 1;
                          } }
      /^  Delay:/         { delay = $NF; sub("x100$", "", delay); }
      /^  Dispose:/       { dispose = $NF; }
      /^  Geometry:/      { size    = substr($NF,0,match($NF, /\+/ )-1 ); }
      /^  Page geometry:/ { canvas  = substr($NF,0,match($NF, /\+/ )-1 );
                            offset  = substr($NF,  match($NF, /\+/ ) );
                          }
      # End of Image Frame -- output collected image meta-data changes
      /^  Tainted: / {
          if ( ! canvas_set ) {
              print "-page " canvas;
              canvas_set = 1;
          }
          if ( time_sync )  printf "# TIME SYNC %d\n", time;
          if ( delay != last_delay ) {
              print "-delay " delay;
              last_delay = delay;
          }
          time+=delay;
          delay=0;   # no delay is given if it remains set to zero
          if ( dispose != last_dispose ) {
              print "-dispose " dispose;
              last_dispose = dispose;
          }
          if ( offset != last_offset ) {
              printf "%-16s ", "-page " offset;
              last_offset = offset;
          } else {
              printf "%-16s ", "";
          }
          printf "BASENAME_'"$fmt.$sfx"'  # %s%s\n",
                  frame++, size, offset;
        }
      END { if ( time_sync ) printf "# LOOP TIME %d\n", time; }
    ' >"$anim_output"
  [ "$VERBOSE" -a -f "$anim_output" ] &&
    echo_n "$(grep -c \\."$sfx"\$ "${name}.anim" | tr -d ' ') frames"

  # split up GIF animation (reseting animation meta-data)
  if [ "$frames" ]; then
    convert "$i" $coalesce \
      +repage -set delay 0 -set dispose none -loop 0 \
      -scene "$init" +adjoin "${name}_${fmt}.${sfx}"

    # Fix the X Pixmap color tables
    # using the AIcons Library "xpm-fix" script
    if [ "$sfx" = 'xpm' ] && cmd_found xpm-fix; then
      [ "$VERBOSE" ] && echo_n " xpm-fixes"
      xpm-fix $(sed -n "s/ *BASENAME_/${name}_/p" ${name}.anim)
    fi
  fi

  [ "$VERBOSE" ] && echo " DONE"

done
