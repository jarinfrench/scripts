## usage avg(column_number)
## This will probably need to be changed
## Taken from https://stackoverflow.com/a/19149931
function avg(col) {
    {sum += $col} END {if (NR > 0) print sum/NR}
}
