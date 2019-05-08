pause mouse keypress

# y or Y respectively
if (MOUSE_KEY == 121 || MOUSE_KEY == 89) {
    confirm = 1
} else { # n or N respectively
    if (MOUSE_KEY == 110 || MOUSE_KEY == 78) {
        confirm = 0
    } else {
        print "Please enter Y|N"
        reread
    }
}
