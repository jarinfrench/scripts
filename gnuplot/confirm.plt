pause mouse keypress

if (MOUSE_KEY != 121 && MOUSE_KEY != 89 && MOUSE_KEY != 110 && MOUSE_KEY != 78) {
    print "Please enter Y|n"
    reread
} else {
    # y or Y respectively
    if (MOUSE_KEY == 121 || MOUSE_KEY == 89) {
        confirm = 1
        } else { # n or N respectively
            confirm = 0
    }
}
