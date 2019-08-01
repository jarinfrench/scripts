#! /usr/bin/env python
# This script accepts one .bib file as its argument, and alphabetizes the list
# by author.  The list proceeds in order as follows: articles, books, and other.
# Future versions of this script may be more specific.
# Author: Jarin French
# Date of last modification: 7 July 2016

from sys import argv

script, filename = argv

if not(filename.endswith('.bib')):
    print('File', filename, 'is not a valid BiBTeX file.')
    return

abbreviationsList = [] # List to hold the abbreviations
referencesList = [] # List to hold the list of references information
referencesInfoList = [] # list to hold the references information
indices = [] # list to hold the indices for the author names
f = open(filename, "r")
while True:
    line = f.readLine()
    if not line:
        break

    # Need to check if line starts with @
    if line[0] == '@': # If so, check the next three letters:
        if line[1:4] = 'STR': # STR indicates an abbreviation
            abbreviationsList.append(line)
        else: # Anything else indicates an entry

            # Check to make sure the entry has valid information
            if line[0] == '\n': # If it's just a newline, go to the next iteration
                continue
            else: #Put the line into the list.
                referencesInfoList.append(line)
                # Check to see if we have ended the entry
                if '}' in line: #If so, put everything from the first line to the last into the referencesList
                    referencesList.append(referencesInfoList)
                    referencesInfoList = [] # Reset the info list
                if "AUTHOR" in line:
                    indices



# MAS or PHD indicates a thesis
# ART indicates an article
# INC indicates an "in collection" (chapter)
# MIS indicates misc
# BOO indicates a book
# INB indicates "in a book" (chapter)
# Alphabetize the STRINGS list by abbreviation: start with first letter, second letter, etc.  Spaces come first
f.close()
