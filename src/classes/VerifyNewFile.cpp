#include <sstream>

#include "VerifyNewFile.h"

VerifyNewFile::VerifyNewFile() : index(0), filebase("none"), extension("none") {}

VerifyNewFile::VerifyNewFile(std::string file) : index(0) {
  int i = file.rfind("."); // find the last period in the filename
  if (i != std::string::npos) {
    setFileBase(file.substr(0,i)); // the basename is everything up to the last dot
    setFileExtension(file.substr(i + 1)); // the extension is everything else
  } else { // otherwise, there is no extension, so we assume only the basename is given
    setFileBase(file);
    setFileExtension("");
  }
}

std::string VerifyNewFile::getFullFilename() const {
  std::stringstream ss;
  if (getIndex() == 0) {
    ss << getFileBase() << "." << getFileExtension();
    return ss.str();
  }
  else {
    ss << getFileBase() << "_" << getIndex() << "." << getFileExtension();
    return ss.str();
  }
}

std::string VerifyNewFile::validNewFile() {
  std::string newfile = getFullFilename();
  while (fileExists(newfile)) {
    setLastFile(newfile);
    increaseIndex();
    newfile = getFullFilename();
  }

  return newfile;
}
