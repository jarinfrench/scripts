#ifndef VERIFY_NEW_FILE_H
#define VERIFY_NEW_FILE_H

#include <string>
#include <sys/stat.h>

class VerifyNewFile {
private:
  int index; // the index for the filename being examined
  std::string filebase; // the file's basename to test against
  std::string extension; // the file's extension

  // the last found existing file (e.g. if files tmp.txt, tmp_1.txt,
  // and tmp_2.txt are found, this would be tmp_2.txt, and the valid new file
  // would be tmp_3.txt)
  std::string last_file;
  void setLastFile(std::string str) {this->last_file = str;}

public:
  VerifyNewFile();
  VerifyNewFile(std::string file);

  int getIndex() const {return this->index;}
  std::string getFileBase() const {return this->filebase;}
  std::string getFileExtension() const {return this->extension;}
  std::string getLastFile() const {return this->last_file;}
  std::string getFullFilename() const;

  void setIndex(int i) {this->index = i;}
  void increaseIndex() {setIndex(this->index + 1);}
  void setFileBase(std::string str) {this->filebase = str;}
  void setFileExtension(std::string str) {this->extension = str;}

  std::string validNewFile();

  inline bool fileExists(const std::string& file) {
    struct stat buffer;
    return (stat (file.c_str(), &buffer) == 0);
  }
};

#endif // VERIFY_NEW_FILE_H
