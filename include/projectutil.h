#ifndef ERIKS_PROJECTUTIL_H
#define ERIKS_PROJECTUTIL_H

#include <string>
#include <regex>

namespace ANALYSIS {
  /* Returns "/home/jeppe/devel/projects/analysistemplate" */
  std::string getProjectRoot() {
    return "/home/jeppe/devel/projects/analysistemplate";
  }

  /* Returns "/home/jeppe" */
  std::string getHomeDir() {
    return "/home/jeppe";
  }

  /* Returns "/home/jeppe/path/to/file.xyz" of input "~/path/to/file.xyz" */
  std::string expandPathWithTilde(const std::string &path) {
    std::string result = path.substr(path.find('~') + 1);
    result.insert(0, getHomeDir());
    return result;
  }

  /* Returns "file.xyz" of input "path/to/file.xyz", both relative and absolute paths can be specified; "~" can also be at start of path */
  std::string getBasename(const std::string &path) {
    const std::regex base_regex(R"(\w*\.\w*(?![\~\.\/]))");
    std::smatch m;
    std::regex_search(path, m, base_regex);
    return m[0];
  }

  /* Returns "file" of input "path/to/file.xyz", both relative and absolute paths can be specified; "~" can also be at start of path */
  std::string getStem(const std::string &path) {
    const std::regex stem_regex(R"((?![\~\.\/])\w*(?=[\.]))");
    std::smatch m;
    std::regex_search(path, m, stem_regex);
    return m[0];
  }

  /* Returns "xyz" of input "path/to/file.xyz", both relative and absolute paths can be specified; "~" can also be at start of path */
  std::string getExtension(const std::string &path) {
    const std::regex ext_regex(R"(\w*(?![\~\.\/\w]))");
    std::smatch m;
    std::regex_search(path, m, ext_regex);
    return m[0];
  }

  /* Returns "path/to/" of input "path/to/file.xyz", both relative and absolute paths can be specified; "~" can also be at start of path but should (probably) be expanded with expandPathWithTilde() */
  std::string getParentPath(const std::string &path) {
    std::string base = getBasename(path);
    return path.substr(0, path.size() - base.size());
  }

}
#endif //ERIKS_PROJECTUTIL_H
