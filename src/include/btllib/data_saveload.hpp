#ifndef BTLLIB_DATA_SAVELOAD_HPP
#define BTLLIB_DATA_SAVELOAD_HPP

#include "status.hpp"
#include "util.hpp"

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <csignal>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

#include <dlfcn.h>
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

namespace btllib {

enum SaveloadOp
{
  READ,
  WRITE,
  APPEND
};

class _Pipeline
{

public:
  _Pipeline() {}

  _Pipeline(FILE* file, pid_t pid_first, pid_t pid_last)
    : file(file)
    , pid_first(pid_first)
    , pid_last(pid_last)
  {}

  FILE* file = nullptr;
  pid_t pid_first = -1;
  pid_t pid_last = -1;
};

inline std::string
get_saveload_cmd(const std::string& path, SaveloadOp op);
inline _Pipeline
run_saveload_cmd(const std::string& cmd, SaveloadOp op);

class DataSource
{
public:
  DataSource(const std::string& source)
  {
    if (source == "-") {
      pipeline.file = stdin;
      pipeline.pid_first = -1;
      pipeline.pid_last = -1;
    } else {
      const auto cmd = get_saveload_cmd(source, READ);
      check_error(cmd.empty(), "Error loading from " + source);
      pipeline = run_saveload_cmd(cmd, READ);
    }
  }

  ~DataSource() { close(); }

  void close()
  {
    if (!closed) {
      if (pipeline.file != stdin) {
        int status;
        kill(pipeline.pid_first, SIGTERM);
        waitpid(pipeline.pid_last, &status, 0);
        std::fclose(pipeline.file);
      }
      closed = true;
    }
  }

  FILE* operator*() { return pipeline.file; }
  FILE* operator->() { return pipeline.file; }
  operator FILE*() { return pipeline.file; }

  _Pipeline pipeline;
  bool closed = false;
};

class DataSink
{
public:
  DataSink(const std::string& sink, bool append)
  {
    if (sink == "-") {
      pipeline.file = stdout;
      pipeline.pid_first = -1;
      pipeline.pid_last = -1;
    } else {
      const auto cmd = get_saveload_cmd(sink, append ? APPEND : WRITE);
      check_error(cmd.empty(), "Error saving to " + sink);
      pipeline = run_saveload_cmd(cmd, append ? APPEND : WRITE);
    }
  }

  ~DataSink() { close(); }

  void close()
  {
    if (!closed) {
      if (pipeline.file != stdout) {
        fclose(pipeline.file);
        int status;
        waitpid(pipeline.pid_last, &status, 0);
      }
      closed = true;
    }
  }

  FILE* operator*() { return pipeline.file; }
  FILE* operator->() { return pipeline.file; }
  operator FILE*() { return pipeline.file; }

  _Pipeline pipeline;
  bool closed = false;
};

/** SIGCHLD handler. Reap child processes and report an error if any
 * fail. */
inline void
sigchld_handler(const int sig)
{
  assert(sig == SIGCHLD);
  (void)sig;

  pid_t pid;
  int status;
  while ((pid = waitpid(-1, &status, WNOHANG)) > 0) {
    if (status != 0) {
      if (WIFEXITED(status)) { // NOLINT
        std::cerr << "PID " << pid << " exited with status "
                  << WEXITSTATUS(status) << std::endl; // NOLINT
      } else if (WIFSIGNALED(status)) {                // NOLINT
        if (WTERMSIG(status) == SIGTERM) {
          return;
        }
        std::cerr << "PID " << pid << " killed by signal "
                  << WTERMSIG(status) // NOLINT
                  << std::endl;
      } else {
        std::cerr << "PID " << pid << " exited with code " << status
                  << std::endl;
      }
      std::exit(EXIT_FAILURE);
    }
  }
  if (pid == -1 && errno != ECHILD) {
    std::perror("waitpid");
    std::exit(EXIT_FAILURE);
  }
}

bool
data_saveload_init();
static const bool data_saveload_initialized = data_saveload_init();

inline bool
data_saveload_init()
{
  (void)data_saveload_initialized;
  struct sigaction action; // NOLINT
  action.sa_handler = sigchld_handler;
  sigemptyset(&action.sa_mask);
  action.sa_flags = SA_RESTART;
  sigaction(SIGCHLD, &action, nullptr);
  return true;
}

inline std::string
get_saveload_cmd(const std::string& path, SaveloadOp op)
{
  struct Datatype
  {
    std::vector<std::string> prefixes;
    std::vector<std::string> suffixes;
    std::vector<std::string> cmds_check_existence;
    std::vector<std::string> read_cmds;
    std::vector<std::string> write_cmds;
    std::vector<std::string> append_cmds;
  };

  // clang-format off
  static const Datatype DATATYPES[]{
    { { "http://", "https://", "ftp://" }, {}, { "which wget" }, { "wget -O-" }, { "" }, { "" } },
    { {}, { ".url" }, { "which wget" }, { "wget -O- -i" }, { "" }, { "" } },
    { {}, { ".ar" }, { "which ar" }, { "ar -p" }, { "" }, { "" } },
    { {}, { ".tar" }, { "which tar" }, { "tar -xOf" }, { "" }, { "" } },
    { {}, { ".tgz" }, { "which tar" }, { "tar -zxOf" }, { "" }, { "" } },
    { {}, { ".gz", ".z" }, { "which pigz", "which gzip" }, { "pigz -dc", "gzip -dc" }, { "pigz >", "gzip >" }, { "pigz >>", "gzip >>" } },
    { {}, { ".bz2" }, { "which bzip2" }, { "bunzip2 -dc" }, { "bzip2 >" }, { "bzip2 >>" } },
    { {}, { ".xz" }, { "which xz" }, { "unxz -dc" }, { "xz -T0 >" }, { "xz -T0 >>" } },
    { {}, { ".7z" }, { "which 7z" }, { "7z -so e" }, { "7z -si a" }, { "7z -si a" } },
    { {}, { ".zip" }, { "which zip" }, { "unzip -p" }, { "" }, { "" } },
    { {}, { ".bam", ".cram" }, { "which samtools" }, { "samtools view -h" }, { "samtools -Sb - >" }, { "samtools -Sb - >>" } },
  };
  // clang-format on
  std::string default_cmd = "cat";
  if (op == WRITE) {
    default_cmd += " >";
  } else if (op == APPEND) {
    default_cmd += " >>";
  }

  std::string path_trimmed = path;
  std::vector<std::string> cmd_layers;
  for (;;) {
    bool found_datatype = false;
    for (const auto& datatype : DATATYPES) {
      size_t trim_start = 0, trim_end = 0;
      bool this_datatype = false;
      for (const auto& prefix : datatype.prefixes) {
        if (starts_with(path_trimmed, prefix)) {
          this_datatype = true;
          trim_start += prefix.size();
          break;
        }
      }
      for (const auto& suffix : datatype.suffixes) {
        if (ends_with(path_trimmed, suffix)) {
          this_datatype = true;
          trim_end += suffix.size();
          break;
        }
      }

      if (this_datatype) {
        found_datatype = true;
        bool found_cmd = false;
        int cmd_idx = 0;
        for (const auto& existence_cmd : datatype.cmds_check_existence) {
          bool good = true;
          auto sub_cmds = split(existence_cmd, "&&");
          std::for_each(sub_cmds.begin(), sub_cmds.end(), trim);
          for (const auto& sub_cmd : sub_cmds) {
            auto args = split(sub_cmd, " ");
            std::for_each(args.begin(), args.end(), trim);

            pid_t pid = fork();
            if (pid == 0) {
              int null_fd = open("/dev/null", O_WRONLY, 0);
              dup2(null_fd, STDOUT_FILENO);
              dup2(null_fd, STDERR_FILENO);
              close(null_fd);

              switch (args.size()) {
                case 1:
                  execlp(args[0].c_str(), args[0].c_str(), NULL);
                // fall through
                case 2:
                  execlp(
                    args[0].c_str(), args[0].c_str(), args[1].c_str(), NULL);
                // fall through
                case 3:
                  execlp(args[0].c_str(),
                         args[0].c_str(),
                         args[1].c_str(),
                         args[2].c_str(),
                         NULL);
                // fall through
                case 4:
                  execlp(args[0].c_str(),
                         args[0].c_str(),
                         args[1].c_str(),
                         args[2].c_str(),
                         args[3].c_str(),
                         NULL);
                // fall through
                default:
                  log_error("Invalid number of arguments supplied to execlp (" +
                            std::to_string(args.size()) + ").");
                  std::exit(EXIT_FAILURE);
              }
              log_error("execlp failed.");
              std::exit(EXIT_FAILURE);
            } else {
              check_error(pid == -1, "Error on fork.");
              int status;
              waitpid(pid, &status, 0);
              if (WIFEXITED(status) && WEXITSTATUS(status) != 0) { // NOLINT
                good = false;
                break;
              }
            }
          }
          if (good) {
            found_cmd = true;
            break;
          }
          cmd_idx++;
        }

        if (found_cmd) {
          std::string cmd;
          switch (op) {
            case READ:
              cmd = datatype.read_cmds[cmd_idx];
              break;
            case WRITE:
              cmd = datatype.write_cmds[cmd_idx];
              break;
            case APPEND:
              cmd = datatype.append_cmds[cmd_idx];
              break;
          }
          if (cmd.empty()) {
            log_warning("Filetype recognized for '" + path +
                        "', but no tool available to work with it.");
          } else {
            cmd_layers.push_back(cmd);
          }
        } else {
          log_warning("Filetype recognized for '" + path +
                      "', but no tool available to work with it.");
        }
        path_trimmed.erase(0, trim_start);
        path_trimmed.erase(path_trimmed.size() - trim_end);
      }
    }
    if (!found_datatype) {
      break;
    }
  }
  if (cmd_layers.empty()) {
    cmd_layers.push_back(default_cmd);
  }
  if (op == WRITE || op == APPEND) {
    std::reverse(cmd_layers.begin(), cmd_layers.end());
  }

  std::string result_cmd;
  for (size_t i = 0; i < cmd_layers.size(); i++) {
    auto& cmd = cmd_layers[i];
    if (op == WRITE || op == APPEND) {
      if (i == cmd_layers.size() - 1) {
        if (cmd.back() == '>') {
          cmd += path;
        } else {
          cmd += " ";
          cmd += path;
        }
      } else {
        if (cmd.back() == '>') {
          while (cmd.back() == '>' || cmd.back() == ' ') {
            cmd.pop_back();
          }
        } else {
          cmd += " -";
        }
      }
    } else {
      if (i == 0) {
        cmd += " ";
        cmd += path;
      } else {
        cmd += " -";
      }
    }
    if (i > 0) {
      result_cmd += " | ";
    }
    result_cmd += cmd;
  }

  return result_cmd;
}

inline _Pipeline
run_saveload_cmd(const std::string& cmd, SaveloadOp op)
{
  static const int READ_END = 0;
  static const int WRITE_END = 1;

  auto individual_cmds = split(cmd, " | ");
  assert(!individual_cmds.empty());
  std::reverse(individual_cmds.begin(), individual_cmds.end());
  std::vector<pid_t> pids;
  std::vector<std::vector<int>> fds;
  int input_fd[2], output_fd[2];
  input_fd[READ_END] = -1;
  input_fd[WRITE_END] = -1;
  output_fd[READ_END] = -1;
  output_fd[WRITE_END] = -1;
  if (op == READ) {
    check_error(pipe2(output_fd, O_CLOEXEC) == -1, "Error opening a pipe.");
  }
  size_t i = 0;
  for (const auto& individual_cmd : individual_cmds) {
    auto args = split(individual_cmd, " ");
    std::for_each(args.begin(), args.end(), trim);

    std::string stdout_to_file;
    decltype(args)::iterator it;
    for (it = args.begin(); it != args.end(); ++it) {
      if (it->front() == '>') {
        stdout_to_file = it->substr(1);
        break;
      }
    }
    if (it != args.end()) {
      args.erase(it);
    }

    if (op == READ) {
      if (i < individual_cmds.size() - 1) {
        check_error(pipe2(input_fd, O_CLOEXEC) == -1, "Error opening a pipe.");
      }
    } else {
      check_error(pipe2(input_fd, O_CLOEXEC) == -1, "Error opening a pipe.");
    }

    pid_t pid = fork();
    if (pid == 0) {
      if (op == READ) {
        dup2(output_fd[WRITE_END], STDOUT_FILENO);
        close(output_fd[READ_END]);
        close(output_fd[WRITE_END]);

        if (i > 0) {
          close(fds.front()[READ_END]);
          close(fds.front()[WRITE_END]);
        }

        if (i < individual_cmds.size() - 1) {
          dup2(input_fd[READ_END], STDIN_FILENO);
          close(input_fd[READ_END]);
          close(input_fd[WRITE_END]);
        }

        switch (args.size()) {
          case 1:
            execlp(args[0].c_str(), args[0].c_str(), NULL);
          // fall through
          case 2:
            execlp(args[0].c_str(), args[0].c_str(), args[1].c_str(), NULL);
          // fall through
          case 3:
            execlp(args[0].c_str(),
                   args[0].c_str(),
                   args[1].c_str(),
                   args[2].c_str(),
                   NULL);
          // fall through
          case 4:
            execlp(args[0].c_str(),
                   args[0].c_str(),
                   args[1].c_str(),
                   args[2].c_str(),
                   args[3].c_str(),
                   NULL);
          // fall through
          default:
            log_error("Invalid number of arguments supplied to execlp (" +
                      std::to_string(args.size()) + ").");
            std::exit(EXIT_FAILURE);
        }
        log_error("execlp failed.");
        std::exit(EXIT_FAILURE);
      } else {
        dup2(input_fd[READ_END], STDIN_FILENO);
        close(input_fd[READ_END]);
        close(input_fd[WRITE_END]);

        if (!stdout_to_file.empty()) {
          int outfd =
            open(stdout_to_file.c_str(),
                 O_WRONLY | O_CREAT | (op == APPEND ? O_APPEND : 0),
                 S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH);
          dup2(outfd, STDOUT_FILENO);
          close(outfd);
        } else if (i > 0) {
          dup2(output_fd[WRITE_END], STDOUT_FILENO);
          close(output_fd[READ_END]);
          close(output_fd[WRITE_END]);
        }

        switch (args.size()) {
          case 1:
            execlp(args[0].c_str(), args[0].c_str(), NULL);
          // fall through
          case 2:
            execlp(args[0].c_str(), args[0].c_str(), args[1].c_str(), NULL);
          // fall through
          case 3:
            execlp(args[0].c_str(),
                   args[0].c_str(),
                   args[1].c_str(),
                   args[2].c_str(),
                   NULL);
          // fall through
          case 4:
            execlp(args[0].c_str(),
                   args[0].c_str(),
                   args[1].c_str(),
                   args[2].c_str(),
                   args[3].c_str(),
                   NULL);
          // fall through
          default:
            log_error("Invalid number of arguments supplied to execlp (" +
                      std::to_string(args.size()) + ").");
            std::exit(EXIT_FAILURE);
        }
        log_error("execlp failed.");
        exit(EXIT_FAILURE);
      }
    }
    check_error(pid == -1, "Error on fork.");
    pids.push_back(pid);
    if (op == READ) {
      fds.push_back({ output_fd[READ_END], output_fd[WRITE_END] });
    } else {
      fds.push_back({ input_fd[READ_END], input_fd[WRITE_END] });
    }
    if (i > 0) {
      close(output_fd[READ_END]);
      close(output_fd[WRITE_END]);
    }
    output_fd[READ_END] = input_fd[READ_END];
    output_fd[WRITE_END] = input_fd[WRITE_END];
    i++;
  }

  if (op == READ) {
    close(fds.front()[WRITE_END]);
    return _Pipeline(
      fdopen(fds.front()[READ_END], "r"), pids.back(), pids.front());
  }
  close(fds.back()[READ_END]);
  return _Pipeline(
    fdopen(fds.back()[WRITE_END], "w"), pids.back(), pids.front());
}

} // namespace btllib

#endif
