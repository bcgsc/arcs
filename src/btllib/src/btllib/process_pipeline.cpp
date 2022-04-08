#include "btllib/process_pipeline.hpp"
#include "btllib/status.hpp"
#include "btllib/util.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <cassert>
#include <cerrno>
#include <csignal>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <vector>

#include <dlfcn.h>     // NOLINT
#include <fcntl.h>     // NOLINT
#include <sys/stat.h>  // NOLINT
#include <sys/types.h> // NOLINT
#include <sys/wait.h>  // NOLINT
#include <unistd.h>    // NOLINT

namespace btllib {

static constexpr int PIPE_READ_END = 0;
static constexpr int PIPE_WRITE_END = 1;
static constexpr int COMM_BUFFER_SIZE = 1024;
static constexpr mode_t OPEN_MODE =
  S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH;

/// @cond HIDDEN_SYMBOLS
enum class PipelineOperation
{
  RUN,
  END
};

class ProcessPipelineInternal;
/// @endcond

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static bool process_spawner_initialized = false;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static int process_spawner_user2spawner_fd[2];
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static int process_spawner_spawner2user_fd[2];
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static std::mutex process_spawner_comm_mutex;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static std::string pipepath_prefix;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
static std::map<PipelineId, ProcessPipelineInternal> pipeline_map;

static PipeId
new_pipe_id()
{
  static PipeId last_pipe_id = 0;
  return last_pipe_id++;
}

static PipelineId
new_pipeline_id()
{
  static PipelineId last_pipeline_id = 0;
  return last_pipeline_id++;
}

static std::string
get_pipepath(const PipeId id)
{
  return pipepath_prefix + "btllib-" + std::to_string(getpid()) + "-" +
         std::to_string(id);
}

static bool
process_spawner_init();
static const bool PROCESS_PIPELINE_INITIALIZER = process_spawner_init();

template<int fd[]>
static bool
read_from(void* buf, const size_t n)
{
  ssize_t so_far = 0, ret;
  while (so_far < ssize_t(n)) {
    ret =
      read(fd[PIPE_READ_END], (uint8_t*)(buf) + so_far, ssize_t(n) - so_far);
    if (ret <= 0) {
      if (ret < 0 && errno == EINTR) {
        continue;
      }
      return false;
    }
    so_far += ret;
  }
  return true;
}

template<int fd[]>
static bool
write_to(const void* buf, size_t n)
{
  ssize_t so_far = 0, ret;
  while (so_far < ssize_t(n)) {
    ret =
      write(fd[PIPE_WRITE_END], (uint8_t*)(buf) + so_far, ssize_t(n) - so_far);
    if (ret <= 0) {
      if (ret < 0 && errno == EINTR) {
        continue;
      }
      return false;
    }
    so_far += ret;
  }
  return true;
}

#define read_from_spawner(buf, count)                                          \
  read_from<process_spawner_spawner2user_fd>(buf, count)
#define write_to_spawner(buf, count)                                           \
  write_to<process_spawner_user2spawner_fd>(buf, count)
#define read_from_user(buf, count)                                             \
  read_from<process_spawner_user2spawner_fd>(buf, count)
#define write_to_user(buf, count)                                              \
  write_to<process_spawner_spawner2user_fd>(buf, count)

/// @cond HIDDEN_SYMBOLS
struct IORedirection
{
  std::string in, out, err;
  bool out_append = false, err_append = false;
};
/// @endcond

static IORedirection
extract_io_redirection(std::string& cmd)
{
  IORedirection redirection;

  auto& in = redirection.in;
  auto& out = redirection.out;
  auto& err = redirection.err;
  auto& out_append = redirection.out_append;
  auto& err_append = redirection.err_append;

  auto args = split(cmd, " ");
  std::for_each(args.begin(), args.end(), (void (*)(decltype(args[0])&))trim);
  args.erase(std::remove_if(args.begin(),
                            args.end(),
                            [](const std::string& arg) { return arg.empty(); }),
             args.end());

  auto in_it = args.end();
  for (auto it = args.begin(); it != args.end(); ++it) {
    if (!it->empty() && it->at(0) == '<') {
      in_it = it;
      if (it->size() == 1) {
        ++it;
        check_error(it == args.end(),
                    "Process pipeline: Missing stdin source.");
        in = *it;
      } else {
        in = it->substr(1);
      }
      break;
    }
  }
  if (in_it != args.end()) {
    if (in_it->size() == 1) {
      args.erase(in_it + 1);
    }
    args.erase(in_it);
  }

  auto out_it = args.end();
  for (auto it = args.begin(); it != args.end(); ++it) {
    if (!it->empty() && it->at(0) == '>') {
      out_it = it;
      if (it->size() == 1) {
        ++it;
        check_error(it == args.end(),
                    "Process pipeline: Missing stdout target.");
        out = *it;
      } else if (it->at(1) == '>') {
        out_append = true;
        if (it->size() == 2) {
          ++it;
          check_error(it == args.end(),
                      "Process pipeline: Missing stdout target.");
          out = *it;
        } else {
          out = it->substr(2);
        }
      } else {
        out = it->substr(1);
      }
      break;
    }
  }
  if (out_it != args.end()) {
    if ((out_it->size() == 1) || (out_it->size() == 2 && (*out_it) == ">>")) {
      args.erase(out_it + 1);
    }
    args.erase(out_it);
  }

  auto err_it = args.end();
  for (auto it = args.begin(); it != args.end(); ++it) {
    if (it->size() >= 2 && it->at(0) == '2' && it->at(1) == '>') {
      err_it = it;
      if (it->size() == 2) {
        ++it;
        check_error(it == args.end(),
                    "Process pipeline: Missing stderr target.");
        err = *it;
      } else if (it->at(2) == '>') {
        err_append = true;
        if (it->size() == 3) {
          ++it;
          check_error(it == args.end(),
                      "Process pipeline: Missing stderr target.");
          err = *it;
        } else {
          err = it->substr(3);
        }
      } else {
        err = it->substr(2);
      }
      break;
    }
  }
  if (err_it != args.end()) {
    if ((err_it->size() == 2) || (err_it->size() == 3 && (*err_it) == "2>>")) {
      args.erase(err_it + 1);
    }
    args.erase(err_it);
  }

  cmd = join(args, " ");
  return redirection;
}

static void
open_comm_pipes(
  const std::string* io_filepaths,
  int* const comm_pipe_fd) // NOLINT(readability-non-const-parameter)
{
  size_t len;
  char confirmation = 0;
  std::string pipepath;
  const int flags[] = { O_RDONLY, O_WRONLY };
  for (size_t i = 0; i < 2; i++) {
    const std::string& filepath = io_filepaths[i];
    int& pipe_fd = comm_pipe_fd[i];
    if (filepath.empty()) {
      pipepath = get_pipepath(new_pipe_id());
      if (access(pipepath.c_str(), F_OK) != -1) {
        unlink(pipepath.c_str());
      }
      mkfifo(pipepath.c_str(), OPEN_MODE);

      len = pipepath.size() + 1;
      check_error(len > COMM_BUFFER_SIZE,
                  "Pipe path length too large for the buffer.");
      check_error(!write_to_user(&len, sizeof(len)) ||
                    !write_to_user(pipepath.c_str(), len),
                  "Process pipeline: Communication failure.");

      if (i > 0) {
        check_error(!read_from_user(&confirmation, sizeof(confirmation)),
                    "Process pipeline: Communication failure.");
      }

      pipe_fd = open(pipepath.c_str(), flags[i] | O_NONBLOCK | O_CLOEXEC);
      check_error(pipe_fd < 0,
                  "Process pipeline: opening comm pipe failed: " +
                    get_strerror());
      check_error(!write_to_user(&confirmation, sizeof(confirmation)),
                  "Process pipeline: Communication failure.");

      if (i == 0) {
        check_error(!read_from_user(&confirmation, sizeof(confirmation)),
                    "Process pipeline: Communication failure.");
      }

      unlink(pipepath.c_str());

      const auto status_flags = fcntl(pipe_fd, F_GETFL);
      check_error(status_flags == -1,
                  "Process pipeline: fcntl error: " + get_strerror());
      check_error(fcntl(pipe_fd, F_SETFL, status_flags & ~O_NONBLOCK) == -1,
                  "Process pipeline: fcntl error: " + get_strerror());
      check_error(!write_to_user(&confirmation, sizeof(confirmation)),
                  "Process pipeline: Communication failure.");
    } else {
      len = 0;
      check_error(!write_to_user(&len, sizeof(len)),
                  "Process pipeline: Communication failure.");
    }
  }
}

static void
redirect_io(const int in_fd, const int out_fd, const int err_fd)
{
  if (in_fd != STDIN_FILENO) {
    check_error(dup2(in_fd, STDIN_FILENO) == -1,
                "Process pipeline: dup2 failed: " + get_strerror());
    check_error(close(in_fd) != 0,
                "Process pipeline: File descriptor close error: " +
                  get_strerror());
  }
  if (out_fd != STDOUT_FILENO) {
    check_error(dup2(out_fd, STDOUT_FILENO) == -1,
                "Process pipeline: dup2 failed: " + get_strerror());
    check_error(close(out_fd) != 0,
                "Process pipeline: File descriptor close error: " +
                  get_strerror());
  }
  if (err_fd != STDERR_FILENO) {
    check_error(dup2(err_fd, STDERR_FILENO) == -1,
                "Process pipeline: dup2 failed: " + get_strerror());
    check_error(close(err_fd) != 0,
                "Process pipeline: File descriptor close error: " +
                  get_strerror());
  }
}

static void
assign_process_cmd(const std::vector<std::string>& args)
{
  char* const* argv = new char*[args.size() + 2];
  ((char*&)(argv[0])) = (char*)(args[0].c_str());
  for (size_t i = 0; i < args.size(); i++) {
    ((char*&)(argv[i + 1])) = (char*)(args[i].c_str());
  }
  ((char*&)(argv[args.size() + 1])) = nullptr;

  execvp(argv[0], argv + 1);

  std::string argv_print = argv[0];
  for (int i = 1; argv[i] != nullptr; i++) {
    argv_print += " " + std::string(argv[i]);
  }
  log_error("exec failed: " + argv_print);
  exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
}

static void
rm_pipes()
{
  const auto last_id = new_pipe_id();
  for (PipeId id = 0; id < last_id; id++) {
    const auto fpath = get_pipepath(id);
    if (access(fpath.c_str(), F_OK) != -1) {
      unlink(fpath.c_str());
    }
  }
}

static bool
check_child_failure(const int status,
                    const pid_t pid,
                    const std::string& cmd = "")
{
  if (status != 0) {
    if (WIFSIGNALED(status) && WTERMSIG(status) == SIGPIPE) {
      return false;
    }
    std::stringstream ss;
    ss << "A helper process has finished unsuccessfully:\n";
    if (!cmd.empty()) {
      ss << "Command: " << cmd << '\n';
    }
    ss << "PID: " << pid << '\n' << "Outcome: ";
    if (WIFEXITED(status)) {                              // NOLINT
      ss << "exited with status " << WEXITSTATUS(status); // NOLINT
    } else if (WIFSIGNALED(status)) {                     // NOLINT
      ss << "killed by signal " << WTERMSIG(status);      // NOLINT
    } else {
      ss << "exited with code " << status;
    }
    log_error(ss.str());
    return true;
  }
  return false;
}

static bool
check_children_failures()
{
  // Checks if any children have failed so the caller can be a disappointed
  // parent.
  int status;
  pid_t pid;
  bool failed = false;
  while ((pid = waitpid(-1, &status, WNOHANG)) > 0) {
    failed = check_child_failure(status, pid) || failed;
  }
  return failed;
}

static void
install_signal_handlers_spawner()
{
  struct sigaction action; // NOLINT
  action.sa_flags = SA_RESTART;
  sigemptyset(&action.sa_mask);

  action.sa_handler = [](const int sig) {
    (void)sig;
    const auto prev_errno = errno;
    if (check_children_failures()) {
      rm_pipes();
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    }
    errno = prev_errno;
  };
  sigaction(SIGCHLD, &action, nullptr);

  action.sa_handler = [](const int sig) {
    (void)sig;
    rm_pipes();
    std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
  };
  sigemptyset(&action.sa_mask);
  sigaction(SIGHUP, &action, nullptr);
  sigaction(SIGQUIT, &action, nullptr);
  sigaction(SIGILL, &action, nullptr);
  sigaction(SIGABRT, &action, nullptr);
  sigaction(SIGBUS, &action, nullptr);
  sigaction(SIGSEGV, &action, nullptr);
  sigaction(SIGPIPE, &action, nullptr);
  sigaction(SIGTERM, &action, nullptr);
}

/// @cond HIDDEN_SYMBOLS
class ProcessPipelineInternal
{

public:
  ProcessPipelineInternal() {}

  void end();

  PipelineId id = 0;
  std::vector<std::pair<std::string, pid_t>> cmds;
  bool ended = false;
};
/// @endcond

void
ProcessPipelineInternal::end()
{
  if (!ended) {
    ended = true;
    int status;
    for (const auto& cmd : cmds) {
      status = 0;
      const auto ret = waitpid(cmd.second, &status, 0);
      check_error((ret == -1) && (errno != ECHILD),
                  "Process pipeline: waitpid failed: " + get_strerror());
      if (ret != -1 && check_child_failure(status, cmd.second, cmd.first)) {
        std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
      }
    }
    if (check_children_failures()) {
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    }
  }
}

static void
set_comm_pipes(
  const std::vector<IORedirection>& redirections,
  const size_t idx,
  const size_t last_idx,
  int* const comm_pipe_fd, // NOLINT(readability-non-const-parameter)
  int& in_fd,
  int& out_fd)
{
  if ((idx == 0) && (redirections[idx].in.empty())) {
    in_fd = comm_pipe_fd[0];
  }
  if ((idx == last_idx) && redirections[idx].out.empty()) {
    out_fd = comm_pipe_fd[1];
  }
}

static void
open_redirection_files(const IORedirection& redirection,
                       int& in_fd,
                       int& out_fd,
                       int& err_fd)
{
  static const int open_read_flags = O_RDONLY;
  static const int open_write_flags = O_WRONLY | O_CREAT;

  if (!redirection.in.empty()) {
    in_fd = open(redirection.in.c_str(), open_read_flags, OPEN_MODE);
  }
  if (!redirection.out.empty()) {

    out_fd =
      open(redirection.out.c_str(),
           open_write_flags | (redirection.out_append ? O_APPEND : O_TRUNC),
           OPEN_MODE);
  }
  if (!redirection.err.empty()) {
    err_fd =
      open(redirection.err.c_str(),
           open_write_flags | (redirection.err_append ? O_APPEND : O_TRUNC),
           OPEN_MODE);
  }
}

static void
run_cmd()
{
  char buf[COMM_BUFFER_SIZE];
  size_t len;

  check_error(!read_from_user(&len, sizeof(len)) || !read_from_user(buf, len),
              "Process pipeline: Communication failure.");

  auto individual_cmds = split(buf, "|");
  check_error(individual_cmds.empty(), "Process spawner: Invalid command.");
  std::for_each(individual_cmds.begin(),
                individual_cmds.end(),
                (void (*)(decltype(individual_cmds[0])&))trim);

  const auto original_cmds = individual_cmds;

  std::vector<IORedirection> redirections;
  size_t idx = 0;
  for (auto& cmd : individual_cmds) {
    const auto redirection = extract_io_redirection(cmd);
    if (idx > 0) {
      check_error(
        !redirection.in.empty(),
        "Process pipeline: Only the first command may redirect stdin.");
    }
    if (idx < individual_cmds.size() - 1) {
      check_error(
        !redirection.out.empty(),
        "Process pipeline: Only the last command may redirect stdout.");
    }
    redirections.push_back(redirection);
    idx++;
  }

  std::string io_filepaths[2] = { redirections[0].in,
                                  redirections[redirections.size() - 1].out };
  int comm_pipe_fd[2] = { -1, -1 };
  open_comm_pipes(io_filepaths, comm_pipe_fd);

  int chainpipe_in_fd[2], chainpipe_out_fd[2];
  chainpipe_in_fd[PIPE_READ_END] = -1;
  chainpipe_in_fd[PIPE_WRITE_END] = -1;
  chainpipe_out_fd[PIPE_READ_END] = -1;
  chainpipe_out_fd[PIPE_WRITE_END] = -1;

  ProcessPipelineInternal pipeline;
  pipeline.id = new_pipeline_id();
  check_error(!write_to_user(&pipeline.id, sizeof(pipeline.id)),
              "Process pipeline: Communication failure.");

  const auto last_idx = ssize_t(individual_cmds.size() - 1);
  for (ssize_t idx = last_idx; idx >= 0; idx--) {
    const auto& individual_cmd = individual_cmds[idx];

    auto args = split(individual_cmd, " ");
    std::for_each(args.begin(), args.end(), (void (*)(decltype(args[0])&))trim);
    args.erase(
      std::remove_if(args.begin(),
                     args.end(),
                     [](const std::string& arg) { return arg.empty(); }),
      args.end());
    std::for_each(args.begin(), args.end(), [](decltype(args[0])& s) {
      if ((startswith(s, "'") && endswith(s, "'")) ||
          (startswith(s, "\"") && endswith(s, "\""))) {
        s.erase(0, 1);
        s.erase(s.size() - 1, 1);
      }
    });

    if (idx > 0) {
      check_error(pipe(chainpipe_in_fd) == -1,
                  "Process pipeline: Error opening a pipe.");
      check_error(fcntl(chainpipe_in_fd[PIPE_READ_END], F_SETFD, FD_CLOEXEC) ==
                    -1,
                  "Process pipeline: fcntl error: " + get_strerror());
      check_error(fcntl(chainpipe_in_fd[PIPE_WRITE_END], F_SETFD, FD_CLOEXEC) ==
                    -1,
                  "Process pipeline: fcntl error: " + get_strerror());
    }

    const pid_t pid = fork();
    if (pid == 0) {
      int in_fd = chainpipe_in_fd[PIPE_READ_END];
      int out_fd = chainpipe_out_fd[PIPE_WRITE_END];
      int err_fd = STDERR_FILENO;

      set_comm_pipes(redirections, idx, last_idx, comm_pipe_fd, in_fd, out_fd);
      open_redirection_files(redirections[idx], in_fd, out_fd, err_fd);

      check_error(in_fd < 0 || out_fd < 0 || err_fd < 0,
                  "Process pipeline: Invalid file descriptor.");
      redirect_io(in_fd, out_fd, err_fd);

      assign_process_cmd(args);
    }
    check_error(pid == -1, "Process pipeline: Error on fork.");

    pipeline.cmds.emplace_back(original_cmds[idx], pid);

    if (idx < last_idx) {
      check_error(close(chainpipe_out_fd[PIPE_READ_END]) != 0,
                  "Process pipeline: Pipe close error: " + get_strerror());
      check_error(close(chainpipe_out_fd[PIPE_WRITE_END]) != 0,
                  "Process pipeline: Pipe close error: " + get_strerror());
    }

    if (idx > 0) {
      chainpipe_out_fd[PIPE_READ_END] = chainpipe_in_fd[PIPE_READ_END];
      chainpipe_out_fd[PIPE_WRITE_END] = chainpipe_in_fd[PIPE_WRITE_END];
    }
  }

  for (const int pipe_fd : comm_pipe_fd) {
    if (pipe_fd != -1) {
      check_error(close(pipe_fd) != 0,
                  "Process pipeline: Pipe close error: " + get_strerror());
    }
  }

  pipeline_map[pipeline.id] = pipeline;
}

static void
end_cmd()
{
  ProcessPipelineInternal pipeline;
  PipelineId pipeline_id;
  char confirmation = 0;

  check_error(!read_from_user(&pipeline_id, sizeof(pipeline_id)),
              "Process pipeline: Communication failure.");
  pipeline = pipeline_map[pipeline_id];
  pipeline.end();
  pipeline_map.erase(pipeline_id);
  check_error(!write_to_user(&confirmation, sizeof(confirmation)),
              "Process pipeline: Communication failure.");
}

static void
process_spawner_operation()
{
  PipelineOperation op;
  for (;;) {
    if (read_from_user(&op, sizeof(op))) {
      switch (op) {
        case PipelineOperation::RUN:
          run_cmd();
          break;
        case PipelineOperation::END:
          end_cmd();
          break;
        default:
          log_error("Pipeline process: Invalid pipeline operation.");
          std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
      }
    } else {
      break;
    }
  }
}

static void
set_pipepath_prefix()
{
  const auto* const tmpdir =
    std::getenv("TMPDIR"); // NOLINT(concurrency-mt-unsafe)
  if (tmpdir != nullptr) {
    pipepath_prefix = tmpdir;
    const auto tmpdirlen = std::strlen(tmpdir);
    if (tmpdirlen > 0 && tmpdir[tmpdirlen - 1] != '/') {
      pipepath_prefix += '/';
    }
  } else {
    struct stat info
    {};
    if ((stat("/tmp", &info) == 0) && ((info.st_mode & S_IFMT) == S_IFDIR)) {
      pipepath_prefix = "/tmp/";
    } else {
      pipepath_prefix = "";
    }
  }
}

static std::array<int, 2>
start_watchdog()
{
  std::array<int, 2> watchdog_pipe{ -1, -1 };
  check_error(pipe(watchdog_pipe.data()) == -1,
              "Process pipeline: Error opening a pipe.");
  (new std::thread([watchdog_pipe]() {
    char dummy;
    if (read(watchdog_pipe[PIPE_READ_END], &dummy, sizeof(dummy)) <= 0) {
      log_error("Process pipeline: Spawner process failed.");
      std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
    }
  }))
    ->detach();
  return watchdog_pipe;
}

static bool
process_spawner_init()
{
  (void)PROCESS_PIPELINE_INITIALIZER;
  if (!process_spawner_initialized) {
    process_spawner_user2spawner_fd[PIPE_READ_END] = -1;
    process_spawner_user2spawner_fd[PIPE_WRITE_END] = -1;
    process_spawner_spawner2user_fd[PIPE_READ_END] = -1;
    process_spawner_spawner2user_fd[PIPE_WRITE_END] = -1;
    check_error(pipe(process_spawner_user2spawner_fd) == -1,
                "Process pipeline: Error opening a pipe.");
    check_error(pipe(process_spawner_spawner2user_fd) == -1,
                "Process pipeline: Error opening a pipe.");

    set_pipepath_prefix();

    const auto watchdog_pipe = start_watchdog();

    const pid_t pid = fork();
    if (pid == 0) {
      // This (and setpgid in parent) are necessary in order to prevent
      // the spawner process from receiving the same signals as the parent
      // (e.g. SIGINT). If the parent receives SIGINT, it can decide how
      // to handle it or whether to ignore it. The spawner will simply
      // exist in the background and die if the parent dies as well.
      check_error(setpgid(0, 0) != 0,
                  "Process pipeline: setpgid failed in spawner process.");

      check_error(close(process_spawner_user2spawner_fd[PIPE_WRITE_END]) != 0,
                  "Process pipeline: Pipe close error: " + get_strerror());
      check_error(close(process_spawner_spawner2user_fd[PIPE_READ_END]) != 0,
                  "Process pipeline: Pipe close error: " + get_strerror());
      check_error(close(watchdog_pipe[PIPE_READ_END]) != 0,
                  "Process pipeline: Pipe close error: " + get_strerror());

      install_signal_handlers_spawner();

      process_spawner_operation();

      if (check_children_failures()) {
        rm_pipes();
        std::exit(EXIT_FAILURE); // NOLINT(concurrency-mt-unsafe)
      }
      std::exit(EXIT_SUCCESS); // NOLINT(concurrency-mt-unsafe)
    }
    check_error(setpgid(pid, pid) < 0 && errno != EACCES,
                "Process pipeline: setpgid failed in parent process.");

    check_error(close(process_spawner_user2spawner_fd[PIPE_READ_END]) != 0,
                "Process pipeline: Pipe close error: " + get_strerror());
    check_error(close(process_spawner_spawner2user_fd[PIPE_WRITE_END]) != 0,
                "Process pipeline: Pipe close error: " + get_strerror());
    check_error(close(watchdog_pipe[PIPE_WRITE_END]) != 0,
                "Process pipeline: Pipe close error: " + get_strerror());

    process_spawner_initialized = true;
  }
  return true;
}

ProcessPipeline::ProcessPipeline(const std::string& cmd)
{
  std::unique_lock<std::mutex> lock(process_spawner_comm_mutex);

  const auto op = PipelineOperation::RUN;
  check_error(!write_to_spawner(&op, sizeof(op)),
              "Process pipeline: Communication failure.");

  size_t len = cmd.size() + 1;

  check_error(len > COMM_BUFFER_SIZE,
              "Process pipeline: Stream path length too large for the buffer.");
  check_error(!write_to_spawner(&len, sizeof(len)) ||
                !write_to_spawner(cmd.c_str(), len),
              "Process pipeline: Communication failure.");

  char buf[COMM_BUFFER_SIZE];
  char confirmation = 0;

  FILE** handles[] = { &in, &out };
  const int flags[] = { O_WRONLY, O_RDONLY };
  const char* modes[] = { "w", "r" };
  for (size_t i = 0; i < 2; i++) {
    check_error(!read_from_spawner(&len, sizeof(len)),
                "Process pipeline: Communication failure.");
    if (len > 0) {
      check_error(!read_from_spawner(buf, len),
                  "Process pipeline: Communication failure.");

      if (i == 0) {
        check_error(!read_from_spawner(&confirmation, sizeof(confirmation)),
                    "Process pipeline: Communication failure.");
      }

      const int pipe_fd = open(buf, flags[i] | O_NONBLOCK);
      check_error(pipe_fd < 0,
                  "Process pipeline: opening comm pipe failed: " +
                    get_strerror());
      check_error(!write_to_spawner(&confirmation, sizeof(confirmation)),
                  "Process pipeline: Communication failure.");

      if (i > 0) {
        check_error(!read_from_spawner(&confirmation, sizeof(confirmation)),
                    "Process pipeline: Communication failure.");
      }

      const auto status_flags = fcntl(pipe_fd, F_GETFL);
      check_error(status_flags == -1,
                  "Process pipeline: fcntl error: " + get_strerror());
      check_error(fcntl(pipe_fd, F_SETFL, status_flags & ~O_NONBLOCK) == -1,
                  "Process pipeline: fcntl error: " + get_strerror());
      check_error(!read_from_spawner(&confirmation, sizeof(confirmation)),
                  "Process pipeline: Communication failure.");

      (*(handles[i])) = fdopen(pipe_fd, modes[i]);
    }
  }

  check_error(!read_from_spawner(&id, sizeof(id)),
              "Process pipeline: Communication failure.");
}

static void
closefile(FILE*& f)
{
  if (f != nullptr) {
    check_error(std::fclose(f) != 0,
                "Process spawner: Error closing file: " + get_strerror());
    f = nullptr;
  }
}

void
ProcessPipeline::close_in()
{
  bool in_closed_expected = false;
  if (in_closed.compare_exchange_strong(in_closed_expected, true)) {
    closefile(in);
  }
}

void
ProcessPipeline::close_out()
{
  bool out_closed_expected = false;
  if (out_closed.compare_exchange_strong(out_closed_expected, true)) {
    closefile(out);
  }
}

void
ProcessPipeline::end()
{
  bool ended_expected = false;
  if (ended.compare_exchange_strong(ended_expected, true)) {
    close_in();
    close_out();

    std::unique_lock<std::mutex> lock(process_spawner_comm_mutex);

    const auto op = PipelineOperation::END;
    check_error(!write_to_spawner(&op, sizeof(op)) ||
                  !write_to_spawner(&id, sizeof(id)),
                "Process pipeline: Communication failure.");

    char confirmation = 0;
    check_error(!read_from_spawner(&confirmation, sizeof(confirmation)),
                "Process pipeline: Communication failure.");
  }
}

} // namespace btllib
