#include "SysUtils.h"
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

// --------------
// Low level
// --------------

#ifdef __cplusplus
extern "C" {
#endif
  static void execute_command(char **argv) {
    pid_t pid = fork();
    if (pid == 0) {
      execvp(argv[0], argv);
      perror(argv[0]);
      _exit(127);
    } else if (pid > 0) {
      int status;
      waitpid(pid, &status, 0);
    } else {
      fprintf(stderr, "Error executing command: %s\n", strerror(errno));
    }
  }
#ifdef __cplusplus
}
#endif

// ------------
// Interface
// ------------

void execute_peak_fit_bin(::std::string mz_b,
			  ::std::string mz_av,
			  ::std::string output_path) {

  char **argv = (char **) malloc(sizeof(char *) * 5); // 4 args + a null terminator
  char c_mz_b[mz_b.size() + 1];
  char c_mz_av[mz_av.size() + 1];
  char c_output_path[output_path.size() + 1];

  strcpy(c_mz_b, mz_b.c_str());
  strcpy(c_mz_av, mz_av.c_str());
  strcpy(c_output_path, output_path.c_str());
  
  argv[0] = ::InferenceAPI::python_executable_path;
  argv[1] = &c_mz_b[0];
  argv[2] = &c_mz_av[0];
  argv[3] = &c_output_path[0];
  argv[4] = NULL;
  
  execute_command(argv);
  free(argv);
}
