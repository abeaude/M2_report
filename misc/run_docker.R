run_container <- function(container_id, name = NULL, ports = NULL, volumes = NULL, workspace = NULL, environment_variables = NULL, verbose = FALSE){
  config <- babelwhale::get_default_config()
  
  ###############################
  #### PREPROCESS PARAMETERS ####
  ###############################
  # determine executable
  processx_command <- config[["docker_cmd"]]#Sys.which(config$backend) %>% unname()
  
  # add safe tempdir to volumes
  # safe_tmp <- dynutils::safe_tempdir("tmp")
  # on.exit(unlink(safe_tmp, recursive = TRUE))
  # volumes <- c(volumes, paste0(fix_windows_path(safe_tmp), ":/tmp2"))
  
  if (config$backend == "docker") {
    volumes <- unlist(purrr::map(volumes, function(x) c("-v", x)))
    ports <- unlist(purrr::map(ports, function(x) c("-p", x)))
  } else if (config$backend == "singularity") {
    volumes <- c("-B", paste0(volumes, collapse = ","))
  }
  
  
  # process workspace
  if (!is.null(workspace)) {
    if (config$backend == "docker") {
      workspace <- c("--workdir", workspace)
    } else if (config$backend == "singularity") {
      workspace <- c("--pwd", workspace)
    }
  }
  
  # add tmpdir to environment variables
  # environment_variables <- c(environment_variables, "TMPDIR=/tmp2")
  
  
  #################################
  #### CREATE DOCKER ARGUMENTS ####
  #################################
  if (config$backend == "docker") {
    # is entrypoint given
    # if (!is.null(command)) {
    #   command <- c("--entrypoint", command, "--rm")
    #   if (debug) {
    #     command <- c(command, "-it")
    #   }
    # }
    
    # give it a name
    if(is.null(name)){
      name <- dynutils::random_time_string("container")
    }

    command <- c("--name", name)
    
    # add -e flags to each environment variable
    env <- unlist(purrr::map(environment_variables, function(x) c("-e", x)))
    
    # do not pass env directly to processx
    processx_env <- NULL
  
    # determine command arguments
    processx_args <- c("run","-d", command, env, workspace,ports, volumes, container_id)
    
    ######################################
    #### CREATE SINGULARITY ARGUMENTS ####
    ######################################
  } else if (config$backend == "singularity") {
    # create tmpdir
    tmpdir <- dynutils::safe_tempdir("singularity_tmpdir")
    on.exit(unlink(tmpdir, force = TRUE, recursive = TRUE))
    
    processx_env <- c(
      set_names(
        environment_variables %>% gsub("^.*=", "", .),
        environment_variables %>% gsub("^(.*)=.*$", "SINGULARITYENV_\\1", .)
      ),
      "SINGULARITY_TMPDIR" = tmpdir,
      "SINGULARITY_CACHEDIR" = config$cache_dir,
      "PATH" = Sys.getenv("PATH") # pass the path along
    )
    
    container <- paste0("docker://", container_id)
    
    # determine command arguments
    processx_args <- c(
      ifelse(is.null(command), "run", "exec"),
      "--containall",
      workspace, volumes, container, command, args
    )
  }
  
  
  #########################
  #### EXECUTE COMMAND ####
  #########################
  # if (debug) {
  #   processx_env_str <- if (length(processx_env) > 0) paste0(names(processx_env), "=", processx_env, collapse = " ") else NULL
  #   command <- paste0(c(processx_env_str, processx_command, processx_args), collapse = " ")
  #   message("Use this command to enter the container: \n", crayon::bold(command))
  #   
  #   processx_args <- processx_args[processx_args != "-it"]
  # }

  # run container
  process <- processx::process$new(command = processx_command,
                        args = processx_args,
                        env = processx_env,
                        echo_cmd = verbose,
                        cleanup_tree = TRUE)
  # process <- processx::run(
  #   command = processx_command,
  #   args = processx_args,
  #   env = processx_env,
  #   echo = verbose,
  #   echo_cmd = verbose,
  #   spinner = TRUE,
  #   error_on_status = FALSE,
  #   cleanup_tree = TRUE
  # )
  
  while(!container_is_running(name)){
    Sys.sleep(5)
  }
  
  return(process)
}

# logs
docker_logs <- function(container_id){
  config <- babelwhale::get_default_config()
  processx::run(config[["docker_cmd"]], c("logs",container_id), echo = FALSE)
}

docker_stop <- function(container_id){
  out <- processx::run(config[["docker_cmd"]], c("stop",container_id), echo = FALSE)
  if(out$status != 0) message("Cannot stop container : ", container_id, "\n", out$stderr)
}

docker_rm <- function(container_id){
  config <- babelwhale::get_default_config()
  if(container_is_running(container_id)){
    message("Container : ", container_id, "is running. It will be stopped before removing it")
    docker_stop(container_id)
  }
  out <- processx::run(config[["docker_cmd"]], c("rm",container_id), echo = FALSE)
  if(out$status != 0) message("Cannot remove container : ", container_id, "\n", out$stderr)
}

# status of container
docker_ps <- function(){
  config <- babelwhale::get_default_config()
  columns <- c("ID","Image", "Command","CreatedAt","Status","Ports","Names")
  format <- paste(paste0("{{.", columns, "}}"), collapse = "\t")
  stdout <- processx::run(config[["docker_cmd"]], c("ps","-a", paste0("--format=", format)), echo = FALSE)$stdout
  
  if (stdout != "") {
    read.delim(text = stdout, sep = "\t", header = FALSE, stringsAsFactors = FALSE) %>% magrittr::set_colnames(columns)
  } else {
    purrr::map(columns, ~ character(0)) %>%
      purrr::set_names(columns) %>%
      tibble::as_tibble()
  }
}

container_is_running <- function(container_id){
  status <- docker_ps() %>% 
    dplyr::filter(Names == container_id) %>% 
    dplyr::pull(Status) %>% 
    stringr::str_detect("Up")
  if(length(status) == 0) return(FALSE)
  return(status)
}

