data_prep_stream <- function(sobj, filename, reductions, path, cell_identities){
  fs::path(path,"stream") %>% fs::dir_create()
  write.table(sobj@assays$RNA@counts, file=fs::path(path,"stream",paste0(filename,'.tsv')), quote=FALSE, sep='\t', col.names = TRUE,row.names = TRUE)
  write.table(sobj[[cell_identities]], file=fs::path(path,"stream",paste0(filename,'_cell_identities.tsv')), quote=FALSE, sep='\t', col.names = FALSE,row.names = TRUE)
  
  fs::path(path,"stream","default") %>% fs::dir_create()
  fs::file_copy("Trajectory/stream/template.ipynb", fs::path(path, "stream","default", paste0(filename,"_stream.ipynb")))
  
  # Preserve colors from Seurat plotting
  clusters <- sobj[[cell_identities]] %>% unique() %>% .[[cell_identities]] 
  colors <- scales::hue_pal()(nlevels(clusters)) %>% magrittr::set_names(levels(clusters)) %>% tibble::enframe() %>% dplyr::filter(name %in% clusters)
  write.table(colors , file = fs::path(path,"stream",paste0(filename,'_cell_colors.tsv')), col.names = FALSE, row.names = FALSE, sep = '\t', quote = FALSE)
}

wrapper_stream <- function(counts,raw_counts,dim_red, millestone_network, progression, grouping, root){
  # inspired by https://github.com/dynverse/ti_elpigraph/blob/master/run.R
  output <-dynwrap::wrap_data(cell_ids = progression[['cell_id']]) %>%
    dynwrap::add_trajectory(milestone_network = millestone_network,
                            progressions = progression) %>% 
    dynwrap::add_expression(counts = raw_counts,
                            expression = counts) %>%
    dynwrap::add_dimred(dimred = dim_red) %>%
    dynwrap::add_grouping(grouping = grouping) %>%
    dynwrap::add_root(root_milestone_id = root) %>%
    dynwrap::add_pseudotime()
  return(output)
}

run_stream <- function(path, filename){
  path %<>% fs::path("stream")
  
  dirs <- fs::dir_ls(path,type = "directory")
  
    if(rlang::is_interactive() & length(dirs) == 0){
        # Use specific environment STREAM
        reticulate::use_condaenv("stream")
        # launch notebook with template
        # copy template from ressources to correct directory 
        fs::path(path,"interactive") %>% fs::dir_create()
        fs::file_copy("Trajectory/stream/template.ipynb", fs::path(path,"interactive", paste0(filename,"_stream.ipynb")))
        p <- processx::process$new("jupyter", c("notebook", fs::path(path,"interactive", paste0(filename,"_stream.ipynb"))), echo = FALSE, wd = fs::path(path,"interactive"))
        # processx::run("jupyter", c("notebook", fs::path(path,"interactive", paste0(filename,"_stream.ipynb"))), echo = FALSE, wd = fs::path(path,"interactive") )
        message("Jupyter notebook is started inside the process ", p$get_pid())
        while(p$is_alive()){
          message("\rPlease do not forget to press the ", crayon::bgWhite(crayon::black("QUIT")),
                  " button to shutdown the jupyter server to continue.\n It might take up to 1 min to detect the shutdown", appendLF = FALSE)
          Sys.sleep(60)
        }
        dirs <- c(dirs,fs::path(path,"interactive"))
    } 
  output <- purrr::map(dirs,process_stream_files) %>% 
    purrr::map(~ rlang::exec(wrapper_stream, !!!.)) %>%
    purrr::set_names(., nm = fs::path_file(names(.)))
  
  return(output)
}

stream_docker <- function(docker_tag,path){
  if(is.null(docker_tag)){
    docker_tag <- "latest"
  }
  if(babelwhale::test_docker_installation()){
    config <- babelwhale::create_config()
    config[['docker_cmd']] <- Sys.which(config$backend) %>% unname()
    babelwhale::set_default_config(config, permanent = FALSE)
    
    images <- babelwhale::list_docker_images() %>% 
      dplyr::filter(Repository == 'abeaude/stream' & Tag == docker_tag)
    if(nrow(images) == 0){
      babelwhale::pull_container(paste0("abeaude/stream:",docker_tag))
    }
    
    # Run container
    stream <- run_container(paste0("abeaude/stream:",docker_tag), name = "stream", ports = "8888:8888", volumes = paste0(path,":/home/jovyan/work"))
    url <- docker_logs("stream")$stderr %>% stringr::str_extract('http://127.0.0.1:8888/\\?token=[:alnum:]*')
    token <- stringr::str_extract(url, '=[:alnum:]*') %>% stringr::str_remove("=") %>% stringr::str_c("Token ",.)
    if(container_is_running("stream")){
      message("The jupyter notebook with stream will open in your browser, please do not forget to press the ",
              crayon::bgWhite(crayon::black("QUIT")), " button to shutdown the jupyter server to continue.")
      Sys.sleep(10)
      utils::browseURL(url)
      copy <- httr::POST(url = "http://127.0.0.1:8888/api/contents/work",
                         config = httr::add_headers("Authorization" = token), 
                         body = '{"copy_from":"/template_stream.ipynb"}')
      
      new_filename <- list(path = stringr::str_c("work/",filename,"_stream_analysis.ipynb"))
      renaming <- httr::PATCH(url = "http://127.0.0.1:8888/api/contents/work/template_stream.ipynb",
                              config = httr::add_headers("Authorization" = token), 
                              body = new_filename, encode = 'json')
      utils::browseURL(paste0("http://127.0.0.1:8888/notebooks/",new_filename))
    } else {
      stop("Container cannot be launched : \n", stream$stderr)
    }
    
      is_running <- container_is_running("stream")
    while(is_running){
      is_running <- container_is_running("stream")
    }
    status_rm <- docker_rm("stream")
    return(TRUE)
  } else {
    stop(crayon::red('\u274C')," No Docker installation detected, please consider installing docker")
  }
}

# ti_stream <- function(params, path){
#   reticulate::use_condaenv(condaenv = "stream")
#   st <- reticulate::import("stream", convert = FALSE)
#   pd <- reticulate::import("pandas", convert = FALSE)
#   mplt <- reticulate::import("matplotlib", convert = FALSE)
#   np <- reticulate::import("numpy", convert = FALSE)
#   nx <- reticulate::import("networkx", convert = FALSE)
#   os <- reticulate::import('os')
#   
#   reticulate:
#   mplt.rcParams.update({'figure.max_open_warning': 0})
# 
# }

guess_stream_dir <- function(path){
  possible_dir <- fs::dir_ls(path,type = "directory") %>% 
    purrr::set_names(.) %>%
    purrr::map_dbl( ~ fs::dir_ls(.,glob = "*.csv") %>% length()) %>%
    purrr::keep(. == 6) %>%
    purrr::imap(function(.x,.y) {
      csv_files <- fs::dir_ls(.y,glob = "*.csv")
      guess_file_name <- fs::path_file(csv_files) %>% Biobase::lcPrefix()
      fs::path_ext_remove(csv_files) %>% fs::path_file() %>% stringr::str_remove(guess_file_name)
    }) %>% 
    purrr::keep(~identical(.,c("counts", "dim_red", "grouping","millestone_network", "progression", "raw_counts"))) %>%
    names()
  
  stream_dir <- possible_dir[1]
  return(stream_dir)
}

process_stream_files <- function(path){
  stream_dir <- guess_stream_dir(path)
  files <- stream_dir %>% 
    fs::dir_ls(glob = "*.csv") 
  
  
  guess_file_name <- fs::path_file(files) %>% Biobase::lcPrefix()
  root <- stringr::str_extract(guess_file_name,'S[0-9]+')
  names(files) <- fs::path_ext_remove(files) %>% fs::path_file() %>% stringr::str_remove(guess_file_name)
  
  col_spec <- list(counts = readr::cols(cell_id = readr::col_character(),.default = readr::col_double()),
                   dim_red = readr::cols(cell_id = readr::col_character(),.default = readr::col_double()),
                   grouping = readr::cols(cell_id = readr::col_character(),label = readr::col_character()) ,
                   millestone_network = readr::cols(from = readr::col_character(),to = readr::col_character(),length = readr::col_double(),directed = readr::col_logical()),
                   progression = readr::cols(cell_id = readr::col_character(),percentage = readr::col_double(),from = readr::col_character(),to = readr::col_character()),
                   raw_counts = readr::cols(cell_id = readr::col_character(),.default = readr::col_double()))
  
  res <- purrr::map2(files, col_spec, ~ readr::read_csv(file = .x, col_types = .y, progress = FALSE)) %>%
    purrr::modify_at(c("counts","raw_counts"), ~ tibble::column_to_rownames(.,var = 'cell_id') %>% as.matrix() %>% as("sparseMatrix") ) %>%
    purrr::modify_at("progression", ~ dplyr::mutate(.,from = stringr::str_remove_all(from,"'"), to = stringr::str_remove_all(to,"'")) %>% dplyr::select(cell_id,percentage,from,to)) %>%
    purrr::modify_at("dim_red", ~ tibble::column_to_rownames(.,var = 'cell_id') %>% magrittr::set_colnames(paste0("comp_",1:ncol(.))) %>% as.matrix()) %>%
    purrr::modify_at('grouping', ~ dplyr::rename(.,group_id = label))
  
  res[["root"]] <- root
  return(res)
}
