#' Create an hmm object
#'
#' @param object A list provided by the fhmm function
#' @param chains number of chains
#' @return An ndp object
#'
hmm <- function(object,chains){


    d <- object$d
    offset <- 5
    if_df = purrr::map_dfr(1:chains,function(y){
                                purrr::map_dfr(1:object$K,function(x){
                                 dplyr::as_tibble(object[[offset+y]]$cluster_intensity[,((x-1)*length(d)+x):(x*length(d)),drop=F ]) %>%
                                 dplyr::mutate(Chain = y, sample_ix = 1:dplyr::n()) %>%
                                 tidyr::gather(dplyr::contains("V"),
                                               key = "column_names",value = "Density") %>%
                                 dplyr::left_join(dplyr::tibble(column_names = paste0("V",1:length(d)),
                                                                Distance = d),
                                                  by = "column_names") %>%
                                 dplyr::select(-column_names) %>%
                                 dplyr::mutate(Intensity_Function = x,
                                               Chain = y)
                                }) } )

	gd_df <- purrr::map_dfr(1:chains,function(y){
								   dplyr::as_tibble(object[[offset+y]]$global_intensity) %>%
								dplyr::mutate(Chain = y, sample_ix = 1 :dplyr::n()) %>%
								tidyr::gather(dplyr::contains("V"),
											  key = "column_names",value = "Global_Density") %>%
								dplyr::left_join(dplyr::tibble(column_names = paste0("V",1:length(d)),
																  Distance = d),
								by = "column_names") %>%
								dplyr::select(-column_names)
								})

    vls <- expand.grid(l = 1:object$L,k = 1:object$K)
    nms <- paste0("k: ", vls$k,", l: ", vls$l)
    mus <- purrr::map(1:chains,function(x){
        mus <- object[[offset+x]]$mu
        colnames(mus) <- paste0("mu ",nms)
        coda::as.mcmc(mus)
        })
    ws <- purrr::map(1:chains,function(x){
        ws <- object[[offset+x]]$w
        colnames(ws) <- paste0("w ", nms)
        coda::as.mcmc(ws)
    })
    pis <- purrr::map(1:chains,function(x){
        pis <- object[[offset+x]]$pi
        colnames(pis) <- paste0("pi ", "k: ",1:object$K)
        coda::as.mcmc(pis)
    })
    pmat <- Reduce("+",lapply(1:chains,function(x) object[[offset+x]]$Cluster_Matrix ))

    num_clusters <- purrr::map(1:chains,function(x) coda::as.mcmc(apply(object[[offset+x]]$Cluster_Assignment,1,function(z) c("Num_Clusters" = length(unique(z))))))

    cluster_assignment  <- purrr::map(1:chains,function(x) object[[offset+x]]$Cluster_Assignment)

    if(ncol(object[[offset+1]]$tau)>1)
        tau_nms <- paste0("tau ",nms)
    else
        tau_nms <- c("tau")

    tau <- purrr::map(1:chains,function(x){
                          taus <- object[[offset+x]]$tau
                          colnames(taus) <- tau_nms
                          return(coda::as.mcmc(taus))})

    out <- list(call = object$call,
                n = object$n,
                J = nrow(object$X),
                if_df = if_df,
				global_density = gd_df,
                pmat = pmat / chains,
				cluster_assignment = cluster_assignment,
                num_clusters = coda::as.mcmc.list(num_clusters),
                pi = coda::as.mcmc.list(pis),
                w = coda::as.mcmc.list(ws),
                mu = coda::as.mcmc.list(mus),
                tau = coda::as.mcmc.list(tau)
    )


    structure(out, class = c("hmm"))
}
