# ---------------------------------------------------------------
# Muestreo estratificado (1..4 clases) -> KML con estilos (My Maps)
# Requisitos: sf, dplyr, xml2
# ---------------------------------------------------------------

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(xml2)
})

CLASS_ORDER  <- c("CONS","DOSE","EDIF","SUEL")
HEX_COLORES  <- c("DOSE"="#428b07", "SUEL"="#dde78e", "EDIF"="#4a4a4a", "CONS"="#d2d2d2")

# -------- utilidades KML --------
hex_to_kml <- function(hex_rgb, alpha = 0.6) {
  hex_rgb <- gsub("#","",hex_rgb)
  r <- as.integer(strtoi(substr(hex_rgb,1,2),16L))
  g <- as.integer(strtoi(substr(hex_rgb,3,4),16L))
  b <- as.integer(strtoi(substr(hex_rgb,5,6),16L))
  a <- as.integer(round(alpha*255))
  sprintf("%02X%02X%02X%02X", a, b, g, r)  # aabbggrr
}
close_ring <- function(m){ if(!all(m[1,1:2]==m[nrow(m),1:2])) m <- rbind(m,m[1,]); m }
add_linear_ring <- function(parent, ring_mat){
  ring_mat <- close_ring(ring_mat)
  coords <- paste0(apply(ring_mat[,1:2,drop=FALSE],1,function(x) sprintf("%.*f,%.*f,0",7,x[1],7,x[2])), collapse=" ")
  lr <- xml_add_child(parent,"LinearRing"); xml_add_child(lr,"coordinates",coords)
}
add_sfg_polygon <- function(parent, sfg_geom){
  cls <- class(sfg_geom)
  if ("MULTIPOLYGON" %in% cls) {
    for (pg in sfg_geom) add_sfg_polygon(parent, structure(pg, class=c("POLYGON","sfg")))
  } else if ("POLYGON" %in% cls) {
    poly <- xml_add_child(parent,"Polygon"); xml_add_child(poly,"tessellate","1")
    outer <- xml_add_child(poly,"outerBoundaryIs"); add_linear_ring(outer, sfg_geom[[1]])
    if (length(sfg_geom)>1) for (i in 2:length(sfg_geom)) { inner <- xml_add_child(poly,"innerBoundaryIs"); add_linear_ring(inner, sfg_geom[[i]]) }
  } else stop("Geometría no soportada para KML: ", paste(cls, collapse="/"))
}

# -------- reparto proporcional (Hamilton) con tope por disponibilidad --------
allocate_counts <- function(avail_named_int, n_total) {
  stopifnot(!is.null(names(avail_named_int)))
  avail <- pmax(0L, as.integer(avail_named_int))
  nm    <- as.character(names(avail)); names(avail) <- nm
  
  if (length(avail)==1L) {
    out <- min(n_total, avail[1]); names(out) <- nm; return(out)
  }
  
  total_av <- sum(avail)
  if (total_av == 0L) return(setNames(integer(length(avail)), nm))
  
  targets <- n_total * (avail / total_av)         # objetivos teóricos
  base    <- pmin(floor(targets), avail)          # suelo + tope
  left    <- n_total - sum(base)
  
  if (left > 0) {
    rema <- targets - floor(targets)             # restos
    order_idx <- order(rema, decreasing = TRUE)
    i <- 1L
    while (left > 0 && any(base < avail)) {
      k <- order_idx[i]
      if (base[k] < avail[k]) { base[k] <- base[k] + 1L; left <- left - 1L }
      i <- i + 1L
      if (i > length(order_idx)) i <- 1L
    }
  }
  setNames(as.integer(base), nm)
}

# -------- función principal --------
sample_hex_to_kml <- function(resol = 11, n = 100, clases = c(1,2,3,4),
                              seed = NULL, alpha_fill = 0.6,
                              alpha_line = 1.0, line_width = 1.2) {
  if (!resol %in% c(11,12)) stop("resol debe ser 11 o 12.")
  
  # normalizar clases pedidas
  if (is.numeric(clases)) {
    if (!all(clases %in% 1:4)) stop("clases numéricas deben estar en 1..4.")
    clases <- CLASS_ORDER[clases]
  } else {
    clases <- toupper(trimws(as.character(clases)))
    if (!all(clases %in% CLASS_ORDER)) stop('clases debe contener solo: "CONS","DOSE","EDIF","SUEL".')
  }
  clases <- unique(clases)
  
  in_file <- sprintf("h3-res-%d-porcentajes-cobertura-clasificado.gpkg", resol)
  if (!file.exists(in_file)) stop("No existe el archivo: ", in_file)
  
  sf_in <- sf::st_read(in_file, quiet = TRUE)
  if (!("clase" %in% names(sf_in))) stop('El archivo no contiene la columna "clase".')
  sf_in$clase <- toupper(trimws(as.character(sf_in$clase)))
  
  # disponibilidad (solo de clases pedidas)
  sf_sel <- dplyr::filter(sf_in, .data$clase %in% clases)
  if (nrow(sf_sel) == 0) {
    dist_all <- sf_in |> st_drop_geometry() |> count(clase, name="n") |> arrange(desc(n))
    stop(paste0("No hay hexágonos en las clases seleccionadas: ",
                paste(clases, collapse=", "),
                "\nDistribución detectada: ",
                paste(paste0(dist_all$clase,"=",dist_all$n), collapse=", ")))
  }
  
  avail_tbl <- sf_sel |> st_drop_geometry() |> count(clase, name="n_avail")
  avail <- setNames(rep(0L, length(clases)), clases)
  if (nrow(avail_tbl) > 0) {
    mt <- match(names(avail), avail_tbl$clase)
    avail[!is.na(mt)] <- as.integer(avail_tbl$n_avail[mt[!is.na(mt)]])
  }
  
  message("Disponibilidad (clases pedidas): ", paste(paste0(names(avail), "=", avail), collapse=", "))
  
  if (!is.null(seed)) set.seed(seed)
  
  
  # --- asignaciones (exactas por proporción si hay disponibilidad) ---
  if (length(clases) == 1L) {
    n_assign <- setNames(min(n, avail[1]), clases)
    
  } else {
    # proporciones entre las clases pedidas
    tot_av <- sum(avail)
    props  <- if (tot_av > 0) (avail / tot_av) else rep(0, length(avail))
    # cuotas teóricas
    targets <- n * props
    base    <- floor(targets)
    left    <- n - sum(base)
    
    # repartir los restos (Hamilton) para cerrar exactamente en n
    if (left > 0) {
      rema <- targets - base
      order_idx <- order(rema, decreasing = TRUE)
      i <- 1L
      while (left > 0) {
        k <- order_idx[i]
        base[k] <- base[k] + 1L
        left <- left - 1L
        i <- i + 1L
        if (i > length(order_idx)) i <- 1L
      }
    }
    quotas <- setNames(as.integer(base), names(avail))
    
    # si todas las cuotas caben en la disponibilidad, usar EXACTAMENTE esas cuotas
    if (all(quotas <= avail)) {
      n_assign <- quotas
    } else {
      # de lo contrario, usar el reparto con tope por disponibilidad (ya implementado)
      n_assign <- allocate_counts(avail, n)
    }
    
    # reindexar exactamente al orden de 'clases' y rellenar con 0
    n_assign <- setNames(n_assign[names(avail)], names(avail))
    n_assign[is.na(n_assign)] <- 0L
    
    if (sum(n_assign) == 0L) {
      dist_all <- sf_in |> st_drop_geometry() |> count(clase, name="n") |> arrange(desc(n))
      stop(paste0(
        "No hay disponibilidad en las clases seleccionadas.\n",
        "Distribución detectada: ",
        paste(paste0(dist_all$clase,"=",dist_all$n), collapse=", ")
      ))
    }
    if (sum(n_assign) < n) message("Aviso: asignados ", sum(n_assign), " de ", n, " (límite por disponibilidad).")
    
    # (opcional) imprime objetivos teóricos y cuotas
    message("Objetivos teóricos: ", paste(paste0(names(avail), "=", sprintf("%.2f", targets)), collapse=", "))
    message("Cuotas exactas (Hamilton): ", paste(paste0(names(quotas), "=", quotas), collapse=", "))
  }
  
  # # asignaciones
  # if (length(clases) == 1L) {
  #   n_assign <- setNames(min(n, avail[1]), clases)
  # } else {
  #   n_assign <- allocate_counts(avail, n)
  #   # reindex estricto + rellenar con 0
  #   n_assign <- setNames(n_assign[names(avail)], names(avail))
  #   n_assign[is.na(n_assign)] <- 0L
  #   # fallback seguro si por algún motivo quedó en 0 pero hay disponibilidad
  #   if (sum(n_assign) == 0L && sum(avail) > 0L) {
  #     message("Asignación proporcional resultó 0; aplico fallback seguro.")
  #     # reparte 1 en ciclo por orden de mayor disponibilidad hasta min(n, sum(avail))
  #     n_target <- min(n, sum(avail))
  #     n_assign[] <- 0L
  #     order_idx <- order(avail, decreasing = TRUE)
  #     i <- 1L; left <- n_target
  #     while (left > 0) {
  #       k <- order_idx[i]
  #       if (n_assign[k] < avail[k]) { n_assign[k] <- n_assign[k] + 1L; left <- left - 1L }
  #       i <- i + 1L; if (i > length(order_idx)) i <- 1L
  #       if (all(n_assign >= avail)) break
  #     }
  #   }
  #   if (sum(n_assign) == 0L) {
  #     dist_all <- sf_in |> st_drop_geometry() |> count(clase, name="n") |> arrange(desc(n))
  #     stop(paste0("No hay disponibilidad en las clases seleccionadas.\nDistribución detectada: ",
  #                 paste(paste0(dist_all$clase,"=",dist_all$n), collapse=", ")))
  #   }
  #   if (sum(n_assign) < n) message("Aviso: asignados ", sum(n_assign), " de ", n, " (límite por disponibilidad).")
  # }
  # 
  # # imprimir objetivos teóricos para transparencia
  # tot_av <- sum(avail)
  # targets <- if (tot_av>0) n * (avail / tot_av) else rep(0, length(avail))
  # message("Objetivos teóricos: ", paste(paste0(names(avail), "=", sprintf("%.2f", targets)), collapse=", "))
  # message("Asignación final: ",  paste(paste0(names(n_assign), "=", n_assign), collapse=", "))
  
  # muestreo por clase
  parts <- list()
  for (cl in clases) {
    k <- unname(n_assign[cl]); if (is.na(k) || k <= 0) next
    pool <- dplyr::filter(sf_sel, .data$clase == cl)
    k_eff <- min(k, nrow(pool))
    idxs <- if (k_eff > 0) sample.int(nrow(pool), size = k_eff, replace = FALSE) else integer(0)
    parts[[cl]] <- pool[idxs, , drop = FALSE]
  }
  if (length(parts) == 0) stop("No se seleccionó ningún hexágono (asignaciones = 0).")
  
  sf_sample <- do.call(rbind, parts)
  
  # geometrías seguras para KML
  sf_sample <- sf::st_make_valid(sf_sample)
  sf_sample <- sf::st_collection_extract(sf_sample, "POLYGON", warn = FALSE)
  sf_sample <- sf_sample[!sf::st_is_empty(sf_sample), , drop = FALSE]
  sf_sample <- suppressWarnings(sf::st_cast(sf_sample, "MULTIPOLYGON", warn = FALSE))
  sf_ll     <- sf::st_transform(sf_sample, 4326)
  
  # ---- KML (carpetas por clase, estilos) ----
  out_kml <- sprintf("muestra-h3-res-%d-%s-n%d.kml", resol, paste(clases, collapse="-"), sum(n_assign))
  kml  <- xml_new_root("kml", xmlns = "http://www.opengis.net/kml/2.2")
  doc  <- xml_add_child(kml, "Document")
  xml_add_child(doc, "name", sprintf("Muestra estratificada (res %d, n=%d)", resol, sum(n_assign)))
  
  darker <- function(hex, f = 0.75) {
    h <- gsub("#","",hex)
    r <- max(0, as.integer(strtoi(substr(h,1,2),16L)) * f)
    g <- max(0, as.integer(strtoi(substr(h,3,4),16L)) * f)
    b <- max(0, as.integer(strtoi(substr(h,5,6),16L)) * f)
    sprintf("#%02X%02X%02X", as.integer(r), as.integer(g), as.integer(b))
  }
  for (cl in clases) {
    style_id <- paste0("style-", cl)
    styn <- xml_add_child(doc, "Style", id = style_id)
    ls <- xml_add_child(styn, "LineStyle")
    xml_add_child(ls, "color", hex_to_kml(darker(HEX_COLORES[cl]), alpha = alpha_line))
    xml_add_child(ls, "width", format(line_width, trim = TRUE))
    ps <- xml_add_child(styn, "PolyStyle")
    xml_add_child(ps, "color", hex_to_kml(HEX_COLORES[cl], alpha = alpha_fill))
    xml_add_child(ps, "fill", "1"); xml_add_child(ps, "outline", "1")
  }
  
  g_all <- st_geometry(sf_ll)
  has_hex_id <- "hex_id" %in% names(sf_ll)
  
  for (cl in clases) {
    folder <- xml_add_child(doc, "Folder")
    xml_add_child(folder, "name", paste0(cl, " (n=", unname(n_assign[cl]), ")"))
    style_url <- paste0("#style-", cl)
    
    rows <- which(sf_ll$clase == cl)
    for (i in rows) {
      pm <- xml_add_child(folder, "Placemark")
      nm <- if (has_hex_id) as.character(sf_ll$hex_id[i]) else paste0(cl, "_", i)
      xml_add_child(pm, "name", nm)
      xml_add_child(pm, "styleUrl", style_url)
      
      ed <- xml_add_child(pm, "ExtendedData")
      add_data <- function(name, value){ d <- xml_add_child(ed, "Data", name = name); xml_add_child(d, "value", as.character(value)) }
      add_data("clase", cl)
      for (p in c("porc_CONS","porc_DOSE","porc_EDIF","porc_SUEL")) if (p %in% names(sf_ll)) add_data(p, round(sf_ll[[p]][i], 6))
      
      add_sfg_polygon(pm, g_all[[i]])
    }
  }
  
  xml2::write_xml(kml, file = out_kml)
  message("Escrito: ", out_kml)
  invisible(sf_ll)
}

# ---- ejemplos de uso ----
# Todas las clases (proporcional a su disponibilidad real):
# sample_hex_to_kml(resol = 11, n = 100, clases = c(2,1,3,4), seed = 2025)

# Dos clases:
# sample_hex_to_kml(resol = 11, n = 40, clases = c("EDIF","SUEL"), seed = 1)

# Una clase (n directo):
# sample_hex_to_kml(resol = 12, n = 25, clases = "DOSE", seed = 7)