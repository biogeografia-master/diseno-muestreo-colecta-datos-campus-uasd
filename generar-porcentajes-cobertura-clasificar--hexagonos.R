# Clasificación de hexágonos por cobertura con filtro: 100% dentro del polígono UASD
# Entradas:
#   - "h3-res-11.gpkg", "h3-res-12.gpkg"
#   - "tipos-cobertura.gpkg" (valores CONS/DOSE/EDIF/SUEL)
#   - "poligono-uasd-32619.gpkg" (polígono de recorte)
# Salidas:
#   - "h3-res-11-porcentajes-cobertura.gpkg", "h3-res-12-porcentajes-cobertura.gpkg"

suppressPackageStartupMessages({
  library(sf)
  library(dplyr)
  library(tidyr)
  library(units)
  library(purrr)
  library(rlang)
})

# ----------------------------- #
# Configuración
# ----------------------------- #
hex_files   <- c("h3-res-11.gpkg", "h3-res-12.gpkg")
cover_file  <- "tipos-cobertura.gpkg"
uasd_file   <- "poligono-uasd-32619.gpkg"
out_files   <- c("h3-res-11-clasificado.gpkg", "h3-res-12-clasificado.gpkg")

cats        <- c("CONS", "DOSE", "EDIF", "SUEL")
THRESH_FULL <- 0.999  # usa 1 si quieres 100% estricto

# ----------------------------- #
# Helpers
# ----------------------------- #

detect_tipo_col <- function(df, categories = cats) {
  for (nm in names(df)) {
    vals <- unique(df[[nm]]) |> as.character()
    if (sum(categories %in% vals) >= 2) return(nm)
  }
  stop("No se encontró columna con valores CONS/DOSE/EDIF/SUEL en 'tipos-cobertura.gpkg'.")
}

choose_projected_crs <- function(g) {
  if (sf::st_is_longlat(g)) {
    message("Geometrías en lon/lat; usando EPSG:32619 para cálculo de áreas.")
    sf::st_crs(32619)
  } else {
    sf::st_crs(g)
  }
}

# Renormaliza los 4 porcentajes a 100.000000 exactos (con ajuste mínimo en la última col)
renormalizar_porcentajes <- function(out, categories = cats, digits = 6) {
  pcols <- paste0("porc_", categories)
  if (!"porc_residuo_hex" %in% names(out)) out$porc_residuo_hex <- 0
  
  df <- st_drop_geometry(out)
  for (nm in pcols) if (!(nm %in% names(df))) df[[nm]] <- 0
  df[pcols] <- lapply(df[pcols], function(x) as.numeric(x))
  
  suma <- rowSums(df[, pcols, drop = FALSE], na.rm = TRUE)
  mult <- ifelse(suma > 0, 100 / suma, 0)
  
  scaled <- sweep(df[, pcols, drop = FALSE], 1, mult, `*`)
  scaled <- lapply(scaled, function(col) round(col, digits))
  scaled <- as.data.frame(scaled)
  
  suma2 <- rowSums(scaled[, pcols, drop = FALSE], na.rm = TRUE)
  err   <- round(100 - suma2, digits)
  last  <- pcols[length(pcols)]
  scaled[[last]] <- round(scaled[[last]] + err, digits)
  
  out[, pcols] <- scaled[pcols]
  out
}

# Calcula porcentajes por hexágono, filtra por 100% dentro UASD y clasifica por categoría dominante
classify_hex <- function(hex_path, cov_path, uasd_path, out_path,
                         categories = cats, full_thresh = THRESH_FULL) {
  
  message("Leyendo hexágonos: ", hex_path)
  hex <- sf::st_read(hex_path, quiet = TRUE) |> sf::st_make_valid()
  
  message("Leyendo coberturas: ", cov_path)
  cov <- sf::st_read(cov_path, quiet = TRUE) |> sf::st_make_valid()
  
  message("Leyendo polígono UASD: ", uasd_path)
  uasd <- sf::st_read(uasd_path, quiet = TRUE) |> sf::st_make_valid()
  uasd_union <- sf::st_union(uasd)  # por si viene multipart
  
  # Detecta columna de tipo
  tipo_col <- detect_tipo_col(cov, categories)
  
  # CRS proyectado
  target_crs <- choose_projected_crs(cov)
  hex  <- sf::st_transform(hex,  target_crs)
  cov  <- sf::st_transform(cov,  target_crs)
  uasd_union <- sf::st_transform(uasd_union, target_crs)
  
  # ID estable de hexágono
  if (!("hex_id" %in% names(hex))) hex$hex_id <- seq_len(nrow(hex))
  
  # --- Filtro: 100% (≈) del área del hex dentro del polígono UASD ---
  message("Filtrando hexágonos 100% dentro del polígono UASD (umbral = ", full_thresh, ")…")
  hex$area_hex <- units::drop_units(sf::st_area(hex))
  
  idx_uasd <- sf::st_intersects(hex, uasd_union)
  any_hits <- lengths(idx_uasd) > 0
  if (!any(any_hits)) stop("Ningún hexágono cae dentro del polígono UASD en: ", hex_path)
  
  hex_sub_for_mask <- hex[any_hits, c("hex_id", "area_hex")]
  inter_uasd <- suppressWarnings(sf::st_intersection(hex_sub_for_mask, uasd_union))
  
  # calcular área de intersección SIN usar '.' en mutate
  inter_uasd$area_uasd <- units::drop_units(sf::st_area(inter_uasd))
  inter_uasd <- inter_uasd |>
    st_drop_geometry() |>
    group_by(hex_id) |>
    summarise(area_uasd = sum(area_uasd, na.rm = TRUE), .groups = "drop")
  
  hex <- hex |>
    left_join(inter_uasd, by = "hex_id") |>
    mutate(area_uasd = ifelse(is.na(area_uasd), 0, area_uasd),
           frac_uasd = ifelse(area_hex > 0, area_uasd / area_hex, 0)) |>
    filter(frac_uasd >= full_thresh) |>
    select(-area_hex, -area_uasd, -frac_uasd)
  
  if (nrow(hex) == 0) stop("Tras el filtro 100% UASD no queda ningún hexágono en: ", hex_path)
  
  # --- Coberturas: solo CONS/DOSE/EDIF/SUEL disueltas ---
  cov4 <- cov |>
    filter(!!sym(tipo_col) %in% categories) |>
    group_by(!!sym(tipo_col)) |>
    summarise(.groups = "drop") |>
    st_make_valid()
  
  # --- Intersecciones con coberturas (solo hex filtrados) ---
  idx_cov <- sf::st_intersects(hex, cov4)
  any_cov <- lengths(idx_cov) > 0
  if (!any(any_cov)) {
    warning("Sin intersecciones hexágonos–coberturas (4 clases) tras el filtro UASD en: ", hex_path)
    out <- hex |>
      mutate(
        porc_CONS = 0, porc_DOSE = 0, porc_EDIF = 0, porc_SUEL = 0,
        porc_residuo_hex = 100,   # todo fuera de las 4 clases
        clase = NA_character_, porc_max = 0
      )
    sf::st_write(out, out_path, delete_dsn = TRUE, quiet = TRUE)
    message("Escrito: ", out_path)
    return(invisible(out))
  }
  
  hex_sub <- hex[any_cov, "hex_id"]
  cov_sub <- cov4[unique(unlist(idx_cov[any_cov])), ]
  
  message("Calculando intersecciones hexágonos–coberturas (4 clases, disueltas)…")
  inter <- suppressWarnings(
    sf::st_intersection(
      dplyr::select(hex_sub, hex_id),
      dplyr::select(cov_sub, !!sym(tipo_col))
    )
  )
  
  if (nrow(inter) == 0) {
    warning("st_intersection devolvió 0 filas (4 clases) en: ", hex_path)
    out <- hex |>
      mutate(
        porc_CONS = 0, porc_DOSE = 0, porc_EDIF = 0, porc_SUEL = 0,
        porc_residuo_hex = 100,
        clase = NA_character_, porc_max = 0
      )
    sf::st_write(out, out_path, delete_dsn = TRUE, quiet = TRUE)
    message("Escrito: ", out_path)
    return(invisible(out))
  }
  
  # Área por pieza de intersección — sin '.' en mutate
  inter$area <- units::drop_units(sf::st_area(inter))
  
  # Áreas por tipo y total conocido (de las 4 clases)
  area_by_type <- inter |>
    st_drop_geometry() |>
    group_by(hex_id, !!sym(tipo_col)) |>
    summarise(area_tipo = sum(area, na.rm = TRUE), .groups = "drop")
  
  area_known <- area_by_type |>
    group_by(hex_id) |>
    summarise(area_known = sum(area_tipo, na.rm = TRUE), .groups = "drop")
  
  # Porcentaje por clase NORMALIZADO a 100 sobre el área conocida (4 clases)
  pct_tbl <- area_by_type |>
    left_join(area_known, by = "hex_id") |>
    mutate(pct = ifelse(area_known > 0, 100 * area_tipo / area_known, 0),
           cat = as.character(!!sym(tipo_col))) |>
    select(hex_id, cat, pct) |>
    tidyr::pivot_wider(names_from = cat,
                       values_from = pct,
                       values_fill = 0,
                       names_prefix = "porc_")
  
  # Asegurar columnas faltantes
  for (nm in paste0("porc_", categories)) {
    if (!(nm %in% names(pct_tbl))) pct_tbl[[nm]] <- 0
  }
  
  # Unir a TODOS los hex filtrados por UASD
  out <- hex |>
    left_join(pct_tbl, by = "hex_id")
  
  # Rellenar NAs con 0 en los 4 porcentajes
  for (nm in paste0("porc_", categories)) {
    if (nm %in% names(out)) out[[nm]][is.na(out[[nm]])] <- 0
  }
  
  # --- Residuo respecto al hexágono completo (diagnóstico) ---
  out$area_hex <- units::drop_units(st_area(out))
  cov4_union <- st_union(cov4)
  hits_cov4 <- st_intersects(out, cov4_union)
  any_hits_cov4 <- lengths(hits_cov4) > 0
  area_cov4 <- numeric(nrow(out))
  if (any(any_hits_cov4)) {
    inter_cov4 <- suppressWarnings(st_intersection(out[any_hits_cov4, "hex_id"], cov4_union))
    inter_cov4$area_cov4 <- units::drop_units(st_area(inter_cov4))
    a_cov4 <- inter_cov4 |>
      st_drop_geometry() |>
      group_by(hex_id) |>
      summarise(area_cov4 = sum(area_cov4, na.rm = TRUE), .groups = "drop")
    area_cov4[match(a_cov4$hex_id, out$hex_id)] <- a_cov4$area_cov4
  }
  out$porc_residuo_hex <- ifelse(out$area_hex > 0, pmax(0, 100 * (out$area_hex - area_cov4) / out$area_hex), 0)
  out$area_hex <- NULL  # no exportar áreas por defecto
  
  # --- Re-normalizar exactamente a 100.000000 (opción B) ---
  out <- renormalizar_porcentajes(out, categories = categories, digits = 6)
  
  # Clasificación por dominante (entre las 4 clases)
  pcols <- paste0("porc_", categories)
  porc_mat <- st_drop_geometry(out)[, pcols, drop = FALSE] |> as.matrix()
  max_idx <- apply(porc_mat, 1, function(v) if (all(is.na(v)) || all(v == 0)) NA_integer_ else which.max(v))
  out$clase    <- ifelse(is.na(max_idx), NA_character_, categories[max_idx])
  out$porc_max <- ifelse(is.na(max_idx), 0, porc_mat[cbind(seq_len(nrow(porc_mat)), max_idx)])
  
  # Escribir
  sf::st_write(out, out_path, delete_dsn = TRUE, quiet = TRUE)
  message("Escrito: ", out_path)
  invisible(out)
}

# ----------------------------- #
# Ejecutar en RStudio
# ----------------------------- #
walk2(hex_files, out_files, ~ classify_hex(.x, cover_file, uasd_file, .y))
message("Listo.")
