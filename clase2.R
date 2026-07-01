install.packages('terra')
install.packages('sf')
install.packages('gstat')
install.packages('dplyr')
install.packages('lwgeom')
library(terra)
library(sf)
library(gstat)
library(dplyr)
library(lwgeom)

# 1. CARGAR DATOS
coleccion <- svc("mapa_rinde.shp")
print(coleccion)

geom_types <- sapply(1:length(coleccion), function(i) geomtype(coleccion[i]))
print(geom_types)

point_index <- which(geom_types %in% c("points", "point", "POINT", "POINTS"))[1]
if (is.na(point_index)) {
  stop("No se encontró una capa de puntos en 'mapa_rinde.shp'. Revisar geometrías mixtas con svc().")
}

puntos_raw <- coleccion[point_index]
lote_raw   <- vect("lote.shp")

# 2. EXTRACCIÓN MANUAL Y SINCRONIZACIÓN
crds_mat <- crds(puntos_raw)
attr_df  <- as.data.frame(puntos_raw)

n_filas <- min(nrow(crds_mat), nrow(attr_df))

df_base <- data.frame(
  x = crds_mat[1:n_filas, 1],
  y = crds_mat[1:n_filas, 2],
  yield = as.numeric(attr_df[1:n_filas, "VRYIELDMAS"])
) #Verificar el nombre de la variable (atributo) del rinde


# 3. LIMPIEZA DE OUTLIERS (Balzarini)
media  <- mean(df_base$yield, na.rm = TRUE)
desvio <- sd(df_base$yield, na.rm = TRUE)

df_limpio <- df_base %>%
  filter(!is.na(yield), 
         yield > 0, 
         yield > (media - 2.5*desvio), 
         yield < (media + 2.5*desvio))

# ---------------------------------------------------------
# BOXPLOT COMPARATIVO (ANTES Y DESPUÉS)
# ---------------------------------------------------------
par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))

boxplot(df_base$yield, 
        main = "Rinde Original\n(Con Outliers)", 
        ylab = "Rinde (unidades)", 
        col = "salmon",
        cex.main = 0.9)

boxplot(df_limpio$yield, 
        main = "Rinde Limpio\n(Criterio Balzarini)", 
        ylab = "Rinde (unidades)", 
        col = "lightblue",
        cex.main = 0.9)

par(mfrow = c(1, 1))
# ---------------------------------------------------------

# 4. CONVERSIÓN A SF Y PROYECCIÓN
puntos_sf <- st_as_sf(df_limpio, coords = c("x", "y"), crs = crs(puntos_raw)) %>% 
  st_transform(32720)

lote_sf <- st_as_sf(lote_raw) %>% st_transform(32720)

# 5. VARIOGRAMA
v_exp <- variogram(yield ~ 1, puntos_sf)
v_fit <- fit.variogram(v_exp, model = vgm(psill = var(puntos_sf$yield), "Sph", range = 150, nugget = 0))

plot(v_exp, model = v_fit, 
     main = "Ajuste del Semivariograma (Modelo Esférico)",
     xlab = "Distancia entre puntos (metros)", 
     ylab = "Semivarianza",
     col = "blue", lwd = 2, pch = 16)

# 6. INTERPOLACIÓN (KRIGING)
grid <- st_make_grid(lote_sf, cellsize = 10, what = "centers") %>%
  st_as_sf() %>% 
  st_filter(lote_sf)

kriging_res <- krige(yield ~ 1, puntos_sf, grid, model = v_fit)
dev.off() 
par(mfrow = c(1, 1))

# 7. RASTERIZACIÓN Y RESULTADO FINAL
res_df <- data.frame(st_coordinates(kriging_res), z = kriging_res$var1.pred)
rinde_raster_m <- rast(res_df, type="xyz", crs="EPSG:32720")

rinde_final_wgs84 <- project(rinde_raster_m, "EPSG:4326")
lote_wgs84 <- project(lote_raw, "EPSG:4326")
rinde_final_wgs84 <- mask(rinde_final_wgs84, lote_wgs84)

# VISUALIZACIÓN DEL MAPA FINAL
colores_rinde <- colorRampPalette(c("red", "yellow", "darkgreen"))(100)

plot(rinde_final_wgs84, 
     col = colores_rinde, 
     main = "Mapa de Rendimiento Final (Interpolado)",
     axes = TRUE)
plot(lote_wgs84, add = TRUE, border = "black", lwd = 2)

writeRaster(rinde_final_wgs84, "rinde_final_ok.tif", overwrite = TRUE)

# ============================================================
# AMBIENTACIÓN: RASTER DE RINDE → POLÍGONOS CONTIGUOS SIN HUECOS
# ============================================================

# --- Paso 1: Definir número de ambientes ---
n_ambientes <- 3  # Podés cambiar a 2, 4, 5, etc.

# --- Paso 2: Clasificar el raster en N zonas por cuantiles ---
valores <- values(rinde_final_wgs84, na.rm = TRUE)
cortes  <- quantile(valores, probs = seq(0, 1, length.out = n_ambientes + 1))

rcl_mat <- cbind(
  cortes[1:n_ambientes],
  cortes[2:(n_ambientes + 1)],
  1:n_ambientes
)

rinde_class <- classify(rinde_final_wgs84, rcl_mat, include.lowest = TRUE)

# --- Paso 3: Suavizar para evitar polígonos fragmentados ---
rinde_suave <- focal(rinde_class, w = 5, fun = "modal", na.rm = TRUE)
rinde_suave <- mask(rinde_suave, lote_wgs84)

# --- Paso 4: Convertir a polígonos y disolver por zona ---
ambientes_vect <- as.polygons(rinde_suave, dissolve = TRUE)
ambientes_sf <- st_as_sf(ambientes_vect)
names(ambientes_sf)[1] <- "ambiente"

# --- Paso 5a: Asegurar contigüidad total (sin huecos) ---
sf_use_s2(FALSE)

ambientes_sf <- st_intersection(ambientes_sf, st_as_sf(lote_wgs84))
ambientes_sf <- st_make_valid(ambientes_sf)
ambientes_sf <- st_buffer(ambientes_sf, dist = 0)

cobertura   <- st_union(ambientes_sf)
cobertura   <- st_make_valid(cobertura)
lote_sf_wgs <- st_as_sf(lote_wgs84)
lote_sf_wgs <- st_make_valid(lote_sf_wgs)
huecos      <- st_difference(lote_sf_wgs, cobertura)

if (nrow(huecos) > 0 & sum(st_area(huecos)) > units::set_units(0, m^2)) {
  for (i in 1:nrow(huecos)) {
    dists <- st_distance(huecos[i, ], ambientes_sf)
    idx_cercano <- which.min(dists)
    ambientes_sf$geometry[idx_cercano] <- st_union(
      ambientes_sf$geometry[idx_cercano], 
      huecos$geometry[i]
    )
  }
}

ambientes_sf <- st_make_valid(ambientes_sf)
sf_use_s2(TRUE)

# --- Paso 5b: Eliminar polígonos menores a 1 ha (versión rápida) ---
ambientes_sf <- st_cast(ambientes_sf, "POLYGON")
ambientes_sf$area_m2 <- as.numeric(st_area(ambientes_sf))
umbral_min <- 10000  # 1 ha

rinde_suave2 <- focal(rinde_class, w = 11, fun = "modal", na.rm = TRUE)
rinde_suave2 <- mask(rinde_suave2, lote_wgs84)

ambientes_vect <- as.polygons(rinde_suave2, dissolve = TRUE)
ambientes_sf <- st_as_sf(ambientes_vect)
names(ambientes_sf)[1] <- "ambiente"

sf_use_s2(FALSE)
ambientes_sf <- st_intersection(ambientes_sf, st_as_sf(lote_wgs84))
ambientes_sf <- st_make_valid(ambientes_sf)
ambientes_sf <- st_buffer(ambientes_sf, dist = 0)

cobertura   <- st_union(ambientes_sf)
cobertura   <- st_make_valid(cobertura)
lote_sf_wgs <- st_make_valid(st_as_sf(lote_wgs84))
huecos      <- st_difference(lote_sf_wgs, cobertura)

if (nrow(huecos) > 0 & sum(st_area(huecos)) > units::set_units(0, m^2)) {
  for (i in 1:nrow(huecos)) {
    dists <- st_distance(huecos[i, ], ambientes_sf)
    idx_cercano <- which.min(dists)
    ambientes_sf$geometry[idx_cercano] <- st_union(
      ambientes_sf$geometry[idx_cercano], 
      huecos$geometry[i]
    )
  }
}

ambientes_sf <- st_make_valid(ambientes_sf)
sf_use_s2(TRUE)

cat("✔ Polígonos suavizados con ventana focal 11x11\n")

# --- Paso 6: Etiquetar ambientes ---
if (n_ambientes == 2) {
  etiquetas <- c("Bajo", "Alto")
} else if (n_ambientes == 3) {
  etiquetas <- c("Bajo", "Medio", "Alto")
} else if (n_ambientes == 4) {
  etiquetas <- c("Bajo", "Medio-Bajo", "Medio-Alto", "Alto")
} else if (n_ambientes == 5) {
  etiquetas <- c("Muy Bajo", "Bajo", "Medio", "Alto", "Muy Alto")
} else {
  etiquetas <- paste("Zona", 1:n_ambientes)
}

ambientes_sf$nombre <- factor(
  ambientes_sf$ambiente,
  levels = 1:n_ambientes,
  labels = etiquetas
)

graphics.off()

ambientes_plot <- st_simplify(ambientes_sf, dTolerance = 0.0002)
ambientes_plot <- st_make_valid(ambientes_plot)

colores_amb <- colorRampPalette(c("red", "gold", "darkgreen"))(n_ambientes)

plot(ambientes_plot["nombre"], 
     col = colores_amb, 
     main = paste("Ambientación -", n_ambientes, "zonas"),
     border = "grey30", lwd = 0.5,
     reset = FALSE)
plot(st_geometry(st_as_sf(lote_wgs84)), add = TRUE, border = "black", lwd = 2)
legend("bottomright", legend = etiquetas, 
       fill = colores_amb, title = "Ambiente")

# --- Paso 7: Visualización ---
colores_amb <- colorRampPalette(c("red", "gold", "darkgreen"))(n_ambientes)

plot(ambientes_sf["nombre"], 
     col = colores_amb, 
     main = paste("Ambientación -", n_ambientes, "zonas"),
     border = "grey30", lwd = 0.5)
plot(st_geometry(st_as_sf(lote_wgs84)), add = TRUE, border = "black", lwd = 2)
legend("bottomright", legend = etiquetas, 
       fill = colores_amb, title = "Ambiente")

# --- Paso 8: Exportar shapefile ---
st_write(ambientes_sf, "ambientes_rinde1.shp", delete_dsn = TRUE)

cat("\n✔ Shapefile exportado: ambientes_rinde1.shp\n")
cat("  Ambientes:", n_ambientes, "\n")
cat("  Áreas por zona:\n")
ambientes_sf$area_ha <- as.numeric(st_area(ambientes_sf)) / 10000
print(ambientes_sf[, c("nombre", "area_ha")] |> st_drop_geometry())

# 8. ANÁLISIS BIVARIADO NDVI
s2_img <- rast("mosaico_lote.tif")
ndvi_raw <- (s2_img[["B8"]] - s2_img[["B4"]]) / (s2_img[["B8"]] + s2_img[["B4"]])

ndvi_resampled <- project(ndvi_raw, rinde_final_wgs84, method = "bilinear")
ndvi_lote <- mask(ndvi_resampled, lote_wgs84)

comparativa_stack <- c(ndvi_lote, rinde_final_wgs84)
names(comparativa_stack) <- c("ndvi", "rinde")
df_correlacion <- as.data.frame(comparativa_stack, na.rm = TRUE)
df_correlacion <- df_correlacion[df_correlacion$rinde > 0, ]

plot(df_correlacion$ndvi, df_correlacion$rinde, 
     pch = 16, cex = 0.6, col = rgb(0.1, 0.4, 0.1, 0.4),
     main = "Correlación Bivariada: NDVI vs RINDE",
     xlab = "Vigor Vegetativo (NDVI Feb)", 
     ylab = "Rendimiento (Marzo)")

modelo_lineal <- lm(rinde ~ ndvi, data = df_correlacion)
abline(modelo_lineal, col = "red", lwd = 3)
r_cuadrado <- summary(modelo_lineal)$r.squared
text(x = min(df_correlacion$ndvi), y = max(df_correlacion$rinde), 
     labels = paste("R² =", round(r_cuadrado, 3)), 
     pos = 4, font = 2, col = "red", cex = 1.2)

# 11. GNDVI
gndvi_raw <- (s2_img[["B8"]] - s2_img[["B3"]]) / (s2_img[["B8"]] + s2_img[["B3"]])
gndvi_resampled <- project(gndvi_raw, rinde_final_wgs84, method = "bilinear")
gndvi_lote <- mask(gndvi_resampled, lote_wgs84)

stack_gndvi <- c(gndvi_lote, rinde_final_wgs84)
names(stack_gndvi) <- c("gndvi", "rinde")
df_gndvi <- as.data.frame(stack_gndvi, na.rm = TRUE)
df_gndvi <- df_gndvi[df_gndvi$rinde > 0, ]

plot(df_gndvi$gndvi, df_gndvi$rinde, 
     pch = 16, cex = 0.6, col = rgb(0.1, 0.1, 0.5, 0.4),
     main = "Correlación Bivariada: GNDVI vs RINDE",
     xlab = "Vigor Vegetativo (GNDVI Feb)", 
     ylab = "Rendimiento (Marzo)")

modelo_gndvi <- lm(rinde ~ gndvi, data = df_gndvi)
abline(modelo_gndvi, col = "darkblue", lwd = 3)
r2_gndvi <- summary(modelo_gndvi)$r.squared
text(x = min(df_gndvi$gndvi), y = max(df_gndvi$rinde), 
     labels = paste("R² GNDVI =", round(r2_gndvi, 3)), 
     pos = 4, font = 2, col = "darkblue", cex = 1.2)


