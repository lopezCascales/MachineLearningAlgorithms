#############################################################################
# PIPELINE DE INTEGRACIÓN MULTIÓMICA - PROGRESIÓN EN ESCLEROSIS MÚLTIPLE
# Maria Teresa Lopez Cascales
# Este script:
#   1. SIMULA 5 archivos de input (clínico, genético, inmunológico,
#      radiológico y de supervivencia) con una estructura realista,
#      y los guarda como CSV en ./data_input/
#   2. Los IMPORTA como si fueran datos reales de tu cohorte
#   3. Preprocesa, integra con SNF, valida estabilidad de clusters,
#      valida clínicamente (KM, log-rank, Cox ajustado) y ajusta un
#      Random Survival Forest sobre los datos integrados
#
# NOTA: los datos son SIMULADOS con una señal artificial insertada a
# propósito (un "subtipo de alto riesgo" con perfil inmuno+genético+
# radiológico coherente) para que el pipeline tenga algo real que
# encontrar. Sustituye la sección 1 por tus propios archivos cuando
# tengas los datos reales - la sección 2 en adelante ya está pensada
# para leer CSVs con esa misma estructura de columnas.
#############################################################################

# ---------------------------------------------------------------------------
# 0. PAQUETES
# ---------------------------------------------------------------------------
paquetes_necesarios <- c(
  "survival", "survminer", "cluster", "fpc", "factoextra",
  "SNFtool", "randomForestSRC", "dplyr", "tidyr"
)

paquetes_faltantes <- paquetes_necesarios[!paquetes_necesarios %in% installed.packages()[, "Package"]]
if (length(paquetes_faltantes) > 0) {
  install.packages(paquetes_faltantes, repos = "https://cloud.r-project.org")
}

invisible(lapply(paquetes_necesarios, library, character.only = TRUE))

set.seed(2026)
dir.create("data_input", showWarnings = FALSE)
dir.create("resultados", showWarnings = FALSE)


#############################################################################
# 1. SIMULACIÓN DE LOS ARCHIVOS DE INPUT
#    (sustituye esto por read.csv() de tus archivos reales cuando los tengas)
#############################################################################

n_pacientes <- 220
id_paciente <- sprintf("MS%03d", 1:n_pacientes)

# --- Subtipo latente "verdadero" (NO estaría en tus datos reales; existe
#     aquí solo para inyectar una señal biológica coherente y que el
#     pipeline tenga estructura real que recuperar) ---
subtipo_latente <- sample(c("alto_riesgo", "bajo_riesgo"), n_pacientes,
                           replace = TRUE, prob = c(0.35, 0.65))

# ---------------------------------------------------------------------------
# 1a. DATOS CLÍNICOS  (clinical_data.csv)
# ---------------------------------------------------------------------------
clinico <- data.frame(
  patient_id          = id_paciente,
  edad                = round(rnorm(n_pacientes, mean = 38, sd = 9)),
  sexo                = sample(c("F", "M"), n_pacientes, replace = TRUE, prob = c(0.7, 0.3)),
  tiempo_evolucion_an  = round(rgamma(n_pacientes, shape = 2, scale = 3), 1),
  edss_basal          = round(pmin(pmax(rnorm(n_pacientes,
                          mean = ifelse(subtipo_latente == "alto_riesgo", 3.2, 2.0), sd = 1), 0), 8), 1),
  tratamiento          = sample(c("interferon", "fingolimod", "natalizumab", "ninguno"),
                                n_pacientes, replace = TRUE, prob = c(0.3, 0.25, 0.25, 0.2)),
  tipo_ms              = sample(c("RRMS", "SPMS", "PPMS"), n_pacientes,
                                replace = TRUE, prob = c(0.75, 0.15, 0.10))
)
# Introducimos algo de missing realista (no todo el mundo tiene EDSS actualizado)
clinico$edss_basal[sample(1:n_pacientes, 8)] <- NA

write.csv(clinico, "data_input/clinical_data.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# 1b. DATOS GENÉTICOS  (genetic_data.csv)
#     Dosis alélica (0/1/2) para un panel reducido de variantes de riesgo
#     conocidas en EM (ej. HLA-DRB1*15:01 es el locus de mayor efecto)
# ---------------------------------------------------------------------------
genetico <- data.frame(
  patient_id = id_paciente,
  HLA_DRB1_1501 = rbinom(n_pacientes, 2, prob = ifelse(subtipo_latente == "alto_riesgo", 0.55, 0.20)),
  IL7R_rs6897932 = rbinom(n_pacientes, 2, prob = 0.30),
  CD58_rs2300747 = rbinom(n_pacientes, 2, prob = 0.25),
  IL2RA_rs2104286 = rbinom(n_pacientes, 2, prob = ifelse(subtipo_latente == "alto_riesgo", 0.40, 0.22)),
  TNFRSF1A_rs1800693 = rbinom(n_pacientes, 2, prob = 0.28),
  EVI5_rs11808092 = rbinom(n_pacientes, 2, prob = 0.20)
)
write.csv(genetico, "data_input/genetic_data.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# 1c. DATOS INMUNOLÓGICOS  (immune_data.csv)
#     Biomarcadores continuos de sangre/LCR
# ---------------------------------------------------------------------------
efecto_riesgo <- ifelse(subtipo_latente == "alto_riesgo", 1, 0)

inmuno <- data.frame(
  patient_id       = id_paciente,
  NfL_pg_ml        = round(rgamma(n_pacientes, shape = 2, scale = 8 + efecto_riesgo * 10), 1),  # neurofilamento de cadena ligera
  CD4_CD8_ratio    = round(rnorm(n_pacientes, mean = 2.2 - efecto_riesgo * 0.5, sd = 0.6), 2),
  IL6_pg_ml        = round(rgamma(n_pacientes, shape = 2, scale = 3 + efecto_riesgo * 2), 2),
  IgG_index        = round(rnorm(n_pacientes, mean = 0.7 + efecto_riesgo * 0.3, sd = 0.25), 2),
  bandas_oligoclonales = rbinom(n_pacientes, 1, prob = 0.85),
  CXCL13_pg_ml     = round(rgamma(n_pacientes, shape = 2, scale = 15 + efecto_riesgo * 20), 1)
)
# missing en un biomarcador menos rutinario (CXCL13, no siempre se pide)
inmuno$CXCL13_pg_ml[sample(1:n_pacientes, 25)] <- NA

write.csv(inmuno, "data_input/immune_data.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# 1d. DATOS RADIOLÓGICOS  (radiology_data.csv)
#     Features ya extraídas de RM (volumetría/conteo de lesiones), es decir,
#     el tipo de output que te daría un pipeline de segmentación por CNN/U-Net
# ---------------------------------------------------------------------------
radiologia <- data.frame(
  patient_id            = id_paciente,
  volumen_lesion_T2_ml  = round(rgamma(n_pacientes, shape = 2, scale = 4 + efecto_riesgo * 5), 2),
  num_lesiones_gad       = rpois(n_pacientes, lambda = 0.5 + efecto_riesgo * 1.2),
  volumen_cerebral_norm  = round(rnorm(n_pacientes, mean = 0.80 - efecto_riesgo * 0.05, sd = 0.04), 3),
  lesiones_medula_espinal = rpois(n_pacientes, lambda = 0.8 + efecto_riesgo * 1.0),
  espesor_capa_fibras_nerviosas_um = round(rnorm(n_pacientes, mean = 92 - efecto_riesgo * 8, sd = 7), 1)
)
write.csv(radiologia, "data_input/radiology_data.csv", row.names = FALSE)

# ---------------------------------------------------------------------------
# 1e. OUTCOME / SUPERVIVENCIA  (outcome_data.csv)
#     tiempo en meses hasta progresión confirmada (o censura), evento 0/1
# ---------------------------------------------------------------------------
riesgo_basal    <- 0.015 + efecto_riesgo * 0.025
tiempo_evento   <- rexp(n_pacientes, rate = riesgo_basal)
tiempo_censura  <- runif(n_pacientes, min = 24, max = 96)  # seguimiento máx. real

outcome <- data.frame(
  patient_id = id_paciente,
  tiempo_meses = round(pmin(tiempo_evento, tiempo_censura), 1),
  evento = as.integer(tiempo_evento <= tiempo_censura)
)
write.csv(outcome, "data_input/outcome_data.csv", row.names = FALSE)

cat("\n>>> Archivos de input simulados creados en ./data_input/\n")
cat("    - clinical_data.csv  (", nrow(clinico), "pacientes,", ncol(clinico) - 1, "variables )\n")
cat("    - genetic_data.csv   (", nrow(genetico), "pacientes,", ncol(genetico) - 1, "variantes )\n")
cat("    - immune_data.csv    (", nrow(inmuno), "pacientes,", ncol(inmuno) - 1, "biomarcadores )\n")
cat("    - radiology_data.csv (", nrow(radiologia), "pacientes,", ncol(radiologia) - 1, "features )\n")
cat("    - outcome_data.csv   (", nrow(outcome), "pacientes, tiempo + evento )\n\n")


#############################################################################
# 2. IMPORTACIÓN  (a partir de aquí, trátalo como si fueran tus datos reales)
#############################################################################

clinico    <- read.csv("data_input/clinical_data.csv", stringsAsFactors = FALSE)
genetico   <- read.csv("data_input/genetic_data.csv",  stringsAsFactors = FALSE)
inmuno     <- read.csv("data_input/immune_data.csv",   stringsAsFactors = FALSE)
radiologia <- read.csv("data_input/radiology_data.csv", stringsAsFactors = FALSE)
outcome    <- read.csv("data_input/outcome_data.csv",  stringsAsFactors = FALSE)

# Vista rápida de completitud por capa (paso obligatorio antes de integrar)
resumen_completitud <- data.frame(
  capa = c("clinico", "genetico", "inmuno", "radiologia", "outcome"),
  n_pacientes = c(nrow(clinico), nrow(genetico), nrow(inmuno), nrow(radiologia), nrow(outcome)),
  pct_missing = c(
    mean(is.na(clinico)) * 100,
    mean(is.na(genetico)) * 100,
    mean(is.na(inmuno)) * 100,
    mean(is.na(radiologia)) * 100,
    mean(is.na(outcome)) * 100
  )
)
print(resumen_completitud)


#############################################################################
# 3. PREPROCESAMIENTO POR CAPA
#############################################################################

imputar_mediana <- function(df) {
  for (col in names(df)) {
    if (is.numeric(df[[col]]) && any(is.na(df[[col]]))) {
      df[[col]][is.na(df[[col]])] <- median(df[[col]], na.rm = TRUE)
    }
  }
  df
}

clinico <- imputar_mediana(clinico)
inmuno  <- imputar_mediana(inmuno)

# Matrices numéricas por capa, escaladas, con patient_id como rownames
preparar_matriz <- function(df, excluir = "patient_id") {
  m <- df[, !(names(df) %in% excluir), drop = FALSE]
  m <- m[, sapply(m, is.numeric), drop = FALSE]   # solo columnas numéricas
  m <- scale(as.matrix(m))
  rownames(m) <- df$patient_id
  m
}

m_clinico    <- preparar_matriz(clinico)
m_genetico   <- preparar_matriz(genetico)
m_inmuno     <- preparar_matriz(inmuno)
m_radiologia <- preparar_matriz(radiologia)

# Aseguramos mismo orden de pacientes en todas las vistas
ids_comunes  <- Reduce(intersect, list(rownames(m_clinico), rownames(m_genetico),
                                        rownames(m_inmuno), rownames(m_radiologia),
                                        outcome$patient_id))

m_clinico    <- m_clinico[ids_comunes, ]
m_genetico   <- m_genetico[ids_comunes, ]
m_inmuno     <- m_inmuno[ids_comunes, ]
m_radiologia <- m_radiologia[ids_comunes, ]
outcome      <- outcome[match(ids_comunes, outcome$patient_id), ]

cat("\n>>> Pacientes con las 4 capas completas tras el emparejamiento:", length(ids_comunes), "\n\n")


#############################################################################
# 4. INTEGRACIÓN MULTIÓMICA CON SNF
#############################################################################

K_vecinos  <- 20
alpha_snf  <- 0.5
T_iter     <- 20

construir_afinidad <- function(m) {
  d <- SNFtool::dist2(m, m)
  SNFtool::affinityMatrix(d, K_vecinos, alpha_snf)
}

afin_clinico    <- construir_afinidad(m_clinico)
afin_genetico   <- construir_afinidad(m_genetico)
afin_inmuno     <- construir_afinidad(m_inmuno)
afin_radiologia <- construir_afinidad(m_radiologia)

red_fusionada <- SNFtool::SNF(
  list(afin_clinico, afin_genetico, afin_inmuno, afin_radiologia),
  K_vecinos, T_iter
)

# Estimar número óptimo de clusters sobre la red fusionada
estimacion_k <- SNFtool::estimateNumberOfClustersGivenGraph(red_fusionada, 2:6)
k_elegido <- estimacion_k[[1]]
cat(">>> Número de clusters sugerido por SNF:", k_elegido, "\n\n")

grupos_snf <- SNFtool::spectralClustering(red_fusionada, K = k_elegido)
names(grupos_snf) <- ids_comunes


#############################################################################
# 5. VALIDACIÓN DE ESTABILIDAD DEL CLUSTERING (bootstrap)
#############################################################################

# clusterboot necesita una función de clustering "envuelta"; usamos k-means
# sobre la concatenación escalada de las 4 vistas como aproximación rápida
# de estabilidad (alternativa: bootstrap re-ejecutando SNF completo en cada
# réplica, más costoso pero más riguroso - recomendado para el análisis final)

m_integrado <- cbind(m_clinico, m_genetico, m_inmuno, m_radiologia)

set.seed(123)
estabilidad <- fpc::clusterboot(
  m_integrado, B = 100,
  clustermethod = fpc::kmeansCBI,
  k = k_elegido, seed = 123
)

cat(">>> Estabilidad de clusters (índice de Jaccard promedio, bootstrap):\n")
print(estabilidad$bootmean)
cat("    (> 0.75 = estable; 0.6-0.75 = razonable; < 0.5 = no reportar)\n\n")


#############################################################################
# 6. VALIDACIÓN CLÍNICA: KAPLAN-MEIER + LOG-RANK
#############################################################################

datos_km <- data.frame(
  patient_id = ids_comunes,
  tiempo     = outcome$tiempo_meses,
  evento     = outcome$evento,
  subtipo    = factor(grupos_snf)
)

ajuste_km <- survival::survfit(Surv(tiempo, evento) ~ subtipo, data = datos_km)

grafico_km <- survminer::ggsurvplot(
  ajuste_km, data = datos_km,
  pval = TRUE, risk.table = TRUE, conf.int = TRUE,
  palette = "jco",
  xlab = "Tiempo (meses)", ylab = "Probabilidad de no progresión",
  legend.title = "Subtipo (SNF)",
  surv.median.line = "hv"
)

pdf("resultados/kaplan_meier_subtipos.pdf", width = 8, height = 7)
print(grafico_km)
dev.off()

cat(">>> Log-rank test global:\n")
print(survival::survdiff(Surv(tiempo, evento) ~ subtipo, data = datos_km))

if (k_elegido > 2) {
  cat("\n>>> Comparaciones por pares (corrección BH):\n")
  print(survminer::pairwise_survdiff(Surv(tiempo, evento) ~ subtipo,
                                      data = datos_km, p.adjust.method = "BH"))
}

cat("\n>>> Mediana de tiempo hasta progresión por subtipo:\n")
print(survminer::surv_median(ajuste_km))


#############################################################################
# 7. COX AJUSTADO POR CONFUSORES CLÍNICOS CONOCIDOS
#############################################################################

datos_cox <- merge(datos_km, clinico, by = "patient_id")

modelo_solo_clinico <- survival::coxph(
  Surv(tiempo, evento) ~ edad + tiempo_evolucion_an + edss_basal + sexo,
  data = datos_cox
)

modelo_con_subtipo <- survival::coxph(
  Surv(tiempo, evento) ~ subtipo + edad + tiempo_evolucion_an + edss_basal + sexo,
  data = datos_cox
)

cat("\n>>> Modelo Cox con subtipo SNF ajustado por confusores clínicos:\n")
print(summary(modelo_con_subtipo))

cat("\n>>> Supuesto de riesgos proporcionales (cox.zph):\n")
print(survival::cox.zph(modelo_con_subtipo))

cat("\n>>> ¿El subtipo mejora el modelo significativamente? (Likelihood Ratio Test)\n")
print(anova(modelo_solo_clinico, modelo_con_subtipo, test = "LRT"))

cat("\n>>> C-index: solo clínico vs. clínico + subtipo multiómico\n")
cat("    Solo clínico:        ", survival::concordance(modelo_solo_clinico)$concordance, "\n")
cat("    Clínico + subtipo:   ", survival::concordance(modelo_con_subtipo)$concordance, "\n\n")


#############################################################################
# 8. MODELO PREDICTIVO SUPERVISADO: RANDOM SURVIVAL FOREST
#    (sobre las 4 capas integradas, sin reducir a subtipo, para comparar
#     rendimiento predictivo puro vs. el enfoque de clustering)
#############################################################################

datos_rsf <- data.frame(
  tiempo = outcome$tiempo_meses,
  evento = outcome$evento,
  m_integrado
)
colnames(datos_rsf) <- make.names(colnames(datos_rsf))

modelo_rsf <- randomForestSRC::rfsrc(
  Surv(tiempo, evento) ~ ., data = datos_rsf,
  ntree = 1000, importance = TRUE, seed = -123
)

cat(">>> Random Survival Forest - resumen del modelo:\n")
print(modelo_rsf)

cat("\n>>> Importancia de variables (top 10):\n")
importancia <- sort(modelo_rsf$importance, decreasing = TRUE)
print(head(importancia, 10))

pdf("resultados/importancia_variables_rsf.pdf", width = 7, height = 6)
randomForestSRC::plot.variable(modelo_rsf, plots.per.page = 4)
dev.off()


#############################################################################
# 9. GUARDAR RESULTADOS CLAVE
#############################################################################

resultados_pacientes <- data.frame(
  patient_id = ids_comunes,
  subtipo_snf = grupos_snf,
  tiempo_meses = outcome$tiempo_meses,
  evento = outcome$evento
)
write.csv(resultados_pacientes, "resultados/subtipos_pacientes.csv", row.names = FALSE)

cat("\n>>> Pipeline completo. Resultados guardados en ./resultados/\n")
cat("    - kaplan_meier_subtipos.pdf\n")
cat("    - importancia_variables_rsf.pdf\n")
cat("    - subtipos_pacientes.csv\n")

##########################################################################################################
