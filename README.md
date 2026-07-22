# MachineLearningalghoritms

Aprendizaje supervisado (cuando hay una variable a predecir, ej. progresión de la enfermedad, conversión a EM progresiva)

Regresión logística / regresión de Cox (para supervivencia/tiempo hasta evento) — muy usados en clínica por su interpretabilidad
Random Forest / Gradient Boosting (XGBoost, LightGBM) — buenos con datos tabulares heterogéneos (clínicos + biomarcadores)
Support Vector Machines (SVM) — clásicos en clasificación de imágenes de resonancia magnética
Redes neuronales / Deep Learning (CNNs) — especialmente para datos radiológicos (segmentación de lesiones, volumetría cerebral)

Aprendizaje no supervisado (para descubrir subgrupos o patrones sin etiquetas previas)

Clustering (k-means, clustering jerárquico, clustering basado en densidad) — para identificar subtipos de pacientes o endotipos
Reducción de dimensionalidad (PCA, UMAP, t-SNE) — esenciales para visualizar y simplificar datos ómicos de alta dimensión
Factor analysis / NMF — para encontrar "firmas" latentes que combinen variables genéticas e inmunológicas

Integración multiómica (el corazón de este puesto, ya que combina múltiples tipos de datos)

Métodos de integración multi-view: 
  MOFA/MOFA2, 
  similarity network fusion (SNF), 
  multi-omics factor analysis
Modelos de grafos / redes (network-based methods) para relacionar biomarcadores genéticos-inmunológicos-clínicos
Modelos bayesianos jerárquicos para combinar fuentes de datos con distinta calidad/ruido

Otros conceptos clave que probablemente se esperen del candidato

Validación cruzada y control de overfitting (crucial con cohortes clínicas pequeñas)
Feature selection / feature engineering en datos ómicos de alta dimensión
Explicabilidad (SHAP, LIME) — importante en investigación clínica para justificar decisiones del modelo

Dado que trabajas en R, muchas de estas técnicas tienen paquetes muy sólidos: 

caret/tidymodels para el flujo general, 
randomForest/xgboost, 
survival/survminer para modelos de supervivencia, 
MOFA2 para integración multiómica, y 
Seurat/limma si hay datos de expresión génica.
