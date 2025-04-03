---
title: "Análise de scRNA-seq e Machine Learning para Classificação de AMD"
author: "Barbara Dalmaso"
date: "`r Sys.Date()`"
output: github_document
---

# 📌 Visão Geral  
Este repositório contém análises de **single-cell RNA sequencing (scRNA-seq)** aplicadas ao estudo da **Degeneração Macular Relacionada à Idade (AMD)**.  
Foram utilizadas técnicas de **redução de dimensionalidade**, **clusterização** e **aprendizado de máquina supervisionado** para explorar padrões biológicos e classificar as células nos diferentes estados fisiológicos e patológicos.  

# 🎲 Banco de dados
Os dados analisados nesse projeto foram coletados do estudo de Kuchroo et al (2023) "Single-cell analysis reveals inflammatory interactions driving macular degeneration"  
O artigo está disponível em: https://www.nature.com/articles/s41467-023-37025-7

# 📊 Metodologia  
As seguintes abordagens foram aplicadas:  

🎯 1️. **Redução de Dimensionalidade:**  
   - **PCA** (Análise de Componentes Principais)  
   - **t-SNE** e **UMAP** (para preservar relações locais e globais)  

🎯 2. **Clusterização Não Supervisionada:**  
   - **K-Means**  
   - **Gaussian Mixture Model (GMM)**  
   - **Leiden**  

🎯 3. **Classificação Supervisionada:**  
   - **Random Forest**  
   - **XGBoost**  
   - **LightGBM**  

A acurácia dos métodos foi avaliada por métricas como **precisão, recall, F1-score** e **matriz de confusão**.  

---

# 📌 Contribuições
Sinta-se à vontade para abrir uma issue ou enviar um pull request se quiser sugerir melhorias!

📧 Para dúvidas ou sugestões, entre em contato via barbdalmaso@gmail.com.