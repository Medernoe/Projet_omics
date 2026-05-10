library(shiny)
library(shinydashboard)
library(waiter)
library(ggplot2)
library(DT)
library(plotly)
library(shinyalert)
source("fonctions.R")



#===== chargement des packages temp pour enrichissiment 
library(clusterProfiler)
library(org.Hs.eg.db) #faire un modele variables pour les diff espece ex : murin > org.MM.eg.db
library(org.Mm.eg.db)
library(org.Dm.eg.db)
library(enrichplot)
library(ggplot2)
library(ggarchery)
