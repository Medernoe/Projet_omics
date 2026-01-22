library(shiny)
library(shinydashboard)
library(waiter)
library(shinycssloaders)
library(ggplot2)
library(DT)
library(plotly)
library(shinyalert)
source("fonctions.R")



#===== chargement des packages temp pour enrichissiment 
library(clusterProfiler)
library(org.Hs.eg.db) #faire un modele variables pour les diff espece ex : murin > org.MM.eg.db
library(enrichplot)
library(ggplot2)
library(ggarchery)
