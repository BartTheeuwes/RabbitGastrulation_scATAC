suppressPackageStartupMessages({
    library(ArchR) 
    library(data.table)
    library(purrr)
    library(parallel)
    library(dplyr)
    library(ggpubr)
    library(BSgenome)
    library(biomaRt)
    library(argparse)
    library(gridExtra)
    library(viridis)
})

# I/O
io = list()
io$basedir='/rds/project/rds-SDzz0CATGms/users/bt392/04_rabbit_scATAC/'
io$archR.directory = file.path(io$basedir, 'ArchR/Project')
io$output.directory <- file.path(io$basedir,"ArchR")

io$sup_figures = file.path(io$basedir,"Figures_sup")
io$main_figures = file.path(io$basedir,"Figures_main")
dir.create(file.path(io$sup_figures), showWarnings = FALSE)
dir.create(file.path(io$main_figures), showWarnings = FALSE)


setwd(io$output.directory)

addArchRThreads(threads = 1) 
addArchRVerbose(verbose = FALSE)

options(repr.plot.width=15, repr.plot.height=5)


# options
opts = list()
opts$samples = c('rabbit_BGRGP1', 
                 'rabbit_BGRGP2',
                 'rabbit_BGRGP3', 
                 'rabbit_BGRGP4', 
                 'rabbit_BGRGP5', 
                 'rabbit_BGRGP6', 
                 'rabbit_BGRGP7', 
                 'rabbit_BGRGP8') 

opts$stages = c('GD7',
                'GD8', 
                'GD9')

opts$stage.colors = c("GD7" = "#D56958",
                  "GD8" = "#6EC280",
                  "GD9" = "#6494D8",
                  "GD9_ExE" = "#BEB44D")

opts$celltype_colours <- c(
                      # Existing mouse colours
                      "Definitive Endoderm" = "#EF4E22",
                      "Gut endoderm" = "#EF5A9D",
                      "Megakaryocytes" = "#ac0404",
                      "Visceral YS endoderm" = "#9d506e",
                      "Parietal YS endoderm" = "#1A1A1A",
                      "Epiblast" = "#635547",
                      "Primitive Streak" = "#DABE99",
                      "Caudal epiblast" = "#9E6762",
                      "PGC" = "#FACB12",
                      "Anterior Primitive Streak" = "#C19F70",
                      "Node"="#153B3D",
                      "Notochord" = "#0F4A9C",
                      "Gut tube" = "#EF5A9D",
                      "Hindgut" = "#F397C0",
                      "Midgut"  = "#FF00B2",
                      "Foregut" = "#FFB7FF",
                      "Pharyngeal endoderm"="#95E1FF",
                      "Thyroid primordium"="#97BAD3",
                      "Nephron progenitors"="#E85639",
                      "Nascent mesoderm" = "#C594BF",
                      "Intermediate mesoderm" = "#139992",
                      "Caudal mesoderm" = "#3F84AA",
                      "Lateral plate mesoderm" = "#F9DFE6",
                      "Limb mesoderm" = "#E35F82",
                      "Forelimb" = "#D02D75",
                      "Presomitic mesoderm"="#5581CA",#"#0000ff",#fc
                      "Somitic mesoderm" = "#005579",
                      "Posterior somitic tissues" = "#5ADBE4",#"#40e0d0",#turquoise
                      "Paraxial mesoderm" = "#8DB5CE",
                      "Cranial mesoderm" = "#456722",#"#006400",#darkgreen
                      "Anterior somitic tissues"= "#D5E839",
                      "Sclerotome" = "#E3CB3A",#"#ffff00",#yellow
                      "Dermomyotome" = "#00BFC4",#"#a52a2a",#brown
                      "Pharyngeal mesoderm" = "#C9EBFB",
                      "Cardiopharyngeal progenitors" = "#556789",
                      "Anterior cardiopharyngeal progenitors"="#683ED8",
                      "Allantois" = "#532C8A",
                      "Mesenchyme" = "#CC7818",
                      "Mesothelium" = "#FF7F9C",
                      "Epicardium"="#F79083",
                      "EPDC" = "#FF487D",
                      "Cardiopharyngeal progenitors FHF"="#D780B0",
                      "Cardiomyocytes FHF 1"="#A64D7E",
                      "Cardiomyocytes FHF 2"="#B51D8D",
                      "Cardiopharyngeal progenitors SHF"="#4B7193",
                      "Cardiomyocytes SHF 1"="#5D70DC",
                      "Cardiomyocytes SHF 2"="#332C6C",
                      "Haematoendothelial progenitors" = "#FBBE92",
                      "Blood progenitors" = "#6C4B4C",
                      "Erythroid" = "#C72228",
                      "Erythroid/Masked"="red",
                      "MEP"="#EF4E22",
                      "EMP"="#7C2A47",
                      "YS endothelium"="#FF891C",
                      "Mesothelium-endothelium/Masked"="blue",#"#907b66",
                      "Allantois endothelium"="#2F4A60",
                      "Embryo proper endothelium"="#90E3BF",
                      "Venous endothelium"="#BD3400",
                      "Endocardium"="#9D0049",
                      "NMPs/Mesoderm-biased" = "#89C1F5",
                      "NMPs" = "#8EC792",
                      "Ectoderm" = "#FF675C",
                      "Optic vesicle" = "#BD7300",
                      "Ventral forebrain progenitors"="#A0B689",
                      "Early dorsal forebrain progenitors"="#0F8073",
                      "Late dorsal forebrain progenitors"="#7A9941",
                      "Midbrain/Hindbrain boundary"="#8AB3B5",
                      "Midbrain progenitors"="#9BF981",
                      "Dorsal midbrain neurons"="#12ED4C",
                      "Ventral hindbrain progenitors"="#7E907A",
                      "Dorsal hindbrain progenitors"="#2C6521",
                      "Hindbrain floor plate"="#BF9DA8",
                      "Hindbrain neural progenitors"="#59B545",
                      "Neural tube"="#233629",
                      "Migratory neural crest"="#4A6798",
                      "Branchial arch neural crest"="#BD84B0",
                      "Frontonasal mesenchyme"="#D3B1B1",
                      "Spinal cord progenitors"="#6B2035",
                      "Dorsal spinal cord progenitors"="#E273D6",
                      "Non-neural ectoderm 1" = "#F7F79E",
                      "Non-neural ectoderm 2" = "#FCFF00",
                      "Non-neural ectoderm 3" = "#FFF335",
                      "Non-neural ectoderm 4" = "#FFD731",
                      "Non-neural ectoderm 5" = "#DBB400",
                      "Placodal ectoderm"="#FF5C00",
                      "Otic placode"="#F1A262",
                      "Otic neural progenitors"="#00B000",
                      "Visceral endoderm" = "#F6BFCB",
                      "ExE endoderm" = "#7F6874",
                      "ExE ectoderm" = "#989898",
                      "Parietal endoderm" = "#1A1A1A",
                      "Blood progenitors 1" = "#f9decf",
                      "Blood progenitors 2" = "#c9a997",
                      "Caudal mesoderm" = "#3F84AA",
                      "Endothelium" = "#ff891c",
                      "Megakariocytes" = "#ac0404",
                      "Caudal neurectoderm" = "#354E23",
                      "Unknown"="#000000",
                      "Other"="#000000",

                      # Additional rabbit colours
                      'Cranial neural crest' = "#BD84B0",
                      "Spinal cord"= "#6B2035",
                      "Differentiating neurons"="#12ED4C",
                      "Floor plate" = "#BF9DA8",
                      "Hypoblast"="#7F6874",
                      'Forebrain'= "#FF4A46",
                      'Midbrain'= "#0000A6",
                      'Hindbrain'="#63FFAC",
                      'Roof plate'= "#772600",
                      'Visceral YS endoderm 2'= "#5B4534",
                      'Amnion 1'= '#549E79',
                        'Amnion'= '#549E79',
                      'Amnion 2'= '#f0e400',
                      'Amnion 3'= '#dea93c',
                      'Trophoblast'= '#989898',
                      'Cytotrophoblast'= '#46c5bf',
                      'Syncytiotrophoblast progenitors'= '#7766db',
                      'Syncytiotrophoblast Progenitors'= '#7766db',
                      "Syncytiotrophoblast" = '#7295e2'
)