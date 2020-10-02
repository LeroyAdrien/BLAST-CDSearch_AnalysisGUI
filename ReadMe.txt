0)Prérequis:


-Python 3.8
Modules:
-Matplotlib 3.2.0
-tkinter 8.6
-scipy 1.4.1pip
-numpy 1.18.2
-itertools 8.2
-OS
-RE
-fnmatch



I) Lancement du programme:


Rentrer dans le dossier "Projet/Scripts"
Executer la commande "python3 main.py <path-vers-dossier-des-données>" (ex: python3 main.py ../donnees/)
/!\ le dossier donnees doit contenir tous les dossiers et fichiers avec des noms proches de ceux d'origine (voir le script fileparsing.py pour plus d'infos sur les noms acceptés)

II) Tour des fonctionnalités du programme:

Dans le menu principal:

        1) Calculer une preview des resultats de blast entre deux organismes
           Dans n'importe quel ordre:
           -cliquer sur "Rechercher selon BLAST"
           -sélectionner deux organismes parmis la liste
        2) Calculer une preview des resultats de CD Search entre deux organismes par Annotation fonctionnelle ou foncion de la protéines
          Dans n'importe quel ordre:
          -cliquer sur "Rechercher selon CD Search
          'Sélectionner COG ou FA (COG=annotation FA=Fonction des protéines,Attention le Calcul selon FA est assez long)
          -Sélectionner deux organismes dans les deux listes
         3) Filtrer en temps réel les résultats de la preview
            -Utiliser les curseurs ou les entrées de texte pour sélectionner un filtre à appliquer sur les résultats (format Xe-y pour E-value BLAST et CD Search)
         4) Accéder à plus d'informations sur les organismes sélectionnées dans les listes
            -Apès avoir choisi un organisme, Cliquer sur 'Organismes analysés'
         5) Exporter les résultats de la preview afin d'y appliquer plus de traitements:
            -Une fois que la preview vous convient, appuyer sur 'Dot Plot'
            
            
Dans le Menu Dot Plot:

        1)Zoomer, exporter, se déplacer sur le graph en plein écran:
          -Utiliser la barre d'outils disponible en bas à gauche de l'écran
        2)Afficher les blocs de synténie, ainsi que des informations statistiques sur leur distribution:
          -Cocher le bouton "Afficher Blocs de synténie (Si la distribution est normale (risque Alpha=0.05 test de shapiro, une distribution normale s'ajoute au graph, à tester en changeant le risque alpha très bas dans le script calculusFiltering.py fonction NormalLaw)
        3) Dans le cas ou des Résultats de BLAST et CD Search sont superposés, afficher l'intersection des coordonnées plus que leur union:
          -Cocher le bouton "Intersection des résultats
        4) Obtenir des P-Value d'enrichissement des différentes fonction dansles blocs de synténie grâce à un test exact de Fisher: 
          -Cliquer sur le bouton enrichissement de synténie
        5) Obtenir des diagrammes de Synténie à partir des résultats de DOT plot, ces diagrammes mettent en évidence les correspondances de chaques protines impliquées dans des blcos de synténie entre le protéome d'un organisme et d'un autre:
          -Cliquer sur Diagramme de Synténie
        6) Mettre en valeur certaines fonctions dans les blocs de syténie en vert
          -/!\ Le bouton blocs de synténie doit être sélectionné pour accéder au bouton chercher dans les blocs
          -Cliquer sur la checkbox 'Chercher une fonction dans les blocs'
          -Sélectionner une fonction dans le menu déroulant  



III)Contenu des différents scripts:

-fileparsing.py: Définit les fonction qui iront chercher les chemins vers les différents fichiers de l'analyse et en extrairont les données sous forme de dictionnaires

-calculusFiltering.py: Définit les fonctions qui opéreront différents traiteemnt ou filtres sur les données des dictionanires obtenus (calcul des coordonnées de BLAST, filtagre des données, blocs de synténies, statistiques sur les blocs de Synténie)

-GUIMainScreen.py: Interface Graphique principale, preview, filtres, accès aux autres interfaces

-GUIFullscreen.py: Interface Graphique qui donne accès aux options de synténie, statistiques sur la synténie, intersection des résultats, zoom, enregistrement de l'image, etc...

-FigurePlotting.py: Fonction matplotlib pour mettre à jour les graph de GUIFullscreen.py plus simplement

-GUIMoreInfo.py: Interface qui va chercher les informations sur les organismes sur lesquels on travaille pour en présenter succintement les caractéristiques

-GUIPopUpEnrichment.py: Interface qui cherche les enrichissements de certaines fonctions dans les blocs de Synteny par un test de Fisher exact

-GUISyntenyDiagram.py Interface qui va afficher les diagrammes de synténie 

-main.py Fonction de lancement du programe

