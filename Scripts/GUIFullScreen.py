import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
import tkinter as tk
from calculusFiltering import *
from FigurePlotting import *
from GUIPopUpEnrichment import *
from GUISyntenyDiagram import *


class interface_fullscreen:

        def __init__(self,coordonneesXBLAST,coordonneesYBLAST,coordonneesXCDSearch,coordonneesYCDSearch,titreBLAST,titreCDSearch,dicoProteome,dicoFun,dicoNoms2COG,dicoCOG,dicoBLAST,dicoCDSearch,GCABLAST,GCACDSearch):
        
                #Certaines varialbes doivent être globales
                self.coordonneesXBLAST=coordonneesXBLAST
                self.coordonneesYBLAST=coordonneesYBLAST
                self.coordonneesXCDSearch=coordonneesXCDSearch
                self.coordonneesYCDSearch=coordonneesYCDSearch
                self.GCABLAST=GCABLAST
                self.GCACDSearch=GCACDSearch
                self.titreCDSearch=titreCDSearch
                self.titreBLAST=titreBLAST
                
                #Initialisation des coordonnés
                self.coordonnesXSyntenyCroissantBLAST=[]
                self.coordonnesYSyntenyCroissantBLAST=[]
                self.coordonnesXSyntenyDecroissantBLAST=[]
                self.coordonnesYSyntenyDecroissantBLAST=[]
                
                self.coordonnesXSyntenyCroissantCDSearch=[]
                self.coordonnesYSyntenyCroissantCDSearch=[]
                self.coordonnesXSyntenyDecroissantCDSearch=[]
                self.coordonnesYSyntenyDecroissantCDSearch=[]
                
                self.coordonnesXSyntenyCroissantIntersection=[]
                self.coordonnesYSyntenyCroissantIntersection=[]
                self.coordonnesXSyntenyDecroissantIntersection=[]
                self.coordonnesYSyntenyDecroissantIntersection=[]
                self.coordonneesXIntersection=[]
                self.coordonneesYIntersection=[]
                
                self.dicoProteome=dicoProteome
                self.dicoFun=dicoFun
                self.dicoNoms2COG=dicoNoms2COG
                self.COG=dicoCOG
                self.dicoBLAST=dicoBLAST
                self.dicoCDSearch=dicoCDSearch
                

                #Couleur
                self.couleurTexte='#D8D8D8'
                couleurBackground='#2c2f33'
                
                self.couleurBLASTSelected='#0057e7'
                self.couleurCDSearchSelected='#a80000'
                self.couleurIntersection='#cc6600'
                
                self.couleurSyntenyBLAST='#1a35ff'
                self.couleurSyntenyCDSearch='#ff1600'
                self.couleurSyntenyIntersection='#a35100'
                self.couleurFonction='#0cf00c'
                
                couleurBLASTHighlight='#a8c6fa'
                couleurCDSearchHighlight='#ff822e'
                
                couleurBorder='#3D4849'
                couleurGeneralSelected='#9baaab'
                couleurGeneralHighlight='#d5dcdc'
                
                self.styleSyntenyCroiss='s'
                self.styleSyntenyDecroiss='x'
                self.styleFonction='x'
                
                
                #Calcul de la synténie selon BLAST
                if len(self.coordonneesXBLAST)>0:
                        self.listeCoordonnesSyntenyCroissantBLAST,self.listeCoordonnesSyntenyDecroissantBLAST,self.statSyntenyBLAST=Syntenie(
                                                                                                                                             self.coordonneesXBLAST,
                                                                                                                                             self.coordonneesYBLAST)
                                                       
                        
                        if len(self.listeCoordonnesSyntenyCroissantBLAST)>0:
                                self.coordonnesXSyntenyCroissantBLAST=self.listeCoordonnesSyntenyCroissantBLAST[:,0]
                                self.coordonnesYSyntenyCroissantBLAST=self.listeCoordonnesSyntenyCroissantBLAST[:,1]
                                #Creation d'un liste de fonctions assocée à la coordonnée
                                self.fonctionsXSyntenyCroissantBLAST,self.fonctionsYSyntenyCroissantBLAST=Pos2FA(self.coordonnesXSyntenyCroissantBLAST,
                                                                                                                 self.coordonnesYSyntenyCroissantBLAST,
                                                                                                                 dicoProteome,
                                                                                                                 GCABLAST[0],
                                                                                                                 GCABLAST[1],
                                                                                                                 dicoNoms2COG,
                                                                                                                 dicoCOG)
                        if len(self.listeCoordonnesSyntenyDecroissantBLAST)>0:
                                self.coordonnesXSyntenyDecroissantBLAST=self.listeCoordonnesSyntenyDecroissantBLAST[:,0]
                                self.coordonnesYSyntenyDecroissantBLAST=self.listeCoordonnesSyntenyDecroissantBLAST[:,1]
                                #Creation d'une liste de fonctionns associée aux coordonnés
                                self.fonctionsXSyntenyDecroissantBLAST,self.fonctionsYSyntenyDecroissantBLAST=Pos2FA(self.coordonnesXSyntenyDecroissantBLAST,
                                                                                                                     self.coordonnesYSyntenyDecroissantBLAST,
                                                                                                                     dicoProteome,
                                                                                                                     GCABLAST[0],
                                                                                                                     GCABLAST[1],
                                                                                                                     dicoNoms2COG,
                                                                                                                     dicoCOG)
                                
                        #Calcul d'infos supplémentaires:
                        #Normalité et statistiques
                        self.normaldistribBLAST=NormalLaw(self.statSyntenyBLAST)
                        #Cover des blocs de synténie 
                        self.syntenyCoverBLAST=SyntenieCover(self.coordonnesXSyntenyCroissantBLAST,
                                                             self.coordonnesYSyntenyCroissantBLAST,
                                                             self.coordonnesXSyntenyDecroissantBLAST,
                                                             self.coordonnesYSyntenyDecroissantBLAST,
                                                             dicoProteome,
                                                             GCABLAST[0],
                                                             GCABLAST[1]
                                                             )
                #Calcul de la synténie selon CD Search
                if len(self.coordonneesXCDSearch)>0:
                        self.listeCoordonnesSyntenyCroissantCDSearch,self.listeCoordonnesSyntenyDecroissantCDSearch,self.statSyntenyCDSearch=Syntenie(
                                                                                                                                                      self.coordonneesXCDSearch,
                                                                                                                                                      self.coordonneesYCDSearch
                                                                                                                                                      )
                        

                        if len(self.listeCoordonnesSyntenyCroissantCDSearch)>0:
                                self.coordonnesXSyntenyCroissantCDSearch=self.listeCoordonnesSyntenyCroissantCDSearch[:,0]
                                self.coordonnesYSyntenyCroissantCDSearch=self.listeCoordonnesSyntenyCroissantCDSearch[:,1]
                                self.fonctionsXSyntenyCroissantCDSearch,self.fonctionsYSyntenyCroissantCDSearch=Pos2FA(self.coordonnesXSyntenyCroissantCDSearch,
                                                                                                                       self.coordonnesYSyntenyCroissantCDSearch,
                                                                                                                       dicoProteome,
                                                                                                                       GCACDSearch[0],
                                                                                                                       GCACDSearch[1],
                                                                                                                       dicoNoms2COG,
                                                                                                                       dicoCOG)
                        if len(self.listeCoordonnesSyntenyDecroissantCDSearch)>0:
                                self.coordonnesXSyntenyDecroissantCDSearch=self.listeCoordonnesSyntenyDecroissantCDSearch[:,0]
                                self.coordonnesYSyntenyDecroissantCDSearch=self.listeCoordonnesSyntenyDecroissantCDSearch[:,1]
                                self.fonctionsXSyntenyDecroissantCDSearch,self.fonctionsYSyntenyDecroissantCDSearch=Pos2FA(self.coordonnesXSyntenyDecroissantCDSearch,
                                                                                                                           self.coordonnesYSyntenyDecroissantCDSearch,
                                                                                                                           dicoProteome,
                                                                                                                           GCACDSearch[0],
                                                                                                                           GCACDSearch[1],
                                                                                                                           dicoNoms2COG,
                                                                                                                           dicoCOG)

                        #Calcul d'infos supplémentaires:
                        #Normalité et statistiques
                        self.normaldistribCDSearch=NormalLaw(self.statSyntenyCDSearch)
                        #Cover des blocs de synténie 
                        self.syntenyCoverCDSearch=SyntenieCover(self.coordonnesXSyntenyCroissantCDSearch,
                                                             self.coordonnesYSyntenyCroissantCDSearch,
                                                             self.coordonnesXSyntenyDecroissantCDSearch,
                                                             self.coordonnesYSyntenyDecroissantCDSearch,
                                                             dicoProteome,
                                                             GCACDSearch[0],
                                                             GCACDSearch[1]
                                                             )
                        
                #Calcul de l'intersection de BLAST et CDSearch'
                if len(self.coordonneesXBLAST)>0 and len(self.coordonneesXCDSearch)>0:
                        variablesAIntersecter=(self.coordonneesXBLAST,
                                               self.coordonneesYBLAST,
                                               self.coordonneesXCDSearch,
                                               self.coordonneesYCDSearch
                                               )
                                               
                        self.coordonneesXIntersection,self.coordonneesYIntersection=Intersection(*variablesAIntersecter)
                        
                        #Calcul de la synténie de l'intersection des deux résultats'
                        if len(self.coordonneesXIntersection)>0:
                                X=self.coordonneesXIntersection
                                Y=self.coordonneesYIntersection
                                self.listeCoordonnesSyntenyCroissantIntersection,self.listeCoordonnesSyntenyDecroissantIntersection,self.statSyntenyIntersection=Syntenie(X,Y)
                                #Calcul d'infos supplémentaires:
                                self.fonctionsXSyntenyCroissantIntersection=np.array([])
                                self.fonctionsYSyntenyCroissantIntersection=np.array([])
                                self.fonctionsXSyntenyDecroissantIntersection=np.array([])
                                self.fonctionsYSyntenyDecroissantIntersection=np.array([])
                                if len(self.listeCoordonnesSyntenyCroissantIntersection)>0:
                                        self.coordonnesXSyntenyCroissantIntersection=self.listeCoordonnesSyntenyCroissantIntersection[:,0]
                                        self.coordonnesYSyntenyCroissantIntersection=self.listeCoordonnesSyntenyCroissantIntersection[:,1]
                                        self.fonctionsXSyntenyCroissantIntersection,self.fonctionsYSyntenyCroissantIntersection=Pos2FA(self.coordonnesXSyntenyCroissantIntersection,
                                                                                                                                       self.coordonnesYSyntenyCroissantIntersection,
                                                                                                                                       dicoProteome,
                                                                                                                                       GCABLAST[0],
                                                                                                                                       GCABLAST[1],
                                                                                                                                       dicoNoms2COG,
                                                                                                                                       dicoCOG)
                                if len(self.listeCoordonnesSyntenyDecroissantIntersection)>0:
                                        self.coordonnesXSyntenyDecroissantIntersection=self.listeCoordonnesSyntenyDecroissantIntersection[:,0]
                                        self.coordonnesYSyntenyDecroissantIntersection=self.listeCoordonnesSyntenyDecroissantIntersection[:,1]
                                        self.fonctionsXSyntenyDecroissantIntersection,self.fonctionsYSyntenyDecroissantIntersection=Pos2FA(self.coordonnesXSyntenyDecroissantIntersection,
                                                                                                                                           self.coordonnesYSyntenyDecroissantIntersection,
                                                                                                                                           dicoProteome,
                                                                                                                                           GCABLAST[0],
                                                                                                                                           GCABLAST[1],
                                                                                                                                           dicoNoms2COG,
                                                                                                                                           dicoCOG)
                                #Calcul d'infos supplémentaires:
                                #Normalité et statistiques
                                self.normaldistribIntersection=(False,[],[],('NaN','Nan','Nan','Nan'))
                                if len(self.statSyntenyIntersection)>8:
                                        self.normaldistribIntersection=NormalLaw(self.statSyntenyIntersection)
                                #Cover des blocs de synténie 
                                self.syntenyCoverIntersection=SyntenieCover(self.coordonnesXSyntenyCroissantIntersection,
                                                                     self.coordonnesYSyntenyCroissantIntersection,
                                                                     self.coordonnesXSyntenyDecroissantIntersection,
                                                                     self.coordonnesYSyntenyDecroissantIntersection,
                                                                     dicoProteome,
                                                                     GCACDSearch[0],
                                                                     GCACDSearch[1]
                                                                     )

                interfaceFullscreen=tk.Tk()
                interfaceFullscreen.title('Plein Ecran')
                interfaceFullscreen.geometry('1300x700')
        #Frames

                #Main Frame
                fenetre=tk.Frame(interfaceFullscreen,bg=couleurBackground)
                fenetre.pack(fill="both",expand=True)

                fenetre.grid_columnconfigure(0,weight=5, uniform="group1")
                fenetre.grid_rowconfigure(0,weight=1, uniform="group2")
                fenetre.grid_rowconfigure(1,weight=8, uniform="group2")
                fenetre.grid_rowconfigure(2,weight=2, uniform="group2")


                #Premiers plans de frames:

                #Frame du haut: prompt du BLAST à comparer
                frameFullScreenHaut=tk.Frame(fenetre,bg=couleurBackground)
                frameFullScreenHaut.grid(row=0,column=0,columnspan=2,sticky="nsew")
                #Paramètre du frame
                frameFullScreenHaut.grid_rowconfigure(0,weight=1, uniform='group1')

                #Frame du milieu: Graph en plein écran 
                frameFullScreenMilieu=tk.Frame(fenetre,bg=couleurBackground)
                frameFullScreenMilieu.grid(row=1,column=0,sticky="nsew")
                #Paramètre du frame
                frameFullScreenMilieu.grid_rowconfigure(0,weight=1)
                frameFullScreenMilieu.grid_columnconfigure(0,weight=1)
                
                frameGraph=tk.Frame(frameFullScreenMilieu,bg=couleurBackground)
                frameGraph.grid(row=0,column=0,sticky="nsew")

                #Frame du bas:Options sur le graph, blocs de synténie et close
                frameFullScreenBas=tk.Frame(fenetre,bg=couleurBackground)
                frameFullScreenBas.grid(row=2,column=0,sticky="nsew")
                #Paramètre du frame
                frameFullScreenBas.grid_rowconfigure(0,weight=1, uniform='group1')
                frameFullScreenBas.grid_rowconfigure(1,weight=1, uniform='group1')
                
                frameFullScreenBas.grid_columnconfigure(0,weight=1)
                frameFullScreenBas.grid_columnconfigure(1,weight=1)
                frameFullScreenBas.grid_columnconfigure(2,weight=1)
                frameFullScreenBas.grid_columnconfigure(3,weight=1)
                frameFullScreenBas.grid_columnconfigure(4,weight=1)
                
                
                #frame de la toolbar:
                frameToolbar=tk.Frame(frameFullScreenBas,bg=couleurBackground)
                frameToolbar.grid(row=0,column=0)


        #Widgets:

                #Figure Plein écran 
                self.fig = Figure(figsize=(5, 4), dpi=100)
                self.fig.patch.set_facecolor(couleurBackground)
                #self.fig.patch.set_facecolor(couleurBackground)
                self.ax1=self.fig.add_subplot(111)
                #self.ax1.set_title("Correspondance des deux protéomes",color=self.couleurTexte)
                self.ax1.set_facecolor(self.couleurTexte)
                #self.ax1.set_facecolor(self.couleurTexte)
                self.xlabel=''
                self.ylabel=''
                
                
                if len(self.coordonneesXBLAST)>0:
                        frameFullScreenHaut.grid_columnconfigure(0,weight=1,uniform='group2')
                        fullScreenLabelBLAST=tk.Label(frameFullScreenHaut,text="BLAST:"+titreBLAST[0]+" vs "+titreBLAST[1])
                        fullScreenLabelBLAST.configure(bg=couleurBackground,fg=self.couleurTexte)
                        fullScreenLabelBLAST.grid(row=0,column=0,sticky="nsew")
                        self.scatterBLAST=self.ax1.scatter(self.coordonneesXBLAST,self.coordonneesYBLAST,c=self.couleurBLASTSelected,s=0.1)
                        self.scatterBLAST.set_label("BLASTp")
                        self.xlabel+="BLAST: "+titreBLAST[0]+" "
                        self.ylabel+="BLAST: "+titreBLAST[1]
                if len(self.coordonneesXCDSearch)>0:
                        frameFullScreenHaut.grid_columnconfigure(1,weight=1,uniform='group2')
                        fullScreenLabelCDSearch=tk.Label(frameFullScreenHaut,text="CDSearch:"+titreCDSearch[0]+' vs '+titreCDSearch[1])
                        fullScreenLabelCDSearch.configure(bg=couleurBackground,fg=self.couleurTexte)
                        fullScreenLabelCDSearch.grid(row=0,column=1,sticky="nsew")
                        self.scatterCDSearch=self.ax1.scatter(self.coordonneesXCDSearch,self.coordonneesYCDSearch,c=self.couleurCDSearchSelected,s=0.1)
                        self.scatterCDSearch.set_label("CD Search")
                        self.xlabel+="CDSearch: "+titreCDSearch[0]
                        self.ylabel+="CDSearch: "+titreCDSearch[1]
                
                self.ax1.set_xlabel(self.xlabel,color=self.couleurTexte)
                self.ax1.set_ylabel(self.ylabel,color=self.couleurTexte)        
                box1 = self.ax1.get_position()
                self.ax1.set_position([box1.x0-(box1.x0*0.5), box1.y0,
                                    box1.width*1.1, box1.height ])

                self.ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=2)
                self.ax1.tick_params(axis='x', colors=self.couleurTexte)
                self.ax1.tick_params(axis='y', colors=self.couleurTexte)
                self.fullScreenCanvas = FigureCanvasTkAgg(self.fig, master=frameGraph)  # A tk.DrawingArea.
                self.fullScreenCanvas.draw()
                self.fullScreenCanvas.get_tk_widget().pack(side='top',fill='both',expand=1)
                
                #Toolbar de la figure
                self.toolbar = NavigationToolbar2Tk(self.fullScreenCanvas,frameToolbar)
                self.toolbar.configure(bg=couleurBackground,highlightbackground=couleurBackground)
                self.toolbar.update()
                self.fullScreenCanvas.get_tk_widget().pack(side='top',expand=1)
                
                
                #Bouton Bloc de synténie
                self.syntenyDisplayState=tk.IntVar(interfaceFullscreen)
                self.syntenyDisplayState.set(0)
                syntenyDisplayCheckButton=tk.Checkbutton(master=frameFullScreenBas, text='Afficher Blocs de synténie',indicatoron=0,variable=self.syntenyDisplayState,fg=self.couleurTexte)
                syntenyDisplayCheckButton.configure(highlightbackground=couleurBackground,activebackground=couleurGeneralHighlight,bg=couleurBackground,pady=10,padx=10)
                syntenyDisplayCheckButton.grid(row=0,column=1)
                self.syntenyDisplayState.trace("w",self.DisplayingSynteny)
                
                #Bouton Diagramme de synténies
                syntenyDiagramButton=tk.Button(frameFullScreenBas,text='Diagramme de Synténie',command=self.NewWindowsDiagram,fg=self.couleurTexte,bg=couleurBackground,pady=8)
                syntenyDiagramButton.configure(highlightbackground=couleurBackground,activebackground=couleurGeneralHighlight)
                syntenyDiagramButton.grid(row=0,column=2)
                #Bouton intersection des résultats
                self.intersectionDisplayState=tk.IntVar(interfaceFullscreen)
                self.intersectionDisplayState.set(0)
                intersectionDisplayCheckButton=tk.Checkbutton(master=frameFullScreenBas,text='Intersection des résultats',indicatoron=0,variable=self.intersectionDisplayState,fg=self.couleurTexte)
                intersectionDisplayCheckButton.configure(highlightbackground=couleurBackground,activebackground=couleurGeneralHighlight,bg=couleurBackground,pady=10,padx=10)
                intersectionDisplayCheckButton.grid(row=0,column=3)
                
                if len(self.coordonneesXBLAST)==0 or len(self.coordonneesXCDSearch)==0:
                        intersectionDisplayCheckButton.configure(state='disabled')
                        
                self.intersectionDisplayState.trace("w",self.DisplayingIntersection)
                
                #Bouton Enrichissement des blocs de synténies
                enrichmentButton=tk.Button(frameFullScreenBas,text='Enrichissement de Synténie',command=self.NewWindowEnrichment,fg=self.couleurTexte,bg=couleurBackground,pady=8,padx=70)
                enrichmentButton.configure(highlightbackground=couleurBackground,activebackground=couleurGeneralHighlight)
                enrichmentButton.grid(row=0,column=4)
                
                if len(self.coordonneesXCDSearch)==0 and len(self.coordonneesXBLAST)==0:
                        enrichmentButton.configure(state='disabled')
                        
                #Checkbox Mise en valeur d'une fonction
                self.fonctionHighlightState=tk.IntVar(interfaceFullscreen)
                self.fonctionHighlightState.set(0)
                self.fonctionHighlightState.trace('w',self.SelectionFonction)
                self.fonctionHighlightCheckButton=tk.Checkbutton(master=frameFullScreenBas,text="Chercher une fonction dans les blocs",variable=self.fonctionHighlightState)
                self.fonctionHighlightCheckButton.configure(highlightbackground=couleurBackground,activebackground=couleurGeneralHighlight,bg=couleurBackground,fg=self.couleurTexte,indicatoron=0,pady=10)
                self.fonctionHighlightCheckButton.grid(row=1,column=0)
                if self.syntenyDisplayState.get()==0:
                        self.fonctionHighlightCheckButton.configure(state='disabled')
                        
                #Liste de la fonction à mettre en valeurs
                self.listeMenuFonctionBackend=[]
                self.listeMenuFonctionFrontend=[]
                
                for i in dicoFun.keys():
                        self.listeMenuFonctionBackend.append(i)
                        self.listeMenuFonctionFrontend.append(dicoFun[i])
                self.choixFonction=tk.StringVar(frameFullScreenBas)
                self.choixFonction.set('')
                self.choixFonction.trace("w",self.SelectionFonction)
                self.menuFonction=tk.OptionMenu(frameFullScreenBas,self.choixFonction,*self.listeMenuFonctionFrontend)
                self.menuFonction['menu'].configure(fg=self.couleurTexte,bg=couleurBackground)
                self.menuFonction.config(fg=self.couleurTexte,bg=couleurBackground,highlightbackground=couleurBackground,activebackground=couleurGeneralHighlight,padx=200)
                self.menuFonction.grid(row=1,column=1,columnspan=3)
                if self.fonctionHighlightState.get()==0 or self.syntenyDisplayState.get()==0:
                        self.menuFonction.configure(state='disabled')
                        
                        
                #Bouton quitter
                quitButton = tk.Button(master=frameFullScreenBas, text="Quitter", command=interfaceFullscreen.destroy,fg=self.couleurTexte)
                quitButton.configure(highlightbackground=couleurBackground,activebackground=couleurGeneralHighlight,bg=couleurBackground,pady=8,padx=160)
                quitButton.grid(row=1,column=4)

                tk.mainloop()
                
        #Creation de la Fenêtre d'analyse d'enrichissement'
        def NewWindowEnrichment(self,*args):
                a=interface_popup_enrichment(self.dicoProteome,
                                             self.dicoFun,
                                             self.dicoNoms2COG,
                                             self.COG,
                                             self.coordonnesXSyntenyCroissantCDSearch,
                                             self.coordonnesYSyntenyCroissantCDSearch,
                                             self.coordonnesXSyntenyDecroissantCDSearch,
                                             self.coordonnesYSyntenyDecroissantCDSearch,
                                             self.GCACDSearch,
                                             self.titreCDSearch,
                                             self.coordonnesXSyntenyCroissantBLAST,
                                             self.coordonnesYSyntenyCroissantBLAST,
                                             self.coordonnesXSyntenyDecroissantBLAST,
                                             self.coordonnesYSyntenyDecroissantBLAST,
                                             self.GCABLAST,
                                             self.titreBLAST
                                             )
        #Creation de la fenêtre du diagramme de Synténie
        def NewWindowsDiagram(self,*args):
                tailleBLAST=(0,0)
                tailleCDSearch=(0,0)
                XCDSearchSyntenyCroissantFiltre=[]
                YCDSearchSyntenyCroissantFiltre=[]
                XCDSearchSyntenyDecroissantFiltre=[]
                YCDSearchSyntenyDecroissantFiltre=[]
                XBLASTSyntenyCroissantFiltre=[]
                YBLASTSyntenyCroissantFiltre=[]
                XBLASTSyntenyDecroissantFiltre=[]
                YBLASTSyntenyDecroissantFiltre=[]
                XIntersectionSyntenyCroissantFiltre=[]
                YIntersectionSyntenyCroissantFiltre=[]
                XIntersectionSyntenyDecroissantFiltre=[]
                YIntersectionSyntenyDecroissantFiltre=[]
                #Filtrage des coordonnées de BLAST
                if len(self.coordonnesXSyntenyCroissantBLAST)>0:
                        #Calcul de la taille du protéome
                        tailleBLAST=(len(self.dicoProteome[self.GCABLAST[0]]),len(self.dicoProteome[self.GCABLAST[1]]))
                        #filtrage des blocs de Synténie par E-value 
                        dicoEvalueBLASTCroissant=EvalueSyntenieBLAST(self.coordonnesXSyntenyCroissantBLAST,
                                                                     self.coordonnesYSyntenyCroissantBLAST,
                                                                     self.GCABLAST[0],
                                                                     self.GCABLAST[1],
                                                                     self.dicoBLAST,
                                                                     self.dicoProteome)
                        XBLASTSyntenyCroissantFiltre,YBLASTSyntenyCroissantFiltre=FilterEvalueSyntenie(dicoEvalueBLASTCroissant)
                        #filtrage des blocs de Synténie par E-value 
                        dicoEvalueBLASTDecroissant=EvalueSyntenieBLAST(self.coordonnesXSyntenyDecroissantBLAST,
                                                                     self.coordonnesYSyntenyDecroissantBLAST,
                                                                     self.GCABLAST[0],
                                                                     self.GCABLAST[1],
                                                                     self.dicoBLAST,
                                                                     self.dicoProteome)
                        XBLASTSyntenyDecroissantFiltre,YBLASTSyntenyDecroissantFiltre=FilterEvalueSyntenie(dicoEvalueBLASTDecroissant)
                #Filtrage des coordonnées de CDSearch
                if len(self.coordonnesXSyntenyCroissantCDSearch)>0:
                        tailleCDSearch=(len(self.dicoProteome[self.GCACDSearch[0]]),len(self.dicoProteome[self.GCACDSearch[1]]))
                        #filtrage des blocs de Synténie par E-value 
                        dicoEvalueCDSearchCroissant=EvalueSyntenieCDSearch(self.coordonnesXSyntenyCroissantCDSearch,
                                                                     self.coordonnesYSyntenyCroissantCDSearch,
                                                                     self.GCACDSearch[0],
                                                                     self.GCACDSearch[1],
                                                                     self.dicoCDSearch,
                                                                     self.dicoProteome,
                                                                     self.dicoNoms2COG)
                        XCDSearchSyntenyCroissantFiltre,YCDSearchSyntenyCroissantFiltre=FilterEvalueSyntenie(dicoEvalueCDSearchCroissant)
                        #filtrage des blocs de Synténie par E-value 
                        dicoEvalueCDSearchDecroissant=EvalueSyntenieCDSearch(self.coordonnesXSyntenyDecroissantCDSearch,
                                                                     self.coordonnesYSyntenyDecroissantCDSearch,
                                                                     self.GCACDSearch[0],
                                                                     self.GCACDSearch[1],
                                                                     self.dicoCDSearch,
                                                                     self.dicoProteome,
                                                                     self.dicoNoms2COG)
                        XCDSearchSyntenyDecroissantFiltre,YCDSearchSyntenyDecroissantFiltre=FilterEvalueSyntenie(dicoEvalueCDSearchDecroissant)
                a=interfaceSyntenyDiagram(
                                          self.dicoProteome,
                                          self.dicoBLAST,
                                          self.dicoCDSearch,
                                          self.titreCDSearch,
                                          self.titreBLAST,
                                          np.array(XCDSearchSyntenyCroissantFiltre),
                                          np.array(YCDSearchSyntenyCroissantFiltre),
                                          np.array(XCDSearchSyntenyDecroissantFiltre),
                                          np.array(YCDSearchSyntenyDecroissantFiltre),
                                          np.array(XBLASTSyntenyCroissantFiltre),
                                          np.array(YBLASTSyntenyCroissantFiltre),
                                          np.array(XBLASTSyntenyDecroissantFiltre),
                                          np.array(YBLASTSyntenyDecroissantFiltre),
                                          tailleBLAST,
                                          tailleCDSearch
                             )
                             

                             
                             
                             
        #Affichage des blocs de synténie
        def DisplayingSynteny(self,*args):          
        
                #Si synténie activée
                if self.syntenyDisplayState.get()==1:
                        #Permettre le choix d'une fonction à montrer'
                        self.menuFonction.configure(state='normal')
                        self.fonctionHighlightCheckButton.configure(state='normal')
                        #Mettre à jour le format des plots
                        self.fig.clear()
                        #self.ax1.clear()
                        self.ax1=self.fig.add_subplot(121)
                        self.ax1.tick_params(axis='x', colors=self.couleurTexte)
                        self.ax1.tick_params(axis='y', colors=self.couleurTexte)
                        self.ax1.set_xlabel(self.xlabel,color=self.couleurTexte)
                        self.ax1.set_ylabel(self.ylabel,color=self.couleurTexte)
                        self.ax1.set_facecolor(self.couleurTexte)
                        
                        self.ax2=self.fig.add_subplot(122)
                        self.ax2.set_facecolor(self.couleurTexte)
                        self.ax2.set_xlabel("Longueur des blocs",color=self.couleurTexte)
                        self.ax2.set_ylabel("Nombre de blocs",color=self.couleurTexte)
                        self.ax2.tick_params(axis='x', colors=self.couleurTexte)
                        self.ax2.tick_params(axis='y', colors=self.couleurTexte)
                    
                        self.ax3=self.ax2.twinx()
                        self.ax3.set_facecolor(self.couleurTexte)
                        self.ax3.tick_params(axis='x', colors=self.couleurTexte)
                        self.ax3.tick_params(axis='y', colors=self.couleurTexte)
                        #Et intersection activée
                        if self.intersectionDisplayState.get()==1:
                        
                                X,Y=self.coordonneesXIntersection,self.coordonneesYIntersection
                                self.scatterIntersection=self.ax1.scatter(X,Y,c=self.couleurIntersection,s=0.2)
                                self.scatterIntersection.set_label('Intersection des deux méthodes de recherche')
                                
                        #afficher la synténie de l'intersection 
                                PlottingSyntenyIntersection(
                                                        self.coordonnesXSyntenyCroissantIntersection,
                                                        self.coordonnesYSyntenyCroissantIntersection,
                                                        self.coordonnesXSyntenyDecroissantIntersection,
                                                        self.coordonnesYSyntenyDecroissantIntersection,
                                                        self.statSyntenyIntersection,
                                                        self.normaldistribIntersection,
                                                        self.syntenyCoverIntersection,
                                                        self.couleurSyntenyIntersection,
                                                        self.styleSyntenyCroiss,
                                                        self.styleSyntenyDecroiss,
                                                        self.ax1,
                                                        self.ax2,
                                                        self.ax3
                                                        )
                        #et l'intersection désactivée'
                        if self.intersectionDisplayState.get()==0:
                                #Réafficher les valeurs de base :
                                if len(self.coordonneesXBLAST)>0:
                                        self.scatterBLAST=self.ax1.scatter(self.coordonneesXBLAST,self.coordonneesYBLAST,c=self.couleurBLASTSelected,s=0.1)
                                        self.scatterBLAST.set_label("BLASTp")
                                if len(self.coordonneesXCDSearch)>0:
                                        self.scatterCDSearch=self.ax1.scatter(self.coordonneesXCDSearch,self.coordonneesYCDSearch,c=self.couleurCDSearchSelected,s=0.1)
                                        self.scatterCDSearch.set_label("CD Search")
                                #Mettre à jour le format des plots
                        #afficher la synténie des deux 
                                #SiBlast est là 
                                if len(self.coordonneesXBLAST)>0:
                                        PlottingSyntenyBLAST(
                                                                self.coordonnesXSyntenyCroissantBLAST,
                                                                self.coordonnesYSyntenyCroissantBLAST,
                                                                self.coordonnesXSyntenyDecroissantBLAST,
                                                                self.coordonnesYSyntenyDecroissantBLAST,
                                                                self.statSyntenyBLAST,
                                                                self.normaldistribBLAST,
                                                                self.syntenyCoverBLAST,
                                                                self.couleurSyntenyBLAST,
                                                                self.styleSyntenyCroiss,
                                                                self.styleSyntenyDecroiss,
                                                                self.ax1,
                                                                self.ax2,
                                                                self.ax3
                                                                )
                                #Si CD Search est là
                                if len(self.coordonneesXCDSearch)>0:
                                        PlottingSyntenyCDSearch(
                                                                self.coordonnesXSyntenyCroissantCDSearch,
                                                                self.coordonnesYSyntenyCroissantCDSearch,
                                                                self.coordonnesXSyntenyDecroissantCDSearch,
                                                                self.coordonnesYSyntenyDecroissantCDSearch,
                                                                self.statSyntenyCDSearch,
                                                                self.normaldistribCDSearch,
                                                                self.syntenyCoverCDSearch,
                                                                self.couleurSyntenyCDSearch,
                                                                self.styleSyntenyCroiss,
                                                                self.styleSyntenyDecroiss,
                                                                self.ax1,
                                                                self.ax2,
                                                                self.ax3
                                                                )
                        self.ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=2)
                        self.ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=2)
                                        
                if self.syntenyDisplayState.get()==0:
                        #Mettre à jour le format des plots
                        self.fig.clear()
                        self.ax1=self.fig.add_subplot(111)
                        self.ax1.set_xlabel(self.xlabel,color=self.couleurTexte)
                        self.ax1.set_ylabel(self.ylabel,color=self.couleurTexte)
                        self.ax1.tick_params(axis='x', colors=self.couleurTexte)
                        self.ax1.tick_params(axis='y', colors=self.couleurTexte)
                        #Remplir à nouveaux le plot n°1

                                
                        
                        
                        #et l'intersection activée
                        if self.intersectionDisplayState.get()==1:
                                X,Y=self.coordonneesXIntersection,self.coordonneesYIntersection
                                self.scatterIntersection=self.ax1.scatter(X,Y,c=self.couleurIntersection,s=0.2)
                                self.scatterIntersection.set_label('Intersection des deux méthodes de recherche')
                        
                        #et l'intersection désactivée
                        if self.intersectionDisplayState.get()==0:
                        #Effacer la synténie des deux 
                                #Si Blast est là
                                if len(self.coordonneesXBLAST)>0:
                                        self.scatterBLAST=self.ax1.scatter(self.coordonneesXBLAST,self.coordonneesYBLAST,c=self.couleurBLASTSelected,s=0.1)
                                        self.scatterBLAST.set_label("BLASTp")
                                #Si CD Search est là
                                if len(self.coordonneesXCDSearch)>0:
                                        self.scatterCDSearch=self.ax1.scatter(self.coordonneesXCDSearch,self.coordonneesYCDSearch,c=self.couleurCDSearchSelected,s=0.1)
                                        self.scatterCDSearch.set_label("CD Search")
                                        
                        self.ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=2)
                
                self.fullScreenCanvas.draw()
        #Affichage de l'intersection entre CDSearch et BLAST'
        def DisplayingIntersection(self,*args):
        
                #Si Intersection activée
                if self.intersectionDisplayState.get()==1:
                        self.ax1.set_xlabel(self.xlabel,color=self.couleurTexte)
                        self.ax1.set_ylabel(self.ylabel,color=self.couleurTexte)
                        #et la synténie activée
                        if self.syntenyDisplayState.get()==1:
                                #Nettoyer le graph
                                self.ax1.clear()
                                self.ax2.clear()
                                self.ax3.clear()
                                #Afficher l'intersection des deux 
                                X,Y=self.coordonneesXIntersection,self.coordonneesYIntersection
                                self.scatterIntersection=self.ax1.scatter(X,Y,c=self.couleurIntersection,s=0.2)
                                self.scatterIntersection.set_label('Intersection des deux méthodes de recherche')
                                #Afficher la synténie de l'intersection'
                                PlottingSyntenyIntersection(
                                                        self.coordonnesXSyntenyCroissantIntersection,
                                                        self.coordonnesYSyntenyCroissantIntersection,
                                                        self.coordonnesXSyntenyDecroissantIntersection,
                                                        self.coordonnesYSyntenyDecroissantIntersection,
                                                        self.statSyntenyIntersection,
                                                        self.normaldistribIntersection,
                                                        self.syntenyCoverIntersection,
                                                        self.couleurSyntenyIntersection,
                                                        self.styleSyntenyCroiss,
                                                        self.styleSyntenyDecroiss,
                                                        self.ax1,
                                                        self.ax2,
                                                        self.ax3
                                                        )

                                self.ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=2)
                                self.ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=2)
                         #et la synténie désactivée
                        if self.syntenyDisplayState.get()==0:
                                #Supprimer les deux 
                                self.ax1.clear()
                                #Afficher l'intersection des deux '
                                X,Y=self.coordonneesXIntersection,self.coordonneesYIntersection
                                self.scatterIntersection=self.ax1.scatter(X,Y,c=self.couleurIntersection,s=0.2)
                                self.scatterIntersection.set_label('Intersection de BLAST et CDSearch')

                #si l'intersection est désactivée '
                if self.intersectionDisplayState.get()==0:
                        #et la synténie est activée
                        if self.syntenyDisplayState.get()==1:
                                #Nettoyer le graph
                                self.ax1.clear()
                                self.ax2.clear()
                                self.ax3.clear()
                                #Afficher les deux 
                                #Afficher Blast
                                X,Y=self.coordonneesXBLAST,self.coordonneesYBLAST
                                self.scatterBLAST=self.ax1.scatter(X,Y,s=0.1,c=self.couleurBLASTSelected)
                                self.scatterBLAST.set_label("BLASTp")
                                #Afficher CDSearch
                                X,Y=self.coordonneesXCDSearch,self.coordonneesYCDSearch
                                self.scatterCDSearch=self.ax1.scatter(X,Y,s=0.1,c=self.couleurCDSearchSelected)
                                self.scatterCDSearch.set_label("CD Search")
                                #Afficher la syntenie des deux 
                                #Afficher la synténie de BLAST
                                PlottingSyntenyBLAST(
                                                        self.coordonnesXSyntenyCroissantBLAST,
                                                        self.coordonnesYSyntenyCroissantBLAST,
                                                        self.coordonnesXSyntenyDecroissantBLAST,
                                                        self.coordonnesYSyntenyDecroissantBLAST,
                                                        self.statSyntenyBLAST,
                                                        self.normaldistribBLAST,
                                                        self.syntenyCoverBLAST,
                                                        self.couleurSyntenyBLAST,
                                                        self.styleSyntenyCroiss,
                                                        self.styleSyntenyDecroiss,
                                                        self.ax1,
                                                        self.ax2,
                                                        self.ax3
                                                        )
                                #Afficher la synténie de CDSearch
                                PlottingSyntenyCDSearch(
                                                        self.coordonnesXSyntenyCroissantCDSearch,
                                                        self.coordonnesYSyntenyCroissantCDSearch,
                                                        self.coordonnesXSyntenyDecroissantCDSearch,
                                                        self.coordonnesYSyntenyDecroissantCDSearch,
                                                        self.statSyntenyCDSearch,
                                                        self.normaldistribCDSearch,
                                                        self.syntenyCoverCDSearch,
                                                        self.couleurSyntenyCDSearch,
                                                        self.styleSyntenyCroiss,
                                                        self.styleSyntenyDecroiss,
                                                        self.ax1,
                                                        self.ax2,
                                                        self.ax3
                                                        )
                                self.ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=2)
                                self.ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=2)
                        #et la synténie est désactivée
                        if self.syntenyDisplayState.get()==0:
                                #Supprimer l'intersection
                                self.ax1.clear()
                                #Afficher les deux 
                                #Afficher Blast
                                self.scatterBLAST=self.ax1.scatter(self.coordonneesXBLAST,self.coordonneesYBLAST,s=0.1,c=self.couleurBLASTSelected)
                                self.scatterBLAST.set_label("BLASTp")
                                #Afficher CDSearch
                                self.scatterCDSearch=self.ax1.scatter(self.coordonneesXCDSearch,self.coordonneesYCDSearch,s=0.1,c=self.couleurCDSearchSelected)
                                self.scatterCDSearch.set_label("CD Search")        
                self.ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=2)
                self.fullScreenCanvas.draw()
        def SelectionFonction(self,*args):
        
                #Si le bouton de sélection de fonction est bien appuyé 
                if self.fonctionHighlightState.get()==1:
                        #On nettoie le graph et on réaffiche la synténie
                        self.ax1.clear()
                        self.DisplayingSynteny()
                        #Récupérer la clé qui correspond au choix de la fonction 
                        if self.choixFonction.get() not in self.listeMenuFonctionFrontend: #On sort si rien n'a été choisi dans la liste'
                                return None
                                
                        choix=self.choixFonction.get()
                        indice=self.listeMenuFonctionFrontend.index(choix)
                        cle=self.listeMenuFonctionBackend[indice]
                        
                        #Si l'intersection n'est pas activée'
                        if self.intersectionDisplayState.get()==0:
                         
                                if len(self.coordonneesXBLAST)>0 and len(self.fonctionsXSyntenyCroissantBLAST)>0:
                                        #Ajouter au plot les points qui correspondent à la fonction sélectionnée dans la liste
                                        XCroissant,YCroissant=FiltreFonctionSyntenie(self.coordonnesXSyntenyCroissantBLAST,
                                                                                     self.coordonnesYSyntenyCroissantBLAST,
                                                                                     self.fonctionsXSyntenyCroissantBLAST,
                                                                                     self.fonctionsYSyntenyCroissantBLAST,
                                                                                     cle)
                                        XDecroissant,YDecroissant=FiltreFonctionSyntenie(self.coordonnesXSyntenyDecroissantBLAST,
                                                                                         self.coordonnesYSyntenyDecroissantBLAST,
                                                                                         self.fonctionsXSyntenyDecroissantBLAST,
                                                                                         self.fonctionsYSyntenyDecroissantBLAST,
                                                                                         cle)
                                        PlottingFonctionSynteny(XCroissant,
                                                                     YCroissant,
                                                                     XDecroissant,
                                                                     YDecroissant,
                                                                     self.couleurFonction,
                                                                     self.styleFonction,
                                                                     self.ax1,
                                                                     self.choixFonction.get()
                                                                )
                                if len(self.coordonneesXCDSearch)>0 and len(self.fonctionsXSyntenyCroissantCDSearch)>0:
                                        #Ajouter au plot les points qui correspondent à la fonction sélectionnée dans la liste
                                        XCroissant,YCroissant=FiltreFonctionSyntenie(self.coordonnesXSyntenyCroissantCDSearch,
                                                                                     self.coordonnesYSyntenyCroissantCDSearch,
                                                                                     self.fonctionsXSyntenyCroissantCDSearch,
                                                                                     self.fonctionsYSyntenyCroissantCDSearch,
                                                                                     cle)
                                        XDecroissant,YDecroissant=FiltreFonctionSyntenie(self.coordonnesXSyntenyDecroissantCDSearch,
                                                                                         self.coordonnesYSyntenyDecroissantCDSearch,
                                                                                         self.fonctionsXSyntenyDecroissantCDSearch,
                                                                                         self.fonctionsYSyntenyDecroissantCDSearch,
                                                                                         cle)
                                        PlottingFonctionSynteny(XCroissant,
                                                                     YCroissant,
                                                                     XDecroissant,
                                                                     YDecroissant,
                                                                     self.couleurFonction,
                                                                     self.styleFonction,
                                                                     self.ax1,
                                                                     self.choixFonction.get()
                                                                )
                        if self.intersectionDisplayState.get()==1:
                                if len(self.fonctionsXSyntenyCroissantIntersection)>0:
                                        XCroissant,YCroissant=FiltreFonctionSyntenie(self.coordonnesXSyntenyCroissantIntersection,
                                                                                     self.coordonnesYSyntenyCroissantIntersection,
                                                                                     self.fonctionsXSyntenyCroissantIntersection,
                                                                                     self.fonctionsYSyntenyCroissantIntersection,
                                                                                     cle)
                                        XDecroissant,YDecroissant=FiltreFonctionSyntenie(self.coordonnesXSyntenyDecroissantIntersection,
                                                                                         self.coordonnesYSyntenyDecroissantIntersection,
                                                                                         self.fonctionsXSyntenyDecroissantIntersection,
                                                                                         self.fonctionsYSyntenyDecroissantIntersection,
                                                                                         cle)
                                        PlottingFonctionSynteny(XCroissant,
                                                                     YCroissant,
                                                                     XDecroissant,
                                                                     YDecroissant,
                                                                     self.couleurFonction,
                                                                     self.styleFonction,
                                                                     self.ax1,
                                                                     self.choixFonction.get()
                                                                )
                                                                                

                if self.fonctionHighlightState.get()==0:
                        #On nettoie le graph et on réaffiche la synténie
                        self.ax1.clear()
                        self.DisplayingSynteny()
                self.ax1.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=2)
                self.ax2.legend(loc='upper center', bbox_to_anchor=(0.5, 1.15),ncol=2)
                self.fullScreenCanvas.draw()               
