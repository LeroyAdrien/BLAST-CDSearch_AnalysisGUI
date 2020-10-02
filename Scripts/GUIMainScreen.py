from GUIFullScreen import *
from GUIMoreInfo import *
from calculusFiltering import *
from fileparsing import *
import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure

class interface_principale:

        def __init__(self,path):
                données=FileExtractor(FileParser(path))
                self.proteome,self.CDSearch,self.BLAST,self.genome,self.COG,self.fun,self.noms2COG=données
                self.listeCoordonnesBLAST=[]
                self.listeCoordonnesCDSearch=[]
                self.listeCoordonnesSyntenyCroissant=[]
                self.listeCoordonnesSyntenyDecroissant=[]
                                #Fenêtre
                interface1=tk.Tk()
                interface1.title('Comparateur de Génomes')
                interface1.geometry('1200x600')

        #Couleur
                couleurTexte='#D8D8D8'
                couleurBackground='#2c2f33'
                
                self.couleurBLASTSelected='#0057e7'
                self.couleurCDSearchSelected='#a80000'
                
                couleurBLASTHighlight='#a8c6fa'
                couleurCDSearchHighlight='#ff822e'
                
                couleurBorder='#3D4849'
                couleurGeneralSelected='#9baaab'
                couleurGeneralHighlight='#d5dcdc'

        #Frames

                #Main Frame
                fenetre=tk.Frame(interface1,bg=couleurBackground)
                fenetre.pack(fill="both",expand=True)

                fenetre.grid_columnconfigure(0,weight=1, uniform="group1")
                fenetre.grid_columnconfigure(1,weight=1, uniform="group1")
                fenetre.grid_rowconfigure(0,weight=3, uniform="group2")
                fenetre.grid_rowconfigure(1,weight=6, uniform="group2")

        #Premiers plans de frames:
                
            #Frame du haut: prompt du BLAST à comparer
                frameHaut=tk.Frame(fenetre,bg=couleurBackground)
                frameHaut.grid(row=0,column=0,columnspan=2,sticky="nsew")
                #Paramètre du frame
                frameHaut.grid_rowconfigure(0,weight=3, uniform='group1')
                frameHaut.grid_rowconfigure(1,weight=3, uniform='group1')
                frameHaut.grid_rowconfigure(2,weight=3, uniform='group1')
                frameHaut.grid_rowconfigure(3,weight=2, uniform='group1')

                frameHaut.grid_columnconfigure(0,weight=1,uniform='group2')
                frameHaut.grid_columnconfigure(1,weight=1,uniform='group2')
                frameHaut.grid_columnconfigure(2,weight=1,uniform='group2')
                frameHaut.grid_columnconfigure(3,weight=1,uniform='group2')

            #Frame de gauche: prompt des options de filtrage des resultats
                frameGauche=tk.Frame(fenetre,bg=couleurBackground)
                frameGauche.grid(row=1,column=0,sticky="nsew")
                #Paramètre du frame
                frameGauche.grid_rowconfigure(0,weight=1, uniform='group1')
                frameGauche.grid_rowconfigure(1,weight=1, uniform='group1')
                frameGauche.grid_rowconfigure(2,weight=1, uniform='group1')
                frameGauche.grid_rowconfigure(3,weight=1, uniform='group1')
                frameGauche.grid_rowconfigure(4,weight=1, uniform='group1')
                #frameGauche.grid_rowconfigure(5,weight=1, uniform='group1')

                frameGauche.grid_columnconfigure(0,weight=2, uniform="group2")
                frameGauche.grid_columnconfigure(1,weight=1, uniform="group2")
                frameGauche.grid_columnconfigure(2,weight=2, uniform="group2")
                frameGauche.grid_columnconfigure(3,weight=1, uniform="group2")

            #Frame de droite: 2 Frames: Preview et GO/CLOSE
                frameDroite=tk.Frame(fenetre,bg=couleurBackground)#,highlightbackground="black",highlightthickness=1)
                frameDroite.grid(row=1,column=1,sticky="nsew")
                #Paramètre du frame
                frameDroite.grid_rowconfigure(0,weight=1, uniform='group1')
                frameDroite.grid_rowconfigure(1,weight=6, uniform='group1')
                frameDroite.grid_rowconfigure(2,weight=1, uniform='group1')
                frameDroite.grid_columnconfigure(0,weight=1)
                frameDroite.grid_columnconfigure(1,weight=1)

        #Choix des organismes et méthodes à comparer

        #Checkboxes des méthodes de recherche:
                #Checkbox de BLAST:
                self.BLASTPreviewState=tk.IntVar()
                self.BLASTPreviewState.set(0)
                BlastPreviewCheckbox=tk.Checkbutton(frameHaut,text="Rechercher selon BLAST     ",variable=self.BLASTPreviewState,indicatoron=0,fg=couleurTexte,bg=couleurBackground)
                BlastPreviewCheckbox.configure(highlightbackground=couleurBackground,selectcolor=self.couleurBLASTSelected,activebackground=couleurBLASTHighlight)
                BlastPreviewCheckbox.grid(row=1,column=0,columnspan=2,sticky='nsew')
                self.BLASTPreviewState.trace("w",self.SelectionBLAST)       
                
                #Checkbox de CD Search
                self.CDSearchPreviewState=tk.IntVar()
                self.CDSearchPreviewState.set(0)
                CDSearchPreviewCheckbox=tk.Checkbutton(frameHaut,text="Rechercher selon CDSearch",variable=self.CDSearchPreviewState,indicatoron=0,fg=couleurTexte,bg=couleurBackground)
                CDSearchPreviewCheckbox.configure(highlightbackground=couleurBackground,selectcolor=self.couleurCDSearchSelected,activebackground=couleurCDSearchHighlight)
                CDSearchPreviewCheckbox.grid(row=2,column=0,columnspan=2,sticky='nsew')
                self.CDSearchPreviewState.trace("w",self.SelectionCDSearchCOG)
                self.CDSearchPreviewState.trace("w",self.SelectionCDSearchFA)
                
                #RadioButton entre CDSearch par COG et FA
                self.CDSearchType=tk.StringVar()
                self.CDSearchType.set('COG')
                CDSearchRadioButton=tk.Radiobutton(frameHaut,variable=self.CDSearchType,text='COG',value='COG',indicatoron=0,fg=couleurTexte,bg=couleurBackground)
                CDSearchRadioButton.configure(highlightbackground=couleurBackground,selectcolor=self.couleurCDSearchSelected,activebackground=couleurCDSearchHighlight)
                CDSearchRadioButton.grid(row=3,column=0,sticky="nsew")
                CDSearchRadioButton=tk.Radiobutton(frameHaut,variable=self.CDSearchType,text='FA',value='FA',indicatoron=0,fg=couleurTexte,bg=couleurBackground)
                CDSearchRadioButton.configure(highlightbackground=couleurBackground,selectcolor=self.couleurCDSearchSelected,activebackground=couleurCDSearchHighlight)
                CDSearchRadioButton.grid(row=3,column=1,sticky="nsew")

                self.CDSearchType.trace("w",self.SelectionCDSearchCOG)
                self.CDSearchType.trace("w",self.SelectionCDSearchFA)

        #Label des menus déroulants:
                listeMenuLabel=tk.Label(frameHaut,text="Choix des organismes à comparer et de l'outil de recherche",font=('Work sans',-18,'bold'),fg=couleurTexte,bg=couleurBackground)
                listeMenuLabel.grid(row=0,column=0,columnspan=4,sticky='nsew')

        #Menu déroulant du blast:
            #Graphique:
                self.listeMenuBLASTBackend=[]
                self.listeMenuBLASTFrontend=[]
                
                for i in self.BLAST.keys():
                    self.listeMenuBLASTBackend.append(i)
                    self.listeMenuBLASTFrontend.append(self.genome[i[0]][0][1:-1]+" vs "+self.genome[i[1]][0][1:-1])
                self.choixBLAST=tk.StringVar(frameHaut)
                self.choixBLAST.set('')
                self.choixBLAST.trace("w",self.SelectionBLAST)
                menuBLAST=tk.OptionMenu(frameHaut,self.choixBLAST,*self.listeMenuBLASTFrontend)
                menuBLAST.grid(row=1,column=2,columnspan=2,sticky='nsew')
                menuBLAST['menu'].config(fg=couleurTexte,bg=couleurBackground)
                menuBLAST.config(fg=couleurTexte,bg=couleurBackground,highlightbackground=couleurBackground,activebackground=couleurBLASTHighlight)

        #Menus déroulants de CD-Search:

                #Organisme1:
                self.listeMenuCDSearchOrg1Frontend=[]
                self.listeMenuCDSearchOrg1Backend=[]
                for i in self.CDSearch.keys():
                        self.listeMenuCDSearchOrg1Frontend.append(self.genome[i][0][1:-1])
                        self.listeMenuCDSearchOrg1Backend.append(i)
                self.choixCDSearchOrga1=tk.StringVar(frameHaut)
                self.choixCDSearchOrga1.set('')
                menuCDSearchOrga1=tk.OptionMenu(frameHaut,self.choixCDSearchOrga1,*self.listeMenuCDSearchOrg1Frontend)
                menuCDSearchOrga1.grid(row=2,column=2,sticky='nsew')
                menuCDSearchOrga1['menu'].config(fg=couleurTexte,bg=couleurBackground)
                menuCDSearchOrga1.config(fg=couleurTexte,bg=couleurBackground,highlightbackground=couleurBackground,activebackground=couleurCDSearchHighlight)

                #Organisme2:
                self.listeMenuCDSearchOrg2Frontend=[]
                self.listeMenuCDSearchOrg2Backend=[]
                for i in self.CDSearch.keys():
                        self.listeMenuCDSearchOrg2Frontend.append(self.genome[i][0][1:-1])
                        self.listeMenuCDSearchOrg2Backend.append(i)
                self.choixCDSearchOrga2=tk.StringVar(frameHaut)
                self.choixCDSearchOrga2.set('')
                menuCDSearchOrga2=tk.OptionMenu(frameHaut,self.choixCDSearchOrga2,*self.listeMenuCDSearchOrg2Frontend)
                menuCDSearchOrga2['menu'].config(fg=couleurTexte,bg=couleurBackground)
                menuCDSearchOrga2.config(fg=couleurTexte,bg=couleurBackground,highlightbackground=couleurBackground,activebackground=couleurCDSearchHighlight)
                menuCDSearchOrga2.grid(row=2,column=3,sticky='nsew')
                
                #Comportement

                self.choixCDSearchOrga1.trace("w",self.SelectionCDSearchCOG)
                self.choixCDSearchOrga2.trace("w",self.SelectionCDSearchCOG)
                self.choixCDSearchOrga1.trace("w",self.SelectionCDSearchFA)
                self.choixCDSearchOrga2.trace("w",self.SelectionCDSearchFA)
                

                #Plus d'infos
                moreInfoButton=tk.Button(frameHaut,text="Organismes analysés",command=self.NewWindowMoreInfo,fg=couleurTexte,bg=couleurBackground)
                moreInfoButton.configure(highlightbackground=couleurBackground,activebackground=couleurGeneralHighlight)
                moreInfoButton.grid(row=3,column=2,columnspan=2,sticky="nsew")

        #Options de filtrage
            #Labels
                #Label explicatif:
                evalueLabel=tk.Label(frameGauche,text='Nom du filtre:',fg=couleurTexte,bg=couleurBackground,font=('Work sans',-15,'bold'))
                evalueLabel.grid(row=0,column=0,sticky='nsew')
                
                #Label Evalue Blast:
                evalueLabel=tk.Label(frameGauche,text='E-Value Blast',fg=couleurTexte,bg=couleurBackground,font=('Work sans',-12,'bold'))
                evalueLabel.grid(row=1,column=0,sticky='nsew')
                #Label Evalue CD-Search:
                evalueLabel=tk.Label(frameGauche,text='E-Value CD-Search',fg=couleurTexte,bg=couleurBackground,font=('Work sans',-12,'bold'))
                evalueLabel.grid(row=4,column=0,sticky='nsew')
                #Label Pourcentage d'identité:
                evalueLabel=tk.Label(frameGauche,text="% d'identité",fg=couleurTexte,bg=couleurBackground,font=('Work sans',-12,'bold'))
                evalueLabel.grid(row=3,column=0,sticky='nsew')
                #Label Couverture du hit:
                evalueLabel=tk.Label(frameGauche,text='Couverture du hit',fg=couleurTexte,bg=couleurBackground,font=('Work sans',-12,'bold'))
                evalueLabel.grid(row=2,column=0,sticky='nsew')

            #Sliders 
                #Label explicatif:
                evalueLabel=tk.Label(frameGauche,text='Force du filtre:',fg=couleurTexte,bg=couleurBackground,font=('Work sans',-15,'bold'))
                evalueLabel.grid(row=0,column=1,columnspan=3,sticky='nsew')

                #SLIDERS DE BLAST
                #Variables des sliders: self.evalueFilterStateBLAST,self.pidentFilterState,self.hitCoverState
                
                #Slider Evalue Blast
                        #Evalue seuil
                self.evalueFilterStateBLAST=tk.IntVar()
                self.evalueFilterStateBLAST.trace('w',self.FilterSelectionBLAST)
                        #Slider Evalue
                evalueBlastSlider=tk.Scale(frameGauche,from_=0, to=-200,orient='horizontal',variable=self.evalueFilterStateBLAST,tickinterval=50)
                evalueBlastSlider.config(bg=couleurBackground,highlightbackground=couleurBackground,highlightcolor=couleurBackground,troughcolor=self.couleurBLASTSelected,fg=couleurTexte)
                evalueBlastSlider.grid(row=1,column=1,columnspan=2,sticky='nsew')

                #Slider pourcentage d'identité
                        #piden seuil
                self.pidentFilterState=tk.IntVar()
                self.pidentFilterState.trace('w',self.FilterSelectionBLAST)
                        #Slider pident
                pidentSlider=tk.Scale(frameGauche,from_=0, to=100,orient='horizontal',variable=self.pidentFilterState,tickinterval=25)
                pidentSlider.config(bg=couleurBackground,highlightbackground=couleurBackground,highlightcolor=couleurBackground,troughcolor=self.couleurBLASTSelected,fg=couleurTexte)
                pidentSlider.grid(row=3,column=1,columnspan=2,sticky="nsew")

                #Slider couverture du hit
                        #couverture seuil
                self.hitCoverFilterState=tk.IntVar()
                self.hitCoverFilterState.trace('w',self.FilterSelectionBLAST)
                        #Slider couverture
                hitCoverSlider=tk.Scale(frameGauche,from_=0, to=100,orient='horizontal',variable=self.hitCoverFilterState,tickinterval=25)
                hitCoverSlider.config(bg=couleurBackground,highlightbackground=couleurBackground,highlightcolor=couleurBackground,troughcolor=self.couleurBLASTSelected,fg=couleurTexte)
                hitCoverSlider.grid(row=2,column=1,columnspan=2,sticky='nsew')

                #SLIDER DE CDSEARCH
                
                #Slider Evalue CD-Search
                        #Evalue seuil
                self.evalueFilterStateCDSearch=tk.IntVar()
                self.evalueFilterStateCDSearch.trace("w",self.FilterSelectionCDSearch)
                evalueCDSearchSlider=tk.Scale(frameGauche,from_=0, to=-200,orient='horizontal',variable=self.evalueFilterStateCDSearch,tickinterval=50)
                evalueCDSearchSlider.config(bg=couleurBackground,highlightbackground=couleurBackground,highlightcolor=couleurBackground,troughcolor=self.couleurCDSearchSelected,fg=couleurTexte)
                evalueCDSearchSlider.grid(row=4,column=1,columnspan=2,sticky='nsew')
                
            #Entry
                #Entry Evalue Blast
                self.evalueFilterStateBLASTEntry=tk.StringVar()
                self.evalueBlastEntry=tk.Entry(frameGauche,width=6,textvariable=self.evalueFilterStateBLASTEntry,fg=couleurBackground,bg=couleurTexte)
                self.evalueBlastEntry.grid(row=1,column=3)
                self.evalueBlastEntry.bind('<Return>',self.AdjustSliderEntryEvalueBLAST)
                #Entry couverture du hit
                self.hitCoverFilterStateEntry=tk.StringVar()
                self.hitCoverBlastEntry=tk.Entry(frameGauche,width=6,textvariable=self.hitCoverFilterStateEntry,fg=couleurBackground,bg=couleurTexte)
                self.hitCoverBlastEntry.grid(row=2,column=3)
                self.hitCoverBlastEntry.bind('<Return>',self.AdjustSliderEntryhitCoverBLAST)
                #Entry pourcentage d'identité
                self.pidentFilterStateEntry=tk.StringVar()
                self.pidentBlastEntry=tk.Entry(frameGauche,width=6,textvariable=self.pidentFilterStateEntry,fg=couleurBackground,bg=couleurTexte)
                self.pidentBlastEntry.grid(row=3,column=3)
                self.pidentBlastEntry.bind('<Return>',self.AdjustSliderEntrypidentBLAST)
                #Entry Evalue CDSearch
                self.evalueFilterStateCDSearchEntry=tk.StringVar()
                self.evalueCDSearchEntry=tk.Entry(frameGauche,width=6,textvariable=self.evalueFilterStateCDSearchEntry,fg=couleurBackground,bg=couleurTexte)
                self.evalueCDSearchEntry.grid(row=4,column=3)
                self.evalueCDSearchEntry.bind('<Return>',self.AdjustSliderEntryEvalueCDSearch)

            #Preview et voir en grand
            
            #Label de la preview:
            #previewLabel
                previewLabel=tk.Label(frameDroite,text='Preview du dot plot',fg=couleurTexte,bg=couleurBackground,font=('Work sans',-15,'bold'))
                previewLabel.grid(row=0,column=0,columnspan=2,sticky='nsew')
            #Preview:

                fig = plt.Figure(figsize=(5, 5), dpi=100)
                plt.ion()
                fig.patch.set_facecolor(couleurBackground)
                #t = np.arange(0, 3, .01)
                self.coordonneesXBLAST=[]
                self.coordonneesYBLAST=[]
                self.coordonneesXCDSearch=[]
                self.coordonneesYCDSearch=[]
                
                self.plot1=fig.add_subplot(111)
                self.plot1.set_facecolor(couleurTexte)
                self.plot1.tick_params(axis='x', colors=couleurTexte)
                self.plot1.tick_params(axis='y', colors=couleurTexte)
                self.scatterBLAST=self.plot1.scatter(self.coordonneesXBLAST,self.coordonneesYBLAST,s=0.1,c=self.couleurBLASTSelected)
                self.scatterCDSearch=self.plot1.scatter(self.coordonneesXCDSearch,self.coordonneesYCDSearch,s=0.1,c=self.couleurCDSearchSelected)
                self.previewCanvas = FigureCanvasTkAgg(fig, master=frameDroite)  # A tk.DrawingArea.
                self.previewCanvas.draw()
                self.previewCanvas.get_tk_widget().grid(row=1,column=0,columnspan=2,sticky="nsew")
                
            #Bouton Agrandir:
                fullsizeButton=tk.Button(frameDroite,text='Dot Plot',command=self.NewWindowFullscreen,fg=couleurTexte,bg=couleurBackground)
                fullsizeButton.configure(highlightbackground=couleurBackground,activebackground=couleurGeneralHighlight)
                fullsizeButton.grid(row=2,column=0)
            #Bouton Quitter
                closeButton=tk.Button(frameDroite,text='Quitter',command=interface1.destroy,fg=couleurTexte,bg=couleurBackground)
                closeButton.configure(highlightbackground=couleurBackground,activebackground=couleurGeneralHighlight)
                closeButton.grid(row=2,column=1)

                    #mainloop

                interface1.mainloop()

        #Fonctions de la classe
        #Creation de l'interface fullScreen'
        def NewWindowFullscreen(self):
                #Envoi des coordonnées de la preview
                coordonneesXBLAST=self.coordonneesXBLAST
                coordonneesYBLAST=self.coordonneesYBLAST
                coordonneesXCDSearch=self.coordonneesXCDSearch
                coordonneesYCDSearch=self.coordonneesYCDSearch
                titreBLAST=''
                titreCDSearch=''
                GCABLAST=''
                GCACDSearch=''
                
                #Envoi des Labels des organismes pour Le titre et GCA pour le calcul de cover
                if len(self.coordonneesXBLAST)!=0:
                        #Label
                        titreBLAST=self.choixBLAST.get().split("vs")
                        #GCA
                        indiceBackend=self.listeMenuBLASTFrontend.index(self.choixBLAST.get())
                        GCABLAST=self.listeMenuBLASTBackend[indiceBackend]
                        
                if len(self.coordonneesXCDSearch)!=0:
                        #Label
                        titreCDSearch=[self.choixCDSearchOrga1.get(),self.choixCDSearchOrga2.get()]
                        OrganismesCDSearch=[]
                        #GCA
                        indiceCleOrga1=self.listeMenuCDSearchOrg1Frontend.index(self.choixCDSearchOrga1.get())
                        indiceCleOrga2=self.listeMenuCDSearchOrg2Frontend.index(self.choixCDSearchOrga2.get())
                        cleOrga1=self.listeMenuCDSearchOrg1Backend[indiceCleOrga1]
                        cleOrga2=self.listeMenuCDSearchOrg2Backend[indiceCleOrga2]
                        GCACDSearch=(cleOrga1,cleOrga2)
                #Creation de l'interface'
                a=interface_fullscreen(
                                       coordonneesXBLAST,
                                       coordonneesYBLAST,
                                       coordonneesXCDSearch,
                                       coordonneesYCDSearch,
                                       titreBLAST,
                                       titreCDSearch,
                                       self.proteome,
                                       self.fun,
                                       self.noms2COG,
                                       self.COG,
                                       self.BLAST,
                                       self.CDSearch,
                                       GCABLAST,
                                       GCACDSearch
                                       )

        #creation de l'interface plus d'infos
        def NewWindowMoreInfo(self):
                infosBLASTOrga1=[]
                infosBLASTOrga2=[]
                infosCDSearchOrga1=[]
                infosCDSearchOrga2=[]
                #On récupère la clé qui correspond au choix dans la liste BLAST
                if self.choixBLAST.get()!='':
                        indiceBackend=self.listeMenuBLASTFrontend.index(self.choixBLAST.get())
                        cleBLAST=self.listeMenuBLASTBackend[indiceBackend]
                        infosBLASTOrga1=[self.genome[cleBLAST[0]][x] for x in [0,3,4,5]]
                        infosBLASTOrga1=self.genome[cleBLAST[0]][1][1:-1].split(';')+infosBLASTOrga1
                        infosBLASTOrga1[3]=infosBLASTOrga1[3][1:-1]
                        infosBLASTOrga2=[self.genome[cleBLAST[1]][x] for x in [0,3,4,5]]
                        infosBLASTOrga2=self.genome[cleBLAST[1]][1][1:-1].split(';')+infosBLASTOrga2
                        infosBLASTOrga2[3]=infosBLASTOrga2[3][1:-1]
                if self.choixCDSearchOrga1.get()!='':
                        indiceCleOrga1=self.listeMenuCDSearchOrg1Frontend.index(self.choixCDSearchOrga1.get())
                        cleOrga1CDSearch=self.listeMenuCDSearchOrg1Backend[indiceCleOrga1]
                        infosCDSearchOrga1=[self.genome[cleOrga1CDSearch][x] for x in [0,3,4,5]]
                        infosCDSearchOrga1=self.genome[cleOrga1CDSearch][1][1:-1].split(';')+infosCDSearchOrga1
                        infosCDSearchOrga1[3]=infosCDSearchOrga1[3][1:-1]
                if self.choixCDSearchOrga2.get()!='':
                        indiceCleOrga2=self.listeMenuCDSearchOrg2Frontend.index(self.choixCDSearchOrga2.get())
                        cleOrga2CDSearch=self.listeMenuCDSearchOrg2Backend[indiceCleOrga2]
                        infosCDSearchOrga2=[self.genome[cleOrga2CDSearch][x] for x in [0,3,4,5]]
                        infosCDSearchOrga2=self.genome[cleOrga2CDSearch][1][1:-1].split(';')+infosCDSearchOrga2
                        infosCDSearchOrga2[3]=infosCDSearchOrga2[3][1:-1]
                b=interface_moreInfo(infosBLASTOrga1,infosBLASTOrga2,infosCDSearchOrga1,infosCDSearchOrga2)

        #Creation du Plot de BLAST
        def SelectionBLAST(self,triggerVariable,blank,triggerMode):
                if self.BLASTPreviewState.get()==0:
                        self.coordonneesXBLAST=[]
                        self.coordonneesYBLAST=[]
                        self.plot1.clear()
                        self.scatterBLAST.remove()
                        self.scatterBLAST=self.plot1.scatter(self.coordonneesXBLAST,self.coordonneesYBLAST,s=0.1,c=self.couleurBLASTSelected)
                        self.scatterCDSearch=self.plot1.scatter(self.coordonneesXCDSearch,self.coordonneesYCDSearch,s=0.1,c=self.couleurCDSearchSelected)
                        self.previewCanvas.draw()
                        
                elif self.BLASTPreviewState.get()==1 and self.choixBLAST.get()!='':
                        indiceBackend=self.listeMenuBLASTFrontend.index(self.choixBLAST.get())
                        cle=self.listeMenuBLASTBackend[indiceBackend]
                        if cle in self.BLAST.keys():
                                self.listeCoordonnesBLAST=ResultatsBlast(self.BLAST,self.proteome,cle) 
                                self.listeCoordonnesBLAST=np.array(self.listeCoordonnesBLAST)
                                self.scatterBLAST.remove()
                                self.coordonneesXBLAST=self.listeCoordonnesBLAST[:,0]
                                self.coordonneesYBLAST=self.listeCoordonnesBLAST[:,1]
                                self.plot1.clear()
                                self.scatterBLAST=self.plot1.scatter(self.coordonneesXBLAST,self.coordonneesYBLAST,s=0.1,c=self.couleurBLASTSelected)
                                self.scatterCDSearch=self.plot1.scatter(self.coordonneesXCDSearch,self.coordonneesYCDSearch,s=0.1,c=self.couleurCDSearchSelected)
                                self.previewCanvas.draw()
                                
        #Creation du plot de CDSearch en COG
        def SelectionCDSearchCOG(self,triggerVariable,blank,triggerMode):
                if self.CDSearchPreviewState.get()==0:
                        self.coordonneesXCDSearch=[]
                        self.coordonneesYCDSearch=[]
                        self.scatterCDSearch.remove()
                        self.plot1.clear()
                        self.scatterBLAST=self.plot1.scatter(self.coordonneesXBLAST,self.coordonneesYBLAST,s=0.1,c=self.couleurBLASTSelected)
                        self.scatterCDSearch=self.plot1.scatter(self.coordonneesXCDSearch,self.coordonneesYCDSearch,s=0.1,c=self.couleurCDSearchSelected)
                        self.previewCanvas.draw()

                elif self.choixCDSearchOrga1.get()!='' and self.choixCDSearchOrga2.get()!='' and self.CDSearchType.get()=='COG' and self.CDSearchPreviewState.get()==1:
                        indiceCleOrga1=self.listeMenuCDSearchOrg1Frontend.index(self.choixCDSearchOrga1.get())
                        indiceCleOrga2=self.listeMenuCDSearchOrg2Frontend.index(self.choixCDSearchOrga2.get())
                        cleOrga1=self.listeMenuCDSearchOrg1Backend[indiceCleOrga1]
                        cleOrga2=self.listeMenuCDSearchOrg2Backend[indiceCleOrga2]
                        if cleOrga1 in self.CDSearch.keys() and cleOrga2 in self.CDSearch.keys():
                                self.listeCoordonnesCDSearch=ResultatsCDSearchCOG(self.CDSearch,self.proteome,cleOrga1,cleOrga2)
                                self.listeCoordonnesCDSearch=np.array(self.listeCoordonnesCDSearch)
                                self.scatterCDSearch.remove()
                                self.coordonneesXCDSearch=self.listeCoordonnesCDSearch[:,0]
                                self.coordonneesYCDSearch=self.listeCoordonnesCDSearch[:,1]
                                self.plot1.clear()
                                self.scatterBLAST=self.plot1.scatter(self.coordonneesXBLAST,self.coordonneesYBLAST,s=0.1,c=self.couleurBLASTSelected)
                                self.scatterCDSearch=self.plot1.scatter(self.coordonneesXCDSearch,self.coordonneesYCDSearch,s=0.1,c=self.couleurCDSearchSelected)
                                self.previewCanvas.draw()
                                
        #Creation du plot de CDSearch en Fonctions de protéines
        def SelectionCDSearchFA(self,triggerVariable,blank,triggerMode):
                if self.CDSearchPreviewState.get()==0:
                        self.coordonneesXCDSearch=[]
                        self.coordonneesYCDSearch=[]
                        self.scatterCDSearch.remove()
                        self.plot1.clear()
                        self.scatterBLAST=self.plot1.scatter(self.coordonneesXBLAST,self.coordonneesYBLAST,s=0.1,c=self.couleurBLASTSelected)
                        self.scatterCDSearch=self.plot1.scatter(self.coordonneesXCDSearch,self.coordonneesYCDSearch,s=0.1,c=self.couleurCDSearchSelected)
                        self.previewCanvas.draw()

                elif self.choixCDSearchOrga1.get()!='' and self.choixCDSearchOrga2.get()!='' and self.CDSearchType.get()=='FA' and self.CDSearchPreviewState.get()==1:
                        indiceCleOrga1=self.listeMenuCDSearchOrg1Frontend.index(self.choixCDSearchOrga1.get())
                        indiceCleOrga2=self.listeMenuCDSearchOrg2Frontend.index(self.choixCDSearchOrga2.get())
                        cleOrga1=self.listeMenuCDSearchOrg1Backend[indiceCleOrga1]
                        cleOrga2=self.listeMenuCDSearchOrg2Backend[indiceCleOrga2]
                        if cleOrga1 in self.CDSearch.keys() and cleOrga2 in self.CDSearch.keys():
                                self.listeCoordonnesCDSearch=ResultatsCDSearchFA(self.CDSearch,self.proteome,self.COG,cleOrga1,cleOrga2)
                                self.listeCoordonnesCDSearch=np.array(self.listeCoordonnesCDSearch)
                                self.scatterCDSearch.remove()
                                self.coordonneesXCDSearch=self.listeCoordonnesCDSearch[:,0]
                                self.coordonneesYCDSearch=self.listeCoordonnesCDSearch[:,1]
                                self.plot1.clear()
                                self.scatterBLAST=self.plot1.scatter(self.coordonneesXBLAST,self.coordonneesYBLAST,s=0.1,c=self.couleurBLASTSelected)
                                self.scatterCDSearch=self.plot1.scatter(self.coordonneesXCDSearch,self.coordonneesYCDSearch,s=0.1,c=self.couleurCDSearchSelected)
                                self.previewCanvas.draw()
                                
        #Filtrage de la selection de BLAST
        def FilterSelectionBLAST(self,triggerVariable,blank,triggerMode):
                if self.BLASTPreviewState.get()==1:
                        self.evalueBlastEntry.delete(0,'end')
                        self.evalueBlastEntry.insert(0,10**self.evalueFilterStateBLAST.get())
                        self.hitCoverBlastEntry.delete(0,'end')
                        self.hitCoverBlastEntry.insert(0,self.hitCoverFilterState.get())
                        self.pidentBlastEntry.delete(0,'end')
                        self.pidentBlastEntry.insert(0,self.pidentFilterState.get())
                        if self.listeCoordonnesBLAST!=[]:
                                coordonnes=self.listeCoordonnesBLAST
                                pident=self.pidentFilterState.get()
                                hitcover=self.hitCoverFilterState.get()
                                evalue=10**self.evalueFilterStateBLAST.get()
                                listeCoordonnesBLASTFiltre=FilterBLAST(coordonnes,pident,hitcover,hitcover,evalue)
                                listeCoordonnesBLASTFiltre=np.array(listeCoordonnesBLASTFiltre)
                                
                                if len(listeCoordonnesBLASTFiltre)>0:
                                        self.coordonneesXBLAST=listeCoordonnesBLASTFiltre[:,0]
                                        self.coordonneesYBLAST=listeCoordonnesBLASTFiltre[:,1]

                                else:
                                        self.coordonneesXBLAST=[]
                                        self.coordonneesYBLAST=[]

                                self.scatterBLAST.remove()
                                self.plot1.clear()
                                self.scatterBLAST=self.plot1.scatter(self.coordonneesXBLAST,self.coordonneesYBLAST,s=0.1,c=self.couleurBLASTSelected)
                                self.scatterCDSearch=self.plot1.scatter(self.coordonneesXCDSearch,self.coordonneesYCDSearch,s=0.1,c=self.couleurCDSearchSelected)
                                self.previewCanvas.draw()
                                
        #Filtrage de la selection de CDSearch
        def FilterSelectionCDSearch(self,triggerVariable,blank,triggerMode):
                if self.CDSearchPreviewState.get()==1:
                        self.evalueCDSearchEntry.delete(0,'end')
                        self.evalueCDSearchEntry.insert(0,10**self.evalueFilterStateCDSearch.get())
                        if self.listeCoordonnesCDSearch!=[]:
                                coordonnes=self.listeCoordonnesCDSearch
                                evalue=10**self.evalueFilterStateCDSearch.get()
                                if self.CDSearchType.get()=='COG':
                                        listeCoordonnesCDSearchFiltre=FilterCDSearchCOG(coordonnes,evalue,evalue)
                                elif self.CDSearchType.get()=='FA':
                                        listeCoordonnesCDSearchFiltre=FilterCDSearchFA(coordonnes,evalue,evalue)
                                else:
                                        pass

                                listeCoordonnesCDSearchFiltre=np.array(listeCoordonnesCDSearchFiltre)
                                if len(listeCoordonnesCDSearchFiltre)>0:
                                        self.coordonneesXCDSearch=listeCoordonnesCDSearchFiltre[:,0]
                                        self.coordonneesYCDSearch=listeCoordonnesCDSearchFiltre[:,1]

                                else:
                                        self.coordonneesXCDSearch=[]
                                        self.coordonneesYCDSearch=[]

                                self.scatterCDSearch.remove()
                                self.plot1.clear()
                                self.scatterBLAST=self.plot1.scatter(self.coordonneesXBLAST,self.coordonneesYBLAST,s=0.1,c=self.couleurBLASTSelected)
                                self.scatterCDSearch=self.plot1.scatter(self.coordonneesXCDSearch,self.coordonneesYCDSearch,s=0.1,c=self.couleurCDSearchSelected)
                                self.previewCanvas.draw()     
                                   
        #Synchronisation des sliders et des entry 
        def AdjustSliderEntryEvalueBLAST(self,*args):
                try: 
                        self.evalueFilterStateBLAST.set(np.log10(float(self.evalueFilterStateBLASTEntry.get())))
                except ValueError:
                        pass
        def AdjustSliderEntrypidentBLAST(self,*args):
                try: 
                        self.pidentFilterState.set(float(self.pidentFilterStateEntry.get()))
                except ValueError:
                        pass
        def AdjustSliderEntryhitCoverBLAST(self,*args):
                try: 
                        self.hitCoverFilterState.set(float(self.hitCoverFilterStateEntry.get()))
                except ValueError:
                        pass
        def AdjustSliderEntryEvalueCDSearch(self,*args):
                try: 
                        self.evalueFilterStateCDSearch.set(np.log10(float(self.evalueFilterStateCDSearchEntry.get())))
                except ValueError:
                        pass


