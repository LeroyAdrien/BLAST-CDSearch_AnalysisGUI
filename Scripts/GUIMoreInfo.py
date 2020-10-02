import tkinter as tk

class interface_moreInfo:

        def __init__(self,infosBLASTOrga1,infosBLASTOrga2,infosCDSearchOrga1,infosCDSearchOrga2):
        
                interfaceMoreInfo=tk.Tk()
                interfaceMoreInfo.title("Plus d'infos")
                interfaceMoreInfo.geometry('1300x600')

                couleurBackground='#2c2f33'
                couleurTexte='#D8D8D8'
                couleursCDSearch=['#360808','#421515',
                                  '#5c0e0e','#7c0404',
                                  '#a80000',couleurBackground,
                                  couleurTexte]

                couleursBLAST=['#022a50','#04396d',
                               '#064887','#0959a7',
                               '#2570b8',couleurBackground,
                               couleurTexte]

               #Frame unique de la fenêtre
                self.mainFrame=tk.Frame(interfaceMoreInfo,bg=couleurBackground)
                self.mainFrame.pack(fill="both",expand=True)

                #Appel de la fonction pour creer les boites
                        #Si Un choix dans blast a été sélectionné
                if infosBLASTOrga1:
                        self.mainFrame.grid_rowconfigure(0,weight=1, uniform="group2")
                        self.mainFrame.grid_columnconfigure(0,weight=1, uniform="group1")
                        self.Label_Layout(infosBLASTOrga1,0,0,couleursBLAST)
                if infosBLASTOrga2:
                        self.mainFrame.grid_rowconfigure(0,weight=1, uniform="group2")
                        self.mainFrame.grid_columnconfigure(1,weight=1, uniform="group1")
                        self.Label_Layout(infosBLASTOrga2,0,1,couleursBLAST)
                if infosCDSearchOrga1:
                        self.mainFrame.grid_columnconfigure(0,weight=1, uniform="group1")
                        self.mainFrame.grid_rowconfigure(1,weight=1, uniform="group2")
                        self.Label_Layout(infosCDSearchOrga1,1,0,couleursCDSearch)
                if infosCDSearchOrga2:
                        self.mainFrame.grid_columnconfigure(1,weight=1, uniform="group1")
                        self.mainFrame.grid_rowconfigure(1,weight=1, uniform="group2")
                        self.Label_Layout(infosCDSearchOrga2,1,1,couleursCDSearch)

                interfaceMoreInfo.mainloop()
        # Creer deux boites, une "matriochka"" pour montrer les embranchements, une autre pour donner les caractérsitiques du génome de l'organisme
        def Label_Layout(self,infos,ROW,COLUMN,color):
                if infos:
                        #Frame de la méthode
                        method=tk.Frame(self.mainFrame,bg='green')
                        method.grid(row=ROW,column=COLUMN,sticky='nsew',padx=5,pady=5)
                        method.grid_rowconfigure(0,weight=2,uniform='group1')
                        method.grid_rowconfigure(1,weight=1,uniform='group1')
                        method.grid_columnconfigure(0,weight=1, uniform="group2")
                        
                        
                        #Frames des espèces 
                        #Règne
                        regne=tk.Frame(method,bg=color[0])
                        regne.grid(row=0,column=0,sticky='nsew')
                        
                        regne.grid_rowconfigure(0,weight=1,uniform='group1')
                        regne.grid_rowconfigure(1,weight=10,uniform='group1')
                        regne.grid_rowconfigure(2,weight=3,uniform='group1')
                        
                        regne.grid_columnconfigure(0,weight=1,uniform='group2')
                        regne.grid_columnconfigure(1,weight=10,uniform='group2')
                        regne.grid_columnconfigure(2,weight=1,uniform='group2')
                        
                        #Embranchement
                        embranchement=tk.Frame(regne,bg=color[1])
                        embranchement.grid(row=1,column=1,sticky='nsew')
                        
                        embranchement.grid_rowconfigure(0,weight=1,uniform='group1')
                        embranchement.grid_rowconfigure(1,weight=10,uniform='group1')
                        embranchement.grid_rowconfigure(2,weight=3,uniform='group1')
                        
                        embranchement.grid_columnconfigure(0,weight=1,uniform='group2')
                        embranchement.grid_columnconfigure(1,weight=10,uniform='group2')
                        embranchement.grid_columnconfigure(2,weight=1,uniform='group2')
                        
                        #Ordre
                        ordre=tk.Frame(embranchement,bg=color[2])
                        ordre.grid(row=1,column=1,sticky='nsew')
                        
                        ordre.grid_rowconfigure(0,weight=1,uniform='group1')
                        ordre.grid_rowconfigure(1,weight=5,uniform='group1')
                        ordre.grid_rowconfigure(2,weight=3,uniform='group1')
                        
                        ordre.grid_columnconfigure(0,weight=1,uniform='group2')
                        ordre.grid_columnconfigure(1,weight=10,uniform='group2')
                        ordre.grid_columnconfigure(2,weight=1,uniform='group2')
                        
                        #Espece
                        espece=tk.Frame(ordre,bg=color[3])
                        espece.grid(row=1,column=1,sticky='nsew')
                        
                        espece.grid_rowconfigure(0,weight=1,uniform='group1')
                        espece.grid_rowconfigure(1,weight=3,uniform='group1')
                        espece.grid_rowconfigure(2,weight=1,uniform='group1')
                        
                        espece.grid_columnconfigure(0,weight=1,uniform='group2')
                        espece.grid_columnconfigure(1,weight=5,uniform='group2')
                        espece.grid_columnconfigure(2,weight=1,uniform='group2')
                        
                        #Labels 
                        
                        regneLabel=tk.Label(regne,text=infos[0],bg=color[0],fg=color[6])
                        regneLabel.grid(column=1,row=2)

                        embranchementLabel=tk.Label(embranchement,text=infos[1],bg=color[1],fg=color[6])
                        embranchementLabel.grid(column=1,row=2)
                        
                        
                        ordreLabel=tk.Label(ordre,text=infos[2],bg=color[2],fg=color[6])
                        ordreLabel.grid(column=1,row=2)
                        

                        especeLabel=tk.Label(espece,text=infos[3],bg=color[3],fg=color[6])
                        especeLabel.grid(column=1,row=1)
                        
                        #Frames de infos sur l'organisme'
                        orga=tk.Frame(method,bg=color[5])
                        orga.grid(row=1,column=0,sticky="nsew")
                        orga.grid_columnconfigure(0,weight=1,uniform='group1')
                        orga.grid_columnconfigure(1,weight=1,uniform='group1')
                        orga.grid_columnconfigure(2,weight=1,uniform='group1')
                        
                        orga.grid_rowconfigure(0,weight=1,uniform='group2')
                        
                        tailleGenomeFrame=tk.Frame(orga,bg=color[4])
                        tailleGenomeFrame.grid(row=0,column=0,sticky='nsew',pady=10,padx=10)
                        
                        tailleGenomeFrame.grid_columnconfigure(0,weight=1,uniform='groupe1')
                        tailleGenomeFrame.grid_rowconfigure(0,weight=1,uniform='group2')
                        tailleGenomeFrame.grid_rowconfigure(1,weight=2,uniform='group2')
                        
                        pourcentageGCFrame=tk.Frame(orga,bg=color[4])
                        pourcentageGCFrame.grid(row=0,column=1,sticky='nsew',pady=10,padx=10)
                        
                        pourcentageGCFrame.grid_columnconfigure(0,weight=1,uniform='groupe1')
                        pourcentageGCFrame.grid_rowconfigure(0,weight=1,uniform='group2')
                        pourcentageGCFrame.grid_rowconfigure(1,weight=2,uniform='group2')
                        
                        nbrCDSFrame=tk.Frame(orga,bg=color[4])
                        nbrCDSFrame.grid(row=0,column=2,sticky='nsew',pady=10,padx=10)
                        
                        nbrCDSFrame.grid_columnconfigure(0,weight=1,uniform='groupe1')
                        nbrCDSFrame.grid_rowconfigure(0,weight=1,uniform='group2')
                        nbrCDSFrame.grid_rowconfigure(1,weight=2,uniform='group2')
                        
                        #Labels
                                #Label Taille génome
                        infoTaille=tk.Label(tailleGenomeFrame,text="Taille en Mb:",bg=color[4],fg=color[6])
                        infoTaille.grid(row=0,column=0,sticky='nsew')
                        
                        taille=tk.Label(tailleGenomeFrame,text=infos[4],bg=color[4],fg=color[6])
                        taille.grid(row=1,column=0,sticky='nsew')
                                #Label Pourcentage GC
                        infoGC=tk.Label(pourcentageGCFrame,text='Pourcentage de GC:',bg=color[4],fg=color[6])
                        infoGC.grid(row=0,column=0,sticky='nsew')
                        
                        GC=tk.Label(pourcentageGCFrame,text=infos[5],bg=color[4],fg=color[6])
                        GC.grid(row=1,column=0,sticky='nsew')
                                #Nombre de CDS
                        infoCDS=tk.Label(nbrCDSFrame,text="Nombre de CDS:",bg=color[4],fg=color[6])
                        infoCDS.grid(row=0,column=0,sticky='nsew')
                        
                        CDS=tk.Label(nbrCDSFrame,text=infos[6],bg=color[4],fg=color[6])
                        CDS.grid(row=1,column=0,sticky='nsew')
