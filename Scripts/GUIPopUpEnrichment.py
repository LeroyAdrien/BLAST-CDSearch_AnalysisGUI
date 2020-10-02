import tkinter as tk
import numpy as np
from calculusFiltering import *
class interface_popup_enrichment:


        def __init__(self,
                     dicoProteome,
                     dicoFun,
                     dicoNoms2COG,
                     dicoCOG,
                     coordonnesXSyntenyCroissantCDSearch,
                     coordonnesYSyntenyCroissantCDSearch,
                     coordonnesXSyntenyDecroissantCDSearch,
                     coordonnesYSyntenyDecroissantCDSearch,
                     GCACDSearch,
                     titreCDSearch,
                     coordonnesXSyntenyCroissantBLAST,
                     coordonnesYSyntenyCroissantBLAST,
                     coordonnesXSyntenyDecroissantBLAST,
                     coordonnesYSyntenyDecroissantBLAST,
                     GCABLAST,
                     titreBLAST
                     ):

                #Creation de l'interface'
                interfacePopUpEnrichment=tk.Tk()
                interfacePopUpEnrichment.title('Enrichissement')
                interfacePopUpEnrichment.geometry('1300x650')
                
                #Couleurs
                self.couleurTexte='#D8D8D8'
                couleurBackground='#2c2f33'
                
                #Frames
                fenetre=tk.Frame(interfacePopUpEnrichment,bg=couleurBackground)
                fenetre.pack(fill='both',expand=True)
                fenetre.grid_columnconfigure(0,weight=1)
                fenetre.grid_rowconfigure(0,weight=1)
                fenetre.grid_rowconfigure(1,weight=8)
                
                #Frame Titre et description
                popUpFrameHaut=tk.Frame(fenetre,bg=couleurBackground)
                popUpFrameHaut.grid(row=0,column=0,sticky='nsew')
                #Frame Valeurs 
                popUpFrameBas=tk.Frame(fenetre,bg=couleurBackground)
                popUpFrameBas.grid(row=1,column=0,sticky='nsew')
                #Frame pour intitulé des fonctions
                popUpFrameBas.grid_columnconfigure(0,weight=3,uniform='group1')
                #Labels Descriptifs
                labelFonction=tk.Label(popUpFrameBas,text='Fonction de la protéine')
                labelFonction.configure(bg=couleurBackground,fg=self.couleurTexte)
                labelFonction.grid(row=0,column=0,rowspan=2,sticky="nsew")
                
                labelPValue=tk.Label(popUpFrameHaut,text="P-value de l'enrichissement selon un test de fisher Exact unilatéral (H1:Synténie>Hors Synténie)")
                labelPValue.configure(bg=couleurBackground,fg=self.couleurTexte)
                labelPValue.grid(row=0,column=0,sticky="nsew")
                
                #Si CD Search est là 
                if len(coordonnesXSyntenyCroissantCDSearch)>0:
                        #Calcule du talbeau de contigence pour CD-Search
                        #Concatenation des coordonnes pour plus de simplicité     
                        coordonnesX=np.concatenate((coordonnesXSyntenyCroissantCDSearch,coordonnesXSyntenyDecroissantCDSearch))
                        coordonnesY=np.concatenate((coordonnesYSyntenyCroissantCDSearch,coordonnesYSyntenyDecroissantCDSearch))
                        
                        #Generation d'une tableau de contigence des fonctions dans les blocs et hors bloc
                        table1,table2,alphabetFun=SyntenieContigence(coordonnesX,coordonnesY,dicoProteome,GCACDSearch[0],GCACDSearch[1],dicoFun,dicoNoms2COG,dicoCOG)
                        #Analffyse des enrichissements des diférentes fonctions
                        PValueOrganisme1CDSearch=FisherResults(table1)
                        PValueOrganisme2CDSearch=FisherResults(table2)
                
                        popUpFrameBas.grid_columnconfigure(1,weight=2,uniform='group1')
                        popUpFrameBas.grid_columnconfigure(2,weight=2,uniform='group1')
                        
                        popUpFrameBas.grid_rowconfigure(0,weight=1,uniform='group2')
                        popUpFrameBas.grid_rowconfigure(1,weight=1,uniform='group2')

                        #Label de la méthode 
                        labelmethode=tk.Label(popUpFrameBas,text='CD-Search')
                        labelmethode.configure(bg='#a80000',fg=self.couleurTexte)
                        labelmethode.grid(row=0,column=1,columnspan=2,sticky='nsew')
                        
                        #Label Organisme1
                        labelOrga1=tk.Label(popUpFrameBas,text=titreCDSearch[0])
                        labelOrga1.configure(bg='#a80000',fg=self.couleurTexte)
                        labelOrga1.grid(row=1,column=1,sticky="nsew")
                        
                        #Label Organisme2
                        labelOrga2=tk.Label(popUpFrameBas,text=titreCDSearch[1])
                        labelOrga2.configure(bg='#a80000',fg=self.couleurTexte)
                        labelOrga2.grid(row=1,column=2,sticky="nsew")
                        #Pour chaque fonction
                        for i in range(len(dicoFun.keys())):
                                #Creer la ligne
                                popUpFrameBas.grid_rowconfigure(i+2,weight=1)
                                #Ajouter l'intitulé
                                labelIntitule=tk.Label(popUpFrameBas,text=dicoFun[alphabetFun[i]])
                                labelIntitule.configure(bg=couleurBackground,fg=self.couleurTexte)
                                labelIntitule.grid(row=i+2,column=0,sticky='nsew')
                                #Ajouter PValue Organisme ff1
                                LabelOrga1=tk.Label(popUpFrameBas,text=PValueOrganisme1CDSearch[i])
                                
                                if PValueOrganisme1CDSearch[i]>0.05:
                                        couleurPValue='#5c0e0e'
                                else:
                                        couleurPValue='#286f28'
                                        
                                LabelOrga1.configure(bg=couleurPValue,fg=self.couleurTexte)
                                LabelOrga1.grid(row=i+2,column=1,sticky="ew")
                                #Ajouter PValue Organisme2
                                LabelOrga2=tk.Label(popUpFrameBas,text=PValueOrganisme2CDSearch[i])
                                
                                if PValueOrganisme2CDSearch[i]>0.05:
                                        couleurPValue='#5c0e0e'
                                else:
                                        couleurPValue='#286f28'
                                        
                                LabelOrga2.configure(bg=couleurPValue,fg=self.couleurTexte)
                                LabelOrga2.grid(row=i+2,column=2,sticky="ew")
                if len(coordonnesXSyntenyCroissantBLAST)>0:
                        #Calcule du talbeau de contigence pour CD-Search
                        #Concatenation des coordonnes pour plus de simplicité     
                        coordonnesX=np.concatenate((coordonnesXSyntenyCroissantBLAST,coordonnesXSyntenyDecroissantBLAST))
                        coordonnesY=np.concatenate((coordonnesYSyntenyCroissantBLAST,coordonnesYSyntenyDecroissantBLAST))
                        
                        #Generation d'une tableau de contigence des fonctions dans les blocs et hors bloc
                        table1,table2,alphabetFun=SyntenieContigence(coordonnesX,coordonnesY,dicoProteome,GCABLAST[0],GCABLAST[1],dicoFun,dicoNoms2COG,dicoCOG)
                        #Analffyse des enrichissements des diférentes fonctions
                        PValueOrganisme1BLAST=FisherResults(table1)
                        PValueOrganisme2BLAST=FisherResults(table2)
                
                        popUpFrameBas.grid_columnconfigure(3,weight=2,uniform='group1')
                        popUpFrameBas.grid_columnconfigure(4,weight=2,uniform='group1')
                        
                        popUpFrameBas.grid_rowconfigure(0,weight=1,uniform='group2')
                        popUpFrameBas.grid_rowconfigure(1,weight=1,uniform='group2')
                        
                        #Label de la méthode 
                        labelmethode=tk.Label(popUpFrameBas,text='BLAST')
                        labelmethode.configure(bg='#0057e7',fg=self.couleurTexte)
                        labelmethode.grid(row=0,column=3,columnspan=2,sticky='nsew')
                        

                        labelOrga1=tk.Label(popUpFrameBas,text=titreBLAST[0])
                        labelOrga1.configure(bg='#0057e7',fg=self.couleurTexte)
                        labelOrga1.grid(row=1,column=3,sticky="nsew")
                        
                        labelOrga2=tk.Label(popUpFrameBas,text=titreBLAST[1])
                        labelOrga2.configure(bg='#0057e7',fg=self.couleurTexte)
                        labelOrga2.grid(row=1,column=4,sticky="nsew")
                        #Pour chaque fonction
                        for i in range(len(dicoFun.keys())):
                                #Creer la ligne
                                popUpFrameBas.grid_rowconfigure(i+2,weight=1)
                                #Ajouter l'intitulé
                                labelIntitule=tk.Label(popUpFrameBas,text=dicoFun[alphabetFun[i]])
                                labelIntitule.configure(bg=couleurBackground,fg=self.couleurTexte)
                                labelIntitule.grid(row=i+2,column=0,sticky='nsew')
                                #Ajouter PValue Organisme ff1
                                LabelOrga1=tk.Label(popUpFrameBas,text=PValueOrganisme1BLAST[i])
                                
                                if PValueOrganisme1BLAST[i]>0.05:
                                        couleurPValue='#5c0e0e'
                                else:
                                        couleurPValue='#286f28'
                                        
                                LabelOrga1.configure(bg=couleurPValue,fg=self.couleurTexte)
                                LabelOrga1.grid(row=i+2,column=3,sticky="ew")
                                #Ajouter PValue Organisme2
                                LabelOrga2=tk.Label(popUpFrameBas,text=PValueOrganisme2BLAST[i])
                                
                                if PValueOrganisme2BLAST[i]>0.05:
                                        couleurPValue='#5c0e0e'
                                else:
                                        couleurPValue='#286f28'
                                        
                                LabelOrga2.configure(bg=couleurPValue,fg=self.couleurTexte)
                                LabelOrga2.grid(row=i+2,column=4,sticky="ew")

                tk.mainloop()
