import matplotlib as plt
import tkinter as tk
from matplotlib.backends.backend_tkagg import (
    FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure
from FigurePlotting import *
    

class interfaceSyntenyDiagram:

        def __init__(
                     self,
                     dicoProteome,
                     dicoBLAST,
                     dicoCDSearch,
                     titreCDSearch,
                     titreBLAST,
                     coordonnesXSyntenyCroissantCDSearch,
                     coordonnesYSyntenyCroissantCDSearch,
                     coordonnesXSyntenyDecroissantCDSearch,
                     coordonnesYSyntenyDecroissantCDSearch,
                     coordonnesXSyntenyCroissantBLAST,
                     coordonnesYSyntenyCroissantBLAST,
                     coordonnesXSyntenyDecroissantBLAST,
                     coordonnesYSyntenyDecroissantBLAST,
                     tailleBLAST,
                     tailleCDSearch
                     
                     ):
                #print(len(coordonnesXSyntenyCroissantBLAST))
                #print(len(coordonnesXSyntenyCroissantCDSearch))
             #Creation de l'interface
                interfaceSyntenyDiagram=tk.Tk()
                interfaceSyntenyDiagram.title("Diagramme de Synténie")
                interfaceSyntenyDiagram.geometry('1300x650')
             
                     #Couleur
                self.couleurTexte='#D8D8D8'
                couleurBackground='#2c2f33'
                couleurBLASTSelected='#0057e7'
                couleurCDSearchSelected='#a80000'
                couleurintersectionSelected='#cc6600'
                
                #Frames principal
                fenetre=tk.Frame(interfaceSyntenyDiagram,bg=couleurBackground)
                fenetre.pack(fill='both',expand=True)
                fenetre.grid_columnconfigure(0,weight=1)
                fenetre.grid_rowconfigure(0,weight=1)
                fenetre.grid_rowconfigure(1,weight=8)
                
                #Frame du diagram
                diagramFrame=tk.Frame(fenetre,bg=couleurBackground)
                diagramFrame.grid(row=1,column=0,sticky='nsew')
                diagramFrame.grid_columnconfigure(0,weight=1)
                diagramFrame.grid_columnconfigure(1,weight=1)
                diagramFrame.grid_columnconfigure(2,weight=1)
                diagramFrame.grid_columnconfigure(3,weight=1)
                
                #Diagramme de BLAST
                if len(coordonnesXSyntenyCroissantBLAST)>0:
                #Mise en place des Frames
                        diagramFrame.grid_rowconfigure(0,weight=4,uniform='group1')
                        diagramFrame.grid_rowconfigure(1,weight=1,uniform='group1')
                        
                        diagramFrameBLAST=tk.Frame(diagramFrame,bg=couleurBackground)
                        diagramFrameBLAST.grid(row=0,column=0,columnspan=4,sticky="nsew")
                        
                        toolbarFrameBLAST=tk.Frame(diagramFrame,bg=couleurBackground)
                        toolbarFrameBLAST.grid(row=1,column=0,sticky='nsew')
                        
                #Mise en place de la figure
                        figBLAST = Figure(figsize=(24, 5), dpi=200)
                        figBLAST.patch.set_facecolor(couleurBackground)
                        self.axBLAST=figBLAST.add_subplot(111)
                        self.axBLAST.set_facecolor(couleurBackground)
                        self.axBLAST.axis("off")
               #Remplissage de la figure
                        self.axBLAST=PlottingSyntenyKariotype(
                                                                coordonnesXSyntenyCroissantBLAST,
                                                                coordonnesYSyntenyCroissantBLAST,
                                                                coordonnesXSyntenyDecroissantBLAST,
                                                                coordonnesYSyntenyDecroissantBLAST,
                                                                titreBLAST[0],
                                                                titreBLAST[1],
                                                                tailleBLAST,
                                                                couleurBLASTSelected,
                                                                'purple',
                                                                'green',
                                                                self.axBLAST
                                                              )
               
               
               
               
                        self.diagramBLASTCanvas = FigureCanvasTkAgg(figBLAST, master=diagramFrameBLAST)  # A tk.DrawingArea.
                        self.diagramBLASTCanvas.draw()
                        self.diagramBLASTCanvas.get_tk_widget().pack(side='top',expand=1)
                
                #Toolbar de la figure
                        self.toolbar = NavigationToolbar2Tk(self.diagramBLASTCanvas,toolbarFrameBLAST)
                        self.toolbar.configure(bg=couleurBackground,highlightbackground=couleurBackground)
                        self.toolbar.update()
                #Diagramme de CDSearch
                #Diagramme de CDSearch
                if len(coordonnesXSyntenyCroissantCDSearch)>0:
                #Mise en place des Frames
                        diagramFrame.grid_rowconfigure(2,weight=4,uniform='group1')
                        diagramFrame.grid_rowconfigure(3,weight=1,uniform='group1')
                        
                        diagramFrameCDSearch=tk.Frame(diagramFrame,bg="blue")
                        diagramFrameCDSearch.grid(row=2,column=0,columnspan=4,sticky="nsew")
                        
                        toolbarFrameCDSearch=tk.Frame(diagramFrame,bg=couleurBackground)
                        toolbarFrameCDSearch.grid(row=3,column=0,sticky='nsew')
                        
                #Mise en place de la figure
                        figCDSearch = Figure(figsize=(24, 5), dpi=200)
                        figCDSearch.patch.set_facecolor(couleurBackground)
                        self.axCDSearch=figCDSearch.add_subplot(111)
                        self.axCDSearch.set_facecolor(couleurBackground)
                        self.axCDSearch.axis("off")
               #Remplissage de la figure
                        self.axCDSearch=PlottingSyntenyKariotype(
                                                                coordonnesXSyntenyCroissantCDSearch,
                                                                coordonnesYSyntenyCroissantCDSearch,
                                                                coordonnesXSyntenyDecroissantCDSearch,
                                                                coordonnesYSyntenyDecroissantCDSearch,
                                                                titreCDSearch[0],
                                                                titreCDSearch[1],
                                                                tailleCDSearch,
                                                                couleurCDSearchSelected,
                                                                'purple',
                                                                'green',
                                                                self.axCDSearch
                                                              )
               
               
               
               
                        self.diagramCDSearchCanvas = FigureCanvasTkAgg(figCDSearch, master=diagramFrameCDSearch)  # A tk.DrawingArea.
                        self.diagramCDSearchCanvas.draw()
                        self.diagramCDSearchCanvas.get_tk_widget().pack(side='top',expand=1)
                
                #Toolbar de la figure
                        self.toolbar = NavigationToolbar2Tk(self.diagramCDSearchCanvas,toolbarFrameCDSearch)
                        self.toolbar.configure(bg=couleurBackground,highlightbackground=couleurBackground)
                        self.toolbar.update()

                
                tk.mainloop()
