#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
import re

from tkinter import *
import subprocess
import sys

from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from matplotlib.backend_bases import key_press_handler
from matplotlib.figure import Figure


# In[3]:


if sys.platform == 'darwin':
    def openFolder(path):
        subprocess.check_call(['open','--', path])
elif sys.platform == 'linux2':
    def openFolder(path):
        subprocess.check_call(['xdg-open', '--', path])
elif sys.platform == 'win32':
    def openFolder(path):
        path2=''
        for i in path:
            if i=='/':
                path2+='\\'
            else:
                path2+=i
        subprocess.check_call(['explorer', path2])


# In[ ]:





# In[4]:


def readBlastp(fichier):
    f=open(fichier,"r")
    dico={}
    for ligne in f:
        if ligne[0]!="#":
            ligne=ligne[4:-1]
            ligne=ligne.split("\t")
            dico[ligne[0],ligne[1]]=ligne[2:]
    f.close()
    return dico

#print (readBlastp("./Blastp/QUERY-GCA_000014865.1_ASM1486v1_translated_cds__DB-GCA_000009985.1_ASM998v1_translated_cds.out"))


# In[5]:


#Création d'un dictionnaire de Blastp
def dicoBlast():
    entriesBlastp=os.listdir('./Blastp')
    dicoBlastp={}
    for fichier in entriesBlastp:
        nom="./Blastp/"+fichier
        motif=re.compile('[A-Za-z]{3}_[0-9]{9}[.]?[0-9]{0,2}')
        cle=motif.findall(fichier)
        cle=cle[0],cle[1]
        dico=readBlastp(nom)
        dicoBlastp[cle]=dico

    return dicoBlastp

#print (dicoBlast())


# In[ ]:





# In[6]:


def readCDS(fichier):
    f=open(fichier,"r")
    dico={}
    for ligne in f:
        if ligne[0:2]=="Q#" and ligne.split('\t')[7][0:3]=='COG':
            ligne=ligne[:-1]
            ligne=ligne.split("\t")
            ligne[0]=ligne[0].split(" ")
            nomprot=ligne[0][2][5:]
            dico[nomprot]=ligne[5],ligne[7]
    f.close()
    return dico
#print (readCDS("./CD-search/QUERY-GCA_000009985.1_ASM998v1_translated_cds__DB-COGv3-16.out"))


# In[7]:


def dicoCds():
    entriesCDS=os.listdir('./CD-search')
    dicoCDS={}
    for fichier in entriesCDS:
        nom="./CD-search/"+fichier
        motif=re.compile('[A-Za-z]{3}_[0-9]{9}[.]?[0-9]{0,2}')
        cle=motif.findall(fichier)
        cle=cle[0]
        dico=readCDS(nom)
        dicoCDS[cle]=dico

    return dicoCDS

#print (dicoCds())


# In[ ]:





# In[8]:


def readProteome(fichier):
    f=open(fichier,"r")
    liste=[]
    seq=''
    nom=''
    sousliste=[]
    for ligne in f:
        if ligne[0]!=">":
            seq+=ligne[:-1]
        else:
            sous_liste=[nom,len(seq)]
            liste.append(sous_liste)
            nom=ligne.split(" ")[0][5:]
            seq=''
    sous_liste=[nom,len(seq)]
    liste.append(sous_liste)
    f.close()
    return liste[1:]

#print(readProteome("./Proteomes/GCA_000009985.1_ASM998v1_translated_cds.faa"))


# In[9]:


def dicoProteome():
    entriesProteome=os.listdir('./Proteomes')
    dicoProteom={}
    for fichier in entriesProteome:
        nom="./Proteomes/"+fichier
        motif=re.compile('[A-Za-z]{3}_[0-9]{9}[.]?[0-9]{0,2}')
        cle=motif.findall(fichier)
        cle=cle[0]
        liste=readProteome(nom)
        dicoProteom[cle]=liste

    return dicoProteom

#print (dicoProteome())


# In[ ]:





# In[10]:


def readCOG(fichier):
    f=open(fichier,"r",encoding='windows-1252')
    dico={}
    for ligne in f:
        if ligne[0]!="#":
            ligne=ligne[:-1]
            ligne=ligne.split("\t")
            dico[ligne[0]]=ligne[1],ligne[2]
    return dico

def dicoCog():
    dicoCOG=readCOG("cognames2003-2014.tab.txt")
    return dicoCOG

#print (dicoCog())


# In[11]:


def readfun(fichier):
    f=open(fichier,"r")
    dico={}
    for ligne in f:
        if ligne[0]!="#":
            ligne=ligne[:-1]
            ligne=ligne.split("\t")
            dico[ligne[0]]=ligne[1]
    return dico

def dicoFUN():
    dicofun=readfun("fun2003-2014.tab.txt")
    return dicofun

#print (dicoFUN())


# In[12]:


def readCSV(fichier):
    f=open(fichier,"r")
    dico={}
    for ligne in f:
        if ligne[0]!="#":
            ligne=ligne[:-1]
            ligne=ligne.split(",")
            dico[ligne[5][1:-1]]=ligne[0][1:-1]
    return dico

def dicoCsv():
    dicoCSV=readCSV("prokaryotes_complete-genomes.csv")
    return dicoCSV

#print (dicoCsv())


# In[ ]:





# In[13]:


#Fonction qui donne les coordonnées selon le test BLAST 

def Resultblast(dicoBlast, dicoProteome, ChoixGCA, evalue_blast, identite, couverture):
        
        listeprot1=[]
        for i in dicoProteome[ChoixGCA[0]]:
            listeprot1.append(i[0])
    
        listeprot2=[]
        for i in dicoProteome[ChoixGCA[1]]:
            listeprot2.append(i[0])
            
        listeResultx=[]
        listeResulty=[]

        for i in dicoBlast[ChoixGCA]:
                
            Couv_subject = (float(dicoBlast[ChoixGCA][i][7]) - float(dicoBlast[ChoixGCA][i][6]) +1) / float(dicoProteome[ChoixGCA[1]][int(listeprot2.index(i[1]))][1])
            Couv_query = (float(dicoBlast[ChoixGCA][i][5]) - float(dicoBlast[ChoixGCA][i][4]) +1) / float(dicoProteome[ChoixGCA[0]][int(listeprot1.index(i[0]))][1])
                
            if (float(dicoBlast[ChoixGCA][i][8]) <= evalue_blast) and (float(dicoBlast[ChoixGCA][i][0]) >= identite) and (Couv_query >= couverture) and (Couv_subject >= couverture):

                listeResultx.append(listeprot1.index(i[0]))
                listeResulty.append(listeprot2.index(i[1]))
    
        return listeResultx,listeResulty


# In[ ]:





# In[14]:


def Resultatcds(dicoCDS, dicoProteome, dicoCOG, ChoixGCA, evalue_cds, type_test):
    listeprot1=[]
    for i in dicoProteome[ChoixGCA[0]]:
        listeprot1.append(i[0])
    
    listeprot2=[]
    for i in dicoProteome[ChoixGCA[1]]:
        listeprot2.append(i[0])
        
    listeResultx=[]
    listeResulty=[]
    dCOG={}
    
    if type_test=="COG":
        for organisme in ChoixGCA:
            dicotmpCOG={}
            for proteine in dicoCDS[organisme]:
                if float(dicoCDS[organisme][proteine][0])<=evalue_cds :
                    if dicoCDS[organisme][proteine][1] in dicotmpCOG:
                        dicotmpCOG[dicoCDS[organisme][proteine][1]].append(proteine)
                    else:
                        dicotmpCOG[dicoCDS[organisme][proteine][1]]=[proteine]
            dCOG[organisme]=dicotmpCOG
    
        for cogname in dCOG[ChoixGCA[0]]:
            if cogname in dCOG[ChoixGCA[1]]:
                for prot1 in dCOG[ChoixGCA[0]][cogname]:
                    for prot2 in dCOG[ChoixGCA[1]][cogname]:
                        listeResultx.append(listeprot1.index(prot1))
                        listeResulty.append(listeprot2.index(prot2))
                        
                        
    if type_test=="ANF":
        dANF={}
        for organisme in ChoixGCA:
            dicotmpCOG={}
            dicotmpANF={}
            for proteine in dicoCDS[organisme]:
                if float(dicoCDS[organisme][proteine][0])<=evalue_cds :
                    if dicoCDS[organisme][proteine][1] in dicotmpCOG:
                        dicotmpCOG[dicoCDS[organisme][proteine][1]].append(proteine)
                    else:
                        dicotmpCOG[dicoCDS[organisme][proteine][1]]=[proteine]
                        
            for code in dicotmpCOG:
                if dicoCOG[code][0] in dicotmpANF:
                    for i in dicotmpCOG[code]:
                        dicotmpANF[dicoCOG[code][0]].append(i)
                else:
                    dicotmpANF[dicoCOG[code][0]]=dicotmpCOG[code]
            dANF[organisme]=dicotmpANF
        
        for anf in dANF[ChoixGCA[0]]:
            if anf in dANF[ChoixGCA[1]]:
                for prot1 in dANF[ChoixGCA[0]][anf]:
                    for prot2 in dANF[ChoixGCA[1]][anf]:
                        listeResultx.append(listeprot1.index(prot1))
                        listeResulty.append(listeprot2.index(prot2))
    
    return listeResultx, listeResulty

Rcds=Resultatcds(dicoCds(), dicoProteome(), dicoCog(), ('GCA_000014865.1', 'GCA_000009985.1'), 1e-4, "COG")

#print (Rcds[0])
#print ('\n')
#print (Rcds[1])
    


# In[ ]:





# In[16]:


class interface:
    
    def __init__(self):
        
        self.Blastp=dicoBlast()
        self.CDS=dicoCds()
        self.proteome=dicoProteome()
        self.COG=dicoCog()
        self.fun=dicoFUN()
        self.CSV=dicoCsv()
        
        self.interface()
        

    #Fonctions utilisées pour dessiner et réinitialiser l'affichage    
    def dessiner(self):
        
        #Valeur par défaut
        self.identite=0
        self.evalue= 1e-4 
        self.Hit=0
        
        #Vérification que les valeurs entrées sont valides
        if self.variable_identite.get()==1 and self.entry_identite.get()!="":
            try :
                float(self.entry_identite.get())
                
            except ValueError:
                message=Label(self.frameFullScreenBas, text="La valeur du pourcentage d'identité n'est pas valable", font=('Helvetica',12), fg='red', padx=10, pady=10)
                message.grid(row=0, column=0, columnspan=3,sticky="nsew")
                return None
            
            if not 0<=float(self.entry_identite.get())<=100:
                message=Label(self.frameFullScreenBas, text="La valeur du pourcentage d'identité doit être comprise entre 0 et 100", font=('Helvetica',12), fg='red', padx=10, pady=10)
                message.grid(row=0, column=0, columnspan=3,sticky="nsew")
                return None
            self.identite= float(self.entry_identite.get())
            
               
        if self.variable_evalue.get()==1 and self.entry_evalue.get()!="":
            try :
                float (self.entry_evalue.get())
            
            except ValueError:
                message=Label(self.frameFullScreenBas, text="La valeur de e_value n'est pas valable", font=('Helvetica',12), fg='red', padx=10, pady=10)
                message.grid(row=0, column=0, columnspan=3,sticky="nsew")
                return None
            
            if not 0<=float(self.entry_evalue.get())<=1e-4:
                message=Label(self.frameFullScreenBas, text="La valeur de e_value doit être entre comprise 0 et 1e-4", font=('Helvetica',12), fg='red', padx=10, pady=10)
                message.grid(row=0, column=0, columnspan=3,sticky="nsew")
                return None
            self.evalue= float(self.entry_evalue.get())
            
        
        if self.variable_Hit.get()==1 and self.entryHit.get()!="":
            try:
                float(self.entryHit.get())
                
            except ValueError:
                message=Label(self.frameFullScreenBas, text="La valeur de couverture du Hit n'est pas valable", font=('Helvetica',12), fg='red', padx=10, pady=10)
                message.grid(row=0, column=0, columnspan=3,sticky="nsew")
                return None
            
            if not 0<=float(self.entryHit.get())<=100:
                message=Label(self.frameFullScreenBas, text="La valeur de couverture de Hit doit être entre comprise 0 et 100", font=('Helvetica',12), fg='red', padx=10, pady=10)
                message.grid(row=0, column=0, columnspan=3,sticky="nsew")
                return None
            
            self.Hit= float(self.entryHit.get())/100
            
        self.noms=self.variableIndividus.get().split(' VS ')
        GCA=self.listeGCA[self.listenom.index(self.noms[0])],self.listeGCA[self.listenom.index(self.noms[1])] 
        
        if self.variableTest.get() != "Choix du Test":
            
            self.type_test=self.variable_radio.get()
            
            if self.variableTest.get() == "Blastp":
                self.Xplot,self.Yplot= Resultblast(self.Blastp, self.proteome, GCA, self.evalue, self.identite, self.Hit)
                
            elif self.variableTest.get() == "CD-search":
                self.Xplot,self.Yplot= Resultatcds(self.CDS, self.proteome, self.COG, GCA, self.evalue, self.type_test)
            
            message=Label(self.frameFullScreenBas,text=" ")
            message.grid(row=0, column=0, columnspan=3,sticky="nsew")
            
            self.taille.remove()
            self.taille=self.fig.add_subplot(111)
            self.scatter.remove()
            try:
                self.scatter_s.remove()
                self.scatter_s=self.taille.scatter([],[], s=0.1, c='r')
            except AttributeError:
                None
            self.scatter=self.taille.scatter(self.Xplot,self.Yplot, s=0.1, c='#056BB3')
            self.canvas.draw() 
            
            self.liste_variables=[self.variableTest.get(), self.noms[0], self.noms[1], self.entry_evalue.get(), self.Hit, self.identite, self.type_test]
            
        else:
            message=Label(self.frameFullScreenBas, text="Choisir un test", font=('Helvetica',12), fg='red', padx=10, pady=10)
            message.grid(row=0, column=0, columnspan=3,sticky="nsew")
            
    
    #Fonction qui permet de faire apparaitre les synténies sur le dotplot et de remplir un fichier txt 
    def synt(self):
        try:
            X=self.Xplot
            Y=self.Yplot
        except AttributeError:
            return None
        
        try:
            os.makedirs("Synténies")
        except OSError:
            if not os.path.isdir("Synténies"):
                raise
        os.chdir("Synténies")
                
        differences={}
        
        for i in range (len(X)):
            diff=int(X[i])- int(Y[i])
            if diff in differences:
                differences[diff].append(X[i])
            else:
                differences[diff]=[X[i],Y[i]]
        
        sommes={}
        for i in range (len(X)):
            s=int(X[i]) + int(Y[i])
            if s in sommes:
                sommes[s].append(X[i])
            else:
                sommes[s]=[X[i]]
                
                
        taille=2
        if self.entry_syntenie.get()!='':
            try:
                taille=int(self.entry_syntenie.get())

            except ValueError:
                message=Label(self.frameFullScreenBas, text="La taille des blocs synténie doit être un entier", font=('Helvetica',12), fg='red', padx=10, pady=10)
                message.grid(row=0, column=0, columnspan=3,sticky="nsew")
                return None
            
        if taille<2:
            message=Label(self.frameFullScreenBas, text="La taille des blocs synténie doit être minimum de 2", font=('Helvetica',12), fg='red', padx=10, pady=10)
            message.grid(row=0, column=0, columnspan=3,sticky="nsew")
            
            
        syntenie_croissante=[]
        for i in differences:
            sorted(differences[i])
            sous_liste=[[differences[i][0],differences[i][0]-i]]
            for j in range(0,len(differences[i])-1):
                if differences[i][j]+1==differences[i][j+1] :
                    sous_liste.append([differences[i][j+1], differences[i][j+1]-i])
                else:
                    if len(sous_liste)>=taille:
                        syntenie_croissante.append(sous_liste)
                    sous_liste=[[differences[i][j+1], differences[i][j+1]-i]]
                    
        syntenie_decroissante=[]
        for i in sommes:
            sorted(sommes[i])
            sous_liste=[[sommes[i][0],i-sommes[i][0]]]
            for j in range(0,len(sommes[i])-1):
                if sommes[i][j]+1==sommes[i][j+1] :
                    sous_liste.append([sommes[i][j+1],i-sommes[i][j+1]])
                else:
                    if len(sous_liste)>=taille:
                        syntenie_decroissante.append(sous_liste)
                    sous_liste=[[sommes[i][j+1],i-sommes[i][j+1]]]
            

        titre=self.liste_variables[0]+"   "+self.liste_variables[1]+" VS "+self.liste_variables[2]+".txt"
        
        ligne_ref='# evalue= '+str(self.liste_variables[3])
            
        if self.liste_variables[0]=="Blastp":
            ligne_ref+=", Couverture du hit: "+str(self.liste_variables[4])+" , % d'identité: "+str(self.liste_variables[5])

            
        elif self.liste_variables[0]=="CD-search":
            ligne_ref+=" ; type de test: "+self.liste_variables[6]
            
        ligne_ref+=" ; taille min des blocs: "+str(taille)
            
        try:
            verif=open(str(titre),'r')
            for ligne in verif:
                if ligne_ref in ligne:
                    verif.close()
                    os.chdir('..')
                    
                    try:
                        openFolder('./Synténies')
                    except:
                        None
                    Xsynt=[]
                    Ysynt=[]
                    for i in syntenie_croissante:
                        for j in i:
                            Xsynt.append(j[0])
                            Ysynt.append(j[1])
                    for i in syntenie_decroissante:
                        for j in i:
                            Xsynt.append(j[0])
                            Ysynt.append(j[1])
                    try:
                        self.scatter_s.remove()
                    except AttributeError:
                        None
                    self.scatter_s=self.taille.scatter(Xsynt,Ysynt, s=0.1, c='r')
                    self.canvas.draw()
                    return None
            verif.close()
            fichier=open(str(titre),'a')
            
        except FileNotFoundError:
            fichier=open(str(titre),'a')
            
        fichier.write(ligne_ref)
        fichier.write("  ")
        fichier.write("\n")
        fichier.write (" 1ere prot X    Dernière prot X    1ere prot Y    Derniere prot Y    Taille du bloc")
        fichier.write("\n")
        fichier.write("\n")

        for i in syntenie_croissante:
            fichier.write(str(i[0][0]))
            fichier.write("\t")
            fichier.write(str(i[-1][0]))
            fichier.write("\t")
            fichier.write(str(i[0][1]))
            fichier.write("\t")
            fichier.write(str(i[-1][1]))
            fichier.write("\t")
            fichier.write(str(len(i)))
            fichier.write("\n")
        for j in syntenie_decroissante:
            fichier.write(str(j[0][0]))
            fichier.write("\t")
            fichier.write(str(j[-1][0]))
            fichier.write("\t")
            fichier.write(str(j[0][1]))
            fichier.write("\t")
            fichier.write(str(j[-1][1]))
            fichier.write("\t")
            fichier.write(str(len(j)))
            fichier.write("\n")

        fichier.write("\n")
        fichier.write("\n")

        fichier.close()
        os.chdir('..')
        
        try:
            openFolder('./Synténies')
        except:
            None
        
        Xsynt=[]
        Ysynt=[]
        for i in syntenie_croissante:
            for j in i:
                Xsynt.append(j[0])
                Ysynt.append(j[1])
        for i in syntenie_decroissante:
            for j in i:
                Xsynt.append(j[0])
                Ysynt.append(j[1])
        try:
            self.scatter_s.remove()
        except AttributeError:
            None
        self.scatter_s=self.taille.scatter(Xsynt,Ysynt, s=0.1, c='r')
        self.canvas.draw()
        
        
            
        

      
        
    def effacer(self):
        self.scatter.remove()
        try:
            self.scatter_s.remove()
        except AttributeError:
            None
        self.scatter=self.taille.scatter([],[], s=0.1)
        self.scatter_s=self.taille.scatter([],[], s=0.1)
        self.canvas.draw()
        
        
 

        
        
    def interface(self):
    
        self.fenetre=Tk()
        self.fenetre.title('Comparateur de Génomes')
        self.fenetre.geometry('1100x500')

        
        #Main Frame
        
        principale=Frame(self.fenetre,bg="white")
        principale.pack(fill="both",expand=True)

        principale.grid_columnconfigure(0,weight=5, uniform="group1")
        principale.grid_columnconfigure(1,weight=5, uniform="group1")

        
        principale.grid_rowconfigure(0,weight=1, uniform="group2")
        principale.grid_rowconfigure(1,weight=12, uniform="group2")
        principale.grid_rowconfigure(2,weight=1, uniform="group2")


        #Création des frames
        
        frameFullScreenHaut=Frame(principale,bg='blue')
        frameFullScreenHaut.grid(row=0, column=0, columnspan=4,sticky="nsew")
        
        self.frameFullScreenBas=Frame(principale)
        self.frameFullScreenBas.grid(row=1,column=0,sticky="nsew")
        
        frameFullScreenImage=Frame(principale,bg='pink')
        frameFullScreenImage.grid(row=1,column=1,sticky="nsew")
        
        frameFullScreenBande=Frame(principale)
        frameFullScreenBande.grid(row=2,column=0, columnspan=2, sticky="nsew")
        
        
        #Paramètre des frames
        
        frameFullScreenHaut.grid_rowconfigure(0,weight=1, uniform='group1')
        frameFullScreenHaut.grid_columnconfigure(0,weight=1,uniform='group2')
        frameFullScreenHaut.grid_columnconfigure(1,weight=1,uniform='group2')
        
        
        self.frameFullScreenBas.grid_columnconfigure(0,weight=1,uniform='group2')
        self.frameFullScreenBas.grid_columnconfigure(1,weight=4,uniform='group2')
        self.frameFullScreenBas.grid_columnconfigure(2,weight=4,uniform='group2')
        
        self.frameFullScreenBas.grid_rowconfigure(0,weight=1, uniform='group1')
        self.frameFullScreenBas.grid_rowconfigure(1,weight=5, uniform='group1')
        self.frameFullScreenBas.grid_rowconfigure(2,weight=5, uniform='group1')
        self.frameFullScreenBas.grid_rowconfigure(3,weight=5, uniform='group1')
        self.frameFullScreenBas.grid_rowconfigure(4,weight=2, uniform='group1')
        self.frameFullScreenBas.grid_rowconfigure(5,weight=5, uniform='group1')
        
        
        
        frameFullScreenBande.grid_columnconfigure(0,weight=1,uniform='group2')
        frameFullScreenBande.grid_columnconfigure(1,weight=1,uniform='group2')
        frameFullScreenBande.grid_columnconfigure(2,weight=1,uniform='group2')
        frameFullScreenBande.grid_columnconfigure(3,weight=1,uniform='group2')
        frameFullScreenBande.grid_rowconfigure(0,weight=1, uniform='group1')
        
        

        
        
        #Menu déroulant pour choisir les individus à comparer
        
        self.listeIndividus=[]
        self.listenom=[]
        self.listeGCA=[]
        i=0
        for cle in self.Blastp.keys():

            self.listeGCA.append(cle[0])
            self.listeGCA.append(cle[1])
            self.listenom.append(self.CSV[cle[0]])
            self.listenom.append(self.CSV[cle[1]])

            self.listeIndividus.append(self.listenom[i]+' VS '+self.listenom[i+1])
            i+=2
        
        self.variableIndividus=StringVar(self.fenetre)
        self.variableIndividus.set(self.listeIndividus[0])
        self.boutonIndividus=OptionMenu(frameFullScreenHaut, self.variableIndividus, *self.listeIndividus)        
        self.boutonIndividus.grid(column=1, row=0, sticky="nsew")    
        
        
        #Fonctions qui permettent de griser ou non les entrys 
        
        def griser_identite():
            if self.variable_identite.get()==1:
                self.entry_identite= Entry(self.frameFullScreenBas)
            else:
                self.entry_identite= Entry(self.frameFullScreenBas, state = DISABLED)
            self.entry_identite.grid(column=2,row=1,sticky="w")
            
        def griser_Hit():
            if self.variable_Hit.get()==1:
                self.entryHit= Entry(self.frameFullScreenBas)
            else:
                self.entryHit= Entry(self.frameFullScreenBas, state = DISABLED)
            self.entryHit.grid(column=2,row=2,sticky="w")  
            
        def griser_evalue():
            if self.variable_evalue.get()==1:
                self.entry_evalue= Entry(self.frameFullScreenBas)
            else:
                self.entry_evalue= Entry(self.frameFullScreenBas, state = DISABLED)
            self.entry_evalue.grid(column=2,row=3,sticky="w")
            
        
        #Bouton checks
        
        self.variable_identite=IntVar(self.frameFullScreenBas)
        self.boutoncheckIdentite=Checkbutton(self.frameFullScreenBas,text="% identité", variable=self.variable_identite, command=griser_identite, state=DISABLED)
        self.boutoncheckIdentite.grid(column=1,row=1,sticky="w")
        self.entry_identite= Entry(self.frameFullScreenBas, state=DISABLED)
        self.entry_identite.grid(column=2,row=1,sticky="w")
        
        self.variable_Hit=IntVar(self.frameFullScreenBas)
        self.boutoncheckHit=Checkbutton(self.frameFullScreenBas,text="% Couverture du hit", variable=self.variable_Hit, command=griser_Hit, state=DISABLED)
        self.boutoncheckHit.grid(column=1,row=2,sticky="w")
        self.entryHit= Entry(self.frameFullScreenBas, state=DISABLED)
        self.entryHit.grid(column=2,row=2,sticky="w")        
                
        self.variable_evalue=IntVar(self.frameFullScreenBas)
        self.boutoncheck_evalue=Checkbutton(self.frameFullScreenBas,text="e-value", variable=self.variable_evalue, command=griser_evalue, state=DISABLED)
        self.boutoncheck_evalue.grid(column=1,row=3,sticky="w")
        self.entry_evalue= Entry(self.frameFullScreenBas, state=DISABLED)
        self.entry_evalue.grid(column=2,row=3,sticky="w") 
        

        ##Boutons pour faire le choix entre annotation fonctionnelle et COG
        etiq_radio=['COG', 'ANF']
        self.variable_radio=StringVar()
        self.variable_radio.set(etiq_radio[0])
        self.COGbutton=Radiobutton(self.frameFullScreenBas,variable=self.variable_radio,text=etiq_radio[0],value=etiq_radio[0],indicatoron=0, state=DISABLED)
        self.COGbutton.grid(row=4,column=1,sticky='n')
        self.AFbutton=Radiobutton(self.frameFullScreenBas,variable=self.variable_radio,text=etiq_radio[1],value=etiq_radio[1],indicatoron=0, state=DISABLED)
        self.AFbutton.grid(row=4,column=2,sticky='n')
        
        
        #Menu déroulant pour choisir le test
        def griser(*args):
            if self.variableTest.get()=="Blastp":
                self.boutoncheckIdentite=Checkbutton(self.frameFullScreenBas,text="% identité", variable=self.variable_identite, command=griser_identite)
                self.boutoncheckIdentite.grid(column=1,row=1,sticky="w")
                
                self.boutoncheckHit=Checkbutton(self.frameFullScreenBas,text="% Couverture du hit", variable=self.variable_Hit, command=griser_Hit)
                self.boutoncheckHit.grid(column=1,row=2,sticky="w")
            
                self.boutoncheck_evalue=Checkbutton(self.frameFullScreenBas,text="e-value", variable=self.variable_evalue, command=griser_evalue)
                self.boutoncheck_evalue.grid(column=1,row=3,sticky="w")
                
                self.COGbutton=Radiobutton(self.frameFullScreenBas,variable=self.variable_radio,text=etiq_radio[0],value=etiq_radio[0],indicatoron=0, state=DISABLED)
                self.COGbutton.grid(row=4,column=1,sticky='n')
                self.AFbutton=Radiobutton(self.frameFullScreenBas,variable=self.variable_radio,text=etiq_radio[1],value=etiq_radio[1],indicatoron=0, state=DISABLED)
                self.AFbutton.grid(row=4,column=2,sticky='n')
                
                
            elif self.variableTest.get()=="CD-search":
                self.variable_identite.set(0)
                self.boutoncheckIdentite=Checkbutton(self.frameFullScreenBas, text="% identité", variable=self.variable_identite, command=griser_identite, state=DISABLED)
                self.boutoncheckIdentite.grid(column=1,row=1,sticky="w")
                self.entry_identite= Entry(self.frameFullScreenBas, state = DISABLED)
                self.entry_identite.grid(column=2,row=1,sticky="w")
                
                self.variable_Hit.set(0)
                self.boutoncheckHit=Checkbutton(self.frameFullScreenBas, text="% Couverture du hit", variable=self.variable_Hit, command=griser_Hit, state=DISABLED)
                self.boutoncheckHit.grid(column=1,row=2,sticky="w")
                self.entryHit= Entry(self.frameFullScreenBas, state=DISABLED)
                self.entryHit.grid(column=2,row=2,sticky="w")  
                
                self.boutoncheck_evalue=Checkbutton(self.frameFullScreenBas,text="e-value", variable=self.variable_evalue, command=griser_evalue)
                self.boutoncheck_evalue.grid(column=1,row=3,sticky="w")
                
                self.COGbutton=Radiobutton(self.frameFullScreenBas,variable=self.variable_radio,text=etiq_radio[0],value=etiq_radio[0],indicatoron=0)
                self.COGbutton.grid(row=4,column=1,sticky='n')                
                self.AFbutton=Radiobutton(self.frameFullScreenBas,variable=self.variable_radio,text=etiq_radio[1],value=etiq_radio[1],indicatoron=0)
                self.AFbutton.grid(row=4,column=2,sticky='n')   
            
            
        listeTest=["Choix du Test","Blastp","CD-search"]
        self.variableTest=StringVar(self.fenetre)
        self.variableTest.set(listeTest[0])
        self.boutonTest=OptionMenu(frameFullScreenHaut, self.variableTest, *listeTest[1:], command=griser)
        self.boutonTest.grid(column=0, row=0, sticky="nsew")
        
        
        def Fermer():
            global fenetre
            self.fenetre.destroy()
            

        #Canvas       
        self.fig = Figure(figsize=(5, 4), dpi=100)
        self.taille=self.fig.add_subplot(111)
        self.scatter=self.taille.scatter([],[], s=0.1)
        self.canvas = FigureCanvasTkAgg(self.fig, master=frameFullScreenImage)  
        self.canvas.draw()
        
        #Barre de gestion du graphique
        toolbar = NavigationToolbar2Tk(self.canvas, frameFullScreenImage)
        toolbar.update()


        def on_key_press(event):
            key_press_handler(event, self.canvas, toolbar)


        self.canvas.mpl_connect("key_press_event", on_key_press)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
    
        
        #Boutons de la Frame tout en bas
        boutonReset=Button(frameFullScreenBande, text="Effacer", command=self.effacer)
        boutonReset.grid (column=2, row=0) 
        boutonAfficher=Button(frameFullScreenBande, text="Afficher", command=self.dessiner)
        boutonAfficher.grid (column=1, row=0)
        boutonQuitter=Button(frameFullScreenBande, text="Quitter", command=Fermer)
        boutonQuitter.grid (column=3, row=0)
        boutonSyntenie=Button (frameFullScreenBande, text="Synténie", command=self.synt)
        boutonSyntenie.grid(column=0, row=0)
        
        self.entry_syntenie= Entry(self.frameFullScreenBas)
        self.entry_syntenie.grid(column=2,row=5, sticky='w') 
        message=Label(self.frameFullScreenBas, text="Taille des blocks de synténie", font=('Helvetica',12), padx=10, pady=10)
        message.grid(row=5, column=1, sticky="")
        
        self.fenetre.mainloop()
        

Final=interface()


# In[ ]:




