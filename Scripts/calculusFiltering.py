# *-*coding:utf-8 -*-
import numpy as np
from fileparsing import *
from itertools import chain
from scipy import stats

#Calcul des Blast et CD-Search

#4 critères d'éaluation du BLAST
def ResultatsBlast(dicoBLAST,dicoProteome,choix):
    listeProteine0=[]
    for i in dicoProteome[choix[0]]:
        listeProteine0.append(i[0])
    listeProteine1=[]
    for i in dicoProteome[choix[1]]:
        listeProteine1.append(i[0])
    listeResultats=[]
    for i in dicoBLAST[choix].keys():
        listeResultats.append([listeProteine0.index(i[0]),
                                listeProteine1.index(i[1]),
                                float(dicoBLAST[choix][i][0]),
                                100*(int(dicoBLAST[choix][i][1])/int(dicoProteome[choix[0]][listeProteine0.index(i[0])][1])),
                                100*(int(dicoBLAST[choix][i][1])/int(dicoProteome[choix[1]][listeProteine1.index(i[1])][1])),
                                float(dicoBLAST[choix][i][8])])
    #$0 position protéome 1
    #$1 position protéome 1
    #$2 pourcentage d'identité
    #$3 cover sur la protéine 1
    #$4 cover sur la protéine 2
    #$5 E-value du match
    return listeResultats
        
def ResultatsCDSearchCOG(dicoCDSearch,dicoProteome,choixOrg1,choixOrg2):
    
    listeProteine0=[]
    for i in dicoProteome[choixOrg1]:
        listeProteine0.append(i[0])
        
    listeProteine1=[]
    for i in dicoProteome[choixOrg2]:
        listeProteine1.append(i[0])
        
    listeResultats=[]
        
    for i in dicoCDSearch[choixOrg1].keys():
        if i in dicoCDSearch[choixOrg2].keys():
                for j in dicoCDSearch[choixOrg1][i]:
                        for k in dicoCDSearch[choixOrg2][i]:
                                proteineOrga1=j[4]
                                proteineOrga2=k[4]
                                resultat=[int(listeProteine0.index(proteineOrga1)),
                                          int(listeProteine1.index(proteineOrga2)),
                                          float(j[2]),
                                          float(k[2])]
                                listeResultats.append(resultat)

    #$0 position protéome 1
    #$1 position protéome 2
    #$2 E-value protéine 1
    #$3 E-value protéine 2
    #$4 Nom du COG qui a matché
    return listeResultats

def ResultatsCDSearchFA(dicoCDSearch,dicoProteome,dicoCOG,choixOrg1,choixOrg2):

    listeProteine0=[]

    for i in dicoProteome[choixOrg1]:
        listeProteine0.append(i[0])

    listeProteine1=[]
    for i in dicoProteome[choixOrg2]:
        listeProteine1.append(i[0])
        
#Conversion du dictionnaire:
    dicoTempOrga1={}
    dicoTempOrga2={}
    
    for i in dicoCDSearch[choixOrg1].keys():
        
        #print(dicoTempOrga1.keys())
        cle=dicoCOG[i][0]
        if len(cle)>1:
                for j in cle:
                        if j in dicoTempOrga1.keys():
                                dicoTempOrga1[j]=dicoTempOrga1[j]+dicoCDSearch[choixOrg1][i]
                        else:
                                dicoTempOrga1[j]=dicoCDSearch[choixOrg1][i]
                                
                                
        else:
                if dicoCOG[i][0] in dicoTempOrga1.keys():
                        #print(dicoCDSearch[choixOrg1][i])
                        dicoTempOrga1[dicoCOG[i][0]]=dicoTempOrga1[dicoCOG[i][0]]+dicoCDSearch[choixOrg1][i]
                else:
                        dicoTempOrga1[dicoCOG[i][0]]=dicoCDSearch[choixOrg1][i]


    for i in dicoCDSearch[choixOrg2].keys():
        cle=dicoCOG[i][0]
        if len(cle)>1:
                for j in cle:
                        if j in dicoTempOrga2.keys():
                                dicoTempOrga2[j]=dicoTempOrga2[j]+dicoCDSearch[choixOrg2][i]
                        else:
                                dicoTempOrga2[j]=dicoCDSearch[choixOrg2][i]
                                
                                
        else:
                if dicoCOG[i][0] in dicoTempOrga2.keys():
                        dicoTempOrga2[dicoCOG[i][0]]=dicoTempOrga2[dicoCOG[i][0]]+dicoCDSearch[choixOrg2][i]
                else:
                        dicoTempOrga2[dicoCOG[i][0]]=dicoCDSearch[choixOrg2][i]

    listeResultats=[]

    for i in dicoTempOrga1.keys():
        if i in dicoTempOrga2.keys():
                for j in dicoTempOrga1[i]:
                        for k in dicoTempOrga2[i]:
                                proteineOrga1=j[4]
                                proteineOrga2=k[4]
                                resultat=[int(listeProteine0.index(proteineOrga1)),
                                          int(listeProteine1.index(proteineOrga2)),
                                          float(j[2]),
                                          float(k[2])]
                                listeResultats.append(resultat)
    #$0 position protéome 1
    #$1 position protéome 2
    #$2 E-value protéine 1
    #$3 E-value protéine 2
    #$4 Nom du COG qui a matché

    return listeResultats

#Fonction de filtrage des résultats de Blast et CD Search

#Filtre les résultats de BLAST en fonction de couverture, pourcentage d'identité et Evalue'
def FilterBLAST(listeResultats,pident,coverProt1,coverProt2,Evalue):
    listeResultatsFiltre=[]
    for i in listeResultats:
        if i[2]>pident and i[3]>coverProt1 and i[4]>coverProt2 and i[5]<Evalue:
            listeResultatsFiltre.append(i)
    return listeResultatsFiltre


#Filtre les resultats de CDSearch COG en fonction des Evalue
def FilterCDSearchCOG(listeResultats,Evalue1,Evalue2):
    listeResultatsFiltre=[]
    for i in listeResultats:
        #print(i[2]<Evalue1,i[3]<Evalue2)
        if i[2]<Evalue1 and i[3]<Evalue2:
            listeResultatsFiltre.append(i)
    return listeResultatsFiltre

#Filtre les resultats de CDSearch FA en fonction des Evalue
def FilterCDSearchFA(listeResultats,Evalue1,Evalue2):
    listeResultatsFiltre=[]
    for i in listeResultats:
        #print(i[2]<Evalue1,i[3]<Evalue2)
        if i[2]<Evalue1 and i[3]<Evalue2:
            listeResultatsFiltre.append(i)
    return listeResultatsFiltre

#Calcul les blocs de synténie à partir des coordonnés
def Syntenie(listeCoordonnesX,listeCoordonnesY):
        dicoDiagonales={}
        listeCoordonnes=tuple(zip(listeCoordonnesX,listeCoordonnesY))
        listeCoordonnes=np.array(listeCoordonnes)
        listeCoordonnes=listeCoordonnes[:,0:2]
        
        #Recherche de synténies normales (diagonales croissantes)
        
        for coordonne in listeCoordonnes:
                #Recherche de diagonale croissante
                diagonaleCroissante=coordonne[0]-coordonne[1]
                #Check de l'existence de la diagonale croissante
                if diagonaleCroissante not in dicoDiagonales.keys():
                        dicoDiagonales[diagonaleCroissante]=[]
                #ajout de la coordonné
                dicoDiagonales[diagonaleCroissante].append((coordonne[0],coordonne[1]))
                
        #Enlever les diagonales qui ne contiennent qu'une seule valeur
        
        cleInutiles=[]
        for cle in dicoDiagonales.keys():
                if len(dicoDiagonales[cle])==1:
                        cleInutiles.append(cle)
        for cle in cleInutiles:
                dicoDiagonales.pop(cle)
                
        coordonnesBlocCroissant=[]
        for cle in dicoDiagonales.keys():
                dicoDiagonales[cle].sort()
                i=0
                while i < (len(dicoDiagonales[cle])-1):
                                
                                if dicoDiagonales[cle][i][1]+1==dicoDiagonales[cle][i+1][1]:
                                        blocTransitoire=[dicoDiagonales[cle][i]]
                                        i+=1
                                        while i < (len(dicoDiagonales[cle])-1) and dicoDiagonales[cle][i][1]+1==dicoDiagonales[cle][i+1][1] and dicoDiagonales[cle][i][0]+1==dicoDiagonales[cle][i+1][0]:
                                                blocTransitoire.append(dicoDiagonales[cle][i])
                                                i+=1
                                        
                                        blocTransitoire.append(dicoDiagonales[cle][i])
                                        #print(blocTransitoire)
                                        coordonnesBlocCroissant.append(blocTransitoire)
                                i+=1
        
        #Recherche de bloc de synténie inversés (bandes décroissantes)
        for coordonne in listeCoordonnes:
                #Recherche de diagonale croissante
                diagonaleDecroissant=coordonne[0]+coordonne[1]
                #Check de l'existence de la diagonale croissante
                if diagonaleDecroissant not in dicoDiagonales.keys():
                        dicoDiagonales[diagonaleDecroissant]=[]
                #ajout de la coordonné
                dicoDiagonales[diagonaleDecroissant].append((coordonne[0],coordonne[1]))
                
        #Enlever les diagonales qui ne contiennent qu'une seule valeur
        cleInutiles=[]
        for cle in dicoDiagonales.keys():
                if len(dicoDiagonales[cle])==1:
                        cleInutiles.append(cle)
        for cle in cleInutiles:
                dicoDiagonales.pop(cle)
        #print(len(dicoDiagonales.keys()))
        
        #Recherche des diagonales 
        coordonnesBlocDecroissant=[]
        
        for cle in dicoDiagonales.keys():
                dicoDiagonales[cle].sort()
                i=0
                while i < (len(dicoDiagonales[cle])-1):
                                
                                if dicoDiagonales[cle][i][1]-1==dicoDiagonales[cle][i+1][1]:
                                        blocTransitoire=[dicoDiagonales[cle][i]]
                                        i+=1
                                        while i < (len(dicoDiagonales[cle])-1) and dicoDiagonales[cle][i][1]-1==dicoDiagonales[cle][i+1][1]:
                                                blocTransitoire.append(dicoDiagonales[cle][i])
                                                i+=1
                                        
                                        blocTransitoire.append(dicoDiagonales[cle][i])
                                        #print(blocTransitoire)
                                        coordonnesBlocDecroissant.append(blocTransitoire)
                                i+=1
                                
        #Récupérer des statistiques sur les blocs de synténie (Distribution des tailles des blocs)
        listeStats=[]
        for bloc in coordonnesBlocCroissant:
                listeStats.append(len(bloc))
                
        coordonnesBlocCroissant=list(chain.from_iterable(coordonnesBlocCroissant))
        coordonnesBlocDecroissant=list(chain.from_iterable(coordonnesBlocDecroissant))
        
        return np.array(coordonnesBlocCroissant),np.array(coordonnesBlocDecroissant),listeStats

#Calcule l'intersection entre CDSearch et BLAST
def Intersection(coordonnesX1,coordonnesY1,coordonnesX2,coordonnesY2):
        coordonnesCurated=[]
        coordonnesCuratedx=[]
        coordonnesCuratedy=[]
        coordonnes1=list(zip(coordonnesX1,coordonnesY1))
        coordonnes2=list(zip(coordonnesX2,coordonnesY2))
        for i in coordonnes1:
                if i in coordonnes2:
                        coordonnesCurated.append(i)
        coordonnesCurated=np.array(coordonnesCurated)
        if len(coordonnesCurated)!=0:
                coordonnesCuratedx=coordonnesCurated[:,0]
                coordonnesCuratedy=coordonnesCurated[:,1]
        return coordonnesCuratedx,coordonnesCuratedy
        
#Vérifie que la liste des blocs est normale ou non et calcule des statistiques dessus
def NormalLaw(liste):
        value=stats.normaltest(liste)
        n=len(liste)
        mean=np.mean(liste)
        if value[1]>0.01:
                std=np.std(liste,ddof=len(np.unique(liste))-1)
                minimum=np.min(liste)
                maximum=np.max(liste)
                X=np.linspace(minimum,maximum,100)
                Y=stats.norm.pdf(X,mean,std)
                infos=(value[1],n,mean,std)
                return True,X,Y,infos
                
        else:
                q1=np.quantile(liste,0.25)
                q3=np.quantile(liste,0.75)
                infos=(value[1],n,mean,(q1,q3))
                
                return False,[],[],infos
                
#Calcul la couverture du protéome en blocs de synténie
def SyntenieCover(XSyntenyCroissant,YSyntenyCroissant,XSyntenyDecroissant,YSyntenyDecroissant,dicoProteome,choixOrg1,choixOrg2):
        Xcover=len( np.unique( np.concatenate( (XSyntenyCroissant,XSyntenyDecroissant) ) ) )/len(dicoProteome[choixOrg1])
        Ycover=len( np.unique( np.concatenate( (YSyntenyCroissant,YSyntenyDecroissant) ) ) )/len(dicoProteome[choixOrg2])
        return Xcover,Ycover
        
#Convertit les positions des protéines en leur Fonction
def Pos2FA(listeCoordonnesX,listeCoordonnesY,dicoProteome,choixOrg1,choixOrg2,dicoNoms2COG,dicoCOG):
        listeProteinesX=[dicoProteome[choixOrg1][int(x)][0] for x in listeCoordonnesX]
        listeProteinesY=[dicoProteome[choixOrg2][int(x)][0] for x in listeCoordonnesY]
        listeCOGX=[]
        listeCOGY=[]
        for i in range(len(listeProteinesX)):
                if listeProteinesX[i] in dicoNoms2COG[choixOrg1].keys() and listeProteinesY[i] in dicoNoms2COG[choixOrg2].keys():
                        listeCOGX.append(dicoNoms2COG[choixOrg1][listeProteinesX[i]])
                        listeCOGY.append(dicoNoms2COG[choixOrg2][listeProteinesY[i]]) 
        listeFonctionX=[]
        listeFonctionY=[]
        for i in range(len(listeCOGX)):
                fonctionX=dicoCOG[listeCOGX[i]][0]
                fonctionY=dicoCOG[listeCOGY[i]][0]
                listeFonctionX.append(list(fonctionX))
                listeFonctionY.append(list(fonctionY))
        return np.array(listeFonctionX),np.array(listeFonctionY)
        
def Pos2FAProteome(listeCoordonnes,dicoProteome,choixOrg,dicoNoms2COG,dicoCOG):
        listeProteines=[dicoProteome[choixOrg][int(x)][0] for x in listeCoordonnes]
        listeCOG=[dicoNoms2COG[choixOrg][x] for x in listeProteines if x in dicoNoms2COG[choixOrg].keys()]
        liste=[]
        for i in listeCOG:
                fonction=dicoCOG[i][0]
                liste.append(list(fonction))
        return np.array(liste)
        
def SyntenieContigence(X,Y,dicoProteome,choixOrg1,choixOrg2,dicoFun,dicoNoms2COG,dicoCOG):
        #coordonnées différentes entre blocs et hors blocs
        Xtemp=np.arange(len(dicoProteome[choixOrg1]))
        Ytemp=np.arange(len(dicoProteome[choixOrg2]))
        Xtemp=np.setdiff1d(Xtemp,X)
        Ytemp=np.setdiff1d(Ytemp,Y)
        #Liste des Fonctions dans les blocs de synténie
        FA1Synteny,FA2Synteny =  Pos2FA(X,Y,dicoProteome,choixOrg1,choixOrg2,dicoNoms2COG,dicoCOG)
        FA1Synteny = list(chain.from_iterable(FA1Synteny))
        FA2Synteny = list(chain.from_iterable(FA2Synteny))
        
        FA1NotSynteny = list(chain.from_iterable(Pos2FAProteome(Xtemp,dicoProteome,choixOrg1,dicoNoms2COG,dicoCOG)))
        FA2NotSynteny = list(chain.from_iterable(Pos2FAProteome(Ytemp,dicoProteome,choixOrg2,dicoNoms2COG,dicoCOG)))
        
        #Fabrication d'un tableau de contingence 
        alphabetFun=list(dicoFun.keys())
        #Creation des tables
        table1=np.zeros((len(dicoFun.keys()),2))
        table2=np.zeros((len(dicoFun.keys()),2))
        #Remplissage du tableau 1
        for i in FA1Synteny:
                table1[alphabetFun.index(i)][0]+=1
        for i in FA1NotSynteny:
                table1[alphabetFun.index(i)][1]+=1
        for i in FA2Synteny:
                table2[alphabetFun.index(i)][0]+=1
        for i in FA2NotSynteny:
                table2[alphabetFun.index(i)][1]+=1

        return table1,table2,alphabetFun

#Créer un tableau de contingence 2x2 et clacul le test de Fisher dessus, (Fonction d'intérêt et le reste)
def FisherResults(table):
        listePValue=[]        
        for i in range(len(table)):
                tupleTest=table[i]
                tupleReste=[0,0]
                for j in range(i):
                        tupleReste+=table[j]
                for j in range(i+1,len(table)):
                        tupleReste+=table[j]
                listePValue.append(stats.fisher_exact([tupleTest,tupleReste],'greater')[1])
        return listePValue

def FiltreFonctionSyntenie(X,Y,fonctionsX,fonctionsY,cle):
        XFiltree=[]
        YFiltree=[]
        for i in range(len(fonctionsX)):
                if cle in fonctionsX[i] and cle in fonctionsY[i]:
                        XFiltree.append(X[i])
                        YFiltree.append(Y[i])
        return XFiltree,YFiltree
# Idéalement j'aurais écrit ces fonctions directement dans Synténie mais n'ayant eu l'idée que plus, tard et mon programme étant déjà assez complexe j'ai préféré faire ce calcul par dessus, quitte à perdre un peu de temps

def EvalueSyntenieBLAST(X,Y,orga1,orga2,dicoBLAST,dicoProteome):
        dicoEvalue={}
        for i in range(len(X)):
                proteine1=dicoProteome[orga1][int(X[i])][0]
                proteine2=dicoProteome[orga2][int(Y[i])][0]
                #print(dicoBLAST[(orga1,orga2)][(proteine1,proteine2)])
                Evalue=dicoBLAST[(orga1,orga2)][(proteine1,proteine2)][8]
                if X[i] not in dicoEvalue.keys():
                        dicoEvalue[X[i]]=[]
                dicoEvalue[X[i]].append((Y[i],Evalue))
        return dicoEvalue
        
def EvalueSyntenieCDSearch(X,Y,orga1,orga2,dicoCDSearch,dicoProteome,dicoNoms2COG):
        dicoEvalue={}
        for i in range(len(X)):
                proteine1=dicoProteome[orga1][int(X[i])][0]
                proteine2=dicoProteome[orga2][int(Y[i])][0]
                
                liste=np.array((dicoCDSearch[orga1][dicoNoms2COG[orga1][proteine1]]))
                indice=(np.where(liste[:,4]==proteine1)[0][0])
                
                EvalueX=liste[indice][2]
                
                liste=np.array((dicoCDSearch[orga2][dicoNoms2COG[orga2][proteine2]]))
                indice=(np.where(liste[:,4]==proteine2)[0][0])
                
                EvalueY=liste[indice][2]
                
                Evalue=EvalueX+EvalueY
                
                if X[i] not in dicoEvalue.keys():
                        dicoEvalue[X[i]]=[]
                dicoEvalue[X[i]].append((Y[i],Evalue))
        return dicoEvalue
        
def FilterEvalueSyntenie(dicoEvalue):
        listeCoordonnesFiltreesX=[]
        listeCoordonnesFiltreesY=[]
        for i in dicoEvalue.keys():
                listeCoordonnesFiltreesX.append(i)
                liste=np.array(dicoEvalue[i])
                Y=dicoEvalue[i][np.argmin(liste[:,1])][0]
                listeCoordonnesFiltreesY.append(Y)
        return listeCoordonnesFiltreesX,listeCoordonnesFiltreesY
