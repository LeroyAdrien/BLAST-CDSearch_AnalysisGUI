# *-*coding:utf-8 -*-
import numpy as np
import os
import fnmatch
import re

def FileParser(path):
    
    matchCDDir=False
    CDDirPattern='*[Cc][Dd]*'
    CDSearchDirPath=''
    
    matchProteomeDir=False
    ProteomeDirPattern='*[Pp][Rr][Oo][Tt][éEe][Oo][Mm][Ee]*'
    ProteomeDirPath=''
    
    matchBLASTDir=False
    BLASTDirPattern='*[Bb][Ll][Aa][Ss][Tt]*'
    BLASTDirPath=''
    
    matchCOGNamesFile=False
    COGNamesFilePattern="*[Cc][Oo][Gg]*.txt"
    COGNamesFilePath=''
    
    matchGenomesFile=False
    GenomeFilePattern="*[Gg][Ee][Nn][Oo][Mm][Ee]*.csv"
    GenomeFilePath=''
    
    matchFunFile=False
    FunFilePattern="*[Ff][Uu][Nn]*.txt"
    FunFilePath=''
    
    listeNoms=np.array(["CDSearch",
                        "Proteome",
                        "BLAST",
                        "COGNames",
                        "Genome",
                        "Fun"
                       ])

    #Recherche des dossiers et fichiers dans le répertoire "fournis"
    with os.scandir(path) as entries:
        for entry in entries:
            
            #Verification des dossiers:
            if entry.is_dir():
                
                #Si le dossier CD Search est là
                if fnmatch.fnmatch(entry,CDDirPattern):
                    CDSearchFilePaths=[]
                    matchCDDir=True
                    
                    #Récupération des fichiers du dossier CDSearch
                    
                    for CDSearchEntry in os.scandir(entry.path):
                        
                        if CDSearchEntry.is_file():
                            CDSearchFilePaths.append(CDSearchEntry.path)
                
                # Si le dossier Proteome est là   
                if fnmatch.fnmatch(entry,ProteomeDirPattern):
                    ProteomeFilePaths=[]
                    matchProteomeDir=True
                    
                    #Récupération des fichiers du dossier Proteome
                    for proteomeEntry in os.scandir(entry.path):
                        
                        if proteomeEntry.is_file():
                            ProteomeFilePaths.append(proteomeEntry.path)
                 
                #Si le dossier BLAST est là
                if fnmatch.fnmatch(entry,BLASTDirPattern):
                    BLASTFilePaths=[]
                    matchBLASTDir=True
                    
                    #Récupération des fichiers du dossier BLAST
                    for BLASTEntry in os.scandir(entry.path):
                        
                        if BLASTEntry.is_file():
                            
                            BLASTFilePaths.append(BLASTEntry.path)
                    
            #Verification de la présence des fichiers COGNames et Genome Complete      
            if entry.is_file():
                
                #Si COGNames est là 
                if fnmatch.fnmatch(entry,COGNamesFilePattern):
                    COGNamesFilePath=entry.path
                    matchCOGNamesFile=True
                    
                #Si Genome est là   
                if fnmatch.fnmatch(entry,GenomeFilePattern):
                    GenomeFilePath=entry.path
                    matchGenomesFile=True
                
                #Si Fun est là
                if fnmatch.fnmatch(entry,FunFilePattern):
                    FunFilePath=entry.path
                    matchFunFile=True
                    
    #Vérification que tout est bien là               
    listeVerif=np.array([matchCDDir,matchProteomeDir,
                         matchBLASTDir,
                         matchCOGNamesFile,
                         matchGenomesFile,
                         matchFunFile
                        ])
    
    
    # Si tout est bien là
    if all(listeVerif):
        listePath=np.array([CDSearchFilePaths,
                            ProteomeFilePaths,
                            BLASTFilePaths,
                            COGNamesFilePath,
                            GenomeFilePath,
                            FunFilePath
                           ])
        #Dictionnaire des chemins des fichiers d'intérêt
        dictPath=dict(zip(listeNoms,listePath))
        
        
        # Si les bon nombres de fichiers sont là 
        
        if len(dictPath["Proteome"])/2 == len(dictPath["CDSearch"])/2 == len(dictPath["BLAST"]):
            return dictPath
        
        # Si il en manque dans un dossier
        else:
            
            #Print le nom du dossier où il en manque
            indice=np.argmin(np.array([len(dictPath["Proteome"])/2,len(dictPath["CDSearch"])/2,len(dictPath["BLAST"])]))
            print("Pas assez de fichiers dans le dossier:",listeNoms[0:3][indice])
            pass
    #Si il manque quelque chose
    else:
        
        #Retour des dossier/fichiers manquants:
        indicesErreur=np.where(listeVerif==False)[0]
        print("Les dossiers ou fichiers suivants sont manquants:",listeNoms[indicesErreur])
        pass

#Fonctions de Lectures des différents types de fichiers

def ReadProteome(nameFile):
    flux=open(nameFile,'r')
    lignes=flux.readlines()
    i=0
    listeSeq=[]
    lenSeqRegex=re.compile('\[location=(complement)?[\(]?([0-9]+)\.+([0-9]+)')

    #Recherche de l'id d'un gène/protéine et stockage dans nomSéquence
    for ligne in lignes:
        if ligne[0]=='>':
            nomSequence=ligne[5:-1]
            match=lenSeqRegex.search(ligne)
            if match!=None:
                longueurSequence=(int(match.group(3))-(int(match.group(2))-1))/3
            
            listeSeq.append([nomSequence.split()[0],longueurSequence])
    flux.close()
    return listeSeq

#Lecture d'un fichier CDSearch et récupération sous forme de listes de COG, et liste d'infos et nom par ligne
def ReadCDSearch(nameFile):
    flux=open(nameFile,'r')
    lignes=flux.readlines()
    listeResultats=[]
    listeCOG=[]
    listeNoms=[]
    for ligne in lignes:
        if ligne[0:2]=="Q#":
            listeLigne=ligne.split('\t')
            if listeLigne[7][0:2]!="cl":
                nomProt=listeLigne[0].split()
                listeInfos=[listeLigne[3],
                            listeLigne[4],
                            listeLigne[5],
                            listeLigne[6],
                            nomProt[2][5:]
                            ]
                listeNoms.append(nomProt[2][5:])
                listeResultats.append(listeInfos)
                listeCOG.append(listeLigne[7][0:])
                
    flux.close()
    return listeCOG,listeResultats,listeNoms

#Lecture d'un fichier BLAST et récupération sous forme de liste de noms, et liste d'infos par ligne
def ReadBlast(nameFile):
    flux=open(nameFile,'r')
    lignes=flux.readlines()
    listeNoms=[]
    listeResultats=[]
    for ligne in lignes:
        if ligne[0]!="#":
            listeLigne=ligne.split('\t')
            nomProt=(listeLigne[0][4:],listeLigne[1])
            listeInfos=listeLigne[2:-1]
            listeInfos.append(listeLigne[-1][:-1])
            
            listeNoms.append(nomProt)
            listeResultats.append(listeInfos)
    flux.close()      
    return listeNoms,listeResultats
    
def ReadGenome(nameFile):
    flux=open(nameFile,'r')
    lignes=flux.readlines()
    listeNoms=[]
    listeResultats=[]
    for ligne in lignes:
        if ligne[0]!="#":
            listeLigne=ligne.split(',')
            nomEsp=listeLigne[5]
            listeInfos=[listeLigne[i] for i in [0,1,2,7,8,12,14,15]]
            
            listeNoms.append(nomEsp[1:-1])
            listeResultats.append(listeInfos)
    flux.close()      
    return listeNoms,listeResultats

def ReadCOGNames(nameFile):
    flux=open(nameFile,'r',encoding='windows-1252')
    lignes=flux.readlines()
    listeNoms=[]
    listeResultats=[]
    for ligne in lignes:
        if ligne[0]!="#":
            listeLigne=ligne.split('\t')
            nomCOG=listeLigne[0]
            listeInfos=[listeLigne[1],listeLigne[2][:-1]]
            
            listeNoms.append(nomCOG)
            listeResultats.append(listeInfos)
    flux.close()      
    return listeNoms,listeResultats

def ReadFun(nameFile):
    flux=open(nameFile,'r')
    lignes=flux.readlines()
    listeNomsCourt=[]
    listeNomsLongs=[]
    for ligne in lignes:
        if ligne[0]!="#":
            listeLigne=ligne.split('\t')
            nomCourt=listeLigne[0]
            nomLong=listeLigne[1][:-1]
            listeNomsCourt.append(nomCourt)
            listeNomsLongs.append(nomLong)
    flux.close()      
    return listeNomsCourt,listeNomsLongs

#Fonction d'extraction de tous les fichiers
def FileExtractor(dicoFiles):
  
    dicoProteome={}
    dicoCDSearch={}
    dicoNoms2COG={}
    dicoBLAST={}
    dicoGenome={}
    dicoCOG={}
    
    #Retrouver le nom de l'organisme dans le nom de fichier
    patternOrganism=re.compile('([A-Za-z]{3}_[0-9]{9}[\.]?[0-9]{0,2})+')
    
    #Extraction du Protéome:
    for file in dicoFiles["Proteome"]:
        listeSeq=ReadProteome(file)
        organismName=patternOrganism.findall(file)
        dicoProteome[organismName[0]]=listeSeq

    #Extraction de CDSearch
    for file in dicoFiles["CDSearch"]:
        listeCOG,listeSeq,listeNoms=ReadCDSearch(file)
        organismName=patternOrganism.findall(file)
        dicoCDSearch[organismName[0]]={}
        for i in range(len(listeCOG)):
                if listeCOG[i] not in dicoCDSearch[organismName[0]].keys():
                        dicoCDSearch[organismName[0]][listeCOG[i]]=[]
                dicoCDSearch[organismName[0]][listeCOG[i]].append(listeSeq[i])
        dicoNoms2COG[organismName[0]]=dict(zip(listeNoms,listeCOG))
                
        

    #Extraction de Blastp
    for file in dicoFiles["BLAST"]:
        listeNoms,listeSeq=ReadBlast(file)
        organismName=patternOrganism.findall(file)
        dicoBLAST[organismName[0],organismName[1]]=dict(zip(listeNoms,listeSeq))

    #Extraction du génome
    listeNoms,listeSeq=ReadGenome(dicoFiles["Genome"])
    dicoGenome=dict(zip(listeNoms,listeSeq))

    #Extraction des COGNames
    listeCOG,listeFonction=ReadCOGNames(dicoFiles["COGNames"])
    dicoCOG=dict(zip(listeCOG,listeFonction))
    
    #Extraction de Fun
    listeCode,listeNom=ReadFun(dicoFiles["Fun"])
    dicoFun=dict(zip(listeCode,listeNom))

    return dicoProteome,dicoCDSearch,dicoBLAST,dicoGenome,dicoCOG,dicoFun,dicoNoms2COG

