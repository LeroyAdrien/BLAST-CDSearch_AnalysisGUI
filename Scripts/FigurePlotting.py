import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import *


#Affichage des blocs de synténie et de l'histogramme de leur statistiques '
def PlottingSyntenyBLAST(
                                coordonnesXSyntenyCroissantBLAST,
                                coordonnesYSyntenyCroissantBLAST,
                                coordonnesXSyntenyDecroissantBLAST,
                                coordonnesYSyntenyDecroissantBLAST,
                                statSyntenyBLAST,
                                normaldistribBLAST,
                                syntenyCoverBLAST,
                                couleurSyntenyBLAST,
                                styleSyntenyCroiss,
                                styleSyntenyDecroiss,
                                ax1 , ax2 , ax3
                                
                            ):
                                
        if len(coordonnesXSyntenyCroissantBLAST)>0:
                X,Y=coordonnesXSyntenyCroissantBLAST,coordonnesYSyntenyCroissantBLAST
                scatterSyntenyCroissantBLAST=ax1.scatter(X,Y,c=couleurSyntenyBLAST,s=3,marker='s')
                scatterSyntenyCroissantBLAST.set_label('Synténie BLAST')
        if len(coordonnesXSyntenyDecroissantBLAST)>0:
                X,Y=coordonnesXSyntenyDecroissantBLAST,coordonnesYSyntenyDecroissantBLAST
                scatterSyntenyDecroissantBLAST=ax1.scatter(X,Y,c=couleurSyntenyBLAST,s=3,marker=styleSyntenyDecroiss)
                #scatterSyntenyDecroissantBLAST.set_label('Synténie inversés BLAST')
        if len(statSyntenyBLAST)>0:
                echelle=np.unique(statSyntenyBLAST)
                histSyntenyBLAST=ax2.hist(statSyntenyBLAST,color=couleurSyntenyBLAST,alpha=0.5,bins=echelle,label="BLAST")
                #Si la distribution est normale 
                if normaldistribBLAST[0]:
                        X=normaldistribBLAST[1]
                        Y=normaldistribBLAST[2]
                        ax3.plot(X,Y,c=couleurSyntenyBLAST,label="Fonction de densité associée")
                        test,n,mean,std=normaldistribBLAST[3]
                        ax3.text(0.20,0.95,"Risque Alpha test de Shapiro: "+str(format(test,'.3g')),transform=ax3.transAxes,color=couleurSyntenyBLAST)
                        ax3.text(0.20,0.90,"Nombre de blocs: "+str(n),transform=ax3.transAxes,color=couleurSyntenyBLAST)
                        ax3.text(0.20,0.85,"Moyenne de longueur des blocs: "+str(round(mean,4)),transform=ax3.transAxes,color=couleurSyntenyBLAST)
                        ax3.text(0.20,0.80,"Ecart type de la longueur des blocs: "+str(round(std,3)),transform=ax3.transAxes,color=couleurSyntenyBLAST)
                        ax3.text(0.20,0.75,"Couverture du protéome X par les blocs: "+str(round(syntenyCoverBLAST[0],3)),transform=ax3.transAxes,color=couleurSyntenyBLAST)
                        ax3.text(0.20,0.70,"Couverture du protéome Y par les blocs: "+str(round(syntenyCoverBLAST[1],3)),transform=ax3.transAxes,color=couleurSyntenyBLAST)
                #Si la distribution n'est pas normale (risque alpha à.05)
                else:
                        test,n,mean,(q1,q3)=normaldistribBLAST[3]
                        ax3.text(0.20,0.95,"Valeur P test de Shapiro: "+str(format(test,'.3g')),transform=ax3.transAxes,color=couleurSyntenyBLAST)
                        ax3.text(0.20,0.90,"Nombre de blocs: "+str(n),transform=ax3.transAxes,color=couleurSyntenyBLAST)
                        ax3.text(0.20,0.85,"Moyenne de longueur des blocs: "+str(round(mean,4)),transform=ax3.transAxes,color=couleurSyntenyBLAST)
                        quantile=str("Premier et dernier quartile :"+str(round(q1,3))+'-'+str(round(q3,3)))
                        ax3.text(0.20,0.80,quantile,transform=ax3.transAxes,color=couleurSyntenyBLAST)
                        ax3.text(0.20,0.75,"Couverture du protéome X par les blocs: "+str(round(syntenyCoverBLAST[0],3)),transform=ax3.transAxes,color=couleurSyntenyBLAST)
                        ax3.text(0.20,0.70,"Couverture du protéome Y par les blocs: "+str(round(syntenyCoverBLAST[1],3)),transform=ax3.transAxes,color=couleurSyntenyBLAST)
                return ax1,ax2,ax3
                
#Affichage des blocs de syntény de CD search et de l'histogramme et leur statistiques
def PlottingSyntenyCDSearch(
                                coordonnesXSyntenyCroissantCDSearch,
                                coordonnesYSyntenyCroissantCDSearch,
                                coordonnesXSyntenyDecroissantCDSearch,
                                coordonnesYSyntenyDecroissantCDSearch,
                                statSyntenyCDSearch,
                                normaldistribCDSearch,
                                syntenyCoverCDSearch,
                                couleurSyntenyCDSearch,
                                styleSyntenyCroiss,
                                styleSyntenyDecroiss,
                                ax1 , ax2 , ax3
                                
                            ):
                                
        if len(coordonnesXSyntenyCroissantCDSearch)>0:
                X,Y=coordonnesXSyntenyCroissantCDSearch,coordonnesYSyntenyCroissantCDSearch
                scatterSyntenyCroissantCDSearch=ax1.scatter(X,Y,c=couleurSyntenyCDSearch,s=3,marker='s')
                scatterSyntenyCroissantCDSearch.set_label('Synténie CDSearch')
        if len(coordonnesXSyntenyDecroissantCDSearch)>0:
                X,Y=coordonnesXSyntenyDecroissantCDSearch,coordonnesYSyntenyDecroissantCDSearch
                scatterSyntenyDecroissantCDSearch=ax1.scatter(X,Y,c=couleurSyntenyCDSearch,s=3,marker=styleSyntenyDecroiss)
                #scatterSyntenyDecroissantCDSearch.set_label('Synténie inversés CDSearch')
        if len(statSyntenyCDSearch)>0:
                echelle=np.unique(statSyntenyCDSearch)
                histSyntenyCDSearch=ax2.hist(statSyntenyCDSearch,color=couleurSyntenyCDSearch,alpha=0.5,bins=echelle,label="CDSearch")
                # Si la distribution est normale
                if normaldistribCDSearch[0]:
                        X=normaldistribCDSearch[1]
                        Y=normaldistribCDSearch[2]
                        ax3.plot(X,Y,c=couleurSyntenyCDSearch,label="Fonction de densité associée")
                        test,n,mean,std=normaldistribCDSearch[3]
                        ax3.text(0.2,0.65,"Risque Alpha test de Shapiro: "+str(format(test,'.3g')),transform=ax3.transAxes,color=couleurSyntenyCDSearch)
                        ax3.text(0.2,0.60,"Nombre de blocs: "+str(n),transform=ax3.transAxes,color=couleurSyntenyCDSearch)
                        ax3.text(0.2,0.55,"Moyenne de longueur des blocs: "+str(round(mean,4)),transform=ax3.transAxes,color=couleurSyntenyCDSearch)
                        ax3.text(0.2,0.50,"Ecart type de la longueur des blocs :"+str(round(std,3)),transform=ax3.transAxes,color=couleurSyntenyCDSearch)
                        ax3.text(0.2,0.45,"Couverture du protéome X par les blocs: "+str(round(syntenyCoverCDSearch[0],3)),transform=ax3.transAxes,color=couleurSyntenyCDSearch)
                        ax3.text(0.2,0.40,"Couverture du protéome Y par les blocs: "+str(round(syntenyCoverCDSearch[1],3)),transform=ax3.transAxes,color=couleurSyntenyCDSearch)
                # Si la distribution n'est pas normale (risque alpha 0.05)
                else:
                        test,n,mean,(q1,q3)=normaldistribCDSearch[3]
                        ax3.text(0.2,0.65,"Valeur P test de Shapiro: "+str(format(test,'.3g')),transform=ax3.transAxes,color=couleurSyntenyCDSearch)
                        ax3.text(0.2,0.60,"Nombre de blocs: "+str(n),transform=ax3.transAxes,color=couleurSyntenyCDSearch)
                        ax3.text(0.2,0.55,"Moyenne de longueur des blocs: "+str(round(mean,4)),transform=ax3.transAxes,color=couleurSyntenyCDSearch)
                        quantile=str("Premier et dernier quartile :"+str(round(q1,3))+'-'+str(round(q3,3)))
                        ax3.text(0.2,0.50,quantile,transform=ax3.transAxes,color=couleurSyntenyCDSearch)
                        ax3.text(0.2,0.45,"Couverture du protéome X par les blocs: "+str(round(syntenyCoverCDSearch[0],3)),transform=ax3.transAxes,color=couleurSyntenyCDSearch)
                        ax3.text(0.2,0.40,"Couverture du protéome Y par les blocs: "+str(round(syntenyCoverCDSearch[1],3)),transform=ax3.transAxes,color=couleurSyntenyCDSearch)
                return ax1,ax2,ax3
                

#Affichage des blocs de syntény de l'intersection de CD search et BLAST et l'histogramme et leur statistiques
def PlottingSyntenyIntersection(
                                coordonnesXSyntenyCroissantIntersection,
                                coordonnesYSyntenyCroissantIntersection,
                                coordonnesXSyntenyDecroissantIntersection,
                                coordonnesYSyntenyDecroissantIntersection,
                                statSyntenyIntersection,
                                normaldistribIntersection,
                                syntenyCoverIntersection,
                                couleurSyntenyIntersection,
                                styleSyntenyCroiss,
                                styleSyntenyDecroiss,
                                ax1 , ax2 , ax3
                                
                            ):
                                
        if len(coordonnesXSyntenyCroissantIntersection)>0:
                X,Y=coordonnesXSyntenyCroissantIntersection,coordonnesYSyntenyCroissantIntersection
                scatterSyntenyCroissantIntersection=ax1.scatter(X,Y,c=couleurSyntenyIntersection,s=3,marker='s')
                scatterSyntenyCroissantIntersection.set_label('Synténie Intersection')
        if len(coordonnesXSyntenyDecroissantIntersection)>0:
                X,Y=coordonnesXSyntenyDecroissantIntersection,coordonnesYSyntenyDecroissantIntersection
                scatterSyntenyDecroissantIntersection=ax1.scatter(X,Y,c=couleurSyntenyIntersection,s=3,marker=styleSyntenyDecroiss)
                #scatterSyntenyDecroissantIntersection.set_label('Synténie inversés Intersection')
        if len(statSyntenyIntersection)>0:
                echelle=np.unique(statSyntenyIntersection)
                histSyntenyIntersection=ax2.hist(statSyntenyIntersection,color=couleurSyntenyIntersection,alpha=0.5,bins=echelle,label="Intersection")
                # Si la distribution est normale
                if normaldistribIntersection[0]:
                        X=normaldistribIntersection[1]
                        Y=normaldistribIntersection[2]
                        ax3.plot(X,Y,c=couleurSyntenyIntersection,label="Fonction de densité associée")
                        test,n,mean,std=normaldistribIntersection[3]
                        ax3.text(0.2,0.95,"Risque Alpha test de Shapiro: "+str(format(test,'.3g')),transform=ax3.transAxes,color=couleurSyntenyIntersection)
                        ax3.text(0.2,0.90,"Nombre de blocs: "+str(n),transform=ax3.transAxes,color=couleurSyntenyIntersection)
                        ax3.text(0.2,0.85,"Moyenne de longueur des blocs: "+str(round(mean,4)),transform=ax3.transAxes,color=couleurSyntenyIntersection)
                        ax3.text(0.2,0.80,"Ecart type de la longueur des blocs :"+str(round(std,3)),transform=ax3.transAxes,color=couleurSyntenyIntersection)
                        ax3.text(0.2,0.75,"Couverture du protéome X par les blocs: "+str(round(syntenyCoverIntersection[0],3)),transform=ax3.transAxes,color=couleurSyntenyIntersection)
                        ax3.text(0.2,0.70,"Couverture du protéome Y par les blocs: "+str(round(syntenyCoverIntersection[1],3)),transform=ax3.transAxes,color=couleurSyntenyIntersection)
                # Si la sitribution n'est pas normale (risque alpha 0.05)
                else:
                        test,n,mean,(q1,q3)=normaldistribIntersection[3]
                        ax3.text(0.2,0.95,"Valeur P test de Shapiro: "+str(format(test,'.3g')),transform=ax3.transAxes,color=couleurSyntenyIntersection)
                        ax3.text(0.2,0.90,"Nombre de blocs: "+str(n),transform=ax3.transAxes,color=couleurSyntenyIntersection)
                        ax3.text(0.2,0.85,"Moyenne de longueur des blocs: "+str(round(mean,4)),transform=ax3.transAxes,color=couleurSyntenyIntersection)
                        quantile=str("Premier et dernier quartile :"+str(round(q1,3))+'-'+str(round(q3,3)))
                        ax3.text(0.2,0.80,quantile,transform=ax3.transAxes,color=couleurSyntenyIntersection)
                        ax3.text(0.2,0.75,"Couverture du protéome X par les blocs: "+str(round(syntenyCoverIntersection[0],3)),transform=ax3.transAxes,color=couleurSyntenyIntersection)
                        ax3.text(0.2,0.70,"Couverture du protéome Y par les blocs: "+str(round(syntenyCoverIntersection[1],3)),transform=ax3.transAxes,color=couleurSyntenyIntersection)

                return ax1,ax2,ax3
                
def PlottingFonctionSynteny(
                                XCroissant,
                                YCroissant,
                                XDecroissant,
                                YDecroissant,
                                couleurFonction,
                                styleFonction,
                                ax1,
                                choixFonction
                                ):
                                
                                        
                                X=XCroissant+XDecroissant
                                Y=YCroissant+YDecroissant
                                scatterFonctionSyntenyBLAST=ax1.scatter(X,Y,c=couleurFonction,s=10,marker='s',label=choixFonction)
                                
                                return ax1

def PlottingSyntenyKariotype(
                                XCroissant,
                                YCroissant,
                                XDecroissant,
                                YDecroissant,
                                organisme1,
                                organisme2,
                                taille,
                                couleurSelected,
                                couleurTraitsCroissants,
                                couleurTraitsDecroissants,
                                ax
                             ):
                                couleurTexte='#D8D8D8'
                                #Création d'un dictionnaire de coordonnés pour accélérer la recherche)'
                                #Xinit=np.array(X)
                                #Yinit=np.array(Y)
                                XCroissant=XCroissant.astype(int)
                                YCroissant=YCroissant.astype(int)
                                dicoCroissant=dict(zip(XCroissant,YCroissant))
                                XDecroissant=XDecroissant.astype(int)
                                YDecroissant=YDecroissant.astype(int)
                                dicoDecroissant=dict(zip(XDecroissant,YDecroissant))
                                
                                #Creation des paramètres de la figure en fonction des données fournies
                                hauteurRectangles=0.03
                                longueurRectangle1=taille[0]/5000
                                longueurRectangle2=taille[1]/5000
                                
                                positionYRectangle1=0.07+hauteurRectangles*1.75
                                positionYRectangle2=positionYRectangle1+0.3
                                
                                departLignesY=positionYRectangle1+hauteurRectangles
                                arriveeLignesY=positionYRectangle2
                                
                                scalingXLignesX=0
                                
                                rectangle1=Rectangle((0,positionYRectangle1), width=longueurRectangle1,height=hauteurRectangles,color=couleurSelected)
                                rectangle2=Rectangle((0,positionYRectangle2), width=longueurRectangle2, height=hauteurRectangles,color=couleurSelected)
                                ax.add_patch(rectangle1)
                                ax.add_patch(rectangle2)
                                
                                xticks=np.arange(0,taille[0],500)
                                yticks=np.arange(0,taille[1],500)
                                proteome1=np.arange(0,taille[0],1)
                                proteome2=np.arange(0,taille[1],1)
                                #TESTS
                                #line=plt.Line2D((0,longueurRectangle1),(positionYRectangle2-hauteurRectangles*0.3,positionYRectangle2-hauteurRectangles*0.3),linewidth=0.05,color='black')
                                #ax.add_line(line)
                                #Affichage des synténies croissantes
                                for i in proteome1:
                                    if i in dicoCroissant.keys():
                                        line=plt.Line2D((i/5000,dicoCroissant[i]/5000),(departLignesY,arriveeLignesY), color='purple',linewidth=0.05)
                                        ax.add_line(line)
                                ax.text(0,0,'Synténie normale',color='purple',fontsize=3.5)
                                ax.text(0,0.03,'Synténie inversée',color='green',fontsize=3.5)
                                #Affichage des synténies décroissantes
                                for i in proteome1:
                                    if i in dicoDecroissant.keys():
                                        line=plt.Line2D((i/5000,dicoDecroissant[i]/5000),(departLignesY,arriveeLignesY), color='green',linewidth=0.05)
                                        ax.add_line(line)
                                #Affichage des ticks de la barre orga1
                                xticks=np.arange(0,taille[0],500)
                                for i in range(len(xticks)):
                                        line=plt.Line2D((xticks[i]/5000,xticks[i]/5000),(positionYRectangle1,positionYRectangle1+hauteurRectangles),linewidth=0.5,color=couleurTexte)
                                        ax.add_line(line)
                                        ax.text(xticks[i]/5000,positionYRectangle1-0.03,str(xticks[i]),color=couleurTexte,ha='center',fontsize=3)
                                #Afficahge des ticks de la barre orga2 
                                yticks=np.arange(0,taille[1],500)
                                for i in range(len(yticks)):
                                        line=plt.Line2D((yticks[i]/5000,yticks[i]/5000),(positionYRectangle2,positionYRectangle2+hauteurRectangles),linewidth=0.5,color=couleurTexte)
                                        ax.add_line(line)
                                        ax.text(yticks[i]/5000,positionYRectangle2+hauteurRectangles*1.8,str(yticks[i]),color=couleurTexte,ha='center',fontsize=3)
                                
                                ax.text(longueurRectangle1/2, positionYRectangle1-0.12,organisme1, fontsize=5,color=couleurTexte,ha='center')
                                
                                ax.text(longueurRectangle2/2, positionYRectangle2+hauteurRectangles*3.2,organisme2, fontsize=5,color=couleurTexte,ha='center')
                                #ax.axis("scaled")
                                ax.axis("off")
                                

