# Analysis

## Extract
- 'Leaf Area' (aux échelles GU, Branch, Tree) 
-  'number of inflorescences' (aux échelles Branch, Tree)

Plant_id 	Branch_id 	GU_id 	Leaf area 2018 	Leaf area 2019 	Total number of FL
p1 	 +B1 	 +C1 	53.924 	
	0 	
p1 	 +B1 	 +C2 	
	221.68 	1 	
p1 	 +B1 	 +C3 	71.98 	
	0 	
… 	… 	… 	… 	… 	… 	

	
	
	
	
	
	

Je fais un recap des derniers échanges, pour clarifier l’"état de l'art":

- (4 mars) PE et FR proposent une liste des variables à extraire des mtgs; un dossier partagé avec la dernière versionne du papier et des FAQs pour l’interprétation des mtgs

- (10 avril) CP envoi une extraction à l’échelle de la plant pour tout les arbres (trees2.csv).
Dans les échanges qui suive on clarifie que:

    la variable trt n'est pas correcte
    Taper et Elongation (calculé pour le tronc uniquement):
    élongation: Total Trunk length / mean_diameter
    mean_diameter (diametre moyen du tronc) : (A1.diameter_b + A3.diameter_a) / 2. #mean(Top trunk  diameter, Base Trunk diameter)
    taper : Base Trunk diameter - Top trunk  diametre/ trunk length
    Pour la signification des colonnes F-M: L : long, S: Short, F: Flowering, V: Vegetative. Par exemple S_shoot_F_2018 : Number of Short Shoot in 2018
    La distinction entre Short et Long est : short = (Length < 5)
    Pour les colonnes N à S :
    nb_V : nb vegetative bud
    nb_F: nb floral buds
    nb_L : nb Latent bud

(je vais intégrer tout ça dans la description du materiel dans l'article, de que les fichiers "extractions" seront "stables").

- (3 juin) PE et FR ont vérifié la bonté des extractions faites par CP.

On se met d'accord que il faudrait extraire:

    le nombre de GU, (done)
	la leaf area, (done)
	le nombre de leaves, flowers et fruits, pour chaque ordre de ramification (A, B, ...) (done)
    les angles et diamètres basal (ramification B)
		min, max median angle
    et de corriger la variable trt
		todo
    on verra à la fin si c'est possible enlever les sommes déjà calculées sur les lignes 2017


Nombre GU


# Updater analyse_p1 dans test
#

# Identifier les angles > 180. : 360. - angle