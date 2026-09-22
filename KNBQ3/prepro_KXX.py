import os
import shutil

def lire_fichier_agrif(nom_fichier):
    """Extrait les données des lignes commentées contenant KAGRIF"""
    structure = {}
    try:
        with open(nom_fichier, 'r') as f:
            lignes = [l.strip() for l in f if l.strip().startswith('#')]

        # Recherche de la ligne KAGRIF
        index_kagrif = next((i for i, l in enumerate(lignes) if "KAGRIF" in l), -1)
        
        if index_kagrif == -1:
            return structure

        # Traitement des lignes suivantes
        for ligne in lignes[index_kagrif + 1:]:
            elements = ligne[1:].strip().split()  # Ignore le '#'
            if len(elements) < 2:
                continue  # Ligne incomplète
            #boucle sur la liste des fichiers .F (tous en x exemplaires)
            cle, *sous_cles = elements
            structure['step3d_fast.F'] = {
                'cle_principale': cle,
                'cles_secondaires': sous_cles
            }
            print(cle," ==> ",sous_cles)
        print(structure)
        return structure

    except FileNotFoundError:
        print(f"ERREUR: Fichier {nom_fichier} introuvable")
        return {}

def traiter_cppdefs(structure, fichier_source='cppdefs.h'):
    """Génère les fichiers modifiés"""
    fichiers_cpp = []
    for fichier, donnees in structure.items():
        for sous_cle in donnees['cles_secondaires']:
            nouveau_fichier = f"cppdefs_{sous_cle}.h"
            
            try:
                # Création du fichier modifié
                with open(fichier_source, 'r') as src, open(nouveau_fichier, 'w') as dest:
                    contenu = src.read().replace(
                        donnees['cle_principale'],
                        sous_cle
                    )
                    dest.write(contenu)

                fichiers_cpp.append(nouveau_fichier)
                print(f"SUCCÈS: {nouveau_fichier} créé avec la clé {sous_cle}")

            except FileNotFoundError:
                print(f"ERREUR: Fichier source {fichier_source} manquant")
    return fichiers_cpp


def creer_fichiers_suffixes(structure, fichiers_crees, suffixe_base='.F'):
    """Crée des copies suffixées des fichiers avec les remplacements spécifiés."""
    nouveaux_fichiers = []

    for nom_fichier, donnees in structure.items():
        print(nom_fichier)
        print(donnees)
        # Préparation du nom de base sans extension .F
        base, extension = os.path.splitext(nom_fichier)
        if extension.upper() == '.F':
            nom_base = base
        else:
            nom_base = nom_fichier

        for sous_cle in donnees['cles_secondaires']:
            # Création du nouveau nom de fichier
            nouveau_nom = f"{nom_base}_{sous_cle}{suffixe_base}"

            try:
                # Copie et modification du fichier
                with open(nom_fichier, 'r') as src, open(nouveau_nom, 'w') as dest:
                    contenu = src.read()

                    # Remplacement 1: cppdefs.h -> cppdefs_[sous_cle].h
                    contenu = contenu.replace('cppdefs.h', f'cppdefs_{sous_cle}.h')

                    # Remplacement 2: nom_base -> nom_base_[sous_cle]
                    contenu = contenu.replace(nom_base, f'{nom_base}_{sous_cle}')

                    dest.write(contenu)

                nouveaux_fichiers.append(nouveau_nom)
                print(f"Fichier suffixé créé : {nouveau_nom}")

            except FileNotFoundError:
                print(f"Erreur: Fichier source {nom_fichier} introuvable")
            except Exception as e:
                print(f"Erreur lors du traitement de {nom_fichier} : {str(e)}")

    return nouveaux_fichiers


def modifier_makefile(structure, fichiers_suffixes, makefile='Makefile'):
    """Modifie le Makefile en remplaçant les noms de fichiers originaux par leurs versions suffixées"""
    fichiers_modifies = []

    try:
        # Création de la sauvegarde
        shutil.copy2(makefile, makefile + '.orig')
        print(f"Backup créé : {makefile}.orig")

        # Lecture du contenu original
        with open(makefile + '.orig', 'r') as f:
            contenu = f.read()

        # Remplacement pour chaque fichier de la structure
        for fichier_original in structure:
            base = os.path.splitext(fichier_original)[0]
            fichiers_associes = [f for f in fichiers_suffixes if f.startswith(base + '_')]

            if fichiers_associes:
                remplacement = ' '.join(fichiers_associes)
                contenu = contenu.replace(fichier_original, remplacement)
                print(f"Remplacé : {fichier_original} → {remplacement}")

        # Écriture du nouveau Makefile
        with open(makefile, 'w') as f:
            f.write(contenu)

        fichiers_modifies.append(makefile)
        return fichiers_modifies

    except FileNotFoundError:
        print("ERREUR : Makefile introuvable")
        return []
    except Exception as e:
        print(f"ERREUR lors de la modification : {str(e)}")
        return []


# Exécution
if __name__ == "__main__":
    # Étape 1: Lecture de la configuration
    config = lire_fichier_agrif("../AGRIF_FixedGrids.in")
    #config['pre_step3d.F'] = config['step3d_fast.F']
    print(config)
    print(config['step3d_fast.F'])
    print(config['step3d_fast.F']['cles_secondaires'])
    #on Compile toujours les 4 noyaux
    config['step3d_fast.F']['cles_secondaires']=['KHCOMP', 'KNBQ', 'KH3D', 'KNHINT']
    #config['pre_step3d.F'] = config['step3d_fast.F']
    
    # Étape 2: Génération des cppdefs
    cppdefs_gen = traiter_cppdefs(config)

    # Étape 3: Création des fichiers suffixés
    fichiers_suffixes = creer_fichiers_suffixes(config, cppdefs_gen)

    # 4. Modification du Makefile
    makefiles = modifier_makefile(config, fichiers_suffixes)

    # Rapport final
    print("\nRapport de génération :")
    print(f"- {len(cppdefs_gen)} fichiers cppdefs créés")
    print(f"- {len(fichiers_suffixes)} fichiers suffixés générés")
    print(f"- {len(makefiles)} Makefile modifié")
