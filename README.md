# PIR

Ce projet a pour but de détermniner les coordonnées d'un point dans l'espace sur une image prise par plusieurs caméras sur Raspberry Pi. Ensuite, nous allons utiliser ces coordonnées pour déterminer la position d'avion dans l'espace pour la détection de contrails.

## Prérequis

- Python 3 (https://www.python.org/downloads/)
- Astropy (https://www.astropy.org/)


## Installation

1. Cloner le projet
2. Configurer l'environnement en ajouter la configuration suivante dans le fichier `.env` à la racine du projet :

```bash
ASTROMETRY_API_KEY=VOTRE_API_KEY_ASTROMETRY
```

Exemple de fichier `.env` :

```bash
ASTROMETRY_API_KEY="123456789"
```

## Utilisation

Pour lancer le projet général, il suffit de lancer le script `__main__.py` avec la commande suivante :

```bash
python3 src/__main__.py
```

Pour lancer un script spécifique, il suffit de lancer le script avec la commande suivante :

```bash
python3 src/chemin/vers/le/script.py
```

### Tests

Lancer les tests avec la commande suivante :

```bash
python3 -m unittest -v
```

Example OK:
```bash
python3 -m unittest -v 
test_convert_coordinates (test.test_coordinates_photo.TestCoordinatesPhoto) ... ok

----------------------------------------------------------------------
Ran 1 test in 0.001s

OK
```

Example KO:
```bash
python3 -m unittest -v 
test_convert_coordinates (test.test_coordinates_photo.TestCoordinatesPhoto) ... FAIL

======================================================================
FAIL: test_convert_coordinates (test.test_coordinates_photo.TestCoordinatesPhoto)
----------------------------------------------------------------------
Traceback (most recent call last):
  File "/home/ronan/Documents/INSA/4A/Modules/PIR/PIR/test/test_coordinates_photo.py", line 17, in test_convert_coordinates
    self.assertEqual(object_location, (35.6895, 139.691))
AssertionError: Tuples differ: (35.6895, 139.6917) != (35.6895, 139.691)

First differing element 1:
139.6917
139.691

- (35.6895, 139.6917)
?                  -

+ (35.6895, 139.691)

----------------------------------------------------------------------
Ran 1 test in 0.001s

FAILED (failures=1)
```

## Fonctionnement

- [ ] Système de prise de photo sur le Raspberry Pi.
- [ ] Calibration de la caméra pour déterminer les coordonnées visibles sur l'image. (Astrometry.net)
- [ ] Établissement d'un lien entre les coordonnées visibles sur l'image et les coordonnées terrestres/GPS (lat, long, alt).
- [ ] Détection d'avions sur les images.
- [ ] Estimation la position des avions dans l'espace pour la détection de contrails.
- [ ] Lien entre ces avions et des vols en cours. (Flightradar24)



## Auteurs

- [Achraf Bensebaa]()
- [Idriss Bensouda]()
- [Ronan Bonnet](https://github.com/BloodFutur)
- [Anna Cazeneuve]()
- [Abderrahman El Ouali]()
- [Miguel Fernadez-Cid Castano]()
- [Luz Vera Morales]()
