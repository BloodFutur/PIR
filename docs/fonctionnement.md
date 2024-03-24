# Fonctionnement pour la calibration

1. Prendre une photo de nuit
2. Mettre la photo en noir et blanc
3. De préférence, crop la photo pour ne garder que le ciel
4. Envoie la photo sur Astrometry.net (À voir à implémenter en local)
5. Récupérer les données de calibration
6. Utiliser les données de calibration pour déterminer les coordonnées visibles sur la photo
7. Utiliser coordonnées du centre de l'image dans le système celeste
8. Convertir les coordonnées en coordonnées terrestres (lat, long, alt)

## Flow de données

input: image
-> traitement de l'image par Astrometry.net qui donne les coordonnées des étoiles visibles sur l'image + une calibration pour convertir les coordonnées: (x, y) -> (ra, dec)
-> utilisation des coordonnées du centre de l'image pour déterminer les coordonnées terrestres (azimuth, altitude)
-> conversion des coordonnées terrestres (azimuth, altitude) en coordonnées terrestres (lat, long, alt)
-> output: coordonnées terrestres/GPS (lat, long, alt)

# Resources

Pour query chez astrometry.net, on peut utiliser astroquery:
https://astroquery.readthedocs.io/en/latest/astrometry_net/astrometry_net.html