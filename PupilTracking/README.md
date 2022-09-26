# Dataset 1.
Zoomed-in pupil videos. Black frames are piezo flyback.
Illumination: difracted 2P laser light
CMad58, 62, 65, 67, 68 (AC axons, mouse disco)
CMad73,74 (AC axons, sound at different loudness)
RD10278, CMad85, CMad86 (LM axons) - mobnet for CMad85 and 86?

## Run DeepLabCut (ResNet)
Model: 'AdAxons_MD75dB_8fittedPts-Camille-2021-07-09'
(Pupil fitted with 8 points)

## Run FitEllipseAndSortData to:
(give the path to the right folder to Matlab: 'setpath_EyeTracking.m')

1) Fit an ellipse to these 8 points. 
DLC accuracy > 0.8
DLC accuracy avg > 0.9

2) Sort the data according to 'trialtimeline'
'trialtimeline' file needs to be in the folder containing the postions (.csv)

## Once all the data is sorted, run 
'AnalyzePerAnimal.m'
'BayesianDecoder_Pupil.m' --> /!\ are there trials full of NaNs? in this case you obviously can't decode

# Dataset 2.
Full-face videos.
Illumination using IR light. Filter on the camera
CMad97-99 (AC axons, Bl6 mice, 2-80 kHz sound, mouse disco)
CMad103,104,106,107,109,110,111 (AC axons, CBA mice, 2-80 kHz sound, mouse disco)
CMad119,120 (TEA axons, Bl6 mice, 2-80 kHz sound, mouse disco)

## Run DeepLabCut (ResNet)
Model: 'FaceVideos_ACaxonsProject-Cam-2022-09-15'
(Pupil fitted with 8 points, corneal reflexion (1 point) and nasal and temporal eye position)

## Run FitEllipseAndSortData to:
(give the path to the right folder to Matlab: 'setpath_EyeTracking.m')

1) Fit an ellipse to these 8 points. 
DLC accuracy > 0.8
DLC accuracy avg > 0.9

2) Sort the data according to 'trialtimeline'
'trialtimeline' file needs to be in the folder containing the postions (.csv)

## Once all the data is sorted, run 
'AnalyzePerAnimal.m'
'BayesianDecoder_Pupil.m' --> /!\ are there trials full of NaNs? in this case you obviously can't decode
