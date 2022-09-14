# Dataset 1.
Zoomed-in pupil videos. Black frames are piezo flyback.
Illumination: difracted 2P laser light
CMad58, 62, 65, 67, 68 (AC axons, mouse disco)
RD10278, CMad85, CMad86 (LM axons) - mobnet for CMad85 nad 86?

## Run DeepLabCut (ResNet)
(Pupil fitted with 8 points)

## Run FitEllipseAndSortData to:
('setpath_EyeTracking.m')

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
CMad (AC axons, Bl6 mice, 2-80 kHz sound, mouse disco)
CMad (AC axons, CBA mice, 2-80 kHz sound, mouse disco)
CMad (TEA axons, Bl6 mice, 2-80 kHz sound, mouse disco)

