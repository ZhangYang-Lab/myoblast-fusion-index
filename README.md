# myoblast-fusion-index
README（Quantification of myoblast fusion index）
Overview
This MATLAB code is designed to calculate the myoblast fusion index from confocal fluorescence images. The fusion index is determined by analyzing two channels: Myosin (indicating myotubes) and DAPI (indicating cell nuclei). The fusion index is calculated by dividing the number of Myosin-positive nuclei by the total number of nuclei in the image.
Features
Data Loading: Load confocal fluorescence images with two channels (Myosin and DAPI).
Image Processing: Identify and segment cell nuclei in both channels.
Fusion Index Calculation: Count Myosin-positive nuclei and total nuclei to compute the fusion index.
Visualization: Visualize the segmentation and counting results.
Requirements
MATLAB R2021a or later
Image Processing Toolbox
Statistics and Machine Learning Toolbox
Usage
Data Preparation: Ensure your confocal images are in a compatible format (e.g., TIFF, PNG) and contain separate channels for Myosin and DAPI.
Loading Data: Use the load Images function to load your imaging data.
Image Segmentation: The segment Nuclei function will process the images to identify and segment cell nuclei.
Fusion Index Calculation: The calculate Fusion Index function will compute the fusion index based on the segmented nuclei.
Visualization: Use the visualize Results function to visualize the segmentation and counting results.
Notes
Ensure that your images are properly pre-processed (e.g., noise reduction, contrast enhancement) before running the analysis.
The segmentation accuracy can significantly affect the fusion index results. Fine-tune the segmentation parameters as needed.
The visualization function helps in verifying the segmentation and counting accuracy.
Contact
For any questions or issues, please contact guohs@szbl.ac.cn.
