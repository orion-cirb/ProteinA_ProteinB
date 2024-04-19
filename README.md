# Proteins_Segmentation 

* **Developed for:** Félix
* **Team:** Rouach
* **Date:** April 2024
* **Software:** Fiji

### Images description

3D images

2 channels: 
  1. Protein A
  2. Protein B (optional)

With each image, a *.roi* or *.zip* file can be provided containing one or multiple ROI(s). Each ROI should be named ..._z... to indicate its z-slice position. Analysis will be performed within it, as well as X slices before and X slices after, with the value of X provided in the dialog box. If no ROI is provided, analysis will be conducted on the entire image.

### Plugin description

* Segment protein A channel using median filtering + thresholding + median filtering
* If available, segment protein B channel using the same approach
* Compute the background noise in each provided channel
* Provide protein volume and background-corrected mean intensity in each provided channel

### Dependencies

* **3DImageSuite** Fiji plugin
* **CLIJ2** Fiji plugin

### Version history

Version 1 released on April 19, 2024.
