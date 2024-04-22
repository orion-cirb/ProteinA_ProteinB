package Proteins_Segmentation_Tools;

import fiji.util.gui.GenericDialogPlus;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.PolygonRoi;
import ij.gui.Roi;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.plugin.filter.Analyzer;
import ij.plugin.frame.RoiManager;
import ij.process.AutoThresholder;
import java.awt.Color;
import java.awt.Font;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DInt;
import mcib3d.image3d.ImageHandler;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import org.apache.commons.io.FilenameUtils;


/**
 * @author Héloïse Monnet
 */
public class Tools {
    
    public final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    private final String helpUrl = "https://github.com/orion-cirb/Proteins_Segmentation";
    
    private final CLIJ2 clij2 = CLIJ2.getInstance();
    
    String[] chNames = {"Protein A: ", "Protein B (optional): "};
    public Calibration cal = new Calibration();
    public double pixVol;
    
    // Number of slices analyzed before and after each ROI
    public int nbSlices = 1;
    
    // Protein A segmentation
    public String protAThMethod = "Default";
    public boolean protAStackHistogram = true;
    
    // Protein B segmentation
    public String protBThMethod = "Default";
    public boolean protBStackHistogram = true;
    
    
    /**
     * Display a message in the ImageJ console and status bar
     */
    public void print(String log) {
        System.out.println(log);
        IJ.showStatus(log);
    }
    
    
    /**
     * Check that needed modules are installed
     */
    public boolean checkInstalledModules() {
        // check install
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.log("3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        try {
            loader.loadClass("net.haesleinhuepf.clij2.CLIJ2");
        } catch (ClassNotFoundException e) {
            IJ.log("CLIJ2 not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    
    /**
     * Find images extension
     */
    public String findImageType(String imagesFolder) {
        String ext = "";
        String[] files = new File(imagesFolder).list();
        for (String name : files) {
            String fileExt = FilenameUtils.getExtension(name);
            switch (fileExt) {
                case "nd" :
                   ext = fileExt;
                   break;
                case "nd2" :
                   ext = fileExt;
                   break;
                case "czi" :
                   ext = fileExt;
                   break;
                case "lif"  :
                    ext = fileExt;
                    break;
                case "ics" :
                    ext = fileExt;
                    break;
                case "ics2" :
                    ext = fileExt;
                    break;
                case "lsm" :
                    ext = fileExt;
                    break;
                case "tif" :
                    ext = fileExt;
                    break;
                case "tiff" :
                    ext = fileExt;
                    break;
            }
        }
        return(ext);
    }

        
    /**
     * Find images in folder
     */
    public ArrayList findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No image found in "+imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt) && !f.startsWith("."))
                images.add(imagesFolder + f);
        }
        Collections.sort(images);
        return(images);
    }
    
    
    /**
     * Find image calibration
     */
    public Calibration findImageCalib(IMetadata meta) {
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
        return(cal);
    }
    
    
    /**
     * Find channels name
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels(String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs+1];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelName(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelName(0, n).toString();
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n).toString();
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = (meta.getChannelFluor(0, n).toString().equals("")) ? Integer.toString(n) : meta.getChannelFluor(0, n).toString();
                break;
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break;    
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    channels[n] = meta.getChannelEmissionWavelength(0, n).value().toString();
                break; 
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        channels[chs] = "None";
        return(channels);     
    }
    
    
    /**
     * Generate dialog box
     */
    public String[] dialog(String[] channels) {
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 60, 0);
        gd.addImage(icon);
        
        gd.addMessage("Channels", new Font("Monospace", Font.BOLD, 12), Color.blue);
        for (int n = 0; n < chNames.length; n++) {
            gd.addChoice(chNames[n], channels, channels[n]);
        }
        
        String[] thMethods = AutoThresholder.getMethods();
        gd.addMessage("Proteins thresholding", new Font("Monospace", Font.BOLD, 12), Color.blue);
        gd.addNumericField("Nb of slices b/a ROI: ", nbSlices, 0);
        
        gd.addMessage("Protein A", new Font("Monospace", Font.PLAIN, 12), Color.blue);
        gd.addChoice("Method: ",thMethods, protAThMethod);
        gd.addCheckbox("Stack histogram", protAStackHistogram);
        
        gd.addMessage("Protein B", new Font("Monospace", Font.PLAIN, 12), Color.blue);
        gd.addChoice("Method: ",thMethods, protBThMethod);
        gd.addCheckbox("Stack histogram", protBStackHistogram);
        
        gd.addMessage("Image calibration", new Font("Monospace", Font.BOLD, 12), Color.blue);
        gd.addNumericField("XY calibration (µm): ", cal.pixelHeight, 3);
        gd.addNumericField("Z calibration (µm): ", cal.pixelDepth, 3);
        gd.addHelp(helpUrl);
        gd.showDialog();
        
        String[] chChoices = new String[chNames.length];
        for (int n = 0; n < chChoices.length; n++) 
            chChoices[n] = gd.getNextChoice();
        
        nbSlices = (int) gd.getNextNumber();
        
        protAThMethod = gd.getNextChoice();
        protAStackHistogram = gd.getNextBoolean();
        
        protBThMethod = gd.getNextChoice();
        protBStackHistogram = gd.getNextBoolean();
        
        cal.pixelHeight = cal.pixelWidth = gd.getNextNumber();
        cal.pixelDepth = gd.getNextNumber();
        pixVol = cal.pixelHeight*cal.pixelWidth*cal.pixelDepth;
        
        if (gd.wasCanceled())
            return(null);
        return(chChoices);
    }
    
    
    /**
     * Flush and close an image
     */
    public void closeImage(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    /**
     * Load ROIs, if any provided
     */
    public List<Roi> loadRois(String roiName, ImagePlus img) {
        List<Roi> rois = new ArrayList<>();
        
        roiName = new File(roiName+".zip").exists() ? roiName+".zip" : roiName+".roi";
        if(new File(roiName).exists()) {
            RoiManager rm = new RoiManager(false);
            rm.runCommand("Open", roiName);
            List<Roi> roisTemp = Arrays.asList(rm.getRoisAsArray());
            
            for(Roi roi: roisTemp) {
                if(roi.getZPosition() > 0) {
                    int z = roi.getZPosition();
                    int zStart = (z - nbSlices < 1)? 1 : z - nbSlices;
                    int zStop = (z + nbSlices > img.getNSlices())? img.getNSlices() : z + nbSlices;
                    int zNb = zStop - zStart + 1;
                    roi.setProperty("zStart", String.valueOf(zStart));
                    roi.setProperty("zStop", String.valueOf(zStop));
                    roi.setProperty("zNb", String.valueOf(zNb));
                    rois.add(roi);
                } else {
                    IJ.showMessage("ERROR", "ROI " + roiName + " will not be analyzed, as it is not associated with a particular slice.");
                }
            }
        }
        
        if(rois.isEmpty()) {
            Roi roi = new Roi(0, 0, img.getWidth(), img.getHeight());
            roi.setName("entire image");
            roi.setPosition((int) 0.5*img.getNSlices());
            roi.setProperty("zStart", "1");
            roi.setProperty("zStop", String.valueOf(img.getNSlices()));
            roi.setProperty("zNb", String.valueOf(img.getNSlices()));
            rois.add(roi);
            System.out.println("WARNING: No ROI file provided, entire image will be analyzed");
        }

        return(rois);
    }
    
    
    /**
     * Compute image background noise:
     * z-project over min intensity + read median intensity
     */
    public double computeBackgroundNoise(ImagePlus img) {
      ImagePlus imgProj = zProject(img, ZProjector.MIN_METHOD);
      double bg = imgProj.getProcessor().getStatistics().median;
      System.out.println("Background noise (median of the min projection) = " + bg);
      closeImage(imgProj);
      return(bg);
    }
    
    
    /**
     * Z-projection a stack
     */
    public ImagePlus zProject(ImagePlus img, int param) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(param);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
       return(zproject.getProjection());
    }
    
       
    /**
     * Segment stack in 2D with median filtering + thresholding
     */
    public ImagePlus segmentation(ImagePlus img, String thMethod, boolean stackHistogram) {
        ImagePlus imgMed = median3DSliceBySlice(img, 2);
        ImagePlus imgTh = threshold(imgMed, thMethod, stackHistogram);
        ImagePlus imgOut = median3DSliceBySlice(imgTh, 2);
        imgOut.setCalibration(cal);
        
        closeImage(imgMed);
        return(imgOut);
    }
        
    
    /**
     * 2D median filtering slice by slice using CLIJ2
     */ 
    public ImagePlus median3DSliceBySlice(ImagePlus img, double sizeXY) {
       ClearCLBuffer imgCL = clij2.push(img); 
       ClearCLBuffer imgCLMed = clij2.create(imgCL);
       clij2.median3DSliceBySliceSphere(imgCL, imgCLMed, sizeXY, sizeXY);
       ImagePlus imgMed = clij2.pull(imgCLMed);
       clij2.release(imgCL);
       clij2.release(imgCLMed);
       return(imgMed);
    }
    
    
    /**
     * Automatic thresholding using CLIJ2
     */
    public ImagePlus threshold(ImagePlus img, String thMed, boolean stackHistogram) {
        ImagePlus imgOut = img.duplicate();
        if(stackHistogram) {
            IJ.setAutoThreshold(imgOut, thMed + " dark stack");
            IJ.run(imgOut, "Convert to Mask", "method=" + thMed + " background=Dark");
        } else {
            IJ.setAutoThreshold(imgOut, thMed + " dark");
            IJ.run(imgOut, "Convert to Mask", "method=" + thMed + " background=Dark calculate");
        }
        return(imgOut);
    }
    
    
    /**
     * Compute ROI volume
     */
    public double getRoiVolume(Roi roi, ImagePlus img) {
        PolygonRoi poly = new PolygonRoi(roi.getFloatPolygon(), Roi.FREEROI);
        poly.setLocation(0, 0);

        img.resetRoi();
        img.setRoi(poly);

        ResultsTable rt = new ResultsTable();
        Analyzer analyzer = new Analyzer(img, Analyzer.AREA, rt);
        analyzer.measure();
        double roiArea = rt.getValue("Area", 0);

        return(roiArea * Integer.valueOf(roi.getProperty("zNb")) * cal.pixelDepth);
    }
    
    
    /**
     * Clear mask outside ROI 
     * Measure volume and mean intensity of mask inside ROI
     */
    public Object3DInt getObjectInsideRoi(ImagePlus mask, Roi roi) {
        ImagePlus maskClear = mask.duplicate();
        maskClear.getProcessor().setColor(Color.BLACK);
        for (int s = 1; s <= mask.getNSlices(); s++) {
            maskClear.setSlice(s);
            if(s >= Integer.valueOf(roi.getProperty("zStart")) && s <= Integer.valueOf(roi.getProperty("zStop")))
                maskClear.getProcessor().fillOutside(roi);
            else
                maskClear.getProcessor().fill();
        }
        
        Object3DInt obj = new Object3DInt(ImageHandler.wrap(maskClear));
        return(obj);
    }
    
    
    /**
     * Draw results
     */
    public void drawResults(ImageHandler resProtA, ImageHandler resProtB, ImagePlus imgProtA, ImagePlus imgProtB, String name) {
        ImagePlus[] imgColors;
        if(imgProtB != null)
            imgColors = new ImagePlus[]{resProtA.getImagePlus(), null, resProtB.getImagePlus(), imgProtA, imgProtB};
        else
            imgColors = new ImagePlus[]{resProtA.getImagePlus(), null, null, imgProtA};
        
        ImagePlus imgObjects = new RGBStackMerge().mergeHyperstacks(imgColors, false);
        imgObjects.setCalibration(cal);
        
        FileSaver ImgObjectsFile = new FileSaver(imgObjects);
        ImgObjectsFile.saveAsTiff(name); 
    }
    
}
