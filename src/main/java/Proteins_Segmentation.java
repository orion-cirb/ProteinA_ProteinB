import Proteins_Segmentation_Tools.Tools;
import ij.*;
import ij.gui.Roi;
import ij.plugin.PlugIn;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.SimpleDateFormat;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import loci.common.DebugTools;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.in.ImporterOptions;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang.ArrayUtils;


/**
* Segment two proteins in their respective channel
* Compute various area and intensity measurements
* Perform analysis in provided ROIs only
* @author Héloïse Monnet
*/
public class Proteins_Segmentation implements PlugIn {

    private Proteins_Segmentation_Tools.Tools tools = new Tools();
       
    public void run(String arg) {
        try {
            if ((!tools.checkInstalledModules())) {
                return;
            }
            
            String imageDir = IJ.getDirectory("Choose images directory");
            if (imageDir == null) {
                IJ.showMessage("", "Plugin canceled");
                return;
            }
            
            // Find images with fileExt extension
            String fileExt = tools.findImageType(imageDir);
            ArrayList<String> imageFiles = tools.findImages(imageDir, fileExt);
            if (imageFiles.isEmpty()) {
                IJ.showMessage("ERROR", "No images found with " + fileExt + " extension");
                return;
            }
            
            // Create OME-XML metadata store of the latest schema version
            DebugTools.setRootLevel("warn");
            ServiceFactory factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            
            // Find image calibration
            tools.findImageCalib(meta);
            
            // Find channel names
            String[] chMeta = tools.findChannels(imageFiles.get(0), meta, reader);
            
            String[] chOrder = tools.dialog(chMeta);
            if (chOrder == null) {
                IJ.showMessage("", "Plugin canceled");
                return;
            } else if(chOrder[0] == "None") {
                IJ.showMessage("ERROR", "Protein A channel not defined.");
                return;
            }
            
            // Create output folder
            String thMethods = (!chOrder[1].equals("None"))? tools.protAThMethod + "_" + tools.protBThMethod : tools.protAThMethod;
            String outDirResults = imageDir + File.separator + "Results_" + thMethods + "_" + new SimpleDateFormat("yyyy-MM-dd_HH-mm-ss").format(new Date()) + File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            
            // Write headers results for results files
            FileWriter fwResults = new FileWriter(outDirResults + "results.csv", false);
            BufferedWriter results = new BufferedWriter(fwResults);
                results.write("Image name\tImage vol (µm3)\tROI name\tROI vol (µm3)\tROI slice position\tROI slices nb\tProtein A bg\tProtein A volume (µm3)\t" +
                        "Protein A bg-corr mean int");
            if(!chOrder[1].equals("None"))
                results.write("\tProtein B bg\tProtein B volume (µm3)\tProtein B bg-corr mean int");
            results.write("\n");
            results.flush();
            
            for (String f: imageFiles) {
                String rootName = FilenameUtils.getBaseName(f);
                tools.print("--- ANALYZING IMAGE " + rootName + " ------");
                reader.setId(f);
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setSplitChannels(true);
                options.setQuiet(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                
                // Open Protein A channel
                tools.print("- Opening Protein A channel -");
                int index = ArrayUtils.indexOf(chMeta, chOrder[0]);
                ImagePlus imgProtA = BF.openImagePlus(options)[index];
                
                // Open Protein B channel, if provided
                ImagePlus imgProtB = null;
                if (!chOrder[1].equals("None")) {
                    tools.print("- Opening Protein B channel -");
                    index = ArrayUtils.indexOf(chMeta, chOrder[1]);
                    imgProtB = BF.openImagePlus(options)[index];
                } else
                    System.out.println("WARNING: No Protein B channel provided");
                
                
                // Load ROIs (if provided)
                tools.print("- Loading ROIs -");
                List<Roi> rois = tools.loadRois(imageDir + File.separator + rootName, imgProtA);
                
                // Analyze Protein A channel
                tools.print("- Analyzing Protein A channel -");
                double bgProtA = tools.computeBackgroundNoise(imgProtA);
                ImagePlus segProtA = tools.segmentation(imgProtA, tools.protAThMethod, tools.protAStackHistogram);
                
                // Analyze Protein B channel
                double bgProtB = 0;
                ImagePlus segProtB = null;
                if(imgProtB != null) {
                    tools.print("- Analyzing Protein B channel -");
                    bgProtB = tools.computeBackgroundNoise(imgProtB);
                    segProtB = tools.segmentation(imgProtB, tools.protBThMethod, tools.protBStackHistogram);
                }
                
                // Save results for each ROI
                tools.print("- Saving results -");
                double imgVol = imgProtA.getWidth() * imgProtA.getHeight() * imgProtA.getNSlices() * tools.pixVol;
                ImageHandler resProtA = ImageHandler.wrap(imgProtA).createSameDimensions();
                ImageHandler resProtB = null;
                if(imgProtB != null)
                    resProtB = ImageHandler.wrap(imgProtB).createSameDimensions();
                for(Roi roi: rois) {
                    double roiVol = tools.getRoiVolume(roi, imgProtA);
                    
                    Object3DInt objProtA = tools.getObjectInsideRoi(segProtA, roi);
                    double volProtA = new MeasureVolume(objProtA).getVolumeUnit();
                    double meanIntProtA = new MeasureIntensity(objProtA, ImageHandler.wrap(imgProtA)).getValueMeasurement(MeasureIntensity.INTENSITY_AVG) - bgProtA;

                    Object3DInt objProtB = null;
                    double volProtB = 0, meanIntProtB = 0;
                    if(imgProtB != null) {
                        objProtB = tools.getObjectInsideRoi(segProtB, roi);
                        volProtB = new MeasureVolume(objProtB).getVolumeUnit();
                        meanIntProtB = new MeasureIntensity(objProtB, ImageHandler.wrap(imgProtB)).getValueMeasurement(MeasureIntensity.INTENSITY_AVG) - bgProtB;
                    }
                    
                    // Write results
                    results.write(rootName+"\t"+imgVol+"\t"+roi.getName()+"\t"+roiVol+"\t"+roi.getZPosition()+"\t"+roi.getProperty("zNb")+"\t"+bgProtA+"\t"+volProtA+"\t"+meanIntProtA);
                    if(imgProtB != null)
                         results.write("\t"+bgProtB+"\t"+volProtB+"\t"+meanIntProtB);
                    results.write("\n");
                    results.flush();
                    
                    // Draw results
                    objProtA.drawObject(resProtA, 255);
                    if(imgProtB != null)
                        objProtB.drawObject(resProtB, 255);
                }
                
                // Draw results
                tools.drawResults(resProtA, resProtB, imgProtA, imgProtB, outDirResults+rootName+".tif");
                
                tools.closeImage(imgProtA);
                tools.closeImage(segProtA);
                resProtA.closeImagePlus();
                if(imgProtB != null) {
                    tools.closeImage(imgProtB);
                    tools.closeImage(segProtB);
                    resProtB.closeImagePlus();
                }
            }
            results.close();
        } catch (IOException | DependencyException | ServiceException | FormatException ex) {
            Logger.getLogger(Proteins_Segmentation.class.getName()).log(Level.SEVERE, null, ex);
        }
        tools.print("All done!");
    }
}
