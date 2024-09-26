
import ij.plugin.filter.PlugInFilter;
import java.util.*;
import java.awt.Color;
import java.awt.Frame;
import java.awt.Font;
import java.io.*;

import ij.*;
import ij.gui.*;
import ij.io.*;
import ij.process.*;
import ij.measure.*;
import ij.text.*;
import ij.plugin.filter.ParticleAnalyzer;
import ij.plugin.filter.BackgroundSubtracter;

/**
 * wrMTrck by Jesper Søndergaard Pedersen (JSP)
 * 
 * Uses ImageJ's particle analyzer to track the movement of
 * multiple objects through a stack. Uses the changes in aspect
 * ratio of elipse fitting to quantify "thrashing" or
 * swimming motion of C. elegans worms
 * 
 * Based on MTrack2 by 2003/06/28 Nico Stuurman
 * https://imagescience.org/meijering/software/mtrackj/
 * 
 * History
 * 
 * Build 090408
 * Implemented dos.newLine() instead of hard-coded "\n"
 * Changed maxAreaChange to a percentage instead of pixel value
 * Build 090707
 * Added summarize feature
 * Added Body-bend display on showPositions
 * Added posibility to make an FPS override
 * Now Works with scale-calibrated images
 * Build 090804
 * Fixed a bug in output of raw data.
 * Build 091115
 * Added x-y point smoothing function to decrease the "noise" movement of
 * immobile objects
 * Build 100201
 * Corrected spelling mistake in input parameter "minTrackLength" (was
 * "minTrackLenght") - this may affect scripts.
 * Build 100411
 * Corrected names in some fields on the dialog
 * Build 100706
 * Option to output plots of thrashing analysis for quality control purposes to
 * rapidly asses if threshold value is good
 * Build 101004
 * Now works on both binary and thresholded movies.
 * Posibility to enable on the fly deflickering and simple background
 * subtraction based on RB15 on first frame.
 * This means that wrMTrck can now track animals in raw movies and virtual
 * stacks.
 * Hitting the ESC key aborts analysis
 * Build 110617
 * Added an extra input field to allow adjustment of label font size
 * Build 110622
 * Added output of bend-time histogram (in frames) for each animal in text file
 * output. Switch on/off histogram by adjusting bendDetect parameter
 * Build 111031
 * Corrected undocumented raw data output with rawData =6
 * With rawData <0 headers are suppressed.
 * (Konstantine Palanski) Summary repors frames with maximum and minimum number
 * of objects and the frame with maximum number of tracked objects
 * 
 * TODO
 * (Allow user to provide a thresholded image-stack instead of requiring
 * binary-stack) - done
 * Better support for multicore CPU to speed up analysis (Not sure how much can
 * be gained)
 * Make output of histograms in more user-friendly format
 * Flexible user defineable point smoothing setting instead of hard coded
 * 5-point smoothing
 * 
 * 
 */
public class wrMTrck_ implements PlugInFilter, Measurements {

	private boolean verbose = IJ.debugMode;
	ImagePlus imp;
	int nParticles;
	static String buildNo = "Build 111031";
	static int fontSize = 20;
	static int minSize = 10;
	static int maxSize = 400;
	static float maxVelocity = 10;
	static float maxAreaChange = 30;
	static int minTrackLength = 50;
	static float minAngle = 2;
	static float binSize = 0;
	static int plotBendTrack = 0;
	static boolean bSaveResultsFile = false;
	static boolean bShowLabels = true;
	static boolean bShowPositions = true;
	static boolean bShowPaths = true;
	static boolean bShowPathLengths = true;
	static boolean bShowSummary = true;
	static boolean bRoundCoord = false;
	static boolean bSmoothing = true;
	static boolean bPlotBendTrack = false;

	static int rawData = 0;
	static int bendType = 2;
	static double fps = 0;
	static int backSub = 0;
	static String threshMode = "Otsu";
	static int maxColumns = 75;

	// KP
	static int kpFrm = 0, kpMax = 0;
	static String row = new String();
	private int nMax = 0;
	private int nMin = 99999;
	private int nMaxFrm = 0;
	private int nMinFrm = 0;
	// EOKP

	private TextWindow tw;
	private String directory, filename, rawFilename;

	private static String prevHdr;
	private String summaryHdr = "File\tnObj\tnObjFrm\tnFrames\tnTracks\ttotLength\tObjFrames\tObjSeconds\tavgSpeed\tavgArea\tavgPerim\tstdSpeed\tstdArea\tstdPerim"; // (KP)
	private double pixelWidth = 1.0, pixelHeight = 1.0;
	boolean done;

	// #region Particle Class TODO:

	public class particle {
		float x;
		float y;
		float sx; // smoothed x-coordinate
		float sy; // smoothed y-coordinate
		int z;
		float area; // inserted by JSP
		float major;
		float minor;
		float ar;
		float angle;
		float circularity;
		float dist; // contains the distance traveled since last frame.
		float perimeter; // length of the perimeter of a object
		int trackNum;
		boolean inTrack = false;
		boolean flag = false;

		public void copy(particle source) {
			this.x = source.x;
			this.y = source.y;
			this.z = source.z; // frame
			this.major = source.major;
			this.minor = source.minor;
			this.ar = source.ar;
			this.area = source.area;
			this.angle = source.angle;
			this.circularity = source.circularity;
			this.dist = source.dist;
			this.perimeter = source.perimeter;
			this.inTrack = source.inTrack;
			this.flag = source.flag;
		}

		public float distance(particle p) {
			return (float) Math.sqrt(sqr(this.x - p.x) + sqr(this.y - p.y));
		}

		public float sdistance(particle p) {
			return (float) Math.sqrt(sqr(this.sx - p.sx) + sqr(this.sy - p.sy));
		}
	}

	public class Track {
		List<particle> particles = new ArrayList<particle>();
		int startFrame;
		int endFrame;
		float distance;

		public int length() {
			return (this.endFrame - this.startFrame);
		}
	}

	public int setup(String arg, ImagePlus imp) {
		this.imp = imp;
		if (IJ.versionLessThan("1.17y"))
			return DONE;
		else
			return DOES_8G + NO_CHANGES;
	}

	// ===================================================== max
	public static double max(double[] t) {
		double maximum = t[0]; // start with the first value
		for (int i = 1; i < t.length; i++) {
			if (t[i] > maximum) {
				maximum = t[i]; // new maximum
			}
		}
		return maximum;
	}// end method max

	public void run(ImageProcessor ip) {

		GenericDialog gd = new GenericDialog("wrMTrck by Jesper S. Pedersen, " + buildNo);
		gd.addNumericField("minSize - Minimum Object Area (pixels^2): ", minSize, 0);
		gd.addNumericField("maxSize - Maximum Object Area (pixels^2): ", maxSize, 0);
		gd.addNumericField("maxVelocity - Maximum Velocity (pixels/frame):", maxVelocity, 0);
		gd.addNumericField("maxAreaChange - Maximum area change (%):", maxAreaChange, 0);
		gd.addNumericField("minTrackLength - Minimum track length (frames):", minTrackLength, 0);
		gd.addNumericField("bendThreshold - Threshold for turn :", minAngle, 1);
		gd.addNumericField("binSize - Size of bin for speed histogram (pixels/frame) (0=disable):", binSize, 1);
		gd.addCheckbox("saveResultsFile - Save Results File:", bSaveResultsFile);
		gd.addCheckbox("showPathLengths - Display Path Lengths:", bShowPathLengths);
		gd.addCheckbox("showLabels - Show Labels:", bShowLabels);
		gd.addCheckbox("showPositions - Show Positions:", bShowPositions);
		gd.addCheckbox("showPaths - Show Paths:", bShowPaths);
		gd.addCheckbox("showSummary - show a summary of tracking", bShowSummary);
		gd.addCheckbox("roundCoord - round off coordinates", bRoundCoord);
		gd.addCheckbox("smoothing - point smoothing", bSmoothing);
		gd.addCheckbox("plotBendTrack - Quality control plots for thrashing analysis", bPlotBendTrack);
		gd.addNumericField("rawData - (0=off,1=XYcord,2=Ellipse,3=AreaPerimDist,4=Ellipse+Circ,5=BendCalc):", rawData, 0);
		gd.addNumericField("bendDetect - (0=Off,1=Angle,2=AspectRatio,3=AR+Histogram):", bendType, 0);
		gd.addNumericField("FPS - frames/s (0=try to load from file):", fps, 0);
		gd.addNumericField("backSub - On-the-fly background subtraction (0=off,1=F1RB15):", backSub, 0);
		gd.addChoice("threshMode - Thresholding method (only if backSub>0)", AutoThresholder.getMethods(), "Otsu");
		gd.addNumericField("fontSize - Size of labeling font:", fontSize, 0);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		minSize = (int) gd.getNextNumber();
		maxSize = (int) gd.getNextNumber();
		maxVelocity = (float) gd.getNextNumber();
		maxAreaChange = (float) gd.getNextNumber();
		minTrackLength = (int) gd.getNextNumber();
		minAngle = (float) gd.getNextNumber();
		binSize = (float) gd.getNextNumber();
		bSaveResultsFile = gd.getNextBoolean();
		bShowPathLengths = gd.getNextBoolean();
		bShowLabels = gd.getNextBoolean();
		bShowPositions = gd.getNextBoolean();
		bShowPaths = gd.getNextBoolean();
		bShowSummary = gd.getNextBoolean();
		bRoundCoord = gd.getNextBoolean();
		bSmoothing = gd.getNextBoolean();
		bPlotBendTrack = gd.getNextBoolean();
		rawData = (int) gd.getNextNumber();
		bendType = (int) gd.getNextNumber();
		fps = (float) gd.getNextNumber();
		backSub = (int) gd.getNextNumber();
		threshMode = gd.getNextChoice();
		fontSize = (int) gd.getNextNumber();
		if (bShowPositions)
			bShowLabels = true;

		/* Extract file name of movie file opened in ImageJ */

		if (rawData < 0) {
			rawData = -rawData;
		} else {
		}
		;
		rawFilename = imp.getTitle();
		int dotPos = rawFilename.lastIndexOf('.');
		if (0 < dotPos && dotPos < rawFilename.length() - 2) {
			rawFilename = rawFilename.substring(0, dotPos);
		}
		if (bSaveResultsFile) {
			SaveDialog sd = new SaveDialog("Save Track Results", rawFilename, ".txt");
			directory = sd.getDirectory();
			filename = sd.getFileName();
		}
		if (done)
			return;

		track(imp, minSize, maxSize, maxVelocity, directory, filename);
	}

	@SuppressWarnings({ "deprecation", "rawtypes", "unchecked" })
	public void track(ImagePlus imp, int minSize, int maxSize, float maxVelocity, String directory, String filename) {
		verbose = true; // Era: For debug
		int nFrames = imp.getStackSize();
		if (nFrames < 2) {
			IJ.showMessage("Tracker", "Stack required");
			return;
		}
		IJ.log("Tracking start.");
		IJ.log("Min particle size: " + minSize);
		IJ.log("Max particle size: " + maxSize);
		IJ.log("Max particle velocity: " + maxVelocity);
		IJ.log("Max area change: " + maxAreaChange);
		if (verbose)
			IJ.log("Verbose.");
		Calibration cal = imp.getCalibration();
		pixelWidth = cal.pixelWidth;
		pixelHeight = cal.pixelHeight;

		// Attempt to extract the frame rate for the movie if user did not provide a fps
		// value

		if (fps <= 0) {
			fps = cal.fps;
		}
		;

		// Get information on the thresholding used

		ImageStack stack = imp.getStack();
		ImageProcessor ip = imp.getProcessor();
		double minThresh = ip.getMinThreshold();
		double maxThresh = ip.getMaxThreshold();

		if (verbose)
			IJ.log("min=" + minThresh + ",max=" + maxThresh);

		// Set options for particle analyser and measurements

		int options = (ParticleAnalyzer.EXCLUDE_EDGE_PARTICLES
				// |ParticleAnalyzer.INCLUDE_HOLES
				| ParticleAnalyzer.SHOW_NONE); // set all PA options false
		int measurements = CENTROID + AREA + ELLIPSE + PERIMETER + CIRCULARITY;

		// Initialize results table
		ResultsTable rt = new ResultsTable();
		rt.reset();

		// create storage for particle positions
		List[] particlesInFrame = new ArrayList[nFrames];
		int trackCount = 0;

		// Code required for deflickering and background subtraction on the fly
		double roiAvg[] = new double[nFrames + 1];
		BackgroundSubtracter bs = new BackgroundSubtracter();

		ip = stack.getProcessor(1);
		ip.resetRoi();
		if (verbose) {
			if (ip.isBinary()) {
				IJ.log("The image is binary");
			} else {
				IJ.log("The image is NOT binary");
			}
		}
		ImageProcessor bip = ip.crop();
		if ((backSub > 0) & !ip.isBinary()) { // subtract background if requested, but not if the image is already binary.
			bs.rollingBallBackground(bip, 15, true, true, false, false, true);
		}

		// record particle positions for each frame in an ArrayList
		IJ.showStatus("Finding all particles...");
		IJ.log("Finding all particles.");

		for (int iFrame = 1; iFrame <= nFrames; iFrame++) {

			// Do frame normalization and background subtraction on the fly

			ip = stack.getProcessor(iFrame);
			if (backSub > 0) {
				ImageStatistics is = ImageStatistics.getStatistics(ip, Measurements.MEAN, imp.getCalibration());
				roiAvg[iFrame] = is.mean;
				ip.resetRoi();
				ip.multiply(roiAvg[1] / roiAvg[iFrame]);
				ip.copyBits(bip, 0, 0, Blitter.DIFFERENCE);
				// bs.rollingBallBackground(ip, 15, false, true, false, false , true);
			}
			if (!(minThresh == -808080.0)) { // NO_THRESHOLD )) {
				ip.setThreshold(minThresh, maxThresh, 0);
				// IJ.log("preset threshold");
			} else if ((backSub > 0) & !ip.isBinary()) {
				AutoThresholder at = new AutoThresholder();
				ip.setThreshold(at.getThreshold(threshMode, ip.getHistogram()), 255, 0);

				minThresh = ip.getMinThreshold();
				maxThresh = ip.getMaxThreshold();
				if (verbose)
					IJ.log("AutoThresholding with " + threshMode + ": min=" + minThresh + ",max=" + maxThresh);
			}

			particlesInFrame[iFrame - 1] = new ArrayList();
			rt.reset();
			ParticleAnalyzer pa = new ParticleAnalyzer(options, measurements, rt, minSize, maxSize);
			if (!pa.analyze(imp, ip))
				return;
			float[] sxRes = rt.getColumn(ResultsTable.X_CENTROID);
			float[] syRes = rt.getColumn(ResultsTable.Y_CENTROID);
			float[] areaRes = rt.getColumn(ResultsTable.AREA);
			float[] angleRes = rt.getColumn(ResultsTable.ANGLE);
			float[] majorRes = rt.getColumn(ResultsTable.MAJOR);
			float[] minorRes = rt.getColumn(ResultsTable.MINOR);
			float[] perimRes = rt.getColumn(ResultsTable.PERIMETER);

			if (sxRes == null)
				;
			// IJ.log("Frame"+iFrame+" contains no objects!");//return;// do not exit! there
			// may be objects on the next slice image!!!
			else {
				// IJ.log("Frame"+iFrame+" contains "+sxRes.length+" objects");
				for (int iPart = 0; iPart < sxRes.length; iPart++) {
					particle particle = new particle();
					if (bRoundCoord) {
						particle.x = Math.round(sxRes[iPart]);// sxRes[iPart]; //
						particle.y = Math.round(syRes[iPart]);// syRes[iPart];//
					} else {
						particle.x = sxRes[iPart];
						particle.y = syRes[iPart];
					}

					particle.area = areaRes[iPart];
					particle.angle = angleRes[iPart];
					particle.major = majorRes[iPart];
					particle.minor = minorRes[iPart];
					if (minorRes[iPart] > 0) {
						particle.ar = majorRes[iPart] / minorRes[iPart];
					}
					;
					particle.perimeter = perimRes[iPart];
					// particle.circularity=circRes[iPart];
					particle.z = iFrame - 1;
					particlesInFrame[iFrame - 1].add(particle);
				}
				if (sxRes.length > nMax) {
					nMax = sxRes.length; // record the maximum number of particles found in one frame.
					nMaxFrm = iFrame; // record the frame at which the maximum number of particles is found. (KP)
				}
				if (sxRes.length < nMin) {
					nMin = sxRes.length; // record the minimum number of particles found in one frame. (KP)
					nMinFrm = iFrame; // record the frame at which the minimum number of particles is found. (KP)
				}
			}
			IJ.showProgress((double) iFrame / nFrames);
			if (IJ.escapePressed()) {
				IJ.beep();
				done = true;
				return;
			}
			// IJ.log((int)iFrame+",");
		}

		// #region Track Assembly

		// TODO: Fix here

		IJ.showStatus("Assembling tracks");
		IJ.log("Assembling tracks");

		List<Track> tracks = new ArrayList<Track>();
		for (int i = 0; i <= (nFrames - 1); i++) {
			IJ.showProgress((double) i / nFrames);

			// Abandon if user has pressed escape
			if (IJ.escapePressed()) {
				IJ.beep();
				done = true;
				return;
			}

			// Ignore frame if it contains no particles (not sure about this)
			if (particlesInFrame[i].size() <= 0)
				continue;

			// IJ.log("Frame " + i + " has " + particlesInFrame[i].size() + " particles");

			for (ListIterator particleIterator = particlesInFrame[i].listIterator(); particleIterator.hasNext();) {
				particle particle = (particle) particleIterator.next();

				// Skip a particle if it is already in a track (not sure about this)
				if (particle.inTrack)
					continue;

				// This must be the beginning of a new track
				Track track = new Track();
				track.startFrame = i;
				track.endFrame = 0;
				trackCount++;

				// Add particle to track
				particle.inTrack = true;
				particle.trackNum = trackCount;
				track.particles.add(particle);

				// search in next frames for more particles to be added to track
				boolean trackExists = true;
				particle prevParticle = particle;
				for (int j = i + 1; j <= (nFrames - 1) && trackExists; j++) {
					List<particle> candidates = new ArrayList<>();
					List<Float> candidateDistances = new ArrayList<>();
					for (ListIterator nextParticleIterator = particlesInFrame[j].listIterator(); nextParticleIterator
							.hasNext();) {
						particle testParticle = (particle) nextParticleIterator.next();
						float distance = testParticle.distance(prevParticle);
						float areaChange = 100 * Math.abs(1 - (testParticle.area / prevParticle.area));
						boolean fitsInPath = (distance / pixelWidth < maxVelocity) && (areaChange < maxAreaChange);
						if (fitsInPath) {
							candidates.add(testParticle);
							candidateDistances.add(distance);
						}
					}
					if (candidates.size() == 0) {
						trackExists = false;
						track.endFrame = j - 1;
						IJ.log("Track " + trackCount + " ended at frame " + i + " after " +
								track.particles.size() + " frames");
					} else {
						particle bestCandidate = getBestCandidate(candidates, candidateDistances);
						bestCandidate.inTrack = true;
						bestCandidate.trackNum = trackCount;
						track.particles.add(bestCandidate);
						prevParticle = bestCandidate;
					}
				}
				if (track.endFrame == 0) {
					track.endFrame = nFrames - 1;
				}
				tracks.add(track);
			}
		}

		// #region Raw tracks
		IJ.showStatus("Generating map of raw tracks...");
		IJ.log("Generating map of raw tracks...");
		if (bShowPaths) {
			ip = new ColorProcessor(imp.getWidth(), imp.getHeight());
			ip.setColor(Color.white);
			ip.fill();
			for (ListIterator<Track> tracksIterator = tracks.listIterator(); tracksIterator.hasNext();) {
				Track track = tracksIterator.next();
				ListIterator<particle> particleIterator = track.particles.listIterator();
				particle oldParticle = particleIterator.next();

				if (track.particles.size() < minTrackLength) {
					ip.setColor("grey");
				} else {
					ip.setColor("blue");
				}

				while (particleIterator.hasNext()) {
					particle newParticle = particleIterator.next();
					ip.moveTo((int) (oldParticle.x / pixelWidth), (int) (oldParticle.y / pixelHeight));
					ip.lineTo((int) (newParticle.x / pixelWidth), (int) (newParticle.y / pixelHeight));
					oldParticle = newParticle;
				}
			}
			new ImagePlus("All Paths", ip).show();
		}

		// #region Track Analysis
		float maxStitchDistance = 25f;
		IJ.log("Stitching tracks with max dist: " + maxStitchDistance);
		List<Track> activeTracks = new ArrayList<>();
		List<Track> lostTracks = new ArrayList<>();
		for (int frame = 0; frame < nFrames; frame++) {
			List<Track> tracksStartedThisFrame = getTracksStartedInFrame(tracks, frame);

			// If this is the first frame, set all the tracks to active
			if (frame == 0) {
				activeTracks.addAll(tracksStartedThisFrame);
				continue;
			}

			// Try to stitch the new tracks onto the lost tracks
			Map<Integer, Track> stitchMap = new HashMap<Integer, Track>();
			for (int i = 0; i < lostTracks.size(); i++) {
				Track lostTrack = lostTracks.get(i);
				particle last = lostTrack.particles.get(lostTrack.particles.size() - 1);
				Track bestMatch = null;
				float bestDistance = maxStitchDistance;
				for (int j = 0; j < tracksStartedThisFrame.size(); j++) {
					Track newTrack = tracksStartedThisFrame.get(j);
					particle first = newTrack.particles.get(0);
					float distance = first.distance(last);
					if (distance < bestDistance) {
						bestMatch = newTrack;
						bestDistance = distance;
					}
				}
				if (bestMatch != null) {
					tracksStartedThisFrame.remove(bestMatch);
					stitchMap.put(i, bestMatch);
				}
			}

			Iterator<Map.Entry<Integer, Track>> stitchIterator = stitchMap.entrySet().iterator();
			List<Track> found = new ArrayList<>();
			while (stitchIterator.hasNext()) {
				Map.Entry<Integer, Track> entry = stitchIterator.next();
				Integer index = entry.getKey();
				Track lostTrack = lostTracks.get(index);
				Track newTrack = entry.getValue();
				Track stitchedTrack = new Track();
				stitchedTrack.startFrame = lostTrack.startFrame;
				stitchedTrack.endFrame = newTrack.endFrame;
				stitchedTrack.particles.addAll(lostTrack.particles);
				stitchedTrack.particles.addAll(newTrack.particles);
				found.add(lostTrack);
				activeTracks.add(stitchedTrack);
			}

			for (Track t : found)
				lostTracks.remove(t);

			// Add any active tracks that end on this frame to Lost
			Iterator<Track> activeIterator = activeTracks.iterator();
			List<Track> removals = new ArrayList<>();
			while (activeIterator.hasNext()) {
				Track track = activeIterator.next();
				if (track.endFrame == frame) {
					removals.add(track);
					lostTracks.add(track);
				}
			}
			for (Track t : removals)
				activeTracks.remove(t);
		}

		List<Track> allTracks = new ArrayList<>();
		allTracks.addAll(activeTracks);
		allTracks.addAll(lostTracks);

		// #region Processed tracks
		IJ.showStatus("Generating map of stitched tracks...");
		IJ.log("Generating map of stitched tracks...");
		if (bShowPaths) {
			List<Integer> smoothingLevels = Arrays.asList(1, 10, 30, 60);
			Iterator<Integer> smoothingIterator = smoothingLevels.iterator();
			while (smoothingIterator.hasNext()) {
				Integer smoothing = smoothingIterator.next();
				IJ.log("Using smoothing level: " + smoothing);
				ip = new ColorProcessor(imp.getWidth(), imp.getHeight());
				ip.setColor(Color.white);
				ip.fill();
				String[] colors = { "red", "blue", "green", "orange" };
				int i = 0;
				for (ListIterator<Track> tracksIterator = allTracks.listIterator(); tracksIterator.hasNext();) {
					Track track = tracksIterator.next();
					track.distance = 0f;
					IJ.log("Track length: " + track.length());
					ListIterator<particle> particleIterator = track.particles.listIterator();
					particle oldParticle = particleIterator.next();
					ip.setColor(colors[i % colors.length]);
					int j = 1; // Start at 1 bc of the initial set
					while (particleIterator.hasNext()) {
						particle newParticle = particleIterator.next();
						if ((j % smoothing) == 0) {
							track.distance += oldParticle.distance(newParticle);
							ip.moveTo((int) (oldParticle.x / pixelWidth), (int) (oldParticle.y / pixelHeight));
							ip.lineTo((int) (newParticle.x / pixelWidth), (int) (newParticle.y / pixelHeight));
							oldParticle = newParticle;
						}
						j++;
					}
					particle first = track.particles.get(0);
					float offset = 10f;
					ip.drawString("Track Distance: " + track.distance, Math.round(first.x + offset),
							Math.round(first.y + offset));
					i++;
				}
				new ImagePlus("Stitch distance: " + maxStitchDistance + ", Smoothing: " + smoothing, ip).show();
			}
		}

		if (bShowLabels) {
			// ImageStack newStack = imp.createEmptyStack();
			// int yWidth = newStack.getWidth();
			// for (int i = 0; i <= (nFrames - 1); i++) {
			// if (IJ.escapePressed()) {
			// IJ.beep();
			// done = true;
			// return;
			// }
			// int iFrame = i + 1;
			// ImageProcessor imgProcessor = stack.getProcessor(iFrame);
			// newStack.addSlice(stack.getSliceLabel(iFrame), imgProcessor.crop());
			// ImageProcessor newImgProcessor = newStack.getProcessor(iFrame);
			// newImgProcessor.setFont(new Font("SansSerif", Font.PLAIN, fontSize));
			// newImgProcessor.setColor(Color.pink);

			// Iterator<Track> trackIterator = allTracks.iterator();
			// while (trackIterator.hasNext()) {
			// Track track = trackIterator.next();

			// }
			// IJ.showProgress((double) iFrame / nFrames);
			// }
			// ImagePlus nimp = new ImagePlus(imp.getTitle() + " labels", newStack);
			// nimp.setCalibration(imp.getCalibration());
			// nimp.show();
			// imp.show();
			// nimp.updateAndDraw();
		}

		// TODO: era end

		// Initiate the writing of an output textfile if a filename is given.

		boolean writefile = false;
		if (filename != null) {
			File outputfile = new File(directory, filename);
			if (!outputfile.canWrite()) {
				try {
					outputfile.createNewFile();
				} catch (IOException e) {
					IJ.showMessage("Error", "Could not create " + directory + filename);
				}
			}
			if (outputfile.canWrite())
				writefile = true;
			else
				IJ.showMessage("Error", "Could not write to " + directory + filename);
		}

		// Calculate length of paths, area of objects and number of body bends

		double[][] angles = new double[trackCount + 1][nFrames];
		double[][] dAngles = new double[trackCount + 1][nFrames];
		double[][] sumDAngles = new double[trackCount + 1][nFrames];
		double[][] bendCounter = new double[trackCount + 1][nFrames];
		double[][] dAAreas = new double[trackCount + 1][nFrames];
		double[] times = new double[nFrames];
		double[][] bendTimes = new double[trackCount + 1][nFrames];
		double[] rawX = new double[nFrames + 2];
		double[] rawY = new double[nFrames + 2];
		double[] smoothX = new double[nFrames + 2];
		double[] smoothY = new double[nFrames + 2];

		if (bShowPathLengths) {
			double[] lengths = new double[trackCount];
			double[] distances = new double[trackCount];
			double[] maxspeeds = new double[trackCount];
			double[] areas = new double[trackCount];
			double[] areaStdev = new double[trackCount];
			double[] perims = new double[trackCount];
			double[] perimsStdev = new double[trackCount];
			double[] bends = new double[trackCount];
			double[] avgX = new double[trackCount];
			double[] avgY = new double[trackCount];
			int[] frames = new int[trackCount];
			int[] firstFrames = new int[trackCount];
			int binNum = 0;

			if (binSize > 0)
				binNum = (int) (maxVelocity / binSize + 1);
			else
				binNum = 100;

			/* double x1, y1, x2, y2, */
			double distance, deltaAngle, oldAngle, sumDAngle, oldDAngle1, temp, temp2, avgArea, sumAreaSq,
					avgPerim, sumPerimSq, sumX, sumY;

			// generate smoothed coordinates to eliminate some digital noise -
			// without coordinate smoothing length of track for immobile particles can
			// actually be quite long
			// for now only 5-point smoothing is possible, but the number of points should
			// be user definable in later releases

			if (bSmoothing) {
				IJ.showStatus("Smoothing tracks...");
				IJ.log("Smoothing tracks...");
				for (ListIterator<Track> iT = tracks.listIterator(); iT.hasNext();) {
					Track track = iT.next();
					List bTrack = track.particles;
					if (bTrack.size() >= minTrackLength) {
						ListIterator jT = bTrack.listIterator();
						particle oldParticle = (particle) jT.next();
						particle oldParticle2 = oldParticle;
						rawX[1] = oldParticle.x;
						rawY[1] = oldParticle.y;
						smoothX[1] = oldParticle.x;
						smoothY[1] = oldParticle.y;
						oldParticle.sx = oldParticle.x;
						oldParticle.sy = oldParticle.y;
						int i = 1;
						int t = 0;
						for (; jT.hasNext();) {
							i++;
							particle newParticle = (particle) jT.next();
							rawX[i] = newParticle.x;
							rawY[i] = newParticle.y;
							newParticle.sx = newParticle.x;
							newParticle.sy = newParticle.y;
							smoothX[i] = newParticle.x;
							smoothY[i] = newParticle.y;
							if (i >= 3) {
								temp = 0;
								temp2 = 0;
								for (t = 0; t < 3; t++) {
									temp += rawX[t + i - 2];
									temp2 += rawY[t + i - 2];
								}
								smoothX[i - 1] = (float) temp / 3;
								smoothY[i - 1] = (float) temp2 / 3;
								oldParticle.sx = (float) temp / 3;
								oldParticle.sy = (float) temp2 / 3;

							}
							if (i >= 5) {
								temp = 0;
								temp2 = 0;
								for (t = 0; t < 5; t++) {
									temp += rawX[t + i - 4];
									temp2 += rawY[t + i - 4];
								}
								smoothX[i - 2] = (float) temp / 5;
								smoothY[i - 2] = (float) temp2 / 5;
								oldParticle2.sx = (float) temp / 5;
								oldParticle2.sy = (float) temp2 / 5;

							} else {
							}
							oldParticle2 = oldParticle;
							oldParticle = newParticle;
						}
					}
				}
			}

			int[][] bins = new int[trackCount][(int) binNum]; // used to store speed histograms
			int[][] bBins = new int[trackCount][101]; // used to store bend-histograms
			int[] FrameTrackCount = new int[nFrames]; // store how many objects are tracked in each frame

			// analyze tracks
			// Discard tracks that are smaller than minTrackLength
			// Calculate length of track
			// Calculate bends
			// Calculate how many animals are tracked on each frame

			IJ.showStatus("Analyzing tracks...");
			IJ.log("Analyzing tracks...");

			int displayTrackNum = 0;
			int bendCount = 0;
			int oldBendFrame = 0;
			for (ListIterator<Track> iT = tracks.listIterator(); iT.hasNext();) {
				Track track = iT.next();
				List bTrack = track.particles;

				if (bTrack.size() < minTrackLength)
					continue;

				displayTrackNum++;
				maxspeeds[displayTrackNum - 1] = 0;
				ListIterator jT = bTrack.listIterator();
				particle oldParticle = (particle) jT.next();
				particle firstParticle = new particle();
				FrameTrackCount[firstParticle.z]++;
				firstParticle.copy(oldParticle);

				avgArea = oldParticle.area;
				avgPerim = oldParticle.perimeter;
				sumPerimSq = 0;
				sumAreaSq = 0;
				sumX = oldParticle.x;
				sumY = oldParticle.y;
				firstFrames[displayTrackNum - 1] = oldParticle.z;
				oldAngle = sumDAngle = oldDAngle1 = 0;

				if (bendType == 1)
					oldAngle = oldParticle.angle;
				if (bendType > 1)
					oldAngle = oldParticle.ar;
				angles[displayTrackNum][0] = oldAngle;
				bendCount = 0;
				int i = 1;
				sumDAngle = oldDAngle1 = 0;
				frames[displayTrackNum - 1] = bTrack.size();

				for (; jT.hasNext();) {
					i++;
					particle newParticle = (particle) jT.next();
					FrameTrackCount[newParticle.z]++; // add one object to FrameCount at the given frame
					if (bSmoothing) {
						distance = newParticle.sdistance(oldParticle);
					} else
						distance = newParticle.distance(oldParticle);

					// Calculate average and standard deviation of object area
					temp = (avgArea + (newParticle.area - avgArea) / i);
					sumAreaSq += (newParticle.area - avgArea) * (newParticle.area - temp);
					avgArea = temp;

					// Calculate average and standard deviation of object perimeter
					temp = (avgPerim + (newParticle.perimeter - avgPerim) / i);
					sumPerimSq += (newParticle.perimeter - avgPerim) * (newParticle.perimeter - temp);
					avgPerim = temp;

					// Calculate the average position of the animals - e.g. useful if a single movie
					// cointains for seprate wells
					sumX += newParticle.x;
					sumY += newParticle.y;

					if (bendType > 0) {

						// detect wobbling of particle or bending of objects

						if (bendType == 1) {
							// use the angle of the ellipse for the fitting
							deltaAngle = newParticle.angle - oldAngle;
							oldAngle = newParticle.angle;

							if (deltaAngle > 100) {
								deltaAngle = deltaAngle - 180;
							}
							; // motion probably transverse 0-180 degree boundery
							if (deltaAngle < -100) {
								deltaAngle = deltaAngle + 180;
							}
							; // ....

						} else if (bendType > 1) {
							// use the aspect ratio of the ellipse for counting body-bends
							deltaAngle = newParticle.ar - oldAngle;
							oldAngle = newParticle.ar;
						} else
							deltaAngle = 0;
						angles[displayTrackNum][newParticle.z] = oldAngle;

						if (oldDAngle1 * deltaAngle <= 0) {
							// A change in sign indicates object started turning the other direction
							if (sumDAngle > minAngle) {
								// if area of last turn angle was over threshold degrees then count as a bend.
								binNum = newParticle.z - oldBendFrame;
								// IJ.log("binNum:"+binNum);
								if (binNum < 1)
									binNum = newParticle.z;
								if (binNum > 100)
									binNum = 100;
								bBins[displayTrackNum - 1][binNum]++;
								dAAreas[displayTrackNum][bendCount] = sumDAngle;
								bendTimes[displayTrackNum][bendCount] = newParticle.z;
								bendCount++;
								oldBendFrame = newParticle.z;
							}
							;
							sumDAngle = deltaAngle; // start summing new motion
						} else {
							sumDAngle += deltaAngle;
						}
						bendCounter[displayTrackNum][newParticle.z] = bendCount;
						sumDAngles[displayTrackNum][newParticle.z] = sumDAngle;

						oldDAngle1 = deltaAngle;

						dAngles[displayTrackNum][newParticle.z] = deltaAngle;
						times[newParticle.z] = newParticle.z;

					}
					// Generate a speed histogram
					if (!(newParticle.flag || oldParticle.flag)) { // only accecpt maximum speeds for non-flagged particles
						if (binSize > 0) {
							int binNr = (int) (distance / pixelWidth / binSize);
							bins[displayTrackNum - 1][binNr]++;
						}
						if (distance > maxspeeds[displayTrackNum - 1]) {
							maxspeeds[displayTrackNum - 1] = distance;
						}
						;
					}
					lengths[displayTrackNum - 1] += distance;
					oldParticle = newParticle;
				}
				distances[displayTrackNum - 1] = oldParticle.distance(firstParticle);// Math.sqrt(sqr(oldParticle.x-firstParticle.x)+sqr(oldParticle.y-firstParticle.y));
				areas[displayTrackNum - 1] = avgArea;
				areaStdev[displayTrackNum - 1] = Math.sqrt(sumAreaSq / (bTrack.size() - 1)); // Calculate Stdev of the areas
				perims[displayTrackNum - 1] = avgPerim;
				perimsStdev[displayTrackNum - 1] = Math.sqrt(sumPerimSq / (bTrack.size() - 1)); // Calculate the Stdev of the
																																												// perimeters
				avgX[displayTrackNum - 1] = sumX / bTrack.size();
				avgY[displayTrackNum - 1] = sumY / bTrack.size();

				if (bendType == 1)
					bends[displayTrackNum - 1] = (float) bendCount;
				else if (bendType > 1)
					bends[displayTrackNum - 1] = (float) bendCount / 2;
			}

			if (writefile) {
				try {
					File outputfile = new File(directory, filename);
					BufferedWriter dos = new BufferedWriter(new FileWriter(outputfile)); // append
					if (bendType > 0) {
						dos.write(
								"Track \tLength\tDistance\t#Frames\t1stFrame\ttime(s)\tMaxSpeed\tAvgArea\tStdArea\tAvgPerim\tStdPerim\tAvgSpeed\tBLPS\tavgX\tavgX\tBends\tBBPS\tTEST");
					} else {
						dos.write(
								"Track \tLength\tDistance\t#Frames\t1stFrame\ttime(s)\tMaxSpeed\tAvgArea\tStdArea\tAvgPerim\tStdPerim\tAvgSpeed\tBLPS\tavgX\tavgX");
					}

					dos.newLine();
					for (int i = 0; i < displayTrackNum; i++) {
						String str = "" + (i + 1) + "\t" + (float) lengths[i] + "\t"
								+ (float) distances[i] + "\t"
								+ (int) frames[i] + "\t"
								+ (int) firstFrames[i] + "\t"
								+ (float) (frames[i] / fps) + "\t"
								+ (float) fps * maxspeeds[i] + "\t"
								+ (float) areas[i] + "\t"
								+ (float) areaStdev[i] + "\t"
								+ (float) perims[i] + "\t"
								+ (float) perimsStdev[i] + "\t"
								+ (float) (fps * lengths[i] / frames[i]) + "\t"
								+ (float) (fps * 2 * lengths[i] / perims[i] / frames[i]) + "\t"
								+ (float) avgX[i] + "\t"
								+ (float) avgY[i];

						if (bendType > 0) {
							str += "\t"
									+ (float) bends[i] + "\t"
									+ (float) fps * bends[i] / frames[i];
						}

						dos.write(str);
						dos.newLine();
					}
					if (binSize > 0) {
						dos.newLine();
						dos.write("Bin\t");
						for (int t = 0; t < displayTrackNum; t++) {
							dos.write("T#" + (int) (t + 1) + "\t");
						}
						;
						dos.newLine();

						for (int i = 0; i < (int) (maxVelocity / binSize); i++) {
							String str = "" + (float) i * binSize + "\t";
							for (int t = 0; t < displayTrackNum; t++) {
								str = str + bins[t][i] + "\t";
							}
							;
							dos.write(str);
							dos.newLine();
						}
					}
					if (bendType > 2) {
						dos.newLine();
						String str = "Bin";
						for (int t = 0; t < displayTrackNum; t++) {
							str = str + ("\tTrack" + (int) (t + 1));
						}
						;
						dos.write(str);
						dos.newLine();
						for (int i = 0; i < 101; i++) {
							str = "" + (float) i;
							for (int t = 0; t < displayTrackNum; t++) {
								str += "\t" + bBins[t][i];
							}
							;
							dos.write(str);
							dos.newLine();
						}

					}
					dos.close();
				} catch (IOException e) {
					if (filename != null)
						IJ.error("An error occurred writing the file. \n \n " + e);
				}

			} else {
				ResultsTable mrt = ResultsTable.getResultsTable();
				mrt.reset();

				for (int i = 0; i < displayTrackNum; i++) {
					mrt.incrementCounter();
					mrt.addValue("Track", (i + 1));
					mrt.addValue("Length", (float) lengths[i]);
					mrt.addValue("Distance", (float) distances[i]);
					mrt.addValue("#Frames", (float) frames[i]);
					mrt.addValue("1stFrame", (float) firstFrames[i]);
					mrt.addValue("Time(s)", (float) frames[i] / fps);
					mrt.addValue("MaxSpeed", (float) fps * maxspeeds[i]);
					mrt.addValue("Area", (float) areas[i]);
					mrt.addValue("sdArea", (float) areaStdev[i]);
					mrt.addValue("Perim", (float) perims[i]);
					mrt.addValue("sdPerim", (float) perimsStdev[i]);
					mrt.addValue("avgSpeed", (float) (fps * lengths[i] / frames[i]));
					mrt.addValue("BLPS", (float) (fps * 2 * lengths[i] / perims[i] / frames[i]));
					mrt.addValue("avgX", (float) avgX[i]);
					mrt.addValue("avgY", (float) avgY[i]);
					if (bendType > 0) {
						mrt.addValue("Bends", (float) bends[i]);
						mrt.addValue("BBPS", (float) fps * bends[i] / frames[i]);
						mrt.addValue("TEST", "t");
					}
				}
				mrt.show("Results");

				if (binSize > 0) {

					// Write histogram for individual track in the movie
					IJ.write("");
					String str = "Bin\t";
					for (int t = 0; t < displayTrackNum; t++) {
						str = str + ("Track" + (int) (t + 1) + "\t");
					}
					;
					str = str + "\n";
					IJ.write(str);

					for (int i = 0; i < (int) (maxVelocity / binSize); i++) {
						str = "" + (float) i * binSize + "\t";
						for (int t = 0; t < displayTrackNum; t++) {
							str = str + bins[t][i] + "\t";
						}
						;
						IJ.write(str + "\n");

					}
				}
				// thrashing histogram
				if (bendType > 2) {
					IJ.write("");
					String str = "#Frames";
					for (int t = 0; t < displayTrackNum; t++) {
						str = str + ("\tTrack" + (int) (t + 1));
					}
					;
					IJ.write(str + "\n");
					for (int i = 0; i < 101; i++) {
						str = "" + (float) i;
						for (int t = 0; t < displayTrackNum; t++) {
							str += "\t" + bBins[t][i];
						}
						;
						IJ.write(str /* +(int)mySum */ + "\n");
					}

				}

			}
			// summarize the speed histogram for the movie

			// Generate summarized output
			// filename,N,nTracks,TotLength,totFrames,TotTime,AvgSpeed,AvgArea,AvgPerim,StdSpeed,StdArea,StdPerim
			if (bShowSummary) {
				float sumLengths = 0;
				float sumFrames = 0;
				avgArea = 0;
				sumAreaSq = 0;
				avgPerim = 0;
				sumPerimSq = 0;
				double avgSpeed = 0;
				double sumSpeedSq = 0;
				double speed = 0;
				double sumBends = 0;
				double avgBBPS = 0;
				double sumBBPSSq = 0;
				double BBPS = 0;

				for (int i = 0; i < displayTrackNum; i++) {
					sumLengths += lengths[i];
					sumFrames += frames[i];
					sumBends += bends[i];

					// Calculate average and standard deviation of object area
					temp = (avgArea + (areas[i] - avgArea) / (i + 1));
					sumAreaSq += (areas[i] - avgArea) * (areas[i] - temp);
					avgArea = temp;

					// Calculate average and standard deviation of object perimeter
					temp = (avgPerim + (perims[i] - avgPerim) / (i + 1));
					sumPerimSq += (perims[i] - avgPerim) * (perims[i] - temp);
					avgPerim = temp;

					// Calculate average and standard deviation of object Speed
					speed = lengths[i] / (frames[i] / fps);
					temp = (avgSpeed + (speed - avgSpeed) / (i + 1));
					sumSpeedSq += (speed - avgSpeed) * (speed - temp);
					avgSpeed = temp;

					// Calculate average and standard deviation of BodyBends per seconds
					BBPS = bends[i] / (frames[i] / fps);
					temp = (avgBBPS + (BBPS - avgBBPS) / (i + 1));
					sumBBPSSq += (BBPS - avgBBPS) * (BBPS - temp);
					avgBBPS = temp;

				}
				String aLine = null;
				summaryHdr = "File\tnObjMax\tnObjMaxFrm\tnObjMin\tnObjMinFrm\tkpMax\tkpMaxFrm\tnFrames\tnTracks\ttotLength\tObjFrames\tObjSeconds\tavgSpeed\tavgArea\tavgPerim\tstdSpeed\tstdArea\tstdPerim";
				if (bendType > 0)
					summaryHdr += "\tBends\tavgBBPS\tstdBBPS";

				aLine = rawFilename + "\t" + (int) nMax // number of objects (N-value)
						+ "\t" + (int) nMaxFrm // frame at nMax (KP)
						+ "\t" + (int) nMin // nMin (KP)
						+ "\t" + (int) nMinFrm // frame at nMin (KP)
						+ "\t" + (int) kpMax // kpMax objects (KP)
						+ "\t" + (int) kpFrm // frame at kpMax objects (KP)
						+ "\t" + (int) nFrames // number of frames
						+ "\t" + (int) displayTrackNum // number of tracks
						+ "\t" + (float) sumLengths // total distance covered by all objects
						+ "\t" + (float) sumFrames // total number of object*frames
						+ "\t" + (float) sumFrames / fps // total obj*seconds
						+ "\t" + (float) avgSpeed // average speed (pixels/seconds)
						+ "\t" + (float) avgArea // average worm area
						+ "\t" + (float) avgPerim // average worm perimeter
						+ "\t" + (float) Math.sqrt(sumSpeedSq / (displayTrackNum - 1))// average speed (pixels/seconds)
						+ "\t" + (float) Math.sqrt(sumAreaSq / (displayTrackNum - 1))// average worm area
						+ "\t" + (float) Math.sqrt(sumPerimSq / (displayTrackNum - 1)) // average worm perimeter
				;
				if (bendType > 0) {
					aLine += "\t" + (float) sumBends // total number of body bends in the movie
							+ "\t" + (float) avgBBPS // average BBPS per track
							+ "\t" + (float) Math.sqrt(sumBBPSSq / (displayTrackNum - 1)) // standard deviation of BBPS per track
					;

				}

				// IJ.log(aLine+"\n");
				Frame frame = WindowManager.getFrame("Summary");
				if (frame != null && (frame instanceof TextWindow) && summaryHdr.equals(prevHdr))
					tw = (TextWindow) frame;

				if (tw == null) {
					String title = "Summary";
					tw = new TextWindow(title, summaryHdr, aLine, 450, 300);
					prevHdr = summaryHdr;
				} else
					tw.append(aLine);

			}
		}

		// Now do the fancy stuff when requested:

		// makes a new stack with objects labeled with track nr
		// optionally also displays centroid position
		IJ.showStatus("Generating movie with labels...");

		if (bShowLabels) {
			String strPart;
			ImageStack newstack = imp.createEmptyStack();
			int yWidth = newstack.getWidth();

			for (int i = 0; i <= (nFrames - 1); i++) {
				if (IJ.escapePressed()) {
					IJ.beep();
					done = true;
					return;
				}
				int iFrame = i + 1;
				ip = stack.getProcessor(iFrame);
				newstack.addSlice(stack.getSliceLabel(iFrame), ip.crop());
				ImageProcessor nip = newstack.getProcessor(iFrame);
				nip.setColor(Color.pink);
				// Font f1 = Font.getFont(Font.FACE_SYSTEM, Font.STYLE_PLAIN, Font.SIZE_LARGE);
				// fontSize
				nip.setFont(new Font("SansSerif", Font.PLAIN, fontSize));
				// nip.boldFont = true;
				// hack to only show tracks longerthan minTrackLength
				int trackNum = 0;
				int displayTrackNum = 0;
				for (ListIterator<Track> iT = tracks.listIterator(); iT.hasNext();) {
					Track track = iT.next();
					List bTrack = track.particles;
					trackNum++;
					if (bTrack.size() >= minTrackLength) {
						displayTrackNum++;
						for (ListIterator k = particlesInFrame[i].listIterator(); k.hasNext();) {
							particle particle = (particle) k.next();
							if (particle.trackNum == trackNum) {
								strPart = "" + displayTrackNum;
								if (bShowPositions) {
									// if (bendType==0) strPart += "= "+(int)particle.x+","+(int)particle.y ; //
									// bend detection not enabled plot coordinates instead
									if (bendType == 1)
										strPart += "= b" + bendCounter[displayTrackNum][i];
									if (bendType > 1)
										strPart += "= b" + bendCounter[displayTrackNum][i] / 2;

									// we could plot a number of different infos for the particles tracked
									// {strPart+="="+/*(int)particle.angle+","*/+(int)particle.area+","+(int)particle.x+","+(int)particle.y;}
									// {strPart+="="+/*(int)particle.angle+","+*/(int)particle.area+","+bendCounter[displayTrackNum][i]/2;
									// };
								}
								// we could do someboundary testing here to place the labels better when we are
								// close to the edge
								nip.moveTo((int) (particle.x / pixelWidth + 0),
										doOffset((int) (particle.y / pixelHeight), yWidth, 5));
								// nip.moveTo(doOffset((int)particle.x,xHeight,5),doOffset((int)particle.y,yWidth,5)
								// );
								nip.drawString(strPart);
								if (bendType == 1) {
									nip.setColor(Color.green);
									double rad = Math.toRadians(particle.angle + 90);
									nip.moveTo((int) (particle.x / pixelWidth + 20 * Math.sin(rad)),
											(int) (particle.y / pixelHeight + 20 * Math.cos(rad)));
									nip.lineTo((int) (particle.x / pixelWidth - 20 * Math.sin(rad)),
											(int) (particle.y / pixelHeight - 20 * Math.cos(rad)));
								}
							}
						}
					}
				}
				IJ.showProgress((double) iFrame / nFrames);
			}
			ImagePlus nimp = new ImagePlus(imp.getTitle() + " labels", newstack);
			nimp.setCalibration(imp.getCalibration());
			nimp.show();
			imp.show();
			nimp.updateAndDraw();
		}
	}

	// Utility functions

	List<Track> getTracksStartedInFrame(List<Track> tracks, int frame) {
		List<Track> out = new ArrayList<>();
		for (Track track : tracks) {
			if (track.startFrame == frame)
				out.add(track);
		}
		return out;
	}

	particle getBestCandidate(List<particle> candidates, List<Float> distances) {
		particle best = null;
		float bestDistance = Float.POSITIVE_INFINITY;
		for (int i = 0; i < candidates.size(); i++) {
			particle candidate = candidates.get(i);
			float distance = distances.get(i);
			if (distance < bestDistance) {
				best = candidate;
				bestDistance = distance;
			}
		}
		return best;
	}

	double sqr(double n) {
		return n * n;
	}

	int doOffset(int center, int maxSize, int displacement) {
		if ((center - displacement) < 2 * displacement) {
			return (center + 4 * displacement);
		} else {
			return (center - displacement);
		}
	}

	@SuppressWarnings("removal")
	double s2d(String s) {
		Double d;
		try {
			d = new Double(s);
		} catch (NumberFormatException e) {
			d = null;
		}
		if (d != null)
			return (d.doubleValue());
		else
			return (0.0);
	}

}
