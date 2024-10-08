
import ij.plugin.filter.PlugInFilter;
import java.util.*;
import java.awt.Color;

import ij.*;
import ij.gui.*;
import ij.process.*;
import ij.measure.*;
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
	static int expectedObjects = 0;
	static float pxPerMillimeter = 7;

	static int rawData = 0;
	static double fps = 0;
	static int backSub = 0;
	static String threshMode = "Otsu";
	static int maxColumns = 75;

	// KP
	static int kpFrm = 0, kpMax = 0;
	static String row = new String();
	private int nMax = 0;
	private int nMin = 99999;
	// EOKP

	private String directory, filename, rawFilename;

	private double pixelWidth = 1.0, pixelHeight = 1.0;
	boolean done;

	// #region Particle Class
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

	// #region Track class
	public class Track {
		List<particle> particles = new ArrayList<particle>();
		int startFrame;
		int endFrame;
		float distance;
		public String color;
		public float meanArea;
		public float startToEnd;
		public float meanPxPerSec;
		public float meanMillimetersPerSec;

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

	public static double max(double[] t) {
		double maximum = t[0]; // start with the first value
		for (int i = 1; i < t.length; i++) {
			if (t[i] > maximum) {
				maximum = t[i]; // new maximum
			}
		}
		return maximum;
	}

	// #region Run function (wrapper)
	public void run(ImageProcessor ip) {
		GenericDialog gd = new GenericDialog("wrMTrck by Jesper S. Pedersen, modified by era");
		gd.addNumericField("minSize - Minimum Object Area (pixels^2): ", minSize, 0);
		gd.addNumericField("maxSize - Maximum Object Area (pixels^2): ", maxSize, 0);
		gd.addNumericField("maxVelocity - Maximum Velocity (pixels/frame):", maxVelocity, 0);
		gd.addNumericField("maxAreaChange - Maximum area change (%):", maxAreaChange, 0);
		gd.addNumericField("minTrackLength - Minimum track length (frames):", minTrackLength, 0);
		gd.addNumericField("rawData - (0=off,1=XYcord,2=Ellipse,3=AreaPerimDist,4=Ellipse+Circ,5=BendCalc):", rawData, 0);
		gd.addNumericField("FPS - frames/s (0=try to load from file):", fps, 0);
		gd.addNumericField("backSub - On-the-fly background subtraction (0=off,1=F1RB15):", backSub, 0);
		gd.addChoice("threshMode - Thresholding method (only if backSub>0)", AutoThresholder.getMethods(), "Otsu");
		gd.addNumericField("fontSize - Size of labeling font:", fontSize, 0);
		gd.addNumericField("Maximum stitch distance", 25);
		gd.addNumericField("# of frames per sample (e.g. 30 = one sample per second at 30fps)", 30);
		gd.addNumericField("# of objects (0 = unknown)", 5);
		gd.addNumericField("pixels per mm", 7);
		gd.showDialog();
		if (gd.wasCanceled())
			return;
		minSize = (int) gd.getNextNumber();
		maxSize = (int) gd.getNextNumber();
		maxVelocity = (float) gd.getNextNumber();
		maxAreaChange = (float) gd.getNextNumber();
		minTrackLength = (int) gd.getNextNumber();
		rawData = (int) gd.getNextNumber();
		fps = (float) gd.getNextNumber();
		backSub = (int) gd.getNextNumber();
		threshMode = gd.getNextChoice();
		fontSize = (int) gd.getNextNumber();
		int maxStitchDistance = (int) gd.getNextNumber();
		int framesPerSample = (int) gd.getNextNumber();
		expectedObjects = (int) gd.getNextNumber();
		pxPerMillimeter = (float) gd.getNextNumber();

		/* Extract file name of movie file opened in ImageJ */
		if (rawData < 0) {
			rawData = -rawData;
		}
		rawFilename = imp.getTitle();
		int dotPos = rawFilename.lastIndexOf('.');
		if (0 < dotPos && dotPos < rawFilename.length() - 2) {
			rawFilename = rawFilename.substring(0, dotPos);
		}
		if (done)
			return;

		track(imp, minSize, maxSize, maxVelocity, directory, filename, maxStitchDistance, framesPerSample);
	}

	// #region Track function (main)
	@SuppressWarnings({ "rawtypes", "unchecked" })
	public void track(ImagePlus imp, int minSize, int maxSize, float maxVelocity, String directory, String filename,
			int maxStitchDistance, int framesPerSample) {
		verbose = true;
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
					particle.x = sxRes[iPart];
					particle.y = syRes[iPart];

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
				}
				if (sxRes.length < nMin) {
					nMin = sxRes.length; // record the minimum number of particles found in one frame. (KP)
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

		// #region Track Stitching
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

		// Remove the shortest tracks until we get to the expected number of objects
		while (expectedObjects > 0 && allTracks.size() > expectedObjects) {
			int min = 99999;
			Integer target = null;
			for (ListIterator<Track> tracksIterator = allTracks.listIterator(); tracksIterator.hasNext();) {
				Track track = tracksIterator.next();
				Integer length = track.particles.size();
				if (length < min) {
					min = length;
					target = tracksIterator.previousIndex();
				}
			}
			if (target != null) {
				allTracks.remove((int) target);
			}
		}

		// #region Track analysis
		double secPerFrame = 1 / fps;
		double secPerSample = secPerFrame * framesPerSample;

		IJ.showStatus("Generating map of stitched tracks...");
		IJ.log("Generating map of stitched tracks...");
		Integer smoothing = framesPerSample;
		IJ.log("Using smoothing level: " + smoothing);
		ip = new ColorProcessor(imp.getWidth(), imp.getHeight());
		ip.setColor(Color.white);
		ip.fill();
		String[] colors = { "red", "blue", "green", "orange", "cyan", "pink" };
		int i = 0;
		for (ListIterator<Track> tracksIterator = allTracks.listIterator(); tracksIterator.hasNext();) {
			Track track = tracksIterator.next();
			track.distance = 0f;
			float totalArea = 0f;
			float totalPxPerSec = 0f;
			IJ.log("Track length: " + track.length());
			ListIterator<particle> particleIterator = track.particles.listIterator();
			particle oldParticle = particleIterator.next();
			String color = colors[i % colors.length];
			track.color = color;
			ip.setColor(color);
			int j = 1; // Start at 1 bc of the initial set
			int lastIndex = 0;
			while (particleIterator.hasNext()) {
				particle newParticle = particleIterator.next();
				if ((j % smoothing) == 0) {
					lastIndex = j;
					totalArea += newParticle.area;
					float distance = oldParticle.distance(newParticle);
					float pxPerSec = distance / (float) secPerSample;
					totalPxPerSec += pxPerSec;
					track.distance += distance;
					ip.moveTo((int) (oldParticle.x / pixelWidth), (int) (oldParticle.y / pixelHeight));
					ip.lineTo((int) (newParticle.x / pixelWidth), (int) (newParticle.y / pixelHeight));
					oldParticle = newParticle;
				}
				j++;
			}
			track.meanArea = totalArea / (track.particles.size() / smoothing);
			track.meanPxPerSec = totalPxPerSec / (track.particles.size() / smoothing);
			track.meanMillimetersPerSec = track.meanPxPerSec / pxPerMillimeter;
			particle first = track.particles.get(0);
			particle last = track.particles.get(lastIndex);

			// Draw start and end points
			ip.drawDot(Math.round(first.x), Math.round(first.y));
			ip.drawDot(Math.round(last.x), Math.round(last.y));

			track.startToEnd = first.distance(last);

			// Label the track
			float labelOffset = 10f;
			ip.drawString("0", Math.round(first.x + labelOffset), Math.round(first.y + labelOffset));
			ip.drawString(track.distance + "", Math.round(last.x + labelOffset), Math.round(last.y + labelOffset));
			i++;
		}
		new ImagePlus("Stitch distance: " + maxStitchDistance + ", Smoothing: " + smoothing, ip).show();

		// #region Output table
		ResultsTable mrt = ResultsTable.getResultsTable();
		mrt.reset();

		for (ListIterator<Track> tracksIterator = allTracks.listIterator(); tracksIterator.hasNext();) {
			mrt.incrementCounter();
			Track track = tracksIterator.next();
			mrt.addValue("Colour", track.color);
			mrt.addValue("Total distance (mm)", track.distance / pxPerMillimeter);
			mrt.addValue("Distance from start to end (mm)", track.startToEnd / pxPerMillimeter);
			mrt.addValue("Mean Speed (mm/s)", track.meanMillimetersPerSec);
			mrt.addValue("Mean Object Area (sqpx)", track.meanArea);
		}

		mrt.show("Results");
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
