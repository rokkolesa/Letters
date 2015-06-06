import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.Set;
import java.util.stream.Collectors;

import edu.stanford.math.plex4.api.Plex4;
import edu.stanford.math.plex4.homology.barcodes.AnnotatedBarcodeCollection;
import edu.stanford.math.plex4.homology.barcodes.BarcodeCollection;
import edu.stanford.math.plex4.homology.barcodes.Interval;
import edu.stanford.math.plex4.homology.chain_basis.Simplex;
import edu.stanford.math.plex4.homology.interfaces.AbstractPersistenceAlgorithm;
import edu.stanford.math.plex4.streams.impl.VietorisRipsStream;

public class Letters
{

	public static void main(String[] args)
	{
		File[] listFiles = new File("./complexes").listFiles();
		for (File file : listFiles)
		{
			String fileName = file.getName();
			if (fileName.endsWith(".out"))
			{
				System.out.println("File " + fileName);
				double[][] points = parseData(file);
				classifyLetter(points);

				System.out.println("---------------------------------------------------------------");
			}
		}
	}

	private static int[] computeHomology(double[][] points)
	{
		int dimension = 3;
		int radius = 12;
		int divisions = 1;

		VietorisRipsStream<double[]> stream = Plex4.createVietorisRipsStream(points, dimension, radius, divisions);
		AbstractPersistenceAlgorithm<Simplex> persistence = Plex4.getModularSimplicialAlgorithm(3, 2);
		BarcodeCollection<Double> circleIntervals = persistence.computeIntervals(stream);

		AnnotatedBarcodeCollection<Double, Object> intervals = circleIntervals.getInfiniteIntervals();
		Iterator<Entry<Integer, List<Interval<Double>>>> intervalIterator = intervals.getIntervalIterator();

		int numOfComponents = intervalIterator.hasNext() ? intervalIterator.next().getValue().size() : 0;
		int numOfCycles = intervalIterator.hasNext() ? intervalIterator.next().getValue().size() : 0;

		return new int[] { numOfComponents, numOfCycles };
	}

	private static void classifyLetter(double[][] points)
	{
		// število rezov
		int cuts = 10;

		for (int angle = 0; angle < 360; angle += 90)
		{
			getHomologies(points, cuts, angle);
		}

		// TODO DARKO - klasifikacija črk glede na to, kaj dobiš iz
		// getHomologies glede na nek kot
	}

	private static int[][] getHomologies(double[][] points, int cuts, int angle)
	{
		double[][] rotatedPoints = rotateAndSortPoints(points, angle);

		// izračunamo reze (izhod je 2D tabela, v vsaki vrstici je število
		// komponent in ciklov za posamezen rez)
		int[][] homologies = cutAndComputeHomologies(rotatedPoints, cuts);
		System.out.println("\tRotation: " + angle + " degrees");
		for (int i = 0; i <= cuts; i++)
		{
			System.out.println("\t\tCut number: " + i);
			System.out.println("\t\tNumber of connected components: " + homologies[i][0]);
			System.out.println("\t\tNumber of cycles: " + homologies[i][1]);
			System.out.println("\t\t-------------------------------------------------------");
		}
		return homologies;
	}

	private static int[][] cutAndComputeHomologies(double[][] points, int cuts)
	{
		double max = -Double.MAX_VALUE;
		double min = Double.MAX_VALUE;
		for (double[] point : points)
		{
			if (point[1] > max)
				max = point[1];
			else if (point[1] < min)
				min = point[1];
		}
		double step = (max - min) / (cuts + 1);

		int[][] cutHomologies = new int[cuts + 1][2];

		for (int i = 1; i <= cuts; i++)
		{
			int cutNumber = 0;
			for (int j = 0; j < points.length; j++)
			{
				if (points[j][1] < max - i * step)
				{
					cutNumber = j;
					break;
				}
			}

			double[][] cutPoints = Arrays.copyOfRange(points, 0, cutNumber);
			cutHomologies[i - 1] = computeHomology(cutPoints);
		}
		cutHomologies[cuts] = computeHomology(points);

		return cutHomologies;
	}

	// ========================================================================================================================
	// =================================================UTILITY_METHODS========================================================
	// ========================================================================================================================

	private static double[][] parseData(File file)
	{
		Set<DoubleArray> pointSet = new HashSet<DoubleArray>();
		try
		{
			for (String line : Files.lines(file.toPath()).collect(Collectors.toList()))
			{
				String str = line.trim();
				if (str.startsWith("#") || str.equals(""))
				{
					continue;
				}
				str = str.replaceAll("\\(", "");
				str = str.replaceAll("\\)", "");
				str = str.replaceAll(",", "");
				String[] split = str.trim().split(" ");

				String[] point = new String[] { split[0], split[1] };
				double[] p = Arrays.stream(point).mapToDouble(string -> Double.parseDouble(string)).toArray();
				pointSet.add(new DoubleArray(p));

				point = new String[] { split[2], split[3] };
				p = Arrays.stream(point).mapToDouble(string -> Double.parseDouble(string)).toArray();
				pointSet.add(new DoubleArray(p));
			}
		}
		catch (IOException | NumberFormatException e)
		{
			throw new RuntimeException(e);
		}
		return pointSet.stream().map(a -> a.getArr()).toArray(size -> new double[size][2]);
	}

	/**
	 * Rotates <b>points</b> according to <b>angle</b> (in degrees) and sorts
	 * them according to their <b>y</b> coordinate from top to bottom. If
	 * <b>angle</b> equals 0, only sorting takes place
	 * 
	 * @param points
	 *            2D array of points (size [n][2])
	 * @param angle
	 *            an angle for which to rotate (in degrees)
	 * @return rotated and sorted points (leaves original <b>points</b> intact)
	 */
	private static double[][] rotateAndSortPoints(double[][] points, double angle)
	{
		double[][] clonePoints = clone2D(points);

		// conversion from degrees -> radians
		angle = angle / 180 * Math.PI;
		if (Double.compare(angle, 0) != 0)
		{
			double cos = Math.cos(angle);
			double sin = Math.sin(angle);
			for (int i = 0; i < clonePoints.length; i++)
			{
				double x = clonePoints[i][0];
				double y = clonePoints[i][1];
				clonePoints[i][0] = x * cos - y * sin;
				clonePoints[i][1] = x * sin + y * cos;
			}
		}
		sortPoints(clonePoints);
		return clonePoints;
	}

	/**
	 * Clones array <b>points</b>, creating a copy of size [n][2]. Leaves
	 * original <b>points</b> array intact.
	 * 
	 * @param points
	 *            the array to copy
	 * @return a copy of <b>points</b>
	 */
	private static double[][] clone2D(double[][] points)
	{
		double[][] clonePoints = new double[points.length][2];
		for (int i = 0; i < points.length; i++)
		{
			clonePoints[i] = Arrays.copyOfRange(points[i], 0, 2);
		}
		return clonePoints;
	}

	/**
	 * Sorts <b>points</b> according to their <i>y</i> coordinate from top to
	 * bottom.
	 * 
	 * @param points
	 *            the array (of size [n][2]) to be sorted
	 */
	private static void sortPoints(double[][] points)
	{
		Arrays.sort(points, (a1, a2) ->
		{
			int compareY = Double.compare(a1[1], a2[1]);
			if (compareY == 0)
			{
				return Double.compare(a1[0], a2[0]);
			}
			return -compareY;
		});
	}

	/**
	 * Prints the given array of size [n][2] to (almost) Mathematica-ready
	 * format.
	 */
	private static void print(double[][] points)
	{
		System.out.print("{");
		for (double[] point : points)
		{
			System.out.println("{" + point[0] + ", " + point[1] + "},");
		}
		System.out.print("}");
	}
}
