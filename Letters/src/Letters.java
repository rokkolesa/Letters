import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;

import edu.stanford.math.plex4.api.Plex4;
import edu.stanford.math.plex4.homology.barcodes.AnnotatedBarcodeCollection;
import edu.stanford.math.plex4.homology.barcodes.BarcodeCollection;
import edu.stanford.math.plex4.homology.barcodes.Interval;
import edu.stanford.math.plex4.homology.chain_basis.Simplex;
import edu.stanford.math.plex4.homology.interfaces.AbstractPersistenceAlgorithm;
import edu.stanford.math.plex4.streams.impl.VietorisRipsStream;

public class Letters
{
	private static double[][] points;

	public static void main(String[] args)
	{
		File[] listFiles = new File("./complexes").listFiles();
		for (File file : listFiles)
		{
			String fileName = file.getName();
			if (fileName.endsWith(".out"))
			{
				System.out.println("File " + fileName);
				ArrayList<double[]> pointsList = parseData(file);
				points = pointsList.toArray(new double[0][0]);
				sortPoints(points);
				classifyLetter(points);

				System.out.println("--------------------------------------------------");
			}
		}
	}

	private static int[] computeHomology(double[][] points)
	{
		VietorisRipsStream<double[]> stream = Plex4.createVietorisRipsStream(points, 3, 11, 1);
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
		int[] homology = computeHomology(points);
		int numOfComponents = homology[0];
		int numOfCycles = homology[1];
		System.out.println("\tNumber of connected components: " + numOfComponents);
		System.out.println("\tNumber of cycles: " + numOfCycles);
		// TODO DARKO
	}

	private static void rotateAndSortPoints(double[][] points, double angle)
	{
		double cos = Math.cos(angle);
		double sin = Math.sin(angle);
		for (int i = 0; i < points.length; i++)
		{
			double x = points[i][0];
			double y = points[i][1];
			points[i][0] = x * cos - y * sin;
			points[i][1] = x * sin + y * cos;
		}
		sortPoints(points);
	}

	private static void sortPoints(double[][] points)
	{
		Arrays.sort(points, (a1, a2) ->
		{
			int compareY = Double.compare(a1[1], a2[1]);
			if (compareY == 0)
			{
				return Double.compare(a1[0], a2[0]);
			}
			return compareY;
		});
	}

	private static ArrayList<double[]> parseData(File file)
	{
		ArrayList<double[]> points = new ArrayList<double[]>();
		try
		{
			Files.lines(file.toPath()).forEach(line ->
			{
				String str = line.trim();
				if (str.startsWith("#") || str.equals(""))
				{
					return;
				}
				str = str.replaceAll("\\(", "");
				str = str.replaceAll("\\)", "");
				str = str.replaceAll(",", "");
				String[] split = str.trim().split(" ");

				String[] point = new String[] { split[0], split[1] };
				double[] p = Arrays.stream(point).mapToDouble(string -> Double.parseDouble(string)).toArray();
				points.add(p);

				point = new String[] { split[2], split[3] };
				p = Arrays.stream(point).mapToDouble(string -> Double.parseDouble(string)).toArray();
				points.add(p);
			});
		}
		catch (IOException | NumberFormatException e)
		{
			throw new RuntimeException(e);
		}
		return points;
	}
}
