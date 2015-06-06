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
		int dimension = 3;
		int radious = 12;
		int divisions = 1;
		
		VietorisRipsStream<double[]> stream = Plex4.createVietorisRipsStream(points, dimension, radious, divisions);
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
		// stevilo rezov
		int cuts = 5;
		
		// rotacije example
		// rotateAndSortPoints(points, 180);
		// rotateAndSortPoints(points, 90);
		
		// izrcunamo reze(izhod je 2D tabela, v vsaki vrstici je stevilo komponent in ciklov za posamezen rez)
		int[][] homologyes = cutAndComputeHomologies(points, cuts);
		for (int i = 0; i <= cuts; i++) {
			System.out.println("\tCut number: " + i);
			System.out.println("\tNumber of connected components: " + homologyes[i][0]);
			System.out.println("\tNumber of cycles: " + homologyes[i][1]);
			System.out.println("\t-------------------------------------------------------");
		}
		
		// klasifikacija
		// TODO DARKO
	}
	
	private static int[][] cutAndComputeHomologies(double[][] points, int cuts) {		
		double max = -Double.MAX_VALUE;
		double min = Double.MAX_VALUE;
		for (double[] point : points) {
			if (point[1] > max)
				max = point[1];
			else
				if (point[1] < min)
					min = point[1];
		}
		
		double step = (max - min) / (cuts + 1);
		
		int[][] cutHomologies = new int[cuts + 1][2];
		
		for (int i = 1; i <= cuts; i++) {
			int cutNumber = 0;
			for (int j = 0; j < points.length; j ++) {
				if (points[j][1] < max - i * step) {
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
			return -compareY;
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
