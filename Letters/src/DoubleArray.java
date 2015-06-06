import java.util.Arrays;

/**
 * Class that acts like double[], but has hashCode(), equals() and toString()
 * methods overriden.
 */
public class DoubleArray
{
	private double[] arr;

	public DoubleArray(double[] arr)
	{
		this.arr = arr;
	}

	public double[] getArr()
	{
		return arr;
	}

	public void setArr(double[] arr)
	{
		this.arr = arr;
	}

	@Override
	public int hashCode()
	{
		final int prime = 31;
		int result = 1;
		result = prime * result + Arrays.hashCode(arr);
		return result;
	}

	@Override
	public boolean equals(Object obj)
	{
		if (this == obj)
		{
			return true;
		}
		if (obj == null)
		{
			return false;
		}
		if (!(obj instanceof DoubleArray))
		{
			return false;
		}
		DoubleArray other = (DoubleArray) obj;
		if (!Arrays.equals(arr, other.arr))
		{
			return false;
		}
		return true;
	}

	@Override
	public String toString()
	{
		return Arrays.toString(arr);
	}

}
